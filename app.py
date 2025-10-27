
import os
import requests
import base64
import json
from io import BytesIO
from typing import Optional, Dict, Any, List

import streamlit as st
from streamlit.components.v1 import html as st_html

from dotenv import load_dotenv
load_dotenv()


DATABASE_URL = os.environ.get("DATABASE_URL")  
USE_DB = bool(DATABASE_URL)

if USE_DB:
    import sqlalchemy as sa
    from sqlalchemy import text
    from sqlalchemy.dialects.postgresql import insert
    engine = sa.create_engine(DATABASE_URL, pool_pre_ping=True)
else:
    engine = None


import pubchempy as pcp
import re
import urllib.parse
CRYSTAL_SERVICE_URL = os.environ.get("CRYSTAL_SERVICE_URL")
try:
    LOGIN_AUTH_TOKEN = None
    try:
        LOGIN_AUTH_TOKEN = st.secrets.get("COURIER_AUTH_TOKEN") if hasattr(st, 'secrets') else None
    except Exception:
        LOGIN_AUTH_TOKEN = None
    if not LOGIN_AUTH_TOKEN:
        LOGIN_AUTH_TOKEN = os.environ.get("LOGIN_AUTH_TOKEN")
except Exception:
    LOGIN_AUTH_TOKEN = os.environ.get("LOGIN_AUTH_TOKEN")

try:
    from streamlit_login_auth_ui.widgets import __login__
    LOGIN_LIB_AVAILABLE = True
except Exception:
    __login__ = None
    LOGIN_LIB_AVAILABLE = False


def normalize_formula(formula: Optional[str]) -> Optional[str]:
    if not formula:
        return formula
    toks = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    if not toks:
        return formula
    elems = []
    for sym, cnt in toks:
        elems.append((sym, cnt))
    alkali = {"Li", "Na", "K", "Rb", "Cs", "Fr"}
    alkaline_earth = {"Be", "Mg", "Ca", "Sr", "Ba", "Ra"}
    metals_first = []
    others = []
    for sym, cnt in elems:
        if sym in alkali or sym in alkaline_earth:
            metals_first.append((sym, cnt))
        else:
            others.append((sym, cnt))
    if metals_first:
        ordered = metals_first + others
    else:
        ordered = elems
    parts = []
    for sym, cnt in ordered:
        parts.append(sym + (cnt if cnt else ""))
    return "".join(parts)


st.set_page_config(page_title="Mid Molecule Thing", layout="wide", initial_sidebar_state="expanded")


def inject_css():
    try:
        with open("assets/styles.css", "r", encoding="utf8") as fh:
            css = fh.read()
            st.components.v1.html(f"<style>{css}</style>", height=0)
    except Exception:
        fallback = """
        <style>
        .molevis-card { background: rgba(255,255,255,0.02); padding: 14px; border-radius: 12px; }
        .molevis-muted { color: rgba(255,255,255,0.6) }
        .molevis-title { font-weight:700; font-size:22px }
        .molevis-small { font-size:13px; color: rgba(255,255,255,0.7) }
        .flat-btn { background: transparent; border: 1px solid rgba(255,255,255,0.06); padding:6px 10px; border-radius:8px; color:inherit }
        </style>
        """
        st.components.v1.html(fallback, height=0)

inject_css()

try:
    qp = st.query_params
    if qp and qp.get("login"):
        st.session_state['show_login_ui'] = True
        st.experimental_set_query_params()
except Exception:
    pass


 
def create_table_if_missing():
        if not engine:
                return
        sql = """
        CREATE EXTENSION IF NOT EXISTS pg_trgm;
        CREATE TABLE IF NOT EXISTS molecules (
            cid bigint PRIMARY KEY,
            pref_name text,
            common_names text[],
            smiles text,
            inchi text,
            formula text,
            mol_weight numeric,
            sdf bytea,
            created_at timestamptz DEFAULT now()
        );
        CREATE INDEX IF NOT EXISTS idx_molecules_prefname_trgm ON molecules USING gin (pref_name gin_trgm_ops);
        """
        with engine.begin() as conn:
                conn.execute(text(sql))

def db_suggest(prefix: str, limit: int = 8) -> List[Dict[str, Any]]:
    if not engine or not prefix:
        return []
    q = text("""
    SELECT cid, pref_name, common_names
    FROM molecules
    WHERE pref_name ILIKE :pfx OR array_to_string(common_names,' ') ILIKE :pfx
    ORDER BY similarity(pref_name, :raw) DESC
    LIMIT :lim
    """)
    pfx = f"{prefix}%"
    with engine.connect() as conn:
        rows = conn.execute(q, {"pfx": pfx, "raw": prefix, "lim": limit}).fetchall()
    return [{"cid": int(r[0]), "name": r[1], "syn": r[2]} for r in rows]

def db_get(cid: int) -> Optional[Dict[str, Any]]:
    if not engine:
        return None
    q = text("SELECT cid, pref_name, common_names, smiles, inchi, formula, mol_weight, sdf FROM molecules WHERE cid=:cid")
    with engine.connect() as conn:
        r = conn.execute(q, {"cid": cid}).first()
    if not r:
        return None
    return {
        "cid": int(r[0]),
        "pref_name": r[1],
        "common_names": r[2],
        "smiles": r[3],
        "inchi": r[4],
        "formula": r[5],
        "mol_weight": float(r[6]) if r[6] is not None else None,
        "sdf": bytes(r[7]) if r[7] is not None else None
    }

def db_insert(mol: Dict[str, Any]) -> None:
    if not engine:
        return
    metadata = sa.MetaData()
    mol_table = sa.Table('molecules', metadata, autoload_with=engine)
    row = {
        "cid": int(mol["cid"]),
        "pref_name": mol.get("pref_name"),
        "common_names": mol.get("common_names") or [],
        "smiles": mol.get("smiles"),
        "inchi": mol.get("inchi"),
        "formula": mol.get("formula"),
        "mol_weight": mol.get("mol_weight"),
        "sdf": mol.get("sdf") if isinstance(mol.get("sdf"), (bytes, bytearray)) else (mol.get("sdf").encode() if mol.get("sdf") else None)
    }
    stmt = insert(mol_table).values(**row).on_conflict_do_update(index_elements=['cid'], set_=row)
    with engine.begin() as conn:
        conn.execute(stmt)


USER_SMILES_FILE = os.path.join(os.path.dirname(__file__), 'user_smiles.json')

def _load_user_smiles() -> Dict[str, List[str]]:
    try:
        if os.path.exists(USER_SMILES_FILE):
            with open(USER_SMILES_FILE, 'r', encoding='utf8') as fh:
                return json.load(fh)
    except Exception:
        pass
    return {}

def _save_user_smiles(store: Dict[str, List[str]]):
    try:
        with open(USER_SMILES_FILE, 'w', encoding='utf8') as fh:
            json.dump(store, fh, indent=2)
    except Exception:
        pass

def add_smile_for_user(username: str, smiles: str) -> bool:
    if not username or not smiles:
        return False
    store = _load_user_smiles()
    lst = store.get(username, [])
    if smiles not in lst:
        lst.append(smiles)
    store[username] = lst
    _save_user_smiles(store)
    return True

def get_smiles_for_user(username: str) -> List[str]:
    if not username:
        return []
    store = _load_user_smiles()
    return store.get(username, [])


if USE_DB:
    try:
        create_table_if_missing()
    except Exception as e:
        st.warning("DB initialization failed: " + str(e))


def pubchem_fetch_by_name(name: str):
    """Return first pubchempy Compound or None."""
    try:
        comps = pcp.get_compounds(name, 'name')
        return comps[0] if comps else None
    except Exception:
        return None

def pubchem_fetch_by_smiles(smiles: str):
    try:
        comps = pcp.get_compounds(smiles, 'smiles')
        return comps[0] if comps else None
    except Exception:
        return None

def pubchem_fetch_by_cid(cid: int):
    try:
        return pcp.Compound.from_cid(cid)
    except Exception:
        return None

def pubchem_sdf(cid: int) -> Optional[str]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    r = requests.get(url, timeout=20)
    return r.text if r.ok and r.text.strip() else None

def pubchem_2d_image(cid: int, width: int = 400, height: int = 300) -> Optional[bytes]:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG?image_size={width}x{height}"
    r = requests.get(url, timeout=20)
    return r.content if r.ok else None

def pubchem_get_synonyms(cid: int) -> List[str]:
    try:
        res = pcp.get_synonyms(cid, 'cid')
        return res[0] if res else []
    except Exception:
        return []


 


def is_probably_ionic(formula: Optional[str]) -> bool:
    if not formula:
        return False
    toks = re.findall(r'([A-Z][a-z]?)', formula)
    if not toks:
        return False
    alkali = {"Li", "Na", "K", "Rb", "Cs", "Fr"}
    alkaline_earth = {"Be", "Mg", "Ca", "Sr", "Ba", "Ra"}
    halogens = {"F", "Cl", "Br", "I", "At"}
    has_metal = any(t in alkali or t in alkaline_earth for t in toks)
    has_halogen = any(t in halogens for t in toks)
    return has_metal and has_halogen


def rdkit_2d_image(smiles: str, size=(400, 300)) -> Optional[bytes]:
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        img = Draw.MolToImage(m, size=size)
        buf = BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    except Exception:
        return None

def rdkit_generate_3d_sdf(smiles: str) -> Optional[str]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        if mol is None:
            return None
        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
        AllChem.UFFOptimizeMolecule(mol)
        mol_block = Chem.MolToMolBlock(mol)
        return mol_block
    except Exception:
        return None


def sdf_to_py3dmol_html(sdf_text: str, width: int = 700, height: int = 480) -> str:
    try:
        import py3Dmol
        view = py3Dmol.view(width=width, height=height)
        view.addModel(sdf_text, 'sdf')
        view.setStyle({'stick': {}})
        view.zoomTo()
        return view.show()
    except Exception:
        
        safe_sdf = sdf_text.replace("\\", "\\\\").replace("\n", "\\n").replace("'", "\\'")
        html = f"""
        <html>
        <head><meta charset="utf-8"><script src="https://3dmol.org/build/3Dmol-min.js"></script></head>
        <body style="margin:0;">
        <div id='gldiv' style='width:{width}px; height:{height}px; position: relative;'></div>
        <script>
        const sdf = '{safe_sdf}';
        let viewer = $3Dmol.createViewer(document.getElementById('gldiv'), {{backgroundColor: '0x0b0b0b'}});
        viewer.addModel(sdf, 'sdf');
        viewer.setStyle({{}}, {{stick:{{}}}});
        viewer.zoomTo();
        viewer.render();
        </script>
        </body>
        </html>
        """
        return html


def render_3dmol_from_text(model_text: str, fmt: str = "sdf", width: int = 700, height: int = 480) -> str:
    safe_text = model_text.replace("\\", "\\\\").replace("\n", "\\n").replace("'", "\\'")
    html = f"""
    <html>
    <head><meta charset="utf-8"><script src="https://3dmol.org/build/3Dmol-min.js"></script></head>
    <body style="margin:0;">
    <div id='gldiv' style='width:{width}px; height:{height}px; position: relative;'></div>
    <script>
    const txt = '{safe_text}';
    let viewer = $3Dmol.createViewer(document.getElementById('gldiv'), {{backgroundColor: '0x0b0b0b'}});
    try {{
        viewer.addModel(txt, '{fmt}');
        viewer.setStyle({{}}, {{stick:{{}}}});
        viewer.zoomTo();
        viewer.render();
    }} catch(e) {{
        document.getElementById('gldiv').innerHTML = '<div style="color:#fff;padding:16px">3Dmol failed to parse the model.</div>';
    }}
    </script>
    </body>
    </html>
    """
    return html


def fetch_cif_from_cod(formula: Optional[str]) -> Optional[str]:
    if not formula:
        return None
    try:
        q = urllib.parse.quote_plus(formula)
        search_url = f"https://www.crystallography.net/cod/search.html?formula={q}"
        r = requests.get(search_url, timeout=15)
        if not r.ok:
            return None
        html = r.text
        # find links to .cif files
        matches = re.findall(r'href="([^"]+\.cif)"', html)
        if not matches:
            return None
        for href in matches:
            # make absolute URL if needed
            if href.startswith("/"):
                url = f"https://www.crystallography.net{href}"
            elif href.startswith("http"):
                url = href
            else:
                url = f"https://www.crystallography.net/cod/{href}"
            try:
                r2 = requests.get(url, timeout=15)
                if r2.ok and r2.text and len(r2.text) > 100:
                    return r2.text
            except Exception:
                continue
    except Exception:
        return None
    return None


def _atomic_number_to_symbol_map():
    return {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
        11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
        19: 'K', 20: 'Ca', 26: 'Fe', 29: 'Cu', 30: 'Zn', 35: 'Br', 53: 'I'
    }


def pubchem3d_to_xyz(cid: int) -> Optional[str]:
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/JSON/?record_type=3d"
        r = requests.get(url, timeout=15)
        if not r.ok:
            return None
        j = r.json()
        candidates = []

        def scan(obj):
            if isinstance(obj, list):
                if obj and all(isinstance(i, dict) for i in obj):
                    good = True
                    for item in obj:
                        if not any(k.lower() in item for k in ('x','y','z','coord')):
                            if not any(isinstance(v, (int, float)) for v in item.values()):
                                good = False
                                break
                    if good:
                        candidates.append(obj)
                for i in obj:
                    scan(i)
            elif isinstance(obj, dict):
                for v in obj.values():
                    scan(v)

        scan(j)

        atom_map = _atomic_number_to_symbol_map()

        def parse_list(lst):
            rows = []
            for item in lst:
                if not isinstance(item, dict):
                    return None
                x = None; y = None; z = None
                for k, v in item.items():
                    kl = k.lower()
                    if kl in ('x','xcoord','x_coord') and isinstance(v, (int, float)):
                        x = v
                    if kl in ('y','ycoord','y_coord') and isinstance(v, (int, float)):
                        y = v
                    if kl in ('z','zcoord','z_coord') and isinstance(v, (int, float)):
                        z = v
                if x is None or y is None or z is None:
                    for v in item.values():
                        if isinstance(v, (list, tuple)) and len(v) >= 3 and all(isinstance(n, (int, float)) for n in v[:3]):
                            x, y, z = float(v[0]), float(v[1]), float(v[2])
                            break

                if x is None or y is None or z is None:
                    return None

                sym = None
                for key in ('element','atom','symbol','label'):
                    if key in item and isinstance(item[key], str):
                        sym = item[key].strip()
                        break
                if not sym:
                    for key in ('atomic_number','element_id','Z'):
                        if key in item:
                            try:
                                zn = int(item[key])
                                sym = atom_map.get(zn, None)
                                break
                            except Exception:
                                continue
                if not sym:
                    sym = 'X'

                rows.append((sym, float(x), float(y), float(z)))
            return rows

        parsed = None
        for cand in candidates:
            parsed = parse_list(cand)
            if parsed:
                break

        if not parsed:
            return None

        lines = [str(len(parsed)), f"PubChem CID {cid} 3D extracted"]
        for sym, x, y, z in parsed:
            lines.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
        return "\n".join(lines)
    except Exception:
        return None


 

st.title("Mid Molecule Thing")

inject_css()

try:
    from streamlit_ketcher import st_ketcher
    KETCHER_AVAILABLE = True
except Exception:
    st_ketcher = None
    KETCHER_AVAILABLE = False

tabs_placeholder, draw_tab = st.tabs(["Search", "Draw"])
with draw_tab:
    if KETCHER_AVAILABLE:
        try:
            try:
                st.markdown(
                    """
                    <style>
                    iframe[title^="streamlit-ketcher"], iframe[title^="streamlit-component"] { background: #f0f0f0 !important; }
                    .ketcher-wrapper { background: #f0f0f0; padding: 6px; border-radius: 6px; }
                    </style>
                    """,
                    unsafe_allow_html=True,
                )
            except Exception:
                pass
            try:
                smiles = st_ketcher(background_color="#f0f0f0")
            except TypeError:
                try:
                    
                    smiles = st_ketcher(theme="light", background="#f0f0f0")
                except TypeError:
                
                    smiles = st_ketcher()

            st.markdown("**SMILES from editor**")
            st.code(smiles or "(empty)")
            if smiles:
                st.session_state["search_query"] = smiles
        except Exception as e:
            st.error("Molecule editor failed to initialize: " + str(e))
    else:
        st.info("Install the official Streamlit molecule editor component `streamlit-ketcher` to enable the full editor. Falling back to a simple SMILES input.")
        fallback_smiles = st.text_area("SMILES", value="", height=360)
        if fallback_smiles:
            st.code(fallback_smiles)
            st.session_state["search_query"] = fallback_smiles

sidebar = st.sidebar
sidebar.title("Controls")

if 'show_login_ui' not in st.session_state:
    st.session_state['show_login_ui'] = False

if st.session_state.get('show_login_ui'):
    if LOGIN_LIB_AVAILABLE:
        try:
            if not LOGIN_AUTH_TOKEN:
                st.warning("Login library available but no Courier auth token configured. Password-reset emails may not work. Set COURIER_AUTH_TOKEN in secrets or LOGIN_AUTH_TOKEN in env.")
            login_obj = __login__(auth_token=LOGIN_AUTH_TOKEN or '', company_name='MidMol', width=200, height=220, logout_button_name='Logout', hide_menu_bool=True, hide_footer_bool=True)
            logged = login_obj.build_login_ui()
            if logged:
            
                st.session_state['LOGGED_IN'] = True
                if not st.session_state.get('account_username'):
                    st.session_state['account_username'] = st.session_state.get('username', '') or st.session_state.get('user', '') or ''
                sidebar.markdown("**Logged in**")
        except Exception as e:
            st.error("Login UI failed to initialize: " + str(e))
    else:
        
        st.error('Login UI library not installed. To enable accounts install `streamlit-login-auth-ui` or set up an alternative auth flow.')

if st.session_state.get('LOGGED_IN'):
    sidebar.markdown("---")
    sidebar.markdown("### Account")
    acct_user = sidebar.text_input('Account username (enter your login username)', value=st.session_state.get('account_username', ''))
    st.session_state['account_username'] = acct_user
    acct_smiles = sidebar.text_area('SMILES to save to account (paste or copy from Draw tab)', value='', key='acct_smiles')
    if sidebar.button('Save SMILES to account'):
        if not acct_user:
            sidebar.error('Please enter your account username (the one you used to sign up).')
        elif not acct_smiles.strip():
            sidebar.error('Please paste a SMILES string to save.')
        else:
            ok = add_smile_for_user(acct_user.strip(), acct_smiles.strip())
            if ok:
                sidebar.success('Saved SMILES to account.')
            else:
                sidebar.error('Failed to save SMILES.')

    saved = get_smiles_for_user(acct_user) if acct_user else []
    if saved:
        sel = sidebar.selectbox('Your saved SMILES', options=['(choose)'] + saved)
        if sel and sel != '(choose)':
            if sidebar.button('Use selected SMILES in Search'):
                st.session_state['search_query'] = sel
                st.session_state['do_search'] = True
                st.experimental_rerun()

def _trigger_search():
    st.session_state["do_search"] = True

if "search_query" not in st.session_state:
    st.session_state["search_query"] = "aspirin"
if "do_search" not in st.session_state:
    st.session_state["do_search"] = False

search_query = sidebar.text_input("Search (IUPAC, common name or SMILES)", key="search_query", on_change=_trigger_search)

selected_suggestion = None


do_search = st.button("Search") or st.session_state.get("do_search", False) or (selected_suggestion is not None)
if do_search:
    st.session_state["do_search"] = False
if do_search:
    
    target_cid = None
    if selected_suggestion:
        target_cid = selected_suggestion["cid"]

    result = None
    
    if target_cid:
        if USE_DB:
            result = db_get(target_cid)
        if result:
            result["source"] = "db"
    
    if not result:
        
        comp = None

        q = search_query.strip()
        if q.isdigit():
            comp = pubchem_fetch_by_cid(int(q))

        if not comp:
            try:
                comp = pubchem_fetch_by_smiles(q)
            except Exception:
                comp = None

        if not comp:
            comp = pubchem_fetch_by_name(q)
        if not comp:
            try:
                comps = pcp.get_compounds(q, 'smiles')
                comp = comps[0] if comps else None
            except Exception:
                comp = None
        if not comp:
            result = None
            local_smiles = q
            try:
                from rdkit import Chem
                m = Chem.MolFromSmiles(local_smiles)
            except Exception:
                m = None

            if m:
                try:
                    from rdkit.Chem import rdMolDescriptors
                    formula = rdMolDescriptors.CalcMolFormula(m)
                except Exception:
                    formula = None
                try:
                    mw = rdMolDescriptors.CalcExactMolWt(m)
                except Exception:
                    mw = None
                try:
                    inchi = Chem.MolToInchi(m)
                except Exception:
                    inchi = None

                img_bytes = rdkit_2d_image(local_smiles)

                sdf_text = rdkit_generate_3d_sdf(local_smiles)

                result = {
                    "cid": None,
                    "pref_name": None,
                    "common_names": [],
                    "smiles": local_smiles,
                    "inchi": inchi,
                    "formula": normalize_formula(formula) if formula else None,
                    "mol_weight": float(mw) if mw is not None else None,
                    "sdf": sdf_text.encode() if sdf_text else None,
                    "source": "local"
                }
            else:
                st.warning("No result found on PubChem for that query.")
                result = None
        else:
            cid_raw = getattr(comp, 'cid', None)
            if cid_raw is None:
                # PubChem returned a Compound-like object but no CID — treat as no result
                st.warning("PubChem returned a compound without a CID; skipping PubChem result.")
                result = None
            else:
                try:
                    cid = int(cid_raw)
                except Exception:
                    st.warning(f"Invalid CID returned by PubChem: {cid_raw}")
                    result = None
                else:
                    cached = db_get(cid) if USE_DB else None
                    if cached:
                        cached["source"] = "db"
                        result = cached
                    else:
                        sdf_text = pubchem_sdf(cid)
                        result = {
                            "cid": cid,
                            "pref_name": comp.iupac_name or (comp.synonyms[0] if comp.synonyms else getattr(comp, 'title', None)),
                            "common_names": comp.synonyms or [],
                            "smiles": comp.smiles or comp.smiles,
                            "inchi": comp.inchi,
                            "formula": normalize_formula(comp.molecular_formula),
                            "mol_weight": comp.molecular_weight,
                            "sdf": sdf_text.encode() if sdf_text else None,
                            "source": "pubchem"
                        }

                        if USE_DB:
                            try:
                                db_insert(result)
                            except Exception as e:
                                st.warning("DB insert failed: " + str(e))

    
    if not result:
        st.info("No result to display.")
    else:
        
        left_col, right_col = st.columns([1, 1])
        with left_col:
            st.markdown(f"### {result.get('pref_name') or 'CID ' + str(result['cid'])}")
            st.write("**CID:**", result["cid"])
            st.write("**Formula:**", result.get("formula"))
            st.write("**Mol. weight:**", result.get("mol_weight"))
            st.write("**SMILES:**", result.get("smiles"))
            st.write("**InChI:**", result.get("inchi"))
            if result.get("common_names"):
                st.write("**Synonyms:**", ", ".join(result.get("common_names")[:12]))
            
            if result.get("sdf"):
                
                sdf_bytes = result["sdf"] if isinstance(result["sdf"], (bytes, bytearray)) else (result["sdf"].encode() if result["sdf"] else b"")
                st.download_button("Download SDF", data=sdf_bytes, file_name=f"{result['cid']}.sdf", mime="chemical/x-mdl-sdfile")
            if result.get("smiles"):
                st.download_button("Download SMILES", data=result.get("smiles"), file_name=f"{result['cid']}.smi", mime="chemical/x-daylight-smiles")
            st.write(f"**Source:** {result.get('source', 'unknown')}")

        with right_col:

            img_bytes = None
            if result.get("smiles"):
                try:
                    img_bytes = rdkit_2d_image(result["smiles"])
                except Exception as e:
                    st.warning(f"RDKit failed to render SMILES: {e}")
                    img_bytes = None

            missing_image = img_bytes is None or (isinstance(img_bytes, (bytes, bytearray)) and len(img_bytes) == 0)
            if missing_image:
                try:
                    if result.get("cid"):
                        img_bytes = pubchem_2d_image(result["cid"], width=480, height=360)
                except Exception as e:
                    st.warning(f"PubChem image fetch failed: {e}")
                    img_bytes = None

            if isinstance(img_bytes, (bytes, bytearray)) and len(img_bytes) > 0:
                st.image(img_bytes, caption="2D structure", use_column_width=True)
            else:
                st.write("2D structure not available.")

            sdf_text = None
            if result.get("sdf"):
                sdf_text = result["sdf"].decode() if isinstance(result["sdf"], (bytes, bytearray)) else result["sdf"]
            else:
                
                try:
                    sdf_text = pubchem_sdf(result["cid"])
                except Exception:
                    sdf_text = None

            
            if not sdf_text and result.get("smiles"):
                rdf = rdkit_generate_3d_sdf(result["smiles"])
                if rdf:
                    sdf_text = rdf

            if sdf_text:
                viewer_html = sdf_to_py3dmol_html(sdf_text, width=700, height=480)
                st_html(viewer_html, height=520)
            else:
                
                formula = result.get("formula")
                if is_probably_ionic(formula):
                    st.write("Detected likely inorganic ionic solid — attempting to find a crystal structure automatically (COD)...")
                    cif_text = None
                    try:
                        with st.spinner("Searching COD for a matching CIF..."):
                            cif_text = fetch_cif_from_cod(result.get('formula') or "")
                    except Exception:
                        cif_text = None

                    if cif_text:
                        try:
                            viewer_html = render_3dmol_from_text(cif_text, fmt="cif", width=700, height=480)
                            st_html(viewer_html, height=520)
                        except Exception as e:
                            st.error("Found a CIF on COD but failed to render it: " + str(e))
                    else:
                        
                        pubchem_xyz = None
                        try:
                            with st.spinner("No CIF found — attempting PubChem 3D fallback..."):
                                pubchem_xyz = pubchem3d_to_xyz(result['cid'])
                        except Exception:
                            pubchem_xyz = None

                        if pubchem_xyz:
                            try:
                                viewer_html = render_3dmol_from_text(pubchem_xyz, fmt="xyz", width=700, height=480)
                                st_html(viewer_html, height=520)
                            except Exception as e:
                                st.error("PubChem 3D data found but failed to render: " + str(e))
                        else:
                            cod_search_url = f"https://www.crystallography.net/cod/search.html?formula={urllib.parse.quote_plus(result.get('formula') or '')}"
                            st.write("No CIF or PubChem 3D coordinates found automatically. You can search COD manually:")
                            st.markdown(f"[Search COD for crystal structures]({cod_search_url})")
                else:
                    st.write("3D model not available.")


st.markdown("---")
st.markdown("Built with PubChem (PUG-REST), RDKit (optional), py3Dmol, Neon/Postgres (optional).")


with st.expander("Quick tips"):
    st.write("""
    - You can search by IUPAC or common name (e.g. `aspirin`, `glucose`).
    - Compounds with ionic formulas (e.g., salts) will trigger an automatic search for crystal structures in the Crystallography Open Database (COD), but it will probably not work.
    """)

if sidebar.button('Account / Login'):
    st.session_state['show_login_ui'] = True



