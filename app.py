
import os
import requests
import base64
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
        body { background: 
        .molevis-card { background: rgba(255,255,255,0.02); padding: 14px; border-radius: 12px; }
        .molevis-muted { color: rgba(255,255,255,0.6) }
        .molevis-title { font-weight:700; font-size:22px }
        .molevis-small { font-size:13px; color: rgba(255,255,255,0.7) }
        .flat-btn { background: transparent; border: 1px solid rgba(255,255,255,0.06); padding:6px 10px; border-radius:8px; color:inherit }
        </style>
        """
        st.components.v1.html(fallback, height=0)

inject_css()


 
def create_table_if_missing():
    """Create molecules table and pg_trgm extension if DB present."""
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
    """Return list of suggestion dicts from DB matching prefix on pref_name or synonyms."""
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
    """Return molecule row from DB by CID (or None)."""
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
    """Insert or update molecule into DB. mol should include keys: cid, pref_name, common_names (list), smiles, inchi, formula, mol_weight, sdf (bytes or str)"""
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
    """Heuristic: return True if formula looks like an inorganic ionic salt.

    We treat presence of alkali/alkaline-earth metals alongside non-metals
    (e.g., halogens) as an indicator.
    """
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
    """Generate a 3D SDF text from SMILES using RDKit; returns MolBlock string or None."""
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
    """Try to use py3Dmol if installed; otherwise return HTML using 3Dmol CDN."""
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
    """Return HTML for a 3Dmol viewer with given model text and format (sdf, cif, xyz, pdb, etc.)."""
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
    """Try to find and download a CIF from the Crystallography Open Database (COD) for the given formula.

    This does a best-effort HTML search on COD and attempts to follow the first .cif link found.
    Returns CIF text or None.
    """
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
    # Minimal mapping for common elements; extend if needed
    return {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
        11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar',
        19: 'K', 20: 'Ca', 26: 'Fe', 29: 'Cu', 30: 'Zn', 35: 'Br', 53: 'I'
    }


def pubchem3d_to_xyz(cid: int) -> Optional[str]:
    """Try to obtain 3D coordinates from PubChem record JSON and convert to XYZ text.

    This is a best-effort parser: it recursively searches the returned JSON for a
    list of atom-like dicts containing numeric x/y/z keys (or parallel arrays).
    Returns an XYZ-formatted string or None.
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/JSON/?record_type=3d"
        r = requests.get(url, timeout=15)
        if not r.ok:
            return None
        j = r.json()

        # recursive search for candidate atom lists
        candidates = []

        def scan(obj):
            if isinstance(obj, list):
                # if list of dicts, check for coords
                if obj and all(isinstance(i, dict) for i in obj):
                    # check if items have numeric x/y/z keys
                    good = True
                    for item in obj:
                        if not any(k.lower() in item for k in ('x','y','z','coord')):
                            # try to find numeric keys
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

        # try to parse candidates
        atom_map = _atomic_number_to_symbol_map()

        def parse_list(lst):
            rows = []
            for item in lst:
                if not isinstance(item, dict):
                    return None
                # try common coordinate keys
                x = None; y = None; z = None
                for k, v in item.items():
                    kl = k.lower()
                    if kl in ('x','xcoord','x_coord') and isinstance(v, (int, float)):
                        x = v
                    if kl in ('y','ycoord','y_coord') and isinstance(v, (int, float)):
                        y = v
                    if kl in ('z','zcoord','z_coord') and isinstance(v, (int, float)):
                        z = v
                # fallback: look for a list/tuple of coords
                if x is None or y is None or z is None:
                    # look for numeric list under some key
                    for v in item.values():
                        if isinstance(v, (list, tuple)) and len(v) >= 3 and all(isinstance(n, (int, float)) for n in v[:3]):
                            x, y, z = float(v[0]), float(v[1]), float(v[2])
                            break

                if x is None or y is None or z is None:
                    # not an atom-like dict
                    return None

                # find element symbol
                sym = None
                for key in ('element','atom','symbol','label'):
                    if key in item and isinstance(item[key], str):
                        sym = item[key].strip()
                        break
                if not sym:
                    # try atomic number
                    for key in ('atomic_number','element_id','Z'):
                        if key in item:
                            try:
                                zn = int(item[key])
                                sym = atom_map.get(zn, None)
                                break
                            except Exception:
                                continue
                if not sym:
                    # fallback to X as 'X'
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

        # build XYZ string
        lines = [str(len(parsed)), f"PubChem CID {cid} 3D extracted"]
        for sym, x, y, z in parsed:
            lines.append(f"{sym} {x:.6f} {y:.6f} {z:.6f}")
        return "\n".join(lines)
    except Exception:
        return None


 

st.title("Mid Molecule Thing")

# Drawing tab: embed a lightweight JS molecule editor (JSME) and lookup PubChem names by SMILES.
draw_html = r"""
<!doctype html>
<html>
    <head>
        <meta charset="utf-8" />
        <!-- JSME editor (hosted) -->
        <script src="https://peter-ertl.com/jsme/JSME_2018-12-16/jsme/jsme.nocache.js"></script>
        <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap" rel="stylesheet">
        <style>
            body{background:#0b0b0b;color:#f7f7f7;font-family:Inter,Arial,Helvetica,sans-serif;margin:12px}
            .panel{padding:18px;border-radius:12px;background:linear-gradient(180deg, #0f0f0f 0%, #0b0b0b 100%);border:1px solid rgba(255,255,255,0.04);max-width:1100px}
            .toolbar{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:8px}
            .btn{background:transparent;border:1px solid rgba(255,255,255,0.08);padding:8px 12px;border-radius:8px;color:inherit;cursor:pointer}
            .btn:active{transform:translateY(1px)}
            pre.box{white-space:pre-wrap;color:#e6e6e6;background:#070707;padding:10px;border-radius:8px;border:1px solid rgba(255,255,255,0.03)}
            .muted{color:rgba(255,255,255,0.6);font-size:13px}
        </style>
    </head>
    <body>
        <div class="panel">
            <div class="toolbar">
                <button class="btn" id="get">Get SMILES & lookup PubChem</button>
                <button class="btn" id="copy">Copy SMILES</button>
                <button class="btn" id="addC">Add carbon (C)</button>
                <button class="btn" id="addDouble">Make last bond double</button>
                <button class="btn" id="addBr">Append Br</button>
                <button class="btn" id="addCl">Append Cl</button>
                <button class="btn" id="addF">Append F</button>
                <div style="flex:1"></div>
                <div class="muted">Tip: Use the built-in editor tools (select, bond, atom) for fine edits. The buttons are convenience helpers.</div>
            </div>
            <div id="applet_container" style="width:100%;height:420px;max-width:900px;margin-bottom:8px"></div>

            <div style="display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:8px">
                <div>
                    <strong>SMILES</strong>
                    <pre id="smiles" class="box">(empty)</pre>
                </div>
                <div>
                    <strong>PubChem names (top results)</strong>
                    <pre id="names" class="box">(none)</pre>
                </div>
            </div>
        </div>
        <script>
            // JSME will call window.jsmeOnLoad when ready
            var jsmeApplet = null;
            function jsmeOnLoad() {
                try{
                    jsmeApplet = new JSME('applet_container','100%','420px');
                }catch(e){ console.error('JSME load failed', e); }
            }

            function safeGetSmiles(){ try{ return jsmeApplet && jsmeApplet.smiles ? jsmeApplet.getSmiles() || '' : (jsmeApplet && jsmeApplet.getSmiles ? jsmeApplet.getSmiles()||'' : ''); }catch(e){return '';}}
            function safeSetSmiles(s){ try{ if(jsmeApplet && (jsmeApplet.setSmiles || jsmeApplet.readString)){
                        if(jsmeApplet.setSmiles) jsmeApplet.setSmiles(s);
                        else if(jsmeApplet.readString) jsmeApplet.readString(s);
                    } else if(jsmeApplet && jsmeApplet.setSmiles){ jsmeApplet.setSmiles(s); } }
                catch(e){ console.warn('setSmiles failed', e); }}

            function updateDisplay(){
                var s = '';
                try{ s = (jsmeApplet && jsmeApplet.getSmiles) ? jsmeApplet.getSmiles() : ''; }catch(e){ s=''; }
                document.getElementById('smiles').textContent = s || '(empty)';
            }

            document.addEventListener('DOMContentLoaded', function(){
                document.getElementById('get').onclick = async function(){
                    if(!jsmeApplet){ alert('Editor not ready yet.'); return; }
                    try{
                        var s = safeGetSmiles();
                        document.getElementById('smiles').textContent = s || '(empty)';
                        if(!s){ document.getElementById('names').textContent = 'No SMILES to lookup.'; return; }
                        var enc = encodeURIComponent(s);
                        var url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' + enc + '/synonyms/JSON';
                        document.getElementById('names').textContent = 'Looking up...';
                        var res = await fetch(url, {method:'GET'});
                        if(!res.ok){ document.getElementById('names').textContent = 'No PubChem entry found.'; return; }
                        var j = await res.json();
                        try{
                            var info = j.InformationList.Information[0];
                            var syns = info.Synonym || [];
                            if(syns.length===0) document.getElementById('names').textContent = 'No synonyms returned by PubChem.';
                            else document.getElementById('names').textContent = syns.slice(0,12).join('\n');
                        }catch(e){ document.getElementById('names').textContent = 'Failed to parse PubChem response.'; }
                    }catch(e){ document.getElementById('names').textContent = 'Editor error: '+String(e); }
                };

                document.getElementById('copy').onclick = function(){
                    var t = document.getElementById('smiles').textContent || '';
                    navigator.clipboard && navigator.clipboard.writeText(t);
                };

                // convenience helpers that operate on SMILES text and attempt to set it back
                function appendFragment(f){
                    try{
                        var s = safeGetSmiles() || '';
                        var ns = s + f;
                        safeSetSmiles(ns);
                        updateDisplay();
                    }catch(e){ console.warn(e); }
                }

                document.getElementById('addC').onclick = function(){ appendFragment('C'); };
                document.getElementById('addBr').onclick = function(){ appendFragment('Br'); };
                document.getElementById('addCl').onclick = function(){ appendFragment('Cl'); };
                document.getElementById('addF').onclick = function(){ appendFragment('F'); };

                document.getElementById('addDouble').onclick = function(){
                    try{
                        var s = safeGetSmiles() || '';
                        if(!s){ return; }
                        // naive: replace last occurrence of two atoms without '=' between them with a '=' inserted
                        // e.g. CC -> C=C, COC -> CO=C (best-effort)
                        var ns = s.replace(/([A-Za-z0-9@\[\]\(\)]+)([A-Za-z0-9@\[\]\(\)]+)$/, '$1=$2');
                        if(ns===s){ // fallback: just append =
                            ns = s + '=';
                        }
                        safeSetSmiles(ns);
                        updateDisplay();
                    }catch(e){ console.warn('make double failed', e); }
                };

                // update display periodically in case user edits in the applet
                setInterval(updateDisplay, 800);
            });
        </script>
    </body>
</html>
"""

tabs_placeholder, draw_tab = st.tabs(["Search", "Draw"])
with draw_tab:
        st.components.v1.html(draw_html, height=640)

sidebar = st.sidebar
sidebar.title("Controls")

def _trigger_search():
    st.session_state["do_search"] = True

if "search_query" not in st.session_state:
    st.session_state["search_query"] = "aspirin"
if "do_search" not in st.session_state:
    st.session_state["do_search"] = False

search_query = sidebar.text_input("Search (IUPAC or common name)", value=st.session_state["search_query"], key="search_query", on_change=_trigger_search)
seed_button = sidebar.button("Seed DB (run small seed)")

if seed_button:
    st.info("Seeding DB with a small selection...")
    try:
        seed_script = os.path.join("db", "populate_pubchem.py")
        if os.path.exists(seed_script):
            import subprocess, sys
            proc = subprocess.run([sys.executable, seed_script], capture_output=True, text=True)
            if proc.returncode == 0:
                st.success("Seed script finished.")
            else:
                st.error(f"Seeding failed (exit code {proc.returncode}).\n\nstdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}")
        else:
            seeds = ["glucose", "aspirin", "acetone", "benzene", "ethanol", "caffeine", "nicotine", "paracetamol", "ibuprofen", "adenine"]
            for s in seeds:
                comp = pubchem_fetch_by_name(s)
                if not comp:
                    continue
                cid = int(comp.cid)
                sdf_text = pubchem_sdf(cid)
                try:
                    db_insert({
                        "cid": cid,
                        "pref_name": comp.iupac_name or (comp.synonyms[0] if comp.synonyms else comp.title),
                        "common_names": comp.synonyms or [],
                        "smiles": comp.isomeric_smiles or comp.smiles,
                        "inchi": comp.inchi,
                        "formula": comp.molecular_formula,
                        "mol_weight": comp.molecular_weight,
                        "sdf": sdf_text.encode() if sdf_text else None
                    })
                except Exception as e:
                    st.warning(f"Insert failed for {s}: {e}")
            st.success("Inline seeding done.")
    except Exception as e:
        st.error("Seeding failed: " + str(e))

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
            comp = pubchem_fetch_by_name(search_query)
        if not comp:
            st.warning("No result found on PubChem for that query.")
            result = None
        else:
            cid = int(comp.cid)
            
            cached = db_get(cid) if USE_DB else None
            if cached:
                cached["source"] = "db"
                result = cached
            else:
                
                sdf_text = pubchem_sdf(cid)
                result = {
                    "cid": cid,
                    "pref_name": comp.iupac_name or (comp.synonyms[0] if comp.synonyms else comp.title),
                    "common_names": comp.synonyms or [],
                    "smiles": comp.isomeric_smiles or comp.smiles,
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
                img_bytes = rdkit_2d_image(result["smiles"])
            if not img_bytes:
                try:
                    img_bytes = pubchem_2d_image(result["cid"], width=480, height=360)
                except Exception:
                    img_bytes = None
            if img_bytes:
                st.image(img_bytes, caption="2D structure", use_container_width=True)
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
                        # Try PubChem 3D JSON fallback -> XYZ
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
    - Search by IUPAC or common name (e.g. `aspirin`, `glucose`).
    - If you expect better 2D/3D rendering and automatic conformer generation, install RDKit in your environment (conda recommended).
    - To enable fast suggestions & caching, set `DATABASE_URL` to a Postgres (Neon) connection string and seed the DB.
    - Tailwind CSS: build `assets/styles.css` and the app will inject it at start-up if present.
    """)

