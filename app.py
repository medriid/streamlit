
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


st.set_page_config(page_title="MoleVis", layout="wide", initial_sidebar_state="expanded")


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


def st_lottie_url(url: str, height: int = 240):
    """Embed a Lottie animation from a URL using the Lottie web component.

    This uses the CDN lottie-player so no extra Python package is required.
    """
    lottie_html = f"""
    <script src="https://unpkg.com/@lottiefiles/lottie-player@latest/dist/lottie-player.js"></script>
    <lottie-player
        src="{url}"
        background="transparent"
        speed="1"
        style="width:100%; height:{height}px;"
        loop
        autoplay>
    </lottie-player>
    """
    # use the existing alias imported above
    st_html(lottie_html, height=height)



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


# small decorative animation near the top
try:
    st_lottie_url("https://assets9.lottiefiles.com/packages/lf20_tfb3estd.json", height=160)
except Exception:
    # if animation fails for any reason, continue without breaking the app
    pass

st.title("MoleVis — Interactive Molecule Explorer")
st.markdown("**Organic-first molecule search** — search by IUPAC or common name. Results show IUPAC, SMILES, InChI, 2D image, interactive 3D viewer, and downloads.")


sidebar = st.sidebar
sidebar.title("Controls")
search_query = sidebar.text_input("Search (IUPAC or common name)", value="aspirin")
sidebar.markdown("**Theme**")
theme_choice = sidebar.selectbox("Choose theme", ["Dark (default)", "Light"], index=0)
seed_button = sidebar.button("Seed DB (run small seed)")


if theme_choice.startswith("Light"):
    
    st.components.v1.html(
        "<style>body{background:#ffffff;color:#0b0b0b} .molevis-card{background:rgba(0,0,0,0.03);color:inherit}</style>",
        height=0
    )
else:
    st.components.v1.html(
        "<style>body{background:#0b0b0b;color:#f7f7f7} .molevis-card{background:rgba(255,255,255,0.02);color:inherit}</style>",
        height=0
    )


suggestions = db_suggest(search_query) if search_query else []
selected_suggestion = None
if suggestions:
    opts = [f"{s['cid']}  —  {s['name']}" for s in suggestions]
    pick = sidebar.selectbox("Suggestions (from cache)", options=["(pick one)"] + opts)
    if pick != "(pick one)":
        idx = opts.index(pick)
        selected_suggestion = suggestions[idx]


if seed_button:
    st.info("Seeding DB with a small selection...")
    try:
        
        seed_script = os.path.join("db", "populate_pubchem.py")
        if os.path.exists(seed_script):
            import subprocess
            subprocess.run(["python", seed_script], check=True)
            st.success("Seed script finished.")
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


do_search = st.button("Search") or (selected_suggestion is not None)
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
                    "formula": comp.molecular_formula,
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

