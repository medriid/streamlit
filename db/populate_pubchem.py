import os
import requests
import re
import psycopg2
from psycopg2.extras import execute_values
from psycopg2 import sql
from tqdm import tqdm
import sys


DB_URL = os.getenv("DATABASE_URL")

if not DB_URL:
    print("ERROR: DATABASE_URL is not set. On Streamlit Cloud add DATABASE_URL as a secret pointing to your Postgres instance.")
    sys.exit(2)
if "localhost" in DB_URL or "user:password" in DB_URL:
    print("ERROR: DATABASE_URL appears to be a local placeholder. Replace it with a real Postgres connection string and add it as a Streamlit Cloud secret.")
    sys.exit(2)


def fetch_cid_for_name(name: str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.requote_uri(name)}/cids/JSON"
    r = requests.get(url, timeout=10)
    if not r.ok:
        return None
    j = r.json()
    ids = j.get("IdentifierList", {}).get("CID", [])
    return int(ids[0]) if ids else None


def fetch_properties(cid: int):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IUPACName,CanonicalSMILES,InChI,MolecularWeight,MolecularFormula/JSON"
    r = requests.get(url, timeout=10)
    if not r.ok:
        return None
    props = r.json().get("PropertyTable", {}).get("Properties", [])
    if not props:
        return None
    p = props[0]
    return {
        "pref_name": p.get("IUPACName") or None,
        "smiles": p.get("CanonicalSMILES") or None,
        "inchi": p.get("InChI") or None,
        "formula": p.get("MolecularFormula") or None,
        "mol_weight": p.get("MolecularWeight") or None,
    }


def fetch_synonyms(cid: int):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
    r = requests.get(url, timeout=10)
    if not r.ok:
        return []
    j = r.json()
    return j.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])


def fetch_sdf(cid: int):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    r = requests.get(url, timeout=20)
    if not r.ok or not r.text.strip():
        return None
    return r.content


 


def seed_database(compound_list):
    """Insert fetched compound data into PostgreSQL database using the schema expected by app.py."""
    conn = psycopg2.connect(DB_URL)
    cur = conn.cursor()

    # Create extension and table matching app.py's schema
    cur.execute("""
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
    """)

    rows = []
    for q in tqdm(compound_list, desc="Fetching molecules"):
        try:
            cid = fetch_cid_for_name(q)
            if not cid:
                print(f"No CID for {q}, skipping")
                continue
            props = fetch_properties(cid) or {}
            syns = fetch_synonyms(cid) or []
            sdf = fetch_sdf(cid)

            rows.append((
                int(cid),
                props.get("pref_name") or q,
                syns,
                props.get("smiles"),
                props.get("inchi"),
                props.get("formula"),
                props.get("mol_weight"),
                psycopg2.Binary(sdf) if sdf else None,
            ))
        except Exception as e:
            print(f"Error processing {q}: {e}")

    if rows:
        insert_sql = """
        INSERT INTO molecules (cid, pref_name, common_names, smiles, inchi, formula, mol_weight, sdf)
        VALUES %s
        ON CONFLICT (cid) DO UPDATE SET
          pref_name = EXCLUDED.pref_name,
          common_names = EXCLUDED.common_names,
          smiles = EXCLUDED.smiles,
          inchi = EXCLUDED.inchi,
          formula = EXCLUDED.formula,
          mol_weight = EXCLUDED.mol_weight,
          sdf = EXCLUDED.sdf;
        """
        execute_values(cur, insert_sql, rows, template=None, page_size=100)
        conn.commit()
    # Avoid Unicode/emoji output which can cause encoding errors on some consoles (Windows cp1252)
    print(f"Inserted/updated {len(rows)} molecules.")

    cur.close()
    conn.close()


if __name__ == "__main__":
    common_compounds = [
        "Water", "Ethanol", "Methane", "Glucose", "Caffeine",
        "Acetone", "Aspirin", "Sodium chloride", "Carbon dioxide",
        "Benzene", "Ammonia", "Sulfuric acid", "Acetic acid"
    ]
    seed_database(common_compounds)
