import os
import requests
import psycopg2
from psycopg2.extras import execute_values
from tqdm import tqdm

# Get DB URL from environment or fallback
DB_URL = os.getenv("DATABASE_URL", "postgresql://user:password@localhost:5432/molevis")

def fetch_pubchem_data(query):
    """Fetch compound data from PubChem using the REST API"""
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    try:
        name_url = f"{base_url}/compound/name/{query}/property/IUPACName,CanonicalSMILES,MolecularWeight,MolecularFormula/JSON"
        resp = requests.get(name_url, timeout=10)
        if resp.status_code != 200:
            return None
        data = resp.json()["PropertyTable"]["Properties"][0]
        return {
            "name": query,
            "iupac_name": data.get("IUPACName"),
            "smiles": data.get("CanonicalSMILES"),
            "molecular_formula": data.get("MolecularFormula"),
            "molecular_weight": data.get("MolecularWeight"),
        }
    except Exception as e:
        print("Error fetching", query, ":", e)
        return None


def seed_database(compound_list):
    """Insert fetched compound data into PostgreSQL database"""
    conn = psycopg2.connect(DB_URL)
    cur = conn.cursor()

    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id SERIAL PRIMARY KEY,
        name TEXT,
        iupac_name TEXT,
        smiles TEXT,
        molecular_formula TEXT,
        molecular_weight FLOAT
    );
    """)

    data = []
    for query in tqdm(compound_list, desc="Fetching molecules"):
        result = fetch_pubchem_data(query)
        if result:
            data.append((
                result["name"],
                result["iupac_name"],
                result["smiles"],
                result["molecular_formula"],
                result["molecular_weight"]
            ))

    if data:
        execute_values(cur, """
        INSERT INTO molecules (name, iupac_name, smiles, molecular_formula, molecular_weight)
        VALUES %s
        ON CONFLICT DO NOTHING;
        """, data)
        conn.commit()
        print(f"✅ Inserted {len(data)} molecules.")

    cur.close()
    conn.close()


if __name__ == "__main__":
    common_compounds = [
        "Water", "Ethanol", "Methane", "Glucose", "Caffeine",
        "Acetone", "Aspirin", "Sodium chloride", "Carbon dioxide",
        "Benzene", "Ammonia", "Sulfuric acid", "Acetic acid"
    ]
    seed_database(common_compounds)
