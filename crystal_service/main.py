from fastapi import FastAPI, HTTPException, Request
from pydantic import BaseModel
from typing import Optional
import uvicorn
import io
from pymatgen.core.structure import Structure

app = FastAPI(title="Crystal Convert Service")


class CIFPayload(BaseModel):
    cif: str


@app.post("/convert/cif")
async def convert_cif(payload: CIFPayload):
    """Convert CIF text to an XYZ string (simple conversion)."""
    cif_text = payload.cif
    try:
        struct = Structure.from_str(cif_text, fmt="cif")
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to parse CIF: {e}")

    # Convert fractional coords to Cartesian and generate XYZ
    cart_coords = struct.cart_coords
    symbols = [site.species_string for site in struct.sites]
    buf = io.StringIO()
    buf.write(f"{len(symbols)}\n")
    buf.write(f"Converted from CIF\n")
    for sym, coord in zip(symbols, cart_coords):
        x, y, z = coord
        buf.write(f"{sym} {x:.6f} {y:.6f} {z:.6f}\n")
    return {"xyz": buf.getvalue()}


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
