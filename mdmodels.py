from pydantic import BaseModel
from typing import Optional, List, Tuple, Union
import requests
import numpy as np

class SolventModel(BaseModel):
    """
    Data model for a solvent and certain of its physical and chemical properties.
    """
    name: str
    molecular_formula: Optional[str]
    molar_mass: Optional[float]
    boiling_point: Optional[float]
    melting_point: Optional[float] 
    density: Optional[float]
    vapor_pressure: Optional[float]
    water_solubility: Optional[Union[float, str]]
    dipole_moment: Optional[float]
    dielectric_constant: Optional[float]
    logP: Optional[float]
    smiles: Optional[str]
    refractive_index: Optional[float]
    flash_point: Optional[float]
    heat_capacity: Optional[float]
    viscosity: Optional[float]
    pKa: Optional[float]
    pKb: Optional[float]
    pKw: Optional[float]

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def get_cid(identifier: str) -> Optional[int]:
    """
    Get the PubChem CID for a compound name or CAS number.

    Args:
        identifier (str): Compound name or CAS number.

    Returns:
        Optional[int]: PubChem CID if found, otherwise None.
    """
    for namespace in ["name", "rn"]:
        url = f"{PUBCHEM_BASE}/compound/{namespace}/{identifier}/cids/JSON"
        r = requests.get(url)
        if r.ok:
            data = r.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            if cids:
                return cids[0]
        else:
            print(f"Error response: {r.status_code}")
    return None

def get_solvent_from_pubchem(identifier: str) -> Optional[SolventModel]:
    """
    Fetches solvent parameters from PubChem and returns a SolventModel.

    Args:
        identifier (str): Compound name or CAS number.

    Returns:
        Optional[SolventModel]: SolventModel instance if found, otherwise None.
    """
    cid = get_cid(identifier)
    if not cid:
        print(f"No CID found for '{identifier}'.")
        return None
    # Only request the most reliable properties
    property_fields = [
        "MolecularFormula",
        "MolecularWeight",
        "XLogP",
        "IsomericSMILES",
        "CanonicalSMILES"
    ]
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{','.join(property_fields)}/JSON"
    r = requests.get(url)
    if not r.ok:
        print(f"Error in property response: {r.text}")
        return None
    props = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
    if not props:
        print(f"No properties found in response: {r.json()}")
        return None
    # Fetch experimental values
    exp = get_experimental_properties_from_pubchem(cid)
    def parse_float(val):
        if val is None:
            return None
        try:
            import re
            m = re.search(r"[-+]?[0-9]*\.?[0-9]+", val)
            return float(m.group()) if m else None
        except Exception:
            return None

    def get_smiles(props):
        return (
            props.get("IsomericSMILES") or
            props.get("CanonicalSMILES") or
            props.get("SMILES") or
            props.get("ConnectivitySMILES")
        )
    return SolventModel(
        name=identifier,
        molecular_formula=props.get("MolecularFormula"),
        molar_mass=props.get("MolecularWeight"),
        boiling_point=parse_float(exp.get("Boiling Point")),
        melting_point=parse_float(exp.get("Melting Point")),
        density=parse_float(exp.get("Density")),
        vapor_pressure=parse_float(exp.get("Vapor Pressure")),
        water_solubility=exp.get("Solubility"),
        dipole_moment=parse_float(exp.get("Dipole Moment")),
        dielectric_constant=parse_float(exp.get("Dielectric Constant")),
        logP=props.get("XLogP"),
        smiles=get_smiles(props),
        refractive_index=parse_float(exp.get("Refractive Index")),
        flash_point=parse_float(exp.get("Flash Point")),
        heat_capacity=parse_float(exp.get("Heat Capacity")),
        viscosity=parse_float(exp.get("Viscosity")),
        pKa=parse_float(exp.get("pKa")),
        pKb=parse_float(exp.get("pKb")),
        pKw=parse_float(exp.get("pKw")),
    )

def get_solvent_list(identifiers: list[str]) -> list[SolventModel]:
    """
    Fetches SolventModel instances for a list of names or CAS numbers.

    Args:
        identifiers (list[str]): List of compound names or CAS numbers.

    Returns:
        list[SolventModel]: List of SolventModel instances.
    """
    result = []
    for ident in identifiers:
        sm = get_solvent_from_pubchem(ident)
        if sm:
            result.append(sm)
    return result

def solvent_similarity(reference: SolventModel, candidates: List[SolventModel], top_n: int = 10, weights: dict = None) -> List[Tuple[SolventModel, float]]:
    """
    Compares a reference solvent with a list and returns the most similar alternatives.

    For each parameter, the distance is calculated as the relative deviation:
        |x_ref - x_alt| / |x_ref|
    If x_ref == 0, the absolute difference is used instead.
    All distances are weighted according to the weights dict, that is set by the user.
    Additionally, dynamic scaling is applied to the weights:
        - If relative deviation <= 10%: weight is doubled
        - If relative deviation <= 20%: weight is multiplied by 1.5
        - Otherwise: weight is unchanged

    The final similarity score is calculated as:
        similarity = max(0, 1 - distance)
    where 1 means perfect match and 0 means no similarity.

    Args:
        reference (SolventModel): The reference solvent.
        candidates (List[SolventModel]): List of candidate solvents.
        top_n (int, optional): Number of top alternatives to return. Defaults to 10.
        weights (dict, optional): Weighting for each parameter. Defaults to None.

    Returns:
        List[Tuple[SolventModel, float]]: List of tuples (SolventModel, similarity score), sorted by score (descending).
    """
    fields = [
        'boiling_point', 'melting_point', 'density', 'vapor_pressure',
        'dipole_moment', 'dielectric_constant', 'water_solubility',
        'refractive_index', 'flash_point', 'heat_capacity', 'viscosity', 'pKa', 'pKb', 'pKw'
    ]
    if weights is None:
        weights = {f: 1 for f in fields}
    results = []
    for cand in candidates:
        dist_sum = 0.0
        weight_sum = 0.0
        for f in fields:
            base_weight = weights.get(f, 0)
            if base_weight > 0:
                ref_val = getattr(reference, f)
                cand_val = getattr(cand, f)
                if isinstance(ref_val, (float, int)) and isinstance(cand_val, (float, int)):
                    if ref_val == 0:
                        rel_diff = abs(ref_val - cand_val)
                    else:
                        rel_diff = abs(ref_val - cand_val) / abs(ref_val)
                    # Dynamic Scaling
                    if rel_diff <= 0.1:
                        dyn_weight = base_weight * 2
                    elif rel_diff <= 0.2:
                        dyn_weight = base_weight * 1.5
                    else:
                        dyn_weight = base_weight
                    dist_sum += dyn_weight * rel_diff
                    weight_sum += dyn_weight
        if weight_sum == 0:
            continue
        avg_dist = dist_sum / weight_sum
        similarity = max(0.0, 1.0 - avg_dist)
        results.append((cand, similarity))
    results.sort(key=lambda x: x[1], reverse=True)
    return results[:top_n]

def get_experimental_properties_from_pubchem(cid: int) -> dict:
    """
    Extracts experimental values (boiling point, melting point, density, etc.) from the PubChem Record endpoint.

    Args:
        cid (int): PubChem Compound ID (CID).

    Returns:
        dict: Dictionary of experimental property names and their values.
    """
    import requests
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    r = requests.get(url)
    if not r.ok:
        print(f"Error in record response: {r.text}")
        return {}
    data = r.json()
    result = {}
    def find_properties(sections, keys):
        for section in sections:
            toc = section.get("TOCHeading", "")
            if toc in keys and "Information" in section:
                for info in section["Information"]:
                    val = info.get("Value", {}).get("StringWithMarkup", [{}])[0].get("String")
                    if val:
                        result[toc] = val
            if "Section" in section:
                find_properties(section["Section"], keys)
    keys = [
        "Boiling Point", "Melting Point", "Density", "Refractive Index", "Solubility", "Vapor Pressure",
        "Dielectric Constant", "Dipole Moment", "Flash Point", "Autoignition Temperature",
        "Heat Capacity", "Viscosity", "pKa", "pKb", "pKw"
    ]
    root_sections = data.get("Record", {}).get("Section", [])
    find_properties(root_sections, keys)
    return result 