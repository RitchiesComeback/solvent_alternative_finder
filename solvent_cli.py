import json
from mdmodels import SolventModel, get_solvent_from_pubchem, solvent_similarity

# Loads the local database
with open("solvents_db.json", encoding="utf-8") as f:
    db = [SolventModel(**entry) for entry in json.load(f)]

ref_name = input("Reference solvent (name or CAS): ").strip()
ref = get_solvent_from_pubchem(ref_name)
if not ref:
    pass
else:
    results = solvent_similarity(ref, db, top_n=10)
    print(f"{'#':<3} {'Name':<20} {'Score':<8} {'Formula':<12}")
    for i, (alt, score) in enumerate(results, 1):
        print(f"{i:<3} {alt.name:<20} {score:.3f}   {alt.molecular_formula or ''}") 