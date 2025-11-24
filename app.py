from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import json
from mdmodels import SolventModel, get_solvent_from_pubchem, solvent_similarity

app = FastAPI()

# Templates and static files
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates")

# Load the local database
with open("solvents_db.json", encoding="utf-8") as f:
    db = [SolventModel(**entry) for entry in json.load(f)]

# Neue PARAMS-Liste ohne die entfernten Felder
PARAMS = [
    'boiling_point', 'melting_point', 'density', 'vapor_pressure',
    'refractive_index', 'flash_point', 'viscosity'
]

def get_attr(obj, key):
    return getattr(obj, key, "")

@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    weights = {p: 0 for p in PARAMS}
    # water_solubility als Checkbox-Status (default: anzeigen)
    return templates.TemplateResponse("index.html", {"request": request, "results": None, "ref_name": "", "weights": weights, "get_attr": get_attr, "show_water_solubility": True})

@app.post("/similarity", response_class=HTMLResponse)
async def similarity(request: Request):
    form = await request.form()
    ref_name = form.get("ref_name", "")
    weights = {}
    for p in PARAMS:
        try:
            weights[p] = int(form.get(f"w_{p}", 1))
        except Exception:
            weights[p] = 1
    # Checkbox für water_solubility
    show_water_solubility = form.get("show_water_solubility") == "on"
    ref = get_solvent_from_pubchem(ref_name)
    if not ref:
        return templates.TemplateResponse(
            "index.html",
            {
                "request": request,
                "results": [],
                "ref_name": ref_name,
                "error": f"Could not find '{ref_name}'.",
                "weights": weights,
                "get_attr": get_attr,
                "show_water_solubility": show_water_solubility
            }
        )
    results = solvent_similarity(ref, db, top_n=10, weights=weights)
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "results": results,
            "ref_name": ref_name,
            "error": None,
            "weights": weights,
            "get_attr": get_attr,
            "show_water_solubility": show_water_solubility
        }
    ) 