# Solvent Alternative Finder

This project provides a tool for finding alternative solvents based on quantitative physicochemical properties. It leverages PubChem data and allows users to compare solvents and identify the most similar alternatives for a given reference solvent.

## Features
- **Automated data collection**: Fetches experimental and calculated solvent properties from PubChem for hundreds of common solvents.
- **Customizable similarity scoring**: Users can select which parameters to consider and assign individual weights to each property.
- **Web interface**: Simple local web app (FastAPI + Bootstrap) for interactive use.
- **Command-line interface**: For quick comparisons in the terminal.
- **Extensible data model**: Supports a wide range of solvent properties (boiling point, density, pKa, viscosity, etc.).

## Main Components
- `mdmodels.py`: Data model and all core functions (data fetching, similarity calculation).
- `solvent_db_builder.py`: Script to build a local database of solvent properties from a list of names.
- `solvent_cli.py`: Command-line tool for similarity search.
- `app.py` + `templates/index.html`: FastAPI web app for interactive solvent comparison.
- `solvents_list.txt`: List of solvent names to build the database.
- `solvents_db.json`: Local database of solvent properties (generated).

## Setup
1. **Install dependencies**
   ```bash
   pip install fastapi uvicorn jinja2 pydantic requests numpy python-multipart
   ```
2. **Build the local database OR use the provided one**
   - Edit `solvents_list.txt` to include your solvents (one per line).
   - Run:
     ```bash
     python solvent_db_builder.py
     ```
   - This will create/update `solvents_db.json`.

## Usage
### Web App
1. Start the server:
   ```bash
   python -m uvicorn app:app --reload
   ```
2. Open your browser at [http://127.0.0.1:8000/](http://127.0.0.1:8000/)
3. Enter a reference solvent (name or CAS), select parameters and weights, and view the top alternatives.

### Command-Line Interface
```bash
python solvent_cli.py
```
- Enter the reference solvent when prompted.
- The top 10 most similar alternatives will be displayed.

## Customization
- You can extend `solvents_list.txt` with more solvents.
- The data model in `mdmodels.py` can be expanded to include additional properties.
- The web interface can be styled or extended as needed.

## Example Properties Used
- Boiling point, melting point, density, vapor pressure
- Dipole moment, dielectric constant, water solubility
- Refractive index, flash point, heat capacity, viscosity
- pKa, pKb, pKw

## Notes
- Not all properties are available for every solvent in PubChem.
- The similarity score is based only on available numeric data for the selected parameters.
- The tool is intended for educational and research purposes.

## Outlook
As the current tool is good for comparing solvents in terms of their physical properties, it might be useful for engineering purposes. But it does not provide enough high quality data to ensure more in-depth chemical behavior (pKa, donors, etc.). This lack is due to the poor accessability of this data in PubChem via its API. This project would improve a lot by the usage of commercially availabe databases (e.g. Reaxys). But at a certain point more advanced tasks couldn't be tackled by the use of a single database (there is (probably) none that is so extensive, that it includes almost any data for a given compound, e.g. IR, Raman, UV/Vis, NMR, EPR, ...). This would need the use of multiple databases, many of them without an API.
The tool is designed in a modular way, so it can be easily extended to incorporate additional data sources (such as commercial or experimental databases) in the future.

## License
MIT License 
