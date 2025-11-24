import json
import time
from mdmodels import get_solvent_from_pubchem

# Read solvent names from a text file (one name per line)
with open("solvents_list.txt", encoding="utf-8") as f:
    names = [line.strip() for line in f if line.strip()]

solvents = []
for i, name in enumerate(names, 1):
    attempt = 0
    while attempt < 3:
        try:
            sm = get_solvent_from_pubchem(name)
            if sm:
                solvents.append(sm)
            else:
                print(f"  -> No data found")
            break
        except Exception as e:
            attempt += 1
            print(f"  -> Error: {e} (Attempt {attempt}/3)")
            time.sleep(3)
    # Save progress after each fetch
    with open("solvents_db.json", "w", encoding="utf-8") as f:
        json.dump([s.model_dump() for s in solvents], f, ensure_ascii=False, indent=2)
    time.sleep(1)  # 1 second break between requests 