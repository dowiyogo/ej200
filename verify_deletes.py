import csv
import random
import os
import subprocess
import sys

delete_log = "/home/reriosto/SHiP/consolidation/delete_log_20261002_022114.csv"
manifest = "/home/reriosto/SHiP/consolidation/msi_root_manifest_v3.csv"

# Load manifest
manifest_hashes = {}
with open(manifest, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        manifest_hashes[row['ruta']] = row['sha256']

# Load delete log and filter
borrados = []
with open(delete_log, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['estado'] == 'borrado':
            borrados.append(row)

# Randomly select 10
selected = random.sample(borrados, min(10, len(borrados)))

print("| ruta | t0_verificado | MSI_Borrado | T0_Existe_Y_Hash_Igual | Notas |")
print("|---|---|---|---|---|")

all_good = True

for row in selected:
    ruta = row['ruta']
    t0_verif = row['t0_verificado']
    
    # 1. Check if ruta exists locally
    msi_borrado = not os.path.exists(ruta)
    
    # 2. Check remote hash
    t0_ok = False
    notas = ""
    
    if not msi_borrado:
        all_good = False
        notas += "Ruta todavia existe en MSI. "
    
    if t0_verif:
        expected_hash = manifest_hashes.get(ruta)
        if not expected_hash:
            notas += "No se encontro hash en manifest para esta ruta. "
            all_good = False
        else:
            # ssh and get sha256sum
            cmd = ["ssh", "-o", "BatchMode=yes", "rrios@t0minidaq", f"sha256sum {t0_verif}"]
            try:
                res = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
                if res.returncode == 0:
                    remote_hash = res.stdout.split()[0]
                    if remote_hash == expected_hash:
                        t0_ok = True
                    else:
                        notas += f"Hash no coincide: esperado {expected_hash}, remoto {remote_hash}. "
                        all_good = False
                else:
                    notas += "Archivo no existe en remoto o error ssh. "
                    all_good = False
            except Exception as e:
                notas += f"Error ssh: {e}. "
                all_good = False
    else:
        notas += "t0_verificado esta vacio. "
        t0_ok = "N/A"
        # If t0_verificado is empty, the requirement "la ruta ... (si no esta vacia)" is implicitly satisfied? 
        # Actually it's N/A. Let's not fail the whole check just because it's empty, unless it shouldn't be empty.

    print(f"| {ruta} | {t0_verif} | {msi_borrado} | {t0_ok} | {notas.strip()} |")

print(f"\nVerdict: {'ALL CHECKS PASSED' if all_good else 'SOME CHECKS FAILED'}")
