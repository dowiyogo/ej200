import csv
import random
import subprocess
import os

csv_file = "/home/reriosto/SHiP/consolidation/msi_root_manifest_v3.csv"

rows_transferir = []
with open(csv_file, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row['decision'] == 'transferir':
            rows_transferir.append(row)

# Ensure at least one from 7, 14, 16
c7 = [r for r in rows_transferir if r['campaign_id'] == '7']
c14 = [r for r in rows_transferir if r['campaign_id'] == '14']
c16 = [r for r in rows_transferir if r['campaign_id'] == '16']

# Pick one from each required campaign
r7 = random.choice(c7)
r14 = random.choice(c14)
r16 = random.choice(c16)

selected = [r7, r14, r16]

remaining_pool = [r for r in rows_transferir if r not in selected]
selected.extend(random.sample(remaining_pool, 12))
random.shuffle(selected)

print("| rel | sha_manifiesto[:12] | sha_MSI[:12] | sha_t0[:12] | coincide |")
print("|---|---|---|---|---|")

for row in selected:
    ruta = row['ruta']
    rel = row['rel']
    destino = row['destino']
    sha_manifiesto = row['sha256'][:12]

    # Calculate local sha256
    local_cmd = ["sha256sum", ruta]
    try:
        local_out = subprocess.check_output(local_cmd, stderr=subprocess.DEVNULL).decode('utf-8')
        sha_msi = local_out.split()[0][:12]
    except Exception as e:
        sha_msi = "ERROR"

    # Calculate remote sha256
    remote_path = f"/home/rrios/ej200/{destino}"
    remote_cmd = ["ssh", "-o", "StrictHostKeyChecking=no", "-o", "BatchMode=yes", "rrios@t0minidaq", f"sha256sum {remote_path}"]
    try:
        remote_out = subprocess.check_output(remote_cmd, stderr=subprocess.DEVNULL).decode('utf-8')
        sha_t0 = remote_out.split()[0][:12]
    except Exception as e:
        sha_t0 = "ERROR"
        
    coincide = "Sí" if (sha_manifiesto == sha_msi and sha_msi == sha_t0 and sha_msi != "ERROR") else "No"
    
    print(f"| {rel} | {sha_manifiesto} | {sha_msi} | {sha_t0} | {coincide} |")
