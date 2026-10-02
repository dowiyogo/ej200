import csv
import subprocess

prior_rels = [
    "ej200_edge_scan/build/photon_hits_run010.root",
    "ej200_endonly/output/endonly_mylar_msi_20260621_092144/photon_hits_x-500mm.root",
    "t0minidaq/validation_skinfix_msi_300ev_31pos/raw/ej204_endtop/ej204_endtop_xm414mm_300ev_a0368c4.root",
    "ej200_edge_scan/results/scan_hi_2026-06-08/photon_hits_run023.root",
    "t0minidaq/validation_skinfix_msi_300ev_31pos/raw/ej230_endtop/ej230_endtop_xm598mm_300ev_ca2f1c3.root",
    "ej200_endonly/output/endonly_mylar_msi_20260621_092144/photon_hits_x-150mm.root",
    "ej200_endonly/out/EXEC_26/raw/rough_R0_ground_s0p30/photon_hits.root",
    "ej200_endonly/output/endonly_mylar_msi_20260620_222723/photon_hits_x-150mm.root",
    "ej200_endonly/output/endonly_mylar_msi_20260621_092144/photon_hits_x-100mm.root",
    "t0minidaq/sslg4/exec08b_window_dip/photon_hits_run_B_x-642mm.root",
    "ej200_endonly/output/endonly_mylar_msi_20260620_222723/photon_hits_x-50mm.root",
    "ej200_endonly/output/endonly_mylar_msi_20260620_222723/photon_hits_x-200mm.root",
    "t0minidaq/validation_skinfix_msi_300ev_31pos/raw/ej230_endtop/ej230_endtop_xp690mm_300ev_ca2f1c3.root",
    "event_display_ej204_ej230/outputs/ej230/xm690/photon_hits_run000.root",
    "t0minidaq/validation_skinfix_msi_300ev_31pos/raw/ej204_endtop/ej204_endtop_xm552mm_300ev_a0368c4.root"
]

with open('/home/reriosto/SHiP/consolidation/msi_root_manifest_v3.csv') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

found_rels = {}
for row in rows:
    if row['rel'] in prior_rels:
        found_rels[row['rel']] = row

def run_cmd(cmd):
    try:
        return subprocess.check_output(cmd, shell=True).decode('utf-8').strip()
    except subprocess.CalledProcessError as e:
        return f"ERROR"

for rel in prior_rels:
    if rel in found_rels:
        row = found_rels[rel]
        ruta = row['ruta']
        destino = row['destino']
        sha_manifest = row['sha256'][:12]
        
        sha_local = run_cmd(f"sha256sum {ruta}")
        sha_local = sha_local.split()[0][:12] if "ERROR" not in sha_local else "ERROR"
        
        # Test exact flawed path requested in prompt
        t0_path_flawed = f"/home/rrios/ej200/data/ej200_campaigns/{destino}"
        sha_t0_flawed = run_cmd(f"ssh -o StrictHostKeyChecking=no rrios@t0minidaq 'sha256sum {t0_path_flawed} 2>/dev/null'")
        sha_t0_flawed = sha_t0_flawed.split()[0][:12] if sha_t0_flawed and "ERROR" not in sha_t0_flawed else "ERROR"

        # Test fixed path
        t0_path_fixed = f"/home/rrios/ej200/{destino}"
        sha_t0_fixed = run_cmd(f"ssh -o StrictHostKeyChecking=no rrios@t0minidaq 'sha256sum {t0_path_fixed} 2>/dev/null'")
        sha_t0_fixed = sha_t0_fixed.split()[0][:12] if sha_t0_fixed and "ERROR" not in sha_t0_fixed else "ERROR"
        
        print(f"rel: {rel}")
        print(f"sha_manifiesto: {sha_manifest}")
        print(f"sha_MSI: {sha_local}")
        print(f"sha_t0_flawed: {sha_t0_flawed}")
        print(f"sha_t0_fixed: {sha_t0_fixed}")
        print("-" * 20)
