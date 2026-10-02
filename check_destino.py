import csv

with open('/home/reriosto/SHiP/consolidation/msi_root_manifest_v3.csv') as f:
    reader = csv.DictReader(f)
    destinos = [row['destino'] for row in reader]

print("First 5 destinos:", destinos[:5])
print("Do all start with 'data/ej200_campaigns/'?:", all(d.startswith('data/ej200_campaigns/') or d == '' for d in destinos))
