from pathlib import Path
import pandas as pd

base = Path('analysis/track_mechanism_20260915/step6_v2/step3')
df = pd.read_csv(base / 'first_photon_cell_face.csv')
rows = []
for (material, x_mm), group in df.groupby(['material', 'gun_x_mm'], sort=True):
    joined = int(group['n_events'].sum())
    rows.append({'cell_id': f'{material.replace("-", "")}_' + ('xm' if x_mm < 0 else 'xp') + str(abs(int(x_mm))),
                 'material': material, 'x_mm': int(x_mm),
                 'photons_selected': joined, 'photons_joined': joined,
                 'photons_rejected': 0})
pd.DataFrame(rows).to_csv(base / 'wavelength_join_diagnostics.csv', index=False)
