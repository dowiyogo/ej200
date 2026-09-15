#!/usr/bin/env python3
"""Momentos del pozo de emisión por cara en EJ200_xp0; sin seleccionar primeros hits."""
import json
import shlex
import sys
import numpy as np
import uproot
from diagnose_decay import OUTPUT, CELL_ID, CHUNK, SCINTILLATION, csv_write, json_write, ratio_cluster_se
from exec46_schema import CAMPAIGN_DIR, TREE_NAME

FACE_LABELS = ['left', 'right', 'top']

def main():
    done=json.loads((CAMPAIGN_DIR/'cells'/CELL_ID/'.DONE').read_text())
    events=done['N_generated']
    counts=np.zeros((3,events))
    sums=np.zeros_like(counts)
    with uproot.open(done['root_path']) as f:
        for a in f[TREE_NAME].iterate(['event_id','face_type','source_type','t_creation_ns'],step_size=CHUNK,library='np'):
            for face in range(3):
                mask=(a['source_type']==SCINTILLATION)&(a['face_type']==face)
                counts[face]+=np.bincount(a['event_id'][mask],minlength=events)
                sums[face]+=np.bincount(a['event_id'][mask],weights=a['t_creation_ns'][mask],minlength=events)
    rows=[]
    for face,label in enumerate(FACE_LABELS):
        if np.any(counts[face]==0): raise RuntimeError('An event has no scintillation hit on this face')
        pooled,pooled_se=ratio_cluster_se(sums[face],counts[face])
        event_mean=sums[face]/counts[face]
        rows.append(dict(pool=label,scintillation_photons=int(counts[face].sum()),
                         pooled_mean_global_creation_ns=pooled,pooled_se_event_ns=pooled_se,
                         equal_event_mean_global_creation_ns=float(event_mean.mean()),
                         equal_event_se_ns=float(event_mean.std(ddof=1)/np.sqrt(events))))
    counts_end=counts[:2].sum(axis=0); sums_end=sums[:2].sum(axis=0)
    pooled,pooled_se=ratio_cluster_se(sums_end,counts_end)
    event_face_mean=(sums[0]/counts[0]+sums[1]/counts[1])/2
    rows.append(dict(pool='END_equal_event_equal_face',scintillation_photons=int(counts_end.sum()),
        pooled_mean_global_creation_ns=pooled,pooled_se_event_ns=pooled_se,
        equal_event_mean_global_creation_ns=float(event_face_mean.mean()),
        equal_event_se_ns=float(event_face_mean.std(ddof=1)/np.sqrt(events))))
    csv_write(OUTPUT/'emission_pool_by_face.csv',rows)
    json_write(OUTPUT/'emission_pool_by_face.meta.json',dict(command=shlex.join([sys.executable,*sys.argv]),
        source_root=done['root_path'],source_root_sha256=done['root_sha256'],events=events,
        source_selection='source_type==1',time_origin_ns=0,rows=rows,
        interpretation='Moment-matched exponential surrogate; not a measured exponential decay lifetime',
        new_events_generated=0,other_cells_scanned=0,first_photon_selection_performed=False))
    print(json.dumps(rows,indent=2))

if __name__=='__main__': main()
