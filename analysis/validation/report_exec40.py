#!/usr/bin/env python3
"""Publish existing EXEC40 results and sidecar figures; never rerun transport."""
import argparse
import csv
import datetime as dt
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uproot
from check_exec40 import sha, dump, ratio, sidecar, N, SEED, REPLICATES

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--previous-golden',type=Path,required=True)
    a=p.parse_args()
    out=a.analysis
    r=json.loads((out/'results.json').read_text())
    repo=Path(__file__).resolve().parents[2]
    meta=r['metadata']
    closure={}
    if r['G_I_1']['status']=='PASS':
        # Explicit post-hoc diagnostic of the observed matching failure. This
        # never substitutes matched detections for the preregistered numerator.
        channels=uproot.open(out/'sipm_event_counts.root')['data'].arrays(library='np')
        rng=np.random.default_rng(SEED)
        weights=np.stack([np.bincount(rng.integers(0,N,N),minlength=N) for _ in range(REPLICATES)]).astype(float)
        for name,face in [('left',0),('right',1),('top',2),('all',None)]:
            mask=np.ones(len(channels['event_id']),dtype=bool) if face is None else channels['face_type']==face
            def by_event(key):
                return np.bincount(channels['event_id'][mask],weights=channels[key][mask],minlength=N)
            incident=by_event('incident');matched=by_event('matched_detected')
            closure[name]=dict(matched_ratio=ratio(matched,incident,weights),
                matched_minus_expected=ratio(matched-by_event('expected_surface_pde_sum'),incident,weights),
                matching_equals_surface_detection=bool(np.array_equal(matched,by_event('surface_detection'))),
                duplicate_detections=int((by_event('detected')-by_event('detected_unique')).sum()))
        dump(out/'posthoc_detection_closure.json',dict(post_hoc=True,primary_result_unchanged=True,arms=closure))
        sidecar(out,'posthoc_detection_closure',dict(arm=list(closure),
            matched_ratio=[v['matched_ratio']['value'] for v in closure.values()],
            matched_ratio_se=[v['matched_ratio']['se'] for v in closure.values()],
            residual=[v['matched_minus_expected']['value'] for v in closure.values()],
            residual_se=[v['matched_minus_expected']['se'] for v in closure.values()]),dict(meta,
                post_hoc=True,primary_result_unchanged=True,source_sidecar_sha256=sha(out/'sipm_event_counts.root')))
    run_link=f"[R40: N=2000, seeds 26092601/8349041]({Path(meta['source_root']).parent/'invocation.meta.json'})"
    golden=json.loads(a.previous_golden.read_text())
    golden['exec40_single_cell']={'scope':'EJ204 x=0 EndTop70 only; no grid extension',
        'previous_golden_path':str(a.previous_golden),'previous_golden_sha256':sha(a.previous_golden),
        'results_path':str(out/'results.json'),'results_sha256':sha(out/'results.json'),
        'posthoc_detection_closure_path':str(out/'posthoc_detection_closure.json'),
        'preregistration_path':str(repo/'analysis/validation/EXEC40_PREREGISTRATION.md'),
        'preregistration_sha256':sha(repo/'analysis/validation/EXEC40_PREREGISTRATION.md'),
        'report_path':str(a.report),**r}
    golden['ready_for_acceptance']=False
    golden['overall_status']='INCOMPLETE_WITH_FAILED_HYPOTHESES_AND_INCIDENT_CLOSURE'
    golden['updated_utc']=dt.datetime.now(dt.timezone.utc).isoformat()
    golden['title']='EXEC38 golden reference with EXEC40 single-cell observations'
    for test in golden['tests']:
        if test['id'] in ('V1','V2','V5'):
            test['status_before_exec40']=test['status']
            test['single_cell_exec40']=r[test['id']]
            test['status']='PARTIAL_SINGLE_CELL_'+r[test['id']]['status']
            test['scope_note']='Original measured/config_ids remain the archived grid; EXEC40 adds one distinct cell only.'
    golden_path=repo/'analysis/validation/GOLDEN_REFERENCE_20260913.json'
    dump(golden_path,golden)
    g1,g2,g3=r['G_I_1'],r['G_I_2'],r['G_I_3']
    lines=['# EXEC40 — Production optical observations and one-cell validation\n',
        '| Gate / hypothesis | Outcome |\n|---|---|',
        f"| G-I.1 exact physical invariance | {g1['status']}; {g1['npe_end_mean']:.7f} ± {g1['npe_end_sem']:.7f} pe/end; difference {g1['difference']:.7f}; {g1['events_with_different_counts']} different events |",
        f"| G-I.2 readable, complete observations | {g2['status']}; "
        + (f"{r['V5']['all']['unmatched_detected']:,} unmatched SD detections |" if 'all' in r['V5'] else 'not evaluated |'),
        f"| G-I.3 output cost | {g3['root_ratio']:.3f}× D1 ROOT, {g3['wall_s']:.2f} s; proposal required: {g3['proposal_required']} |",
        *[f"| {name} | {r[name]['status']} |" for name in ('V1','V2','V5')],
        f'\nAll new measurements below refer to {run_link}. Only one cell generated transport events. '
        'No grid, physics retuning, injected electronics, BLUE, deck edits, push or main merge. '
        '**`ready_for_acceptance=false`; the grid requires separate approval.**\n',
        '## Registered contract and provenance\n',
        f"[Preregistration]({repo}/analysis/validation/EXEC40_PREREGISTRATION.md), "
        f"SHA256 `{meta['preregistration_sha256']}`, was frozen before the successful transport run. "
        'G-I.1 requires exact zero difference and identical event counts, otherwise abort. '
        'V1 requires unity within 3 event-bootstrap SE; V2 requires overlap with 0.36–0.40 and '
        'a p≥0.01 cosine goodness of fit; V5 requires agreement with 0.619397 within 3 SE. '
        'Bootstrap: 500 paired generated-event resamples, seed 40091301; reported ± values are one SE.\n',
        f"Worktree `{repo}`, branch `diag/exec40-20260913`, base `420addf`, rollback tag "
        '`pre-exec40-20260913` created before editing. Main remains unchanged. '
        'The run metadata preserves the pre-regression source commit, exact working patch, individual '
        'C++ source hashes, macro/material-resource hashes and binary hash. The final regression '
        'commit packages that source state with tests, documentation and this golden update.\n',
        'Configuration: EndTop, N_TOP=70, EJ-204/OPSC-101, vertical mu− 1 GeV, x=0 mm, '
        'N=2000, seeds 26092601 8349041, 4 workers, eventModulo 1, jitter=0 ns, diagnostics OFF. '
        'EXEC33 had already measured exactly equal Npe/end at 1/4/12/24 workers; this is explicitly '
        'recorded in the new cell metadata.\n',
        'Exact successful command (working directory is the parent of the ROOT below):\n',
        f"```bash\n{meta['simulation_command']}\n```\n",
        f"ROOT: `{meta['source_root']}`\n\nSHA256: `{meta['source_root_sha256']}`.\n",
        f"Read-only analysis command:\n\n```bash\nOPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "
        f"/home/rrios/exec35_20260912/venv/bin/python {meta['analysis_command']}\n```\n",
        '## G-I: validation and cost\n',
        f"D1 ROOT counts were reconstructed and verified against its archived CSV. New L/R/T totals: "
        f"{g1['totals']}; event ledger equals the new hit table: {g1['event_ledger_matches_hits']}. "
        f"[All 2000 event comparisons]({out}/event_regression.csv).\n",
        'Builds: diagnostics OFF and ON compile. Existing SSLG4 material and readout geometry tests '
        'pass without beam events. Four non-transport regression tests check exact agreement, '
        'compensating left/right errors hidden by the mean, a TOP-only change and an inconsistent '
        'event ledger. The ROOT checker tests field population, independent detection matching '
        'and first-photon uniqueness in the actual validation cell. No additional transport smoke tests ran.\n',
        f"ROOT: {g3['root_bytes']:,} bytes versus D1 {g3['D1_root_bytes']:,} bytes. "
        f"Wall: {g3['wall_s']:.2f} s versus D1 {g3['D1_wall_s']:.2f} s ({g3['wall_ratio']:.3f}×). "
        'D1 used one worker with diagnostics ON, so that wall ratio is not a clean instrumentation '
        f"overhead estimate. Archived EXEC33 S2 used four workers with diagnostics OFF: 358.38 s; "
        f"the new ratio is {g3['matching_S2_wall_ratio']:.3f}×. Builds on other available cores overlapped "
        'part of this run; this is one wall-time observation, not a controlled timing benchmark.\n',
        '| Table | Compressed branch bytes |\n|---|---:|',
        *[f'| {name} | {size:,} |' for name,size in g3['tree_compressed_bytes'].items()],
        f"\n[Per-branch storage sidecars]({out}/root_branch_sizes.csv). Linear 10000-event projection "
        f"is {5*g3['root_bytes']/1e9:.3f} GB, compared with {5*g3['D1_root_bytes']/1e9:.3f} GB for D1.\n",
        'Production choice: deposited energy, first encounters and accepted incident counts interpret '
        'physical results, so all three are available with diagnostics OFF. Existing detailed bounce '
        'and terminal censuses remain optional. '
        + ('The >2× threshold is exceeded. **Proposal only:** move first-encounter `energy_eV`, '
           '`normal_norm`, `normal_valid`, `normal_orientation_valid`, and toolkit `boundary_status` '
           'behind the diagnostic flag, retaining the portable outcome, signed cosine, identity, '
           'source and volume pair in production. The normal checks are debugging evidence and '
           'the per-photon energy is not required for the requested first-angle observable. '
        'Removing verification flags alone will not halve the file. '
           'A further reviewed option is volume dictionaries and reduced precision for angle/energy. '
           '**No proposed storage change has been applied.**'
           if g3['proposal_required'] else 'The >2× threshold is not exceeded.')+'\n',
        'A first initialization attempt in `exec40_20260913/cell/` failed because the relative '
        '`sslg4` resource link was missing. Its fatal material-property lookup explicitly had no '
        'track/step, before event generation; log and 204-byte incomplete ROOT are preserved. '
        'The runner now checks the four required material files against D1, creates the local link '
        'and refuses an existing output directory. The only completed transport cell is `cell_validated`.\n']
    if g1['status']!='PASS':
        lines += ['## Abort\n','The exact invariance gate failed. V1/V2/V5 were not evaluated; no retuning or additional cell was run.']
        a.report.write_text('\n'.join(lines)+'\n')
        return
    v1,v2,v5=r['V1'],r['V2'],r['V5']
    primary=v2['primary_scint_exiting'];gof=v2['cosine_gof']
    lines += ['## V1 — generated scintillation photons per deposited energy\n',
        f"Primary ratio: **{v1['value']:.8f} ± {v1['se']:.8f}**, expected 1; "
        f"**{v1['status']}**. Produced scintillation photons: {v1['produced_scint']:,}; "
        f"all optical photons created in BarLV: {v1['produced_optical']:,}. "
        f"Total deposited energy: {v1['total_edep_MeV']:.8f} MeV, optical component "
        f"{v1['optical_edep_MeV']:.8f} MeV, nonionizing component {v1['nonionizing_edep_MeV']:.8f} MeV. "
        f"The denominator uses all-particle `GetTotalEnergyDeposit()` exactly as preregistered. {run_link}.\n",
        f"Declared secondary ratio after subtracting only the separately recorded optical deposit: "
        f"{v1['secondary_nonoptical_denominator']['value']:.8f} ± {v1['secondary_nonoptical_denominator']['se']:.8f}. "
        'This is not substituted for the primary result. An all-particle energy sum includes energy '
        'deposited by the optical photons themselves; the two recorded denominators distinguish '
        'that contribution from the scintillation source response. The secondary value does not '
        'change the registered primary V1 outcome. '
        f"[Event-energy sidecars]({out}/event_observables.csv).\n",
        '## V2 — first surface encounter and signed incidence angle\n',
        f"Primary first-encounter escape fraction: **{primary['value']:.8f} ± {primary['se']:.8f}**, "
        f"{int(primary['numerator']):,}/{int(primary['denominator']):,}; expected interval 0.36–0.40; "
        f"**{v2['escape_status']}**. Escape means transmission/refraction out of BarLV to air/world; "
        f"sensor absorption/detection does not count as escape. {run_link}.\n",
        '| Population | Escape fraction ± SE |\n|---|---:|',
        *[f"| {name} | {v2[name]['value']:.8f} ± {v2[name]['se']:.8f} |" for name in
          ('primary_scint_exiting','scint_bar_air','all_exiting')],
        f"\nCosine null `p(mu)=2mu`, 20 fixed bins: event-bootstrap covariance statistic "
        f"χ²={gof['chi_square_event_bootstrap']:.4g}, rank={gof['rank']}, p={gof['p_value']:.6g}; "
        f"**{gof['status']}** at p≥0.01. Photon Pearson χ²={gof['pearson_chi_square_photons']:.4g} "
        f"(19 df, p={gof['pearson_p_value']:.6g}) is secondary because photons share events. "
        'Both printed zero p-values are floating-point underflow (p<1e-300), not exact zero probability. '
        'The bootstrap statistic is an asymptotic covariance-based GOF, not an exact finite-sample p-value.\n',
        f"Normal/identity checks: `{json.dumps(v2['observation_checks'],sort_keys=True)}`. "
        'Cosines use the pre-step direction and an outward-from-pre-volume normal with no absolute '
        'value or fitted sign. For daughter sensors the normal comes from the actual sensor solid, '
        'not the bar outer face. ±1e-5 mm solid probes check its orientation.\n',
        'The registered isotropic-flux law applies to equilibrium flux through a surface. First '
        'hits of a localized source in a finite rectangular bar are a different conditional '
        'population. The literal requested result is retained, but a failure here does not by '
        'itself establish a faulty isotropic-emission generator. No parameters or population '
        'definitions were changed to obtain a more favorable result.\n',
        f"![First-encounter cosine]({out}/first_encounter_cosine.png)\n",
        '## V5 — accepted incident flux versus detected photons\n',
        'Incidents are recorded independently in stepping, detections independently in the '
        'existing sensitive detector. Accepted terminal absorption/detection at the sensitive '
        'R=0 surface counts as reaching it; reflection without entry is excluded. Matches are '
        'joined at event end, so callback order cannot bias the denominator.\n',
        f"**G-I.2 fails the preregistered matching check:** {v5['all']['unmatched_detected']:,} "
        'SD detections have no accepted surface incident under the declared criterion. The '
        'counts are readable, and all 172000 channel/event rows exist, but completeness of the '
        'physical incident/detection accounting is not established. This may expose a pre-existing '
        'SD path or an unobserved class of arrivals; these records do not identify the missing '
        'history. It is not evidence that the added instrument altered transport (G-I.1 passed). '
        'The numerical comparisons below remain FAIL, while V5 is overall NOT VALIDATED. '
        'Do not reinterpret these ratios as a validated PDE measurement.\n',
        '| Arm | Detected / incident | Ratio ± SE | Emission-PDE test | Incident-spectrum surface PDE ± SE | Unmatched |\n|---|---:|---:|---|---:|---:|',
        *[f"| {name} | {int(v5[name]['numerator']):,}/{int(v5[name]['denominator']):,} | "
          f"{v5[name]['value']:.8f} ± {v5[name]['se']:.8f} | {v5[name]['status']} | "
          f"{v5[name]['incident_spectrum_surface_pde']['value']:.8f} ± {v5[name]['incident_spectrum_surface_pde']['se']:.8f} | "
          f"{v5[name]['unmatched_detected']} |" for name in ('left','right','top','all')],
        f"\nEvery primary row is compared with the fixed emitted-spectrum EJ204 reference **0.619397**, "
        f"within 3 paired-event SE. Aggregate V5: **{v5['status']}**. {run_link}.\n",
        '| Arm | Surface Detection | Surface Absorption | Transmitted | Missing surface PDE | Scintillation incident fraction |\n|---|---:|---:|---:|---:|---:|',
        *[f"| {name} | {v5[name]['surface_detected']:,} | {v5[name]['surface_absorbed']:,} | "
          f"{v5[name]['transmitted']:,} | {v5[name]['unknown_surface_pde']:,} | {v5[name]['incident_scint_fraction']:.6f} |"
          for name in ('left','right','top','all')],
        '\nIncident-spectrum expectations read the actual directed border EFFICIENCY with energy '
        'interpolation, independent of whether that photon was detected. The legacy emitted-spectrum '
        'prediction integrates the wavelength-interpolated PDE. Differences can therefore reflect '
        'transport spectral selection, non-scintillation photons, interface paths or interpolation '
        'conventions; primary agreement/disagreement alone does not isolate one cause. No sensor '
        'curve or transport setting is corrected.\n',
        '**Post-hoc closure diagnostic**, added after observing the unmatched population; '
        'not an alternative primary estimator:\n',
        '| Arm | Matched detections / accepted incidents ± SE | (Matched − surface-PDE expectation) / incidents ± SE |\n|---|---:|---:|',
        *[f"| {name} | {closure[name]['matched_ratio']['value']:.8f} ± {closure[name]['matched_ratio']['se']:.8f} | "
          f"{closure[name]['matched_minus_expected']['value']:.8f} ± {closure[name]['matched_minus_expected']['se']:.8f} |"
          for name in ('left','right','top','all')],
        '\nMatched detections equal recorded boundary Detection counts for every event in each arm; '
        'there are no duplicate SD detections. The surplus SD count is therefore an observed '
        'unmatched population, not a duplicate callback total. These statements are from the '
        f"[closure sidecars]({out}/posthoc_detection_closure.meta.json), with the same run, N and seeds. "
        'No unmatched detection is discarded from the primary numerator and no incident is '
        'created from a detection to close the gap. Further trajectory evidence would require '
        'separate work; no additional cell was run.\n',
        f"![Per-channel detection fraction]({out}/sipm_channel_summary.png)\n",
        '## Deliverables and final gate\n',
        f"[Output contract and commands]({repo}/analysis/validation/EXEC40_README.md), "
        f"[machine results]({out}/results.json), [golden reference]({golden_path}). "
        'The golden retains every prior grid result, explicitly scopes the new values to this '
        'single EJ204 center cell, and remains `ready_for_acceptance=false`.\n',
        'CSV, ROOT and metadata sidecars accompany the event ledger, exact regression, first-angle '
        'histogram and event counts, first-volume-pair audit, sensor/event ledger, channel summary, '
        'branch-size audit and validation-status table. Figures share the identically named '
        'CSV/ROOT/meta sidecars; artifact hashes are in `ARTIFACTS.json`. Original simulation ROOT '
        'and execution metadata remain intact.\n',
        '**STOP: no 21-cell grid has been launched. Separate authorization is required.**\n']
    # Figures use only already computed sidecars; no new physical fits or selections.
    h=np.genfromtxt(out/'first_encounter_cosine.csv',delimiter=',',names=True)
    fig,ax=plt.subplots(figsize=(7,4))
    centers=(h['bin_low']+h['bin_high'])/2
    ax.errorbar(centers,h['observed_probability'],yerr=h['probability_se'],fmt='o',label='First scintillation encounter')
    ax.plot(centers,h['expected_probability'],label='Isotropic flux: p(μ)=2μ')
    ax.set(xlabel='Signed cosine of incidence',ylabel='Probability per 0.05 bin',title='EXEC40: EJ204, EndTop70, x=0, N=2000')
    ax.legend();fig.tight_layout();fig.savefig(out/'first_encounter_cosine.png',dpi=180);plt.close(fig)
    h=np.genfromtxt(out/'sipm_channel_summary.csv',delimiter=',',names=True)
    fig,ax=plt.subplots(figsize=(8,4))
    ax.errorbar(h['global_id'],h['detection_ratio'],yerr=h['ratio_se'],fmt='.',label='SD detections / accepted incidents (closure fails) ± SE')
    ax.axhline(.619397,color='C1',label='Emission reference 0.619397')
    ax.axvline(7.5,color='gray',alpha=.5);ax.axvline(15.5,color='gray',alpha=.5)
    ax.set(xlabel='Global SiPM ID (0–7 left, 8–15 right, 16–85 TOP)',ylabel='Detection fraction',title='EXEC40: EJ204, EndTop70, x=0, N=2000')
    ax.legend();fig.tight_layout();fig.savefig(out/'sipm_channel_summary.png',dpi=180);plt.close(fig)
    for stem in ('first_encounter_cosine','sipm_channel_summary'):
        path=out/f'{stem}.meta.json'; m=json.loads(path.read_text())
        m['figure_sha256']=sha(out/f'{stem}.png');m['figure_script_sha256']=sha(Path(__file__))
        dump(path,m)
    a.report.write_text('\n'.join(lines)+'\n')
    artifacts={str(path):dict(sha256=sha(path),bytes=path.stat().st_size) for path in out.iterdir()
        if path.is_file() and path.name!='ARTIFACTS.json'}
    artifacts[str(a.report)]=dict(sha256=sha(a.report),bytes=a.report.stat().st_size)
    artifacts[str(golden_path)]=dict(sha256=sha(golden_path),bytes=golden_path.stat().st_size)
    dump(out/'ARTIFACTS.json',artifacts)
    print(a.report)
    print(golden_path)

if __name__=='__main__':
    main()
