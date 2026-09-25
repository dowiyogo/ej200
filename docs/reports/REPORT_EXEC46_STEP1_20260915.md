# EXEC_46 Step 1 — inventory, schema, and integrity checkpoint

Date: 2026-09-15

Branch: `diag/exec46-track-mechanism-20260915`

Input: `/home/rrios/exec46_20260915/full_grid`

## Verdict

**PASS.** All 21 cells and 1,216,222,416 detected-photon rows were read. All
21 ROOT hashes match `.DONE`; each file has the expected 23-branch schema and
10,000 events. The key `(event_id, track_id)` is unique. The reference END
yields, clock identity, physical-integrity constraints, source labels, and the
approved material-specific decay gate all pass. No simulation was run and
Step 2 has not started.

The gate selects `source_type == 1` and
`abs(x_creation_mm-gun_x_mm) < 0.1 mm`. This is the local coordinate at every
gun position and equals the literal `abs(x_creation_mm) < 0.1 mm` at x=0. The
decay constant is read independently from OPSC-100/101/106; the untruncated
mean excess is tested above 3, 4, and 5 decay constants. The 3% tolerance is a
declared analysis assumption and each cut must contain at least 200 events.

## Exact run

```bash
/usr/bin/time -f 'WALL=%e MAXRSS_KB=%M EXIT=%x' \
  python3 analysis/track_mechanism_20260915/analyze_step1.py \
  --processes 4 --sha-file /tmp/exec46_step1_sha256.txt
```

Final start/end: `2026-09-15T20:04:43Z` / `20:09:33Z`; wall time 290.17 s
(4.84 min). Every cell uses N=10,000 and seeds `[26092601, 8349041]`.
Total ROOT size is 257,862,178,782 bytes. Exact paths, sizes, mtimes, hashes,
counts, and per-cell metrics are in
`analysis/track_mechanism_20260915/exec46_inventory.csv`.

## Photon budget

Every event has at least one hit on both END faces. Values below are detected
PE per event.

| Cell | photons | Npe L | Npe R | Npe END | Npe TOP |
|---|---:|---:|---:|---:|---:|
| EJ200_xm650 | 67,720,057 | 2345.3235 | 243.6909 | 2589.0144 | 4182.9913 |
| EJ200_xm500 | 63,997,268 | 1379.9831 | 277.5710 | 1657.5541 | 4742.1727 |
| EJ200_xm200 | 63,348,852 | 743.1841 | 396.9631 | 1140.1472 | 5194.7380 |
| EJ200_xp0 | 62,596,038 | 528.5215 | 528.2879 | 1056.8094 | 5202.7944 |
| EJ200_xp200 | 63,261,381 | 396.4248 | 741.9157 | 1138.3405 | 5187.7976 |
| EJ200_xp500 | 63,932,058 | 277.5965 | 1378.6384 | 1656.2349 | 4736.9709 |
| EJ200_xp650 | 68,056,391 | 244.9545 | 2357.5863 | 2602.5408 | 4203.0983 |
| EJ204_xm650 | 63,710,913 | 2405.7183 | 140.7303 | 2546.4486 | 3824.6427 |
| EJ204_xm500 | 58,285,995 | 1303.6846 | 169.1917 | 1472.8763 | 4355.7232 |
| EJ204_xm200 | 56,379,503 | 609.1654 | 273.4047 | 882.5701 | 4755.3802 |
| EJ204_xp0 | 55,640,955 | 398.1965 | 398.0877 | 796.2842 | 4767.8113 |
| EJ204_xp200 | 56,378,895 | 273.5438 | 608.7005 | 882.2443 | 4755.6452 |
| EJ204_xp500 | 58,336,461 | 169.6250 | 1304.9172 | 1474.5422 | 4359.1039 |
| EJ204_xp650 | 63,448,319 | 140.3753 | 2395.2265 | 2535.6018 | 3809.2301 |
| EJ230_xm650 | 54,935,742 | 2158.1344 | 91.4689 | 2249.6033 | 3243.9709 |
| EJ230_xm500 | 49,643,803 | 1124.1163 | 113.8181 | 1237.9344 | 3726.4459 |
| EJ230_xm200 | 47,447,287 | 486.4227 | 197.8695 | 684.2922 | 4060.4365 |
| EJ230_xp0 | 46,946,914 | 304.1811 | 304.3238 | 608.5049 | 4086.1865 |
| EJ230_xp200 | 47,421,723 | 197.7986 | 486.4402 | 684.2388 | 4057.9335 |
| EJ230_xp500 | 49,733,221 | 114.0377 | 1126.6588 | 1240.6965 | 3732.6256 |
| EJ230_xp650 | 55,000,640 | 91.4549 | 2160.5723 | 2252.0272 | 3248.0368 |

The nine references at x=0 and ±650 mm reproduce to 0.01 PE precision.

## Clock gate on every cell

For every row, `time_ns-t_detection_ns` is bit-identically zero.

| Material | x tested [mm] | hits tested | max abs delta [ns] | RMS [ns] | nonzero fraction |
|---|---|---:|---:|---:|---:|
| EJ-200 | -650,-500,-200,0,200,500,650 | 452,912,045 | 0 | 0 | 0 |
| EJ-204 | -650,-500,-200,0,200,500,650 | 412,181,041 | 0 | 0 | 0 |
| EJ-230 | -650,-500,-200,0,200,500,650 | 351,129,330 | 0 | 0 | 0 |

Per-cell clock rows:

| Cell | hits tested | max abs delta [ns] | RMS [ns] | nonzero fraction |
|---|---:|---:|---:|---:|
| EJ200_xm650 | 67,720,057 | 0 | 0 | 0 |
| EJ200_xm500 | 63,997,268 | 0 | 0 | 0 |
| EJ200_xm200 | 63,348,852 | 0 | 0 | 0 |
| EJ200_xp0 | 62,596,038 | 0 | 0 | 0 |
| EJ200_xp200 | 63,261,381 | 0 | 0 | 0 |
| EJ200_xp500 | 63,932,058 | 0 | 0 | 0 |
| EJ200_xp650 | 68,056,391 | 0 | 0 | 0 |
| EJ204_xm650 | 63,710,913 | 0 | 0 | 0 |
| EJ204_xm500 | 58,285,995 | 0 | 0 | 0 |
| EJ204_xm200 | 56,379,503 | 0 | 0 | 0 |
| EJ204_xp0 | 55,640,955 | 0 | 0 | 0 |
| EJ204_xp200 | 56,378,895 | 0 | 0 | 0 |
| EJ204_xp500 | 58,336,461 | 0 | 0 | 0 |
| EJ204_xp650 | 63,448,319 | 0 | 0 | 0 |
| EJ230_xm650 | 54,935,742 | 0 | 0 | 0 |
| EJ230_xm500 | 49,643,803 | 0 | 0 | 0 |
| EJ230_xm200 | 47,447,287 | 0 | 0 | 0 |
| EJ230_xp0 | 46,946,914 | 0 | 0 | 0 |
| EJ230_xp200 | 47,421,723 | 0 | 0 | 0 |
| EJ230_xp500 | 49,733,221 | 0 | 0 | 0 |
| EJ230_xp650 | 55,000,640 | 0 | 0 | 0 |

All violation counts are zero for creation after detection, path shorter than
the microscopic chord, negative boundary count, apparent speed above group
speed, gun-position mismatch, duplicate key, invalid source, and invalid sensor
map. The configured
RINDEX is 1.58 for every material, giving 189.742062 mm/ns.

## Decay gate

Configured decay/rise values are 2.1/0.9 ns (EJ-200), 1.8/0.7 ns (EJ-204),
and 1.5/0.5 ns (EJ-230). All 63 local tests pass. The worst absolute relative
differences are 0.280% (`EJ200_xp0`, 5 tau), 0.259% (`EJ204_xp200`, 5 tau),
and 0.235% (`EJ230_xp200`, 5 tau), below 3%. Full per-cell estimates and
event-cluster errors are in `exec46_tau_diagnostics.csv`.

The complete-pool scans at 3, 4, 5, 8, 12, and 20 tau remain informational.
The 3.098019 ns full END-pool mean is rejected as an exponential benchmark and
is not used.

## Source census

All-hit Cherenkov fractions are only 3.1–3.6%, but near-END first-photon
fractions are 74.10%/73.56% for EJ-200, 62.90%/63.42% for EJ-204, and
51.77%/52.33% for EJ-230 (near left/near right). At x=0 the left/right values
are 4.91%/5.15%, 4.37%/4.30%, and 3.76%/3.76%. Exact results at all 21 cells
are in `exec46_source_census.csv`; primary and other fractions are zero. This
shows why scintillation order statistics and Cherenkov composition must be
analyzed as a coupled mixture.

## Independent order-statistics check

The requested `N^-1/2` exponent is verified. For the actual Geant4 11.4
rejection sampler (`G4Scintillation.cc:557-570`),

```text
f_G4(t) = (tau_r+tau_d)/tau_d^2 * exp(-t/tau_d) * [1-exp(-t/tau_r)]
E[min_N] -> tau_d*sqrt(pi*tau_r/[2*(tau_r+tau_d)])*N^-1/2.
```

Its prefactors are 1.441584, 1.193745, and 0.939986 ns for EJ-200/204/230.
The standard-biexponential expression in the approval gives 1.723022,
1.406842, and 1.085402 ns and is retained as a separate reference curve. The
derivation and numerical integration are in
`exec46_order_stat_asymptotic.{csv,json}`. The empirical `t_creation_ns`
distribution remains the accepted prediction for Step 4.

The detailed canonical report is
`/home/rrios/ej200/analysis/track_mechanism_20260915/REPORT_STEP1.md`.
No figure was required or produced. No push, merge, reset, rebase, deck edit,
or simulation was performed. Step 2 awaits explicit checkpoint approval.
