# EXEC_46 Step 1 — inventory, schema, and integrity checkpoint

Date: 2026-09-15

Branch: `diag/exec46-track-mechanism-20260915`

Input: `/home/rrios/exec46_20260915/full_grid`

## Verdict

**PASS.** All 21 cells and 1,216,222,416 detected-photon rows were read. The
schema, ROOT hashes, event counts, reference END yields, clock identity,
physical-integrity constraints, source labels, and revised material-specific
decay gate all pass. No simulation was run and Step 2 has not started.

The approved decay gate selects `source_type == 1` and
`abs(x_creation_mm-gun_x_mm) < 0.1 mm`. Subtracting `gun_x_mm` applies the
local definition at every position and equals `abs(x_creation_mm) < 0.1 mm`
at x=0. The configured decay constant is read independently from
OPSC-100/101/106. The untruncated mean excess is evaluated above 3, 4, and 5
decay constants. The tolerance of 3% is a declared analysis assumption; every
cut must contain at least 200 events.

## Reproducibility and inventory

```bash
/usr/bin/time -f 'WALL=%e MAXRSS_KB=%M EXIT=%x' \
  python3 analysis/track_mechanism_20260915/analyze_step1.py \
  --processes 4 --sha-file /tmp/exec46_step1_sha256.txt
```

The final read-only scan ran from `2026-09-15T20:04:43Z` through `20:09:33Z`:
290.17 s (4.84 min). Every cell has N=10,000 and seeds
`[26092601, 8349041]`. Total input size is 257,862,178,782 bytes. All 21
recomputed ROOT hashes match `.DONE`; paths, sizes, mtimes, hashes, counts, and
metrics are in `exec46_inventory.csv`.

Every file has the expected 23-branch `sipm_hits` schema. The key
`(event_id, track_id)` is unique. END-left is `global_id` 0–7, END-right 8–15,
and TOP 16–85. A TOP description as IDs 16–35 is incomplete.

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

Every event has at least one hit at each END. The nine END-yield references at
x=0 and ±650 mm reproduce to their quoted 0.01 PE precision.

## Clock and physical-integrity gates

The clock test passes separately in every cell:

| Material | x values [mm] | hits | max abs delta [ns] | RMS [ns] | nonzero fraction |
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

The per-cell hit counts in the preceding table are the independent row counts
used for each zero result. The baseline and decomposition clocks are therefore
bit-identical across the grid.

All eight aggregate violation counts are zero: creation after detection, path
shorter than microscopic chord, negative boundary count, apparent speed above
configured group speed, gun-position mismatch, duplicate key, invalid source
type, and invalid sensor map. All three RINDEX tables are constant at 1.58,
giving a Geant4 11.4 group
speed of 189.742062 mm/ns.

## Revised decay gate

| Cell | tau decay/rise [ns] | fit 3 tau | fit 4 tau | fit 5 tau | max abs rel [%] |
|---|---:|---:|---:|---:|---:|
| EJ200_xm650 | 2.1/0.9 | 2.100167 | 2.100496 | 2.101505 | 0.072 |
| EJ200_xm500 | 2.1/0.9 | 2.100000 | 2.101999 | 2.101672 | 0.095 |
| EJ200_xm200 | 2.1/0.9 | 2.100215 | 2.100277 | 2.100781 | 0.037 |
| EJ200_xp0 | 2.1/0.9 | 2.102404 | 2.103800 | 2.105881 | 0.280 |
| EJ200_xp200 | 2.1/0.9 | 2.099582 | 2.100614 | 2.097814 | 0.104 |
| EJ200_xp500 | 2.1/0.9 | 2.101833 | 2.102314 | 2.100459 | 0.110 |
| EJ200_xp650 | 2.1/0.9 | 2.100575 | 2.098388 | 2.095481 | 0.215 |
| EJ204_xm650 | 1.8/0.7 | 1.798907 | 1.797293 | 1.797650 | 0.150 |
| EJ204_xm500 | 1.8/0.7 | 1.799980 | 1.799611 | 1.796020 | 0.221 |
| EJ204_xm200 | 1.8/0.7 | 1.801433 | 1.803463 | 1.802306 | 0.192 |
| EJ204_xp0 | 1.8/0.7 | 1.800824 | 1.799688 | 1.800572 | 0.046 |
| EJ204_xp200 | 1.8/0.7 | 1.800596 | 1.801800 | 1.804661 | 0.259 |
| EJ204_xp500 | 1.8/0.7 | 1.800505 | 1.800778 | 1.801743 | 0.097 |
| EJ204_xp650 | 1.8/0.7 | 1.799618 | 1.799405 | 1.797233 | 0.154 |
| EJ230_xm650 | 1.5/0.5 | 1.499851 | 1.499031 | 1.499091 | 0.065 |
| EJ230_xm500 | 1.5/0.5 | 1.501120 | 1.503069 | 1.502976 | 0.205 |
| EJ230_xm200 | 1.5/0.5 | 1.500863 | 1.497905 | 1.496892 | 0.207 |
| EJ230_xp0 | 1.5/0.5 | 1.500736 | 1.500922 | 1.498058 | 0.129 |
| EJ230_xp200 | 1.5/0.5 | 1.500295 | 1.500205 | 1.503520 | 0.235 |
| EJ230_xp500 | 1.5/0.5 | 1.500965 | 1.500604 | 1.500607 | 0.064 |
| EJ230_xp650 | 1.5/0.5 | 1.500085 | 1.499008 | 1.497896 | 0.140 |

All 63 local tests pass. Exact event-cluster errors and the full-population
informational scans at 3, 4, 5, 8, 12, and 20 decay constants are stored in
`exec46_tau_diagnostics.csv`. The 3.098019 ns full END-pool mean is rejected as
an exponential benchmark and is not used.

## Source census

| Cell | scint all [%] | Cher all [%] | first Cher L [%] | first Cher R [%] |
|---|---:|---:|---:|---:|
| EJ200_xm650 | 96.566 | 3.434 | 74.10 | 2.81 |
| EJ200_xm500 | 96.452 | 3.548 | 8.44 | 3.39 |
| EJ200_xm200 | 96.463 | 3.537 | 5.87 | 4.65 |
| EJ200_xp0 | 96.441 | 3.559 | 4.91 | 5.15 |
| EJ200_xp200 | 96.461 | 3.539 | 4.57 | 5.98 |
| EJ200_xp500 | 96.445 | 3.555 | 3.33 | 8.90 |
| EJ200_xp650 | 96.575 | 3.425 | 3.16 | 73.56 |
| EJ204_xm650 | 96.885 | 3.115 | 62.90 | 2.20 |
| EJ204_xm500 | 96.792 | 3.208 | 7.44 | 2.62 |
| EJ204_xm200 | 96.802 | 3.198 | 4.89 | 3.43 |
| EJ204_xp0 | 96.798 | 3.202 | 4.37 | 4.30 |
| EJ204_xp200 | 96.810 | 3.190 | 3.65 | 4.98 |
| EJ204_xp500 | 96.793 | 3.207 | 2.86 | 7.64 |
| EJ204_xp650 | 96.870 | 3.130 | 2.20 | 63.42 |
| EJ230_xm650 | 96.657 | 3.343 | 51.77 | 1.75 |
| EJ230_xm500 | 96.584 | 3.416 | 7.18 | 2.06 |
| EJ230_xm200 | 96.601 | 3.399 | 4.28 | 3.15 |
| EJ230_xp0 | 96.607 | 3.393 | 3.76 | 3.76 |
| EJ230_xp200 | 96.603 | 3.397 | 3.07 | 4.75 |
| EJ230_xp500 | 96.586 | 3.414 | 1.91 | 7.02 |
| EJ230_xp650 | 96.661 | 3.339 | 1.84 | 52.33 |

Primary and other fractions are zero. Despite only 3.1–3.6% Cherenkov hits in
the total population, Cherenkov photons dominate the first detected photon at
the near END. This directly supports treating mechanisms 4 and 5 as a coupled
mixture in later steps.

## Independent Step 4 derivation

The requested `N^-1/2` exponent is correct. The actual Geant4 11.4 sampler
(`G4Scintillation.cc:557-570`) has

```text
f_G4(t) = (tau_r+tau_d)/tau_d^2 * exp(-t/tau_d) * [1-exp(-t/tau_r)]
E[min_N] -> tau_d*sqrt(pi*tau_r/[2*(tau_r+tau_d)])*N^-1/2.
```

Its prefactors are 1.441584, 1.193745, and 0.939986 ns for EJ-200/204/230.
The standard-biexponential expression in the approval gives 1.723022,
1.406842, and 1.085402 ns and will remain a separate reference curve. The
derivation and numerical survival-integral check are in
`verify_order_stat_asymptotic.py` and
`exec46_order_stat_asymptotic.{csv,json}`. The accepted empirical
`t_creation_ns` distribution remains the Step 4 prediction.

No figure was required or produced, so figure sidecars do not apply. No push,
merge, reset, rebase, deck edit, or simulation was performed.
