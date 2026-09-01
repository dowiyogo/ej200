# Cross-Branch Physics Coherence

Legend: yes = coherent with canonical value; no = inconsistent; stale = branch
is behind the current `main` physics baseline and should inherit it by rebase.

| Parameter | Canonical | `main` | `feature/edge-scan-and-readout-grouping` | `feature/sipm-electronics-response` |
|---|---|---|---|---|
| EJ-200 rise time | 0.9 ns (`EJ200.pdf`) | **no: absent** | **no: absent** | **no: absent for EJ-200** |
| EJ-200 decay | 2.1 ns (`EJ200.pdf`) | yes: 2.1 ns | yes: 2.1 ns | yes: 2.1 ns |
| EJ-200 yield | 10000 ph/MeV (`EJ200.pdf`) | yes | yes | yes |
| EJ-200 attenuation length | 3.8 m (`EJ200.pdf`) | yes | yes | yes |
| EJ-200 RINDEX | 1.58 (`EJ200.pdf`) | yes | yes | yes |
| EJ-200 emission peak | 425 nm (`EJ200.pdf`) | yes | yes | yes |
| Reflector model | `main` baseline: explicit panels, R=0.98 validated fallback | yes | **stale: Mylar segmented wrap** | **stale: Mylar segmented wrap** |
| Mylar active volume | absent in current baseline | yes: absent | **no: present** | **no: present** |
| PDE comment SiPM | S13360-6050PE (Betancourt 2020 / W4) | **no: says S13360-6025** | **no: says S13360-6025** | **no: says S13360-6025 legacy** |
| PDE peak for prototype curve | ~0.405 @ 460 nm | yes | yes | yes for S13360 curve |
| Default sensor | S13360-6050PE | yes by fixed table, though mislabeled | yes by fixed table, though mislabeled | **no: Broadcom default** |
| Jitter model | SPTR 106 ps + FastIC+ 10 ps in quadrature | **no: 20 ps monolithic** | **no: 20 ps monolithic** | **no: SPTR/electronics split exists but defaults are 150 ps / 30 ps and electronics not applied** |
| Jitter UI | `/sipm/sptrSigma`, `/sipm/electronicsJitter`, legacy `/sipm/jitterSigma` alias | **no: only `/sipm/jitterSigma`** | **no: only `/sipm/jitterSigma`** | partial: `/sipm/sptrSigma`, `/sipm/electronicsSigma`; no legacy alias |
| Bar half X | 700 mm | yes | yes | yes |
| Bar half Y | 30 mm | yes | yes | yes |
| Bar half Z | 5 mm | yes | yes | yes |
| End SiPM active size | 6 x 6 mm2 | yes | yes | yes |
| End SiPM count | 8 per side | yes | yes | yes |
| Top SiPM default | 20 at 70 mm pitch | yes | yes | yes |
| Primary particle | `mu-` | yes | yes | yes |
| Primary energy | 1 GeV | yes | yes | yes |
| Primary position/direction | `(0,0,+60 mm)`, `(0,0,-1)` | yes | yes | yes |
| `/muon/gunX` | present | yes | yes | yes |
| `gun_x_mm` ntuple | present | yes | yes | yes |
| Canonical ntuple branches | 12 canonical branches | yes | yes | additive diagnostics only absent | additive but reordered indices after `time_ns` |
| EJ-230 implementation | future separate feature, not baseline | absent | absent | present and default; **scope drift** |
| Diagnostic scripts | retained in repo, generated outputs ignored | yes | stale/not present until rebase | present from old diag branch if rebased |

## Compact Findings

1. P0: all branches are missing `SCINTILLATIONRISETIME1 = 0.9 ns` for EJ-200.
2. P0/P1: all branches misidentify the prototype PDE table as S13360-6025
   rather than S13360-6050PE.
3. P1: timing smear is inconsistent in all branches. `main` and edge have
   monolithic 20 ps; electronics has a split but wrong defaults and does not
   use electronics jitter in the timestamp.
4. P1: both feature branches are geometrically stale relative to `main`:
   they still use active Mylar wrap volumes and lack the final explicit
   reflector panel baseline.
5. P1/P2: `feature/sipm-electronics-response` has valuable electronics work,
   but it changes default material and sensor to EJ-230 + Broadcom, which is
   not the canonical EJ-200 / S13360-6050PE baseline.

