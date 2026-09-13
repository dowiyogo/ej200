# EXEC36 test-beam comparison: scope and unresolved configuration

The historical numbers and experimental statements below are supplied by René
in the EXEC36 task attachment, not newly measured or independently verified:
`/home/rrios/.codex/attachments/67ebef39-e056-449c-8014-d680670f1a80/pasted-text.txt`.
They remain unchanged. This analysis implements a same-end group difference;
matching that observable does not establish matching hardware response.

| Reference quantity | Supplied test-beam value / scope |
| --- | --- |
| Raw same-end difference | sigma(T1−T2) approximately 125 ps |
| Same-end difference divided by sqrt(2) | approximately 88.4 ps; assumes equal independent group contributions |
| Opposite-end difference | Not measured as that test-beam reference; the older simulated DeltaT_LR is distinct |
| Isolated channel | 200–250 ps; the task requests a single **physical** SiPM check while SiPMs per FastIC input remain unresolved |
| Material and dimensions | EJ-204, 1400×60×10 mm; primary comparison must use EJ-204 |

The single-channel experimental wording is an additional scope limitation:
under V2 a FastIC input receives two physical SiPMs, whereas Block 1 deliberately
tests one physical SiPM as requested. No interpretation of that ambiguity is
selected. V1 has two physical SiPMs per output group; V2 has four. Their choice
is pending experimental confirmation. Approximately sqrt(2) is the supplied
expected width scale for changing independent sensor counts, not an enforced
simulation correction or a verified law for threshold crossings.

The explicit provisional mapping is in `exec36_config.json`: input pairs
0+2 -> CH1 and 4+6 -> CH5; V2 physical groups match the preserved
`GroupIndex` (`upstream/related/420addf/analysis/congruent_sum4_timing.C:216–226`).
V1 retains local IDs 0,2 and 4,6, leaving 1,3,5,7 unused. Right adds 8.
`src/DetectorConstruction.cc:407–419` establishes the local-ID ordering along Y
and negative/positive-X ends; no measured cable map is inferred from geometry.

## Three discrepancies, recorded without correction

| Issue | Test-beam statement supplied in task | Simulation / local documentary evidence | Consequence and status |
| --- | --- | --- | --- |
| Sensor / PDE | “FBK NUV-MT 14M”, OV=10 V; meeting wording is “SiPMs Broadcom (FBK NUV-MT 14M)” | `data/sipm/AFBR-S4N66P024M_pde.txt:1–4` explicitly says Broadcom AFBR-S4N66P024M, curve at 12 V above breakdown; `ELECTRONICS_PROVENANCE.md` records related device/OV conflicts | Manufacturer/device/OV identification remains unresolved. This directly affects the photon budget predicted by the simulation; no PDE is changed |
| Time reference | Older CH8 FastIC trigger included 20–50 ps paddle jitter; recent Constanza analyses use detector channel differences | New same-end observable uses T1−T2 without CH8; isolated-SiPM control uses gun t=0 | Common-reference cancellation is part of the observable definition; no trigger jitter is injected or subtracted numerically |
| Quantization | 24.4 ps TDC LSB; 24.4/sqrt(12) approximately 7.0 ps in quadrature for each timestamp | No quantization in the preserved END primitive; `congruent_sum4_timing.C:50–51` has a different 24 ps metadata-only label. Earlier electronics docs also retain a 25 ps label | Record all versions, do not unify or inject quantization. The per-timestamp term is not automatically the raw two-channel difference contribution |

## Reflector nomenclature — not a configuration discrepancy

René's subsequent clarification supersedes the original task's wrap description:
the physical reflector is, and always was, **Vikuiti 3M ESR**. “Mylar” is an
inherited code name: volumes are `Vikuiti<panel>LV/PV` and the factory is
`CreateMylarReflector`. **R=0.98 is the correct specification of the real
reflector**, retained at `src/DetectorConstruction.cc:313`. This is a naming
clarification, not a simulation/test-beam material discrepancy.

Commit `522a1d0816e606d548f7042d992f71136c3450d4` documents the historical
divergence between volume names and the then-current model in comments. Its
description of that older implementation is not a statement about today's
reflectivity. The physical-material identification above follows René's explicit
clarification. No reflectivity, physics, data or numerical analysis is changed;
the earlier wrap discrepancy and proposed Mylar reflectivity comparison are withdrawn.

## Sensitivity study is not electronics injection

Rise/fall 0.5/5 ns and threshold 4 PE have unverified physical calibration
(`ELECTRONICS_PROVENANCE.md`, END rows). The task supplies a bank DCR threshold
scan Th=35–40 DAC, corresponding to 2–4 PE. EXEC36 varies thresholds 2,3,4 PE
and rise/fall 0.5/5 versus 2/3 ns on EJ204_xp0 only. The peak-normalized
difference-of-exponentials shape is retained; the 2/3 ns pair is the bank-derived
pair recorded in `analysis/timing/pulse_models.py:49–63`, still awaiting formal
full-chain calibration. No SPTR, electronics jitter or TDC quantization is added.
No parameter is selected because it approaches the test-beam number.

The gun reference uses the native arrival times relative to t=0. The installed
Geant4 source `source/event/src/G4ParticleGun.cc:66,219` initializes particle_time
to 0 and supplies it to the primary vertex. `src/PrimaryGeneratorAction.cc` does
not override it; valid campaign macros contain no `/gun/time` change. A constant
time-origin shift would not alter widths, but no time-origin fitting is applied.

The 21 cells use the same seeds 26092601 and 8349041; paired points within a
material are correlated. Right-end results are a symmetry check, with opposite
end at mirrored x as the geometry-matched comparison; left versus right at the
same nonzero x need not be equal. No cross-cell fit, BLUE, arm average or
experimental calibration follows these calculations.
