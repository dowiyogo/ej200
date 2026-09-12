# EXEC35 electronics provenance register

Documentary audit, 2026-09-12. No electronics was injected and no timing width
was recalculated with jitter. All live versions remain separate. `established`
means the stated definition/value and its scope are explicit in the cited local
record, **not** a new experimental validation or confirmation that it applies to
this detector. `unverified` means the physical provenance or transfer is missing;
`contradictory` means incompatible claims or, for SPTR, a source that fails to
distinguish FWHM from sigma. Published measurements below are what the repository
records; the original papers were not independently revalidated in this audit.

Paths below are relative to the worktree root. Short source names expand as follows:

| Key | Exact source |
| --- | --- |
| P | `analysis/timing/pulse_models.py` |
| W | `analysis/timing/sipm_waveform_dcfd.py` |
| S | `docs/branch_diagnosis/SPTR_PROVENANCE.md` |
| E | `docs/branch_diagnosis/ELECTRONICS_PARAMETERS.md` |
| X17 | `analysis/sigma_t/upstream/analysis_core/exec17_corrections.py` |
| X18 | `analysis/sigma_t/upstream/analysis_core/exec18_main.py` |
| X19 | `analysis/sigma_t/upstream/analysis_core/exec19_main.py` |
| C | `analysis/sigma_t/upstream/related/420addf/analysis/congruent_sum4_timing.C` |

## SPTR and operating conditions

| Parameter / live version | Value and width convention | Device / operating point | Evidence (file:line) | Status and scope |
| --- | --- | --- | --- | --- |
| EXEC17 SPTR_ANALOG_PS | 106 ps; used as sigma in quadrature, original FWHM/sigma not established | Comment attributes OPSC-101 SiPM model, no OV; later catalog attributes another device | X17:38,404–407; P:149–161 | **contradictory**: executable sigma usage does not establish measurement units or Broadcom provenance |
| EXEC18 SPTR_PS | 106 ps; sigma usage, unresolved measured convention | No OV in constant; catalog says Hamamatsu S13360-3050CS, 3–4.5 V OV | X18:38,580–587; P:149–161 | **contradictory** |
| EXEC19 SPTR_PS | 106 ps; sigma usage, unresolved measured convention | Same unidentified transfer to Broadcom | X19:34,217,269–270; P:149–161 | **contradictory** |
| Collaboration catalog | `sigma_ps=106`; notes also say 100–120 ps **FWHM** and ~106 ps **RMS** | Hamamatsu S13360-3050CS; prose 3–4.5 V OV, `ov_v=None`; not detector Broadcom | P:148–161; E:121–129 | **contradictory**: these conventions are explicitly incompatible; no choice made |
| Lee intrinsic FWHM | 137 ±4 ps **FWHM** | AFBR-S4N66P014M, 6×6 mm²; bias 48 V, OV approximately 15.5 V; LNHF front end | S:24–40; P:70,77,105–111 | **established** as local record of that measurement; transfer to 024M at another OV unverified |
| Lee detector FWHM | 172 ps **FWHM**, includes Lee electronics | Same 014M and OV; includes Lee front end, not FastIC+ | S:30–40; P:78,113–119 | **established** within stated scope; not a SiPM-only jitter to combine again with electronics |
| Lee intrinsic equivalent sigma | 137/2.355, approximately 58.2 ps **sigma** | Same device and OV, Gaussian-equivalent conversion | P:22,81,108; S:46–50 | **established** algebraic convention conditional on Gaussian conversion |
| Lee detector equivalent sigma | 172/2.355, approximately 73.0 ps **sigma** | Same device and OV; electronics still included | P:22,82,116; S:46–51 | **established** algebraic convention; does not remove electronics |
| Small-area comparison in SPTR record | 45 ±1 ps **FWHM** | Same cited NUV-MT family, 2×2 mm², bias 48 V; separate OV not given | S:61–65 | **unverified** for detector use; units explicit, device size and OV transfer unavailable |
| Intermediate-area comparison | 55 ±1 ps **FWHM** | Approximately 4×4 mm², bias 48 V; separate OV not given | S:61–66 | **unverified** for detector use; not the 6×6 mm² device |
| `banco_propio_low` (7be59bd) | 85 ps **sigma/RMS**, lower derived bound | 014M + FastIC+, catalog OV 10 V; bias 38 V and threshold code 35 asserted | P:121–135; E:49 | **unverified**: explicitly derived, formal event-level EOS measurement pending; operating-voltage conflict below |
| `banco_propio_high` (7be59bd) | 92 ps **sigma/RMS**, upper derived bound | 014M + FastIC+, same claimed OV 10 V | P:137–146; E:50 | **unverified**: derived range, not an independent SPTR measurement |
| Bank central OLS value | 88.8, approximately 89 ps **RMS** from slope 7877 ps² | Same bank; three intensity points, some midpoint inputs | P:97–102; E:53–56 | **unverified**: model-derived from summary points; no new OLS calculation here |
| Formal target SPTR | No value established | AFBR-S4N66P024M at OV 10 V | E:51 | **unverified**, pending measurement; no default invented |
| Old 200 ps value mentioned by waveform CLI | 200 ps, FWHM versus sigma **not specified in that attribution** | SiPM of Lv et al.; CLI says not AFBR-S4N66P024M; exact device and OV not given there | W:143 | **contradictory** width convention; provenance/device/OV additionally unverified; not used |
| Old standalone ~100 ps assertion | ~100 ps, convention unspecified | Attributed to detector without documentary support | S:9–12 | **contradictory** width convention; source records removal, not a live calibration |
| Detector/device identity transfer | 014M described as optically identical to 024M; package difference only | Cross-device equivalence claim, not a same-device SPTR measurement | S:24–27; E:4 | **unverified** as justification for transferring SPTR/calibration |
| Published operating point | Bias 48 V, breakdown approximately 32.5 V, OV approximately 15.5 V | 014M, laser 408 nm/~15 ps, approximately 20 C, LNHF approximately 2 GHz | S:29–33 | **established** local record of cited conditions |
| Bank OV/bias | OV 10 V and bias approximately 38 V versus breakdown 32.5 V and bank bias approximately 42.5 V | Bank device attribution 014M/024M | E:4–5,106–108; P:123,126 | **contradictory**: voltage relations disagree within the documentation; no inferred replacement |
| SPTR operating-limit wording | “overvoltage below approximately 40 V” beside OV 15.5 V and maximum OV 16 V | Published setup description | S:73–77 | **contradictory** voltage terminology; likely bias/OV confusion is not silently repaired |
| PDE operating point | Curve at 12 V above breakdown | AFBR-S4N66P024M simulation PDE | S:104–109; `data/sipm/AFBR-S4N66P024M_pde.txt:3` | **established** curve annotation; consistency with a 10 V electronics model remains unverified |

The /2.355 conversion is a Gaussian FWHM-to-sigma conversion, not a universal
identity for asymmetric distributions (S:46). S:85–98 explicitly warns that
`SPTR/sqrt(k*N_active)` is an independent-mean formula and is not established
for a kth order statistic. X17:405, X18:586 and X19:269 retain their separate
mean-model expressions; none is adopted for EXEC35 TOP. No alternative
propagation law is fabricated.

## Pulse shape and electronics versions

| Parameter / version | Value and units | Evidence (file:line) | Status and interpretation |
| --- | --- | --- | --- |
| EXEC17 FastIC constant | 10 ps, used as sigma in quadrature | X17:39,406–407 | **unverified**: “from SHiP electronics” does not give a measurement or operating point |
| EXEC18 FastIC constant | 10 ps, sigma usage | X18:39,582,587 | **unverified**, same provenance gap |
| EXEC19 FastIC constant | 10 ps, sigma usage | X19:35,217,270 | **unverified**, no proof that this is total FastIC+ jitter |
| Peña-Rodríguez shortened rise | tau_r=2 ns | P:25–35; E:24 | **established** local attribution to another shortened readout setup; transfer to detector unverified |
| Peña-Rodríguez shortened fall | tau_f=3 ns | P:25–35; E:25 | **established** same setup attribution; not intrinsic SiPM fall |
| Broadcom intrinsic rise | None / unmeasured, unit ns | P:37–46; E:26 | **unverified**, requires measurement, no substituted 2 ns |
| Broadcom intrinsic fall | tau_f=55 ns | P:39–42; E:27 | **established** local DS105 attribution; cannot mix with another setup's rise constant |
| FastIC bank rise (7be59bd) | tau_r=2 ns | P:49–63; E:28,35–37 | **unverified** for full chain: scope claims 1 GHz, 6.4 GS/s oscilloscope and shortened response, expressly requires formal FastIC+ confirmation |
| FastIC bank fall (7be59bd) | tau_f=3 ns | P:49–63; E:29,35–37 | **unverified** for full chain, same pending confirmation |
| END SUM4 rise | 0.5 ns | C:46 | **unverified** physical calibration; **implementation fixed** for this analysis |
| END SUM4 fall | 5 ns | C:47 | **unverified** physical calibration; separate live model from 2/3 ns |
| END SUM4 shape | Peak-normalized difference `exp(-t/5)-exp(-t/0.5)` | C:155–167 | **established** executable definition; different from waveform dCFD mathematical shape |
| END leading-edge threshold | 4 PE amplitude units, provisional | C:49,170 onward | **unverified** physical threshold calibration; no threshold changed here |
| END TDC label | 24 ps LSB | C:51 | **unverified** hardware provenance; metadata only in this pipeline; not a quantization applied to EXEC35 arrivals |
| Register-era dCFD shape | `A*(1-exp(-t/tau_r))*exp(-t/tau_f)` | P:9–12; W:275–284 | **established** executable shape, cannot equate its rise constant to a different kernel without calibration |
| dCFD transit jitter interface | Required `transit_sigma_ps`; Gaussian **sigma**, converted to ns for per-arrival jitter | W:132–145,292,300–303 | **established** software units; selected physical SPTR remains subject to the register above |
| dCFD fraction | 0.14 of peak | W:120–123; E:89–96 | **unverified** optimum; no source calibration supplied |
| Waveform time step | 10 ps | W:126–129 | **established** numerical default, not a measured hardware sampling/jitter parameter |
| Waveform SPE amplitude spread | 0.1 PE sigma; Gaussian then clipped at zero | W:168–171,301–302 | **unverified** detector calibration; amplitude fluctuation, not time SPTR |
| Kernel tail length | 8*tau_f (at least 5*dt) | W:189–192,281–283 | **established** numerical truncation default; no physical-tail validation claimed |
| dCFD electronics interface versus use | Named ps Gaussian time jitter; converted ns value added to waveform **amplitude** bins | W:174–185,294,319–320 | **contradictory** dimensional semantics: implemented sample-amplitude noise does not match advertised time jitter; document only, no code repair or run |
| FastIC total electronics bound | <=18 ps sigma, bank-derived | E:54–55,71; W:178–181 | **unverified** as measured total jitter; derived upper bound does not establish the separate 10 ps constant |
| FastIC TDC label | 25 ps LSB | E:69; W:184 | **established** local hardware attribution; distinct from END's 24 ps label, no unification |
| TDC quantization contribution | 25/sqrt(12), approximately 7.1 ps sigma | E:70; W:184 | **established** uniform-quantization derivation, not total electronics resolution |
| CLI registry synchronization | Help says FastIC pulse values pending while registry supplies 2/3 ns | W:105,153,163 versus P:50–63 | **contradictory** live help/state; registry itself still requires physical confirmation |

## Historical snapshots remain identifiable

`9040783d1dc04176094db355ae9a1e0a0ebed8ef` introduced the separated registry.
In **that snapshot** `analysis/timing/pulse_models.py:25–35` contains the
shortened 2/3 ns pair, `:37–46` contains unknown-rise/55 ns intrinsic response,
and `:49–56` has both FastIC constants **None**. `:22,70–75` records explicit
FWHM conversion and the 137/172 ps measurements. The snapshot does not contain
the later bank 85/92 ps entries. Inspect reproducibly with:

```bash
git show 9040783:analysis/timing/pulse_models.py | nl -ba
git show 7be59bd:analysis/timing/pulse_models.py | nl -ba
```

`7be59bda5d1c6fc1e75d9ba0ae81bb7abc9333f8` supplies the bank 2/3 ns pair,
85/92 ps **sigma** bounds and explicit 106 ps ambiguity. Its relevant registry
lines match P:49–63,121–161 in the audited tree. A name such as
`fastic_measured` does not override its explicit full-chain confirmation caveat.
The END 0.5/5 ns model and EXEC17/18/19 post-hoc quadratures are still separate
live paths; the dCFD script is another path with the dimensional conflict above.

## Open questions and recommendation pending René

1. Obtain the original collaboration 106 ps measurement/convention from Gerardo:
   FWHM or sigma, device, OV, temperature, optical setup, and included electronics.
2. Establish the actual detector device/lot and measured breakdown voltage;
   resolve bank bias 38/42.5 V, OV 10 V, and the 12 V PDE operating point.
3. Measure the single-photon arrival distribution and pulse shape through the
   **full detector FastIC+ chain** at that operating point. Archive event-level
   input, pulse fits, uncertainties, reference-laser contribution, and separation
   of sensor SPTR from electronics. Validate 85/92 ps derivation rather than
   treating its range as a formal measurement.
4. Establish total electronics time jitter and distinguish it from amplitude
   noise, TDC quantization and common-mode contributions. Resolve W:319–320's
   dimensional mismatch before any electronics study. Measure threshold/fraction
   behavior for END leading edge and dCFD; neither 4 PE nor 14% is calibrated here.
5. Derive/validate propagation through the actual TOP order statistic. Do not
   import independent-mean sqrt(kN) scaling or double-count detector-level SPTR
   that already includes a front end.

**Recommendation, pending René's explicit decision:** use the apparatus-specific
`fastic_measured` registry entry as the future canonical *record structure*,
after full-chain calibration and resolution of units/device/OV conflicts. Its
present 2/3 ns and 85/92 ps entries remain provisional; this recommendation does
not activate them, replace the END model, choose a numeric SPTR, or unify any
live implementation. Until that decision and evidence exist, retain the current
separate intrinsic TOP and END reporting.
