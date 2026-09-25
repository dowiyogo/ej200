# EXEC_07 SiPM model audit

Date: 2026-06-10

## Active model and data

- Active/default model: `AFBR-S4N66P024M`.
- Selector: `/sipm/model AFBR-S4N66P024M`
  (`src/DetectorConstruction.cc:92-96`).
- Data file: `data/sipm/AFBR-S4N66P024M_pde.txt`.
- Data-file resolver and parser: `src/SiPMModel.cc:36-79`.
- Surface loading: `src/Materials.cc:222-248`.
- SD metadata lookup from the same curve: `src/SiPMSD.cc:117-126`.

The curve is digitized from Broadcom AFBR-S4N66P024M-DS105, Figure 6, at the
datasheet's typical 12 V overvoltage. It contains 36 points from 250 to 900 nm
and includes the datasheet cross-check value 63% at 420 nm.

The curve is installed on every End and Top SiPM optical surface as
`EFFICIENCY(lambda)`. `REFLECTIVITY(lambda)` is zero. Geant4's optical boundary
process applies PDE once and invokes the sensitive detector; the SD does not
perform a second Bernoulli trial.

`tests/readout_config_check.cc:41-57` checks the surface MPT, 63% efficiency at
420 nm, and zero reflectivity.

## Excluded correlated noise

The PDE file does not include optical crosstalk or afterpulsing. Broadcom
DS105 reports approximately 23% total correlated noise/crosstalk and less than
1% afterpulsing under its stated conditions; neither is modeled in EXEC_07.

## Passivation/coupling approximation

The product passivation is a clear epoxy mold compound, not silicone resin.
DS105 does not specify a refractive index for that epoxy. The existing SiO2
coupling approximation at n=1.58 is therefore retained and explicitly marked
for future refinement.

The physical product is a 2 x 1 array with 7 mm element pitch. In this
simulation, each placed SiPM represents one 6 x 6 mm2 element; the physical
package pitch is not applied to the simulated element positions.
