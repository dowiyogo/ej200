# EXEC_07 EndTop geometry and wrapping audit

Date: 2026-06-10

## Channel map

The authoritative counts are in `include/DetectorConstruction.hh:37-39`:

- End left: IDs 0-7.
- End right: IDs 8-15.
- Top left group: IDs 16-50, x=-692 to -12 mm.
- Top right group: IDs 51-85, x=+12 to +692 mm.
- Total: 86 channels.

The fixed Top position function is `src/DetectorConstruction.cc:45-48`. It
keeps 20 mm pitch within each group and the deliberate 24 mm central pitch.
Top placement and individual Bar-to-coupling border surfaces are in
`src/DetectorConstruction.cc:347-364`.

The EXEC_06 End placement loop remains geometrically unchanged at
`src/DetectorConstruction.cc:322-345`.

## Wrapping implementation

This branch uses explicit physical reflector-panel volumes carrying
`BarSkinReflector` optical border surfaces, rather than a bar skin surface.
Therefore the physical-volume route was used for the Top windows.

The `+Y` panel starts as a complete face and receives 70 exact 6 x 6 mm2
subtractions at the Top centers in `src/DetectorConstruction.cc:289-300`. The
resulting panel is placed and connected to the bar by an optical border surface
in `src/DetectorConstruction.cc:301-308`.

Each Top coupling volume is a daughter inside `BarLV` and begins before the
outer bar boundary. Its individual `BarPV -> TopSiPMPV` border surface is
therefore encountered at the active coupling interface and is not blocked by
the outer `+Y` reflector panel.

For EndTop, the two End faces have no reflector panels. All other faces are
wrapped, except the 70 exact Top windows.

## Guardrails

- `tests/readout_config_check.cc:71-97` checks copy numbers 0-85, uniqueness,
  exact Top x/z positions, and the 24 mm central gap.
- `tests/readout_config_check.cc:60-68` checks that Top/EndTop use a
  `G4SubtractionSolid` `+Y` reflector.
- `tests/export_endtop_gdml.cc` exports the EndTop geometry.
- `tests/check_endtop_gdml.py` parses that GDML and checks all 70 Top positions
  within 0.01 mm.
- `tests/check_endtop_balance.py` rejects open-face-level leakage in a 50-event
  center smoke.

## Tedlar note

The hardware's exterior black Tedlar layer is deliberately not modeled. Its
purpose is light tightness, not internal reflection. The non-reflected fraction
`1-R` is already removed by the Mylar-like optical surface, and the Tedlar
material budget is negligible for the muon.
