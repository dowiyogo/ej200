# SYNC_STATUS.md — Sincronización de clones ej200

Fecha: 2026-09-01  
Origin: `git@github.com:dowiyogo/ej200.git`  
Clones: HOST (`/home/rrios/ej200`), W1 (`/home/reriosto/SHiP/ej200`), W2 (`/home/reriosto/SHiP/ej230`)

---

## S1 — Estado previo a la sincronización

| Clone | Archivos ` M` (tracked) | Stash | Worktrees activos |
|-------|------------------------|-------|------------------|
| HOST  | ninguno ✓              | `stash@{0}`: DATA_AUDIT +57 líneas (sesión actual) ⚠ | 1 (feat/bar-end-vikuiti) |
| W1    | ninguno ✓              | ninguno ✓ | 4 (ver §4) |
| W2    | ninguno ✓              | ninguno ✓ | 2 (ver §4) |

**Stash HOST**: el stash contiene el trabajo en curso de `docs/branch-diagnosis-2026-08-31`
(DATA_AUDIT +57 líneas guardado en esta sesión para cambio de rama seguro).
No hay ` M` en el working tree. El fetch y el ff-merge no tocan el stash.
Se procedió con la sincronización; el stash queda intacto.

---

## S2 — Resultado del fetch

```
HOST:  ya al día (sin output)
W1:    origin/feat/bar-end-vikuiti   b007750 → 1bfb827 (actualizado)
       origin/diag/phase7-delta-2026-08-31  (nueva)  → a18e863
       origin/docs/branch-diagnosis-2026-08-31 (nueva) → e67d57f
W2:    origin/feat/bar-end-vikuiti   b007750 → 1bfb827 (actualizado)
       origin/diag/phase7-delta-2026-08-31  (nueva)  → a18e863
       origin/docs/branch-diagnosis-2026-08-31 (nueva) → e67d57f
```

Los tres clones ven los mismos 18 refs remotas con los mismos hashes.

---

## S3 — Fast-forward de ramas locales

### Ramas fast-forwardeadas

| Clone | Rama | Hash anterior | Hash nuevo | Método |
|-------|------|--------------|-----------|--------|
| HOST  | `diag/phase7-delta-2026-08-31` | `6a29b09` | `a18e863` | `fetch origin B:B` |
| HOST  | `feat/endtop-sslg4`            | `4cb2c23` | `5576687` | `fetch origin B:B` |
| W1    | `feat/ej230-bar-tir-only`      | `a95f610` | `b281aea` | `fetch origin B:B` |

### Ramas bloqueadas

| Clone | Rama | Estado | Causa |
|-------|------|--------|-------|
| W1 | `feat/bar-end-vikuiti` | **BLOCKED** (ff no realizado) | `merge --ff-only` falló: archivos no rastreados en el worktree serían sobreescritos por el commit `50ae425` que añadió `analysis/presentation_v4/`, `talks/napkin_first_principles/` al índice. El local está en `b007750`; origin en `1bfb827`. |

**Acción requerida (W1):** Desde el worktree `/mnt/d/SHiP/ej200` (o equivalente en Windows/WSL2), mover o eliminar los archivos en conflicto listados abajo y ejecutar:
```bash
git merge --ff-only origin/feat/bar-end-vikuiti
```
Archivos en conflicto:
- `analysis/presentation_v4/scripts/analysis_v4.py`
- `analysis/presentation_v4/scripts/fig_materials.C`
- `analysis/presentation_v4/scripts/fig_materials.py`
- `analysis/presentation_v4/tables/summary_numbers.json`
- `analysis/presentation_v4/talk_v4.tex`
- `presentations/optim_2026-08-16/talk.tex`
- `talks/napkin_first_principles/fig_gen.C`
- `talks/napkin_first_principles/fig_gen.py`
- `talks/napkin_first_principles/napkin.py`
- `talks/napkin_first_principles/tex/talk.tex`
- `talks/napkin_first_principles/tex/talk_v3.tex`
- `talks/napkin_first_principles/values/convention.json`
- `talks/napkin_first_principles/values/materials.yaml`
- `talks/napkin_first_principles/values/napkin_macros.tex`

### Ramas ya en sincronía (sin acción)

Estas ramas locales coincidían con origin antes o después del fetch:

HOST: `exp/pair-scan-2026-06-11`, `feat/bar-end-vikuiti`, `feat/bar-vikuiti`,
`feat/ej204-bar-tir-only`, `feat/ej230-bar-tir-only`, `main`,
`wip/host-stash-endtop-junio`, `wip/host-uncommitted-2026-08-31`

W1: `exp/pair-scan-2026-06-11`, `feat/bar-vikuiti`, `feat/ej204-bar-tir-only`,
`feat/ej204-event-display-tracks`, `feat/ej228-cylinder`, `feat/ej228-tir-only`,
`feat/ej230-sslg4`, `feat/endonly-mylar`, `feat/endtop-sslg4`, `main`

W2: `feat/ej230-sslg4`, `main`

### Ramas locales ADELANTADAS al origin (push pendiente — no se tocan)

| Clone | Rama | Local | Origin | Estado |
|-------|------|-------|--------|--------|
| HOST  | `docs/branch-diagnosis-2026-08-31` | `2a9645f` | `e67d57f` | Local ahead — necesita push |
| W2    | `feat/ej230-endonly-mylar`          | `04a8047` | `98af271` | Local ahead — necesita push |

---

## S4 — Verificación de consistencia de refs remotas

`ls-remote --heads origin | wc -l` en los tres: **18** ✓

### Tabla clone × rama remota (18 ramas)

Todas las celdas son iguales en los tres clones — sin discrepancias.

| Rama remota | HOST | W1 | W2 |
|-------------|------|----|----|
| `origin/diag/phase7-delta-2026-08-31` | `a18e863` | `a18e863` | `a18e863` |
| `origin/docs/branch-diagnosis-2026-08-31` | `e67d57f` | `e67d57f` | `e67d57f` |
| `origin/exp/pair-scan-2026-06-11` | `f19c093` | `f19c093` | `f19c093` |
| `origin/feat/bar-end-vikuiti` | `1bfb827` | `1bfb827` | `1bfb827` |
| `origin/feat/bar-vikuiti` | `219fbe3` | `219fbe3` | `219fbe3` |
| `origin/feat/ej204-bar-tir-only` | `09f8b18` | `09f8b18` | `09f8b18` |
| `origin/feat/ej204-event-display-tracks` | `47a9a4f` | `47a9a4f` | `47a9a4f` |
| `origin/feat/ej228-cylinder` | `66b674c` | `66b674c` | `66b674c` |
| `origin/feat/ej228-tir-only` | `0006919` | `0006919` | `0006919` |
| `origin/feat/ej230-bar-tir-only` | `b281aea` | `b281aea` | `b281aea` |
| `origin/feat/ej230-endonly-mylar` | `98af271` | `98af271` | `98af271` |
| `origin/feat/ej230-sslg4` | `5b93f4c` | `5b93f4c` | `5b93f4c` |
| `origin/feat/endonly-mylar` | `fb3749d` | `fb3749d` | `fb3749d` |
| `origin/feat/endtop-sslg4` | `5576687` | `5576687` | `5576687` |
| `origin/feature/sipm-electronics-response` | `bb4cf6a` | `bb4cf6a` | `bb4cf6a` |
| `origin/main` | `84e902c` | `84e902c` | `84e902c` |
| `origin/wip/host-stash-endtop-junio` | `9710f34` | `9710f34` | `9710f34` |
| `origin/wip/host-uncommitted-2026-08-31` | `d2c8d4c` | `d2c8d4c` | `d2c8d4c` |

**Conclusión S4**: Las 18 refs remotas son idénticas en los tres clones. ✓

---

## S4 — Ramas locales por clone y worktrees

Las ramas locales difieren entre clones como es esperado (cada clone tiene su propio
conjunto de ramas de trabajo y worktrees distintos).

### HOST — ramas locales

| Rama | Hash | Worktree |
|------|------|----------|
| `diag/phase7-delta-2026-08-31` | `a18e863` | — |
| `docs/branch-diagnosis-2026-08-31` | `2a9645f` | — |
| `exp/pair-scan-2026-06-11` | `f19c093` | — |
| `feat/bar-end-vikuiti` | `1bfb827` | `/home/rrios/ej200` ← rama actual |
| `feat/bar-vikuiti` | `219fbe3` | — |
| `feat/ej204-bar-tir-only` | `09f8b18` | — |
| `feat/ej230-bar-tir-only` | `b281aea` | — |
| `feat/endtop-sslg4` | `5576687` | — |
| `main` | `84e902c` | — |
| `wip/host-stash-endtop-junio` | `9710f34` | — |
| `wip/host-uncommitted-2026-08-31` | `d2c8d4c` | — |

### W1 (ej200) — ramas locales

| Rama | Hash local | Hash origin | Estado | Worktree |
|------|-----------|------------|--------|----------|
| `exp/pair-scan-2026-06-11` | `f19c093` | `f19c093` | ✓ | — |
| `feat/bar-end-vikuiti` | `b007750` | `1bfb827` | **BLOCKED** | `/mnt/d/SHiP/ej200` ← rama actual |
| `feat/bar-vikuiti` | `219fbe3` | `219fbe3` | ✓ | — |
| `feat/ej204-bar-tir-only` | `09f8b18` | `09f8b18` | ✓ | — |
| `feat/ej204-event-display-tracks` | `47a9a4f` | `47a9a4f` | ✓ | `/mnt/d/SHiP/ej200_event_display` |
| `feat/ej228-cylinder` | `66b674c` | `66b674c` | ✓ | — |
| `feat/ej228-tir-only` | `0006919` | `0006919` | ✓ | — |
| `feat/ej230-bar-tir-only` | `b281aea` | `b281aea` | ✓ (ff'd) | — |
| `feat/ej230-sslg4` | `5b93f4c` | `5b93f4c` | ✓ | — |
| `feat/endonly-mylar` | `fb3749d` | `fb3749d` | ✓ | `/mnt/d/SHiP/ej200_endonly` |
| `feat/endtop-sslg4` | `5576687` | `5576687` | ✓ | — |
| `main` | `84e902c` | `84e902c` | ✓ | `/mnt/d/SHiP/ej200_edge_scan` |

**Nota W1**: Los worktrees adicionales viven en rutas Windows (`/mnt/d/SHiP/`) no accesibles desde
el entorno SSH (`bash --noprofile --norc`). Sus estados no pudieron verificarse directamente.

### W2 (ej230) — ramas locales

| Rama | Hash local | Hash origin | Estado | Worktree |
|------|-----------|------------|--------|----------|
| `feat/ej230-endonly-mylar` | `04a8047` | `98af271` | local ahead | `/mnt/d/SHiP/ej230_endonly_mylar` |
| `feat/ej230-sslg4` | `5b93f4c` | `5b93f4c` | ✓ | `/mnt/d/SHiP/ej230` ← rama actual |
| `main` | `84e902c` | `84e902c` | ✓ | — |

---

## S4 — Tags

Los tags difieren entre clones. Origin no publica tags (confirmado: los tres clones
acumularon tags locales distintos en diferentes sesiones de trabajo).

| Clone | Tags presentes |
|-------|---------------|
| HOST  | 17 (checkpoint/*, diag-*, exec07-09, physics-baseline-v1, v-exec21..v27) |
| W1    | 30 (incluyendo EXEC_22b-pre-EXEC_28-pre, pre-exec13-*, pre-exec16-*) |
| W2    | 22 (incluyendo checkpoint/pre-exec14b..e, pre-exec13-230-checkpoint) |

Diferencia esperada — los tags son locales y no forman parte del sync remoto.

---

## S5 — Estado de los working trees después de la sincronización

| Clone | Archivos ` M` | Estado |
|-------|--------------|--------|
| HOST  | ninguno | ✓ limpio |
| W1    | ninguno | ✓ limpio |
| W2    | ninguno | ✓ limpio |

Ningún working tree quedó modificado por las operaciones de sincronización.

---

## Pendientes para el usuario

| # | Clone | Acción | Detalle |
|---|-------|--------|---------|
| 1 | **W1** | Resolver conflicto de archivos no rastreados | Mover/eliminar `analysis/presentation_v4/` y `talks/napkin_first_principles/` del worktree `/mnt/d/SHiP/ej200` y ejecutar `git merge --ff-only origin/feat/bar-end-vikuiti` |
| 2 | **HOST** | Push `docs/branch-diagnosis-2026-08-31` | Local (`2a9645f`) adelantado a origin (`e67d57f`) por el commit DATA_AUDIT |
| 3 | **W2** | Push `feat/ej230-endonly-mylar` | Local (`04a8047`) adelantado a origin (`98af271`) |
| 4 | **HOST** | `stash@{0}` (DATA_AUDIT +57 líneas) | Requiere commit en `docs/branch-diagnosis-2026-08-31` (pendiente de aprobación) |

---

*Sincronización read-only + fast-forward. Sin merges, sin pushes, sin resets.*
*Ningún worktree quedó modificado.*
