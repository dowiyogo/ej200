# Preparación y Verificación de Hard Links para las 9 Celdas `.root` Físicamente Duplicadas (`2026-09-25`)

- **Fecha de verificación:** `2026-09-25T21:35:00+02:00`
- **Fuente de referencia:** [`docs/reports/AUDITORIA_REPO_20260925.md`](file:///home/rrios/ej200/docs/reports/AUDITORIA_REPO_20260925.md#L530-L541) (§6.1, filas 3 y 4)
- **Estado de ejecución:** **SOLO PREPARACIÓN (0 archivos modificados, borrados o enlazados)**.

---

## 1. Hallazgo Crítico de la Verificación `sha256sum` en Vivo (Regla 2 — Casos Dudosos)

Al ejecutar hoy `sha256sum` y `stat` / `stat -f` sobre los **18 archivos (9 pares)** se confirmó lo siguiente:

1. **Mismo filesystem confirmado (`9/9` pares):** Los 18 archivos residen en el mismo dispositivo y sistema de archivos (`st_dev = 64770`, `fsid = fd0200000000`, `xfs` montado en `/home` sobre `/dev/mapper/almalinux_t0minidaq-home`). Todos tienen `nlink = 1`.
2. **Integridad individual intacta (`18/18` archivos):** El `sha256sum` calculado hoy para cada uno de los 18 archivos coincide al **100%** con el `root_sha256` registrado el día de su simulación en su respectivo archivo `.DONE` / `.DONE.json`. Ningún archivo ha cambiado ni se ha corrompido desde el reporte original.
3. **Diferencia byte-a-byte entre canónico y duplicado (`0/9` pares con mismo `sha256sum`):**
   - En [`AUDITORIA_REPO_20260925.md` (§6.1, líneas 532–540)](file:///home/rrios/ej200/docs/reports/AUDITORIA_REPO_20260925.md#L532-L540), estos 9 pares fueron clasificados como **duplicados físicos de simulación Geant4 (`Tipo: Física (100% mismos hits)`)**, **no** como copias exactas byte-a-byte (`Exacta SHA-256`).
   - Cada par proviene de dos corridas independientes de Geant4 con **4 hilos** (`/run/numberOfThreads 4`, `/run/eventModulo 1`) ejecutadas con el mismo binario (`4967ec8`), el mismo `run.mac`, las mismas tablas ópticas `sslg4`, `N = 10,000` eventos y las mismas semillas (`/random/setSeeds 26092601 8349041`). Por ello, ambos archivos contienen **exactamente los mismos fotones y conteos por cara (`root_entries`, `left_total`, `right_total`, `top_total`) hasta el último fotón**, pero difieren en el `TUUID`/timestamp de la cabecera `TFile` y en el orden en que los 4 hilos volcaron sus *baskets* comprimidos a disco (diferencia de tamaño de `2.1 KB` a `34.0 KB` y distinto `sha256sum`).
4. **Efecto colateral sobre compuertas de SHA-256 si se enlazan sin actualizar `.DONE`:**
   - Tanto [`analysis/sigma_t/orchestration/detached_grid.py`](file:///home/rrios/ej200/analysis/sigma_t/orchestration/detached_grid.py#L224-L235) (`done_record()`) como [`analysis/track_mechanism_20260915/analyze_step1.py`](file:///home/rrios/ej200/analysis/track_mechanism_20260915/analyze_step1.py#L400-L421) (`hash_matches_done = digest == cell["done"]["root_sha256"]`) verifican que el tamaño y el `sha256sum` de `photon_hits_run000.root` coincidan con los registrados en `.DONE`.
   - Si se reemplaza el archivo duplicado por un hard link al canónico sin actualizar simultáneamente `root_sha256` y `root_size_bytes` en el `.DONE` / `.DONE.json` del duplicado, la compuerta `hash_matches_done` fallará al re-ejecutar `analyze_step1.py` o `detached_grid.py status`.

---

## 2. Tabla de Verificación de los 9 Pares (`sha256sum` + `stat -f` ejecutados hoy)

| Par # | Archivo canónico (ruta completa) | Archivo duplicado (ruta completa) | Tamaño (`canónico` / `duplicado`) | `sha256` (primeros 16 caracteres: `canónico` / `duplicado`) | Verificado hoy (`sí/no`) |
| :-: | :--- | :--- | :--- | :--- | :--- |
| **1** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xm650/attempts/df6fd49e2d9440c3976f61da0e8cb37e/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm650/attempts/6e0b254358c444daa17306a784154b74/photon_hits_run000.root` | `11,697,290,259 B` (`10.89 GiB`) /<br>`11,697,301,088 B` (`10.89 GiB`) | `1c0a37dbca3de1b7` /<br>`e0de53732d3f2837` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `54,935,742`) |
| **2** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xm500/attempts/0f7971bbc02045b4b8da1a7e22197481/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm500/attempts/73c9aea974224a7bb9fd72f557d8cf28/photon_hits_run000.root` | `10,969,343,243 B` (`10.22 GiB`) /<br>`10,969,320,333 B` (`10.22 GiB`) | `82028fd7e7d3f118` /<br>`81888c7c81442f61` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `49,643,803`) |
| **3** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xm200/attempts/3c5abb8811eb46edb86c9718c1f927e0/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm200/attempts/1e2483b4fd644358bff06d8bbea5c710/photon_hits_run000.root` | `10,660,076,544 B` (`9.93 GiB`) /<br>`10,660,087,595 B` (`9.93 GiB`) | `256c94446962cba3` /<br>`e039297e5381c128` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `47,447,287`) |
| **4** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp0/attempts/8e8a9b5e95f149ed8ad3f8ae6d04dbdb/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp0/attempts/13b525e236174c7698e9d9490eba6b96/photon_hits_run000.root` | `10,684,705,272 B` (`9.95 GiB`) /<br>`10,684,688,107 B` (`9.95 GiB`) | `90695d61bb72f6fc` /<br>`ca77bd6f03854021` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `46,946,914`) |
| **5** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp200/attempts/7bb21318515847d78241e239ef9d37ec/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp200/attempts/9b53b9bdda0e4c2d863e9a965c0a792b/photon_hits_run000.root` | `10,655,057,420 B` (`9.92 GiB`) /<br>`10,655,062,962 B` (`9.92 GiB`) | `6f8188fd40d6142a` /<br>`bddf5b0eefc7f921` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `47,421,723`) |
| **6** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp500/attempts/2dbfc26293a14b269acf70b7bf85a859/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp500/attempts/c9a2ca3b1a2f4d8eb77aab2a44123c35/photon_hits_run000.root` | `10,990,320,129 B` (`10.24 GiB`) /<br>`10,990,317,669 B` (`10.24 GiB`) | `64394dbd75938fc5` /<br>`425b21bf9491f45c` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `49,733,221`) |
| **7** | `/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp650/attempts/9eaf2cbaff094f568538a535a44a11f2/photon_hits_run000.root` | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp650/attempts/8e2c787777d84fba9fae2023e92d5538/photon_hits_run000.root` | `11,717,670,095 B` (`10.91 GiB`) /<br>`11,717,704,090 B` (`10.91 GiB`) | `4292559dd9fdfb4d` /<br>`4409de815b9c737d` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `55,000,640`) |
| **8** | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ200_xm650/attempts/397938f5676a479ca27baf2729058772/photon_hits_run000.root` | `/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/photon_hits_run000.root` | `13,721,650,643 B` (`12.78 GiB`) /<br>`13,721,646,767 B` (`12.78 GiB`) | `4acec57272a60327` /<br>`d81cbe71b133d584` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `68,068,729`) |
| **9** | `/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ204_xm650/attempts/0c742a195a4d4eaf98d5ab87b7493f7b/photon_hits_run000.root` | `/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/photon_hits_run000.root` | `13,314,248,448 B` (`12.40 GiB`) /<br>`13,314,246,330 B` (`12.40 GiB`) | `05504aecb07c5255` /<br>`675210e81616bdf9` | **Sí** (`xfs` mismo FS; **NO idéntico byte-a-byte**, `100%` idéntico en hits: `64,963,360`) |

---

## 3. Bloque de Comandos Bash

### 3.1 Lista Estricta Bajo la Regla 2 (`sha256sum` byte-a-byte idéntico)

Aplicando estrictamente la Regla 2 del protocolo (*"Si algún par YA NO es idéntico [por `sha256sum` byte a byte], exclúyelo de la lista final y repórtalo como caso dudoso — no generes el comando para ese par"*), **los 9 pares quedan excluidos de ejecución automática porque son duplicados físicos multihilo (distinto `sha256sum` de contenedor ROOT), no copias byte-a-byte**.

### 3.2 Referencia de Rutas Reales para los 9 Pares de Duplicados Físicos (Solo si el Usuario Decide Enlazar Duplicados Físicos)

> [!WARNING]
> Ejecutar los siguientes comandos libera **`97.24 GiB` (`104.41 GB`)** reemplazando cada segunda corrida Geant4 por un hard link a la primera corrida físicamente equivalente. Sin embargo, dado que sus `sha256sum` de contenedor ROOT difieren, tras ejecutarlos sería necesario actualizar `root_sha256` y `root_size_bytes` en los archivos `.DONE` / `state.json` / `.DONE.json` de las rutas duplicadas si se desea volver a correr `analyze_step1.py` o `detached_grid.py status`.

```bash
rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm650/attempts/6e0b254358c444daa17306a784154b74/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xm650/attempts/df6fd49e2d9440c3976f61da0e8cb37e/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm650/attempts/6e0b254358c444daa17306a784154b74/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm500/attempts/73c9aea974224a7bb9fd72f557d8cf28/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xm500/attempts/0f7971bbc02045b4b8da1a7e22197481/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm500/attempts/73c9aea974224a7bb9fd72f557d8cf28/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm200/attempts/1e2483b4fd644358bff06d8bbea5c710/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xm200/attempts/3c5abb8811eb46edb86c9718c1f927e0/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xm200/attempts/1e2483b4fd644358bff06d8bbea5c710/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp0/attempts/13b525e236174c7698e9d9490eba6b96/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp0/attempts/8e8a9b5e95f149ed8ad3f8ae6d04dbdb/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp0/attempts/13b525e236174c7698e9d9490eba6b96/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp200/attempts/9b53b9bdda0e4c2d863e9a965c0a792b/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp200/attempts/7bb21318515847d78241e239ef9d37ec/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp200/attempts/9b53b9bdda0e4c2d863e9a965c0a792b/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp500/attempts/c9a2ca3b1a2f4d8eb77aab2a44123c35/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp500/attempts/2dbfc26293a14b269acf70b7bf85a859/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp500/attempts/c9a2ca3b1a2f4d8eb77aab2a44123c35/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp650/attempts/8e2c787777d84fba9fae2023e92d5538/photon_hits_run000.root'
ln '/home/rrios/exec46_20260915/full_grid/cells/EJ230_xp650/attempts/9eaf2cbaff094f568538a535a44a11f2/photon_hits_run000.root' '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ230_xp650/attempts/8e2c787777d84fba9fae2023e92d5538/photon_hits_run000.root'

rm '/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/photon_hits_run000.root'
ln '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ200_xm650/attempts/397938f5676a479ca27baf2729058772/photon_hits_run000.root' '/home/rrios/exec46_20260915/f4_bc408_sensitivity/visible_current_3800mm/photon_hits_run000.root'

rm '/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/photon_hits_run000.root'
ln '/home/rrios/exec46_20260916/full_grid_bc408_bc404_3800_v2/cells/EJ204_xm650/attempts/0c742a195a4d4eaf98d5ab87b7493f7b/photon_hits_run000.root' '/home/rrios/exec46_20260916/i2_bc404_validation/EJ204_xm650/photon_hits_run000.root'

# Total que se liberará si se ejecutan los 9 pares completos: 97.24 GiB (104,410,374,941 bytes / 104.41 GB)
```
