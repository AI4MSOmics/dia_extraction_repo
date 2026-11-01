# Troubleshooting

## Extraction (71)

| Issue | Cause | Fix |
|-------|-------|-----|
| `RawFileReaderAdapter.FileFactory` throws | SelectInstrument race | Retry 3× with 0.5s sleep (built in) |
| `GetCentroidStream` "No method matches" | numpy.int64 scan number | Cast to `int(scan)` before call |
| `RunHeader2.LastScan` not found | API name changed | Use `RunHeaderEx.LastSpectrum` |
| All rows status=no_scan | Window regex mismatch | Check filename CV suffix pattern |
| abs_drt_med > 1 min | Wrong RT column or unit | Verify search table RT is in minutes |

## Profiles (73)

| Issue | Cause | Fix |
|-------|-------|-----|
| Heatmap shows color across entire RT | Dense low-m/z window background | Three-tier removal already built in; check ion floor is 0.10 |
| Correlation profile flat for all peptides | apex scan not found (S14 RT out of range) | Check `rt_s14` column values are per-sample (not global) |
| profiles.npz empty | `TGT`/`ROWSX` env vars not set for extra chain | Set both before running |

## R rendering (78)

| Issue | Cause | Fix |
|-------|-------|-----|
| `setwd` fails | Chinese characters in path | Copy risoPSV/ to ASCII path |
| `invalid input found on input connection` | Git Bash LC_CTYPE=C.UTF-8 | `env -u LC_CTYPE -u LC_ALL -u LANG Rscript` |
| `grey60` color not found | R color name not in matplotlib | R uses it natively; matplotlib version uses `#999999` |
| Panel title clipped at top | oma/mar line budget exceeded | Increase `oma[3]` or decrease `mtext line` |
| `pdf family="Arial"` warning | Not a standard R pdf family | Use `family="Helvetica"` |

## Composition (79)

| Issue | Cause | Fix |
|-------|-------|-----|
| `FzErrorSystem: cannot remove file` | PDF viewer holds lock | Close viewer before re-running |
| Heatmap same size after rect change | insert_image preserves aspect | Match `figsize` aspect to rect aspect exactly |
| Text overflows page edge | Info line too long | `fit_text` auto-shrinks; check max_w parameter |
| Heatmap blank (all white) | All rows below quant threshold | Check that profile data exists; check row_quant_threshold |
| Blank rows missing quantity annotation | `gk.qty` is NaN | Check search table has quantity column |

## Table 1 (81)

| Issue | Cause | Fix |
|-------|-------|-----|
| Many `no_spectrum` placeholders | Peptides not in S14 or extraction | Run extra chain (82→71→73→84) first |
| Page order wrong | Table 1 xlsx structure changed | Verify `parse_table1()` column indices |
| Extra chain pages blank | `profiles_extra` or `topblock_extra` missing | Run 73 and 78 for extra chain |

## Predictions (86)

| Issue | Cause | Fix |
|-------|-------|-----|
| `peptdeep` rc≠0 | Missing model or config | Check `apd_settings_missing.yaml` paths |
| pkl key mismatch | Synthetic pid vs Spectronat-style pid | Add alias: `pred[(ds, our_pid)] = pred[(ds, apd_pid)]` |
| SFTFTLK z=1 has no z=1 entry | AlphaPeptDeep input built with page charge | Use z=2 intensity (same sequence+mods = same fragment intensities) |

## Audit (85)

| Issue | Cause | Fix |
|-------|-------|-----|
| `get_text` can't find heatmap title | Title is inside PNG (rasterized) | Normal — use recalculated values for spot-check |
| `extract_image` md5 mismatch | PDF re-encodes images | Use structural check (bbox + pixel width) instead |
| A9 false positive | Audit formula out of sync with rendering | Update audit A9 to match 79's threshold chain |
