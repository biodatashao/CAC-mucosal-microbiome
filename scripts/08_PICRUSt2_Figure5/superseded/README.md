# Superseded scripts

These three scripts were part of an earlier attempt at the pathway-level LEfSe
analysis. **None of their outputs is read by any script in the final pipeline**,
and none contributed to Figure 5 or Supplementary Figure 2.

| File | Why it is superseded |
|---|---|
| `superseded_prepare_LEfSe_input_progression5_altnaming.R` | Alternative LEfSe input preparation writing files without the `oldstyle_currentdata` suffix. The `.res` file it leads to is read only by `superseded_parse_LEfSe_progression5_v2.R`. |
| `superseded_parse_LEfSe_progression5_v1.R` | Standalone parser writing `FINAL_LEfSe_clean_markers_MetaCyc_pathways_oldstyle_currentdata_*.csv`. The downstream script `08_04` re-parses the `.res` file itself and reads a differently named file, so this output is never consumed. |
| `superseded_parse_LEfSe_progression5_v2.R` | Parser for the alternative naming branch; writes `FINAL_v2_*` files that no downstream script reads. |

They are retained for transparency about the analysis history. Delete this
folder if you prefer the repository to contain only the final pipeline.
