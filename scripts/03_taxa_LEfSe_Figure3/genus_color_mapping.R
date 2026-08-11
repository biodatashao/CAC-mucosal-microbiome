############################################################
## genus_color_mapping.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
##
## Shared genus-to-colour mapping, sourced by the Figure 3A and
## Figure 3B composition panels so that both use identical colours.
############################################################

GENUS_COLORS <- c(
  "UCG_005"                         = "#A92535",
  "Enterococcus"                    = "#CC6C91",
  "Blautia"                         = "#827976",
  "Romboutsia"                      = "#4C9996",
  "Lactobacillus"                   = "#91BDDB",
  "Ligilactobacillus"               = "#83B8B4",
  "Christensenellaceae_R_7_group"   = "#ADA5A1",
  "Lactococcus"                     = "#9A735D",
  "Collinsella"                     = "#F5B970",
  "Lachnospiraceae_NK4A136_group"   = "#F58B24",
  "Dubosiella"                      = "#C89AC1",
  "Clostridium"                     = "#AD7DA7",
  "Faecalibacterium"                = "#FB9695",
  "Escherichia_Shigella"            = "#E85A5E",
  "Ruminococcus"                    = "#F2CB57",
  "Rikenellaceae_RC9_gut_group"     = "#B99B27",
  "Bacteroides"                     = "#8DCC76",
  "Mediterraneibacter"              = "#58A64D",
  "Extibacter"                      = "#557EAA",
  "Streptococcus"                   = "#6F63A6",
  
  ## Genera appearing only in the paired CA/nonCA panel
  "Weissella"                       = "#476F9F",
  "Ammoniphilus"                    = "#72BD68",
  "UCG_009"                         = "#C2A42C"
)

missing_genus_color <- function(genera, color_map = GENUS_COLORS) {
  setdiff(as.character(genera), names(color_map))
}