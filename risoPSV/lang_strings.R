# lang_strings.R — 由 lang.py 生成的当前语言标签（78 号 R 脚本 source 本文件）
# 在运行 78 之前由 Python 端写入：`python lang.py --r > risoPSV/lang_strings.R`
# 或直接 source（默认 en）。
DIA_LANG <- Sys.getenv("DIA_LANG", "en")
if (DIA_LANG == "" || !startsWith(tolower(DIA_LANG), "zh")) {
  L <- function(s) s   # en: 透传
} else {
  .zh <- new.env(parent = emptyenv())
  .zh[["Intensity (%)"]]            <- "强度 (%)"
  .zh[["Observed spectrum"]]        <- "实测谱图"
  .zh[["AlphaPeptDeep prediction (inverted)"]] <- "AlphaPeptDeep 理论谱（倒置）"
  .zh[["AlphaPeptDeep prediction unavailable for this precursor"]] <-
    "该前体无 AlphaPeptDeep 预测"
  .zh[["Mass error of matched ions (observed spectrum)"]] <-
    "匹配离子质量差（实测谱）"
  L <- function(s) {
    v <- .zh[[s]]
    if (is.null(v)) s else v
  }
}
