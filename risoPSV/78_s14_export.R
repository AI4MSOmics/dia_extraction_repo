# 78_s14_export.R —— S14 证据页（isoPSV 引擎原版排版，每肽一页前五区块）
# 区块顺序（用户指定）：信息行(3) -> 离子梯 -> 原始谱图 -> 预测谱图(倒置) -> 质量差分布
# 引擎（all_ions / draw_ladder / 谱面板 / ppm 面板 / find_spec / ION_COL）自
# isoPSV single_spectrum_export.R 原样拷贝；倒置预测面板为其绘图原语的扩展。
# 输入：in/s14_spectra.tsv + in/best.mgf（76 号脚本生成）
# 输出：out/s14_topblock.pdf（每肽一页，8.27 x 7.9 in，供 79 号与热图合成 A4）
# 用法（环境变量）：IN_TSV / OUT_PDF / SEQ（可选，只渲一条）
suppressWarnings(suppressMessages({
  source("compat.R", encoding = "UTF-8")
  source("ReadMGF_script.R", encoding = "UTF-8")
  source("chimeric_check.R", encoding = "UTF-8")
}))
setwd("F:/IPX0018588000/work/risoPSV")
IN_TSV <- Sys.getenv("IN_TSV", "F:/IPX0018588000/work/risoPSV/in/s14_spectra.tsv")
OUT_PDF <- Sys.getenv("OUT_PDF", "F:/IPX0018588000/work/risoPSV/out/s14_topblock.pdf")
SEQ <- Sys.getenv("SEQ", "")
MGF_DIR <- dirname(IN_TSV)

aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
aa_tab <- amino_acid_table

pf <- read.delim(IN_TSV, sep = "\t", check.names = FALSE, na.strings = character(0),
                 quote = "", fill = TRUE, fileEncoding = "UTF-8")
if (nzchar(SEQ)) pf <- pf[pf$Sequence == SEQ, , drop = FALSE]
cat("rows:", nrow(pf), "\n")

# 单文件小 MGF：直接用库的 parse_lines_vectorized（不走并行集群）
MGF_FILE <- file.path(MGF_DIR, "best.mgf")
specs <- parse_lines_vectorized(readLines(MGF_FILE, warn = FALSE))
spec_titles <- vapply(specs, function(o) o@title, character(1))
cat("mgf spectra:", length(specs), "\n")
find_spec <- function(fn) {
  li <- match(fn, spec_titles)
  if (is.na(li)) NULL else specs[[li]]
}

# ---- isoPSV 引擎（原样拷贝自 single_spectrum_export.R）----
ION_COL <- c(a = "#004000", b = "#008000", c = "#CA1D52", x = "#650E29",
             y = "#CA1D52", z = "#943AA4", M = "#808080", I = "#FF8C00", O = "#9ACD32")
ION_PRI <- c(b = 1, y = 2, a = 3, c = 4, x = 5, z = 6, M = 7, I = 8, O = 9)
CO_M <- 27.994915; NH3_M <- 17.026549; H2O_M <- 18.010565; H_M <- 1.007276

all_ions <- function(seq, mod, mz, ints, tol_ppm = 20) {
  for (md in unique(c(mod, ""))) {
    bip <- NULL
    tt <- try(suppressWarnings(invisible(capture.output(
      bip <- by_ions_position(seq, md, "Unlabel", aa, aa_tab)))), silent = TRUE)
    if (inherits(tt, "try-error") || is.null(bip) ||
        length(bip) == 1 || is.null(bip$ions_table)) next
    t <- bip$ions_table; n <- t$Number
    theo <- data.frame(
      mz  = c(t$b_1_charge, t$b_1_charge - CO_M, t$b_1_charge + NH3_M,
              t$b_2_charge, t$b_2_charge - CO_M/2, t$b_2_charge + NH3_M/2,
              t$y_1_charge, t$y_1_charge + CO_M, t$y_1_charge - NH3_M,
              t$y_2_charge, t$y_2_charge + CO_M/2, t$y_2_charge - NH3_M/2),
      lab = c(paste0("b", n, "+"), paste0("a", n, "+"), paste0("c", n, "+"),
              paste0("b", n, "++"), paste0("a", n, "++"), paste0("c", n, "++"),
              paste0("y", n, "+"), paste0("x", n, "+"), paste0("z", n, "+"),
              paste0("y", n, "++"), paste0("x", n, "++"), paste0("z", n, "++")),
      ion = c(rep(c("b", "a", "c"), each = length(n)),
              rep(c("b", "a", "c"), each = length(n)),
              rep(c("y", "x", "z"), each = length(n)),
              rep(c("y", "x", "z"), each = length(n))),
      num = c(rep(n, 3), rep(n, 3), rep(n, 3), rep(n, 3)),
      chg = c(rep(1L, 3 * length(n)), rep(2L, 3 * length(n)),
              rep(1L, 3 * length(n)), rep(2L, 3 * length(n))),
      stringsAsFactors = FALSE)
    if (is.finite(bip$pep_mass)) {
      pm <- bip$pep_mass
      theo <- rbind(theo, data.frame(
        mz = c(pm + H_M, pm + H_M - H2O_M, pm + H_M - NH3_M,
               (pm + 2 * H_M) / 2, (pm + 2 * H_M - H2O_M) / 2, (pm + 2 * H_M - NH3_M) / 2),
        lab = c("M+", "M-H2O+", "M-NH3+", "M++", "M-H2O++", "M-NH3++"),
        ion = "M", num = NA_integer_, chg = c(1L, 1L, 1L, 2L, 2L, 2L),
        stringsAsFactors = FALSE))
    }
    theo <- theo[is.finite(theo$mz), ]
    lo <- theo$mz * (1 - tol_ppm * 1e-6); hi <- theo$mz * (1 + tol_ppm * 1e-6)
    i_lo <- findInterval(lo, mz); i_hi <- findInterval(hi, mz)
    peak_i <- rep(NA_integer_, nrow(theo))
    for (j in seq_len(nrow(theo))) {
      if (i_hi[j] <= i_lo[j]) next
      cand <- (i_lo[j] + 1L):i_hi[j]
      if (length(cand)) peak_i[j] <- cand[which.max(ints[cand])]
    }
    theo$peak <- peak_i
    return(theo)
  }
  NULL
}

draw_ladder <- function(sq, ions, mod = "", chg_v = 0L) {
  n <- nchar(sq)
  if (n < 2L) return(invisible())
  u <- par("usr")
  x0 <- u[1] + 0.015 * (u[2] - u[1]); x1 <- u[2] - 0.015 * (u[2] - u[1])
  u_per_in <- abs(diff(grconvertX(c(0, 1), from = "inches", to = "user")))
  slot_u <- min(0.16 * u_per_in, (x1 - x0) / n)
  xc0 <- (x0 + x1) / 2 - n * slot_u / 2
  xs <- xc0 + (seq_len(n) - 0.5) * slot_u
  gap <- slot_u; ymid <- 1.30; v_h <- 0.045
  site_by <- b_sites <- y_sites <- integer(0)
  if (!is.null(ions)) {
    hit <- ions[!is.na(ions$peak), ]
    if (nrow(hit)) {
      hb <- hit[hit$ion == "b", ]
      if (nrow(hb)) b_sites <- sort(unique(as.integer(hb$num)))
      hy <- hit[hit$ion == "y", ]
      if (nrow(hy)) y_sites <- sort(unique(n - as.integer(hy$num)))
      site_by <- sort(unique(c(b_sites, y_sites)))
    }
  }
  op <- par(xpd = TRUE); on.exit(par(op), add = TRUE)
  slot_in <- slot_u / u_per_in
  cexs <- min(8 / 12, max(6 / 12,
           0.9 * slot_in / strwidth("M", cex = 1, units = "inches")))
  cexn <- min(6 / 12, max(5 / 12,
           0.9 * slot_in / strwidth("b88", cex = 1, units = "inches")))
  if (length(site_by)) {
    vx <- (xs[site_by] + xs[site_by + 1L]) / 2
    segments(vx, ymid - v_h, vx, ymid + v_h, col = "black", lwd = 1.33)
  }
  if (length(b_sites)) {
    vxb <- (xs[b_sites] + xs[b_sites + 1L]) / 2
    segments(vxb, ymid - v_h, vxb - gap * 0.5, ymid - v_h, col = "black", lwd = 1.33)
    text(xs[b_sites], ymid - 0.085,
         as.expression(lapply(b_sites, function(k) bquote(b[.(k)]))),
         cex = cexn, col = ION_COL[["b"]], font = 2)
  }
  if (length(y_sites)) {
    vxy <- (xs[y_sites] + xs[y_sites + 1L]) / 2
    segments(vxy, ymid + v_h, vxy + gap * 0.5, ymid + v_h, col = "black", lwd = 1.33)
    text(xs[y_sites + 1L], ymid + 0.085,
         as.expression(lapply(y_sites, function(k) bquote(y[.(n - k)]))),
         cex = cexn, col = ION_COL[["y"]], font = 2)
  }
  if (is.finite(chg_v) && chg_v > 0L)
    text(xc0 - gap * 1.5, ymid, sprintf("%d+ ", chg_v),
         cex = cexs, font = 2, adj = 1, col = "#CA1D52")
  lab_col <- rep("black", n); lab_font <- rep(1, n)
  if (!is.null(mod) && nzchar(mod) && mod != "-") {
    pm <- suppressWarnings(as.integer(sub(",.*$", "",
                 strsplit(mod, ";", fixed = TRUE)[[1]])))
    pm <- pm[!is.na(pm) & pm >= 1L & pm <= n]
    lab_col[pm] <- ION_COL[["I"]]; lab_font[pm] <- 2
  }
  text(xs, ymid, strsplit(sq, "")[[1]], cex = cexs, font = lab_font, col = lab_col)
}

# 观测谱面板（isoPSV 原版；mar 顶留给梯）
draw_obs_panel <- function(spec_obj, th, chg_v, ann, head) {
  par(mar = c(0.5, 2.2, 3.9, 0.3))
  mz <- as.numeric(spec_obj@spectrum[[1]]); ints <- as.numeric(spec_obj@spectrum[[2]])
  o <- order(mz, na.last = NA); mz <- mz[o]; ints <- ints[o]
  ints_n <- ints / max(ints[is.finite(ints) & ints > 0])
  xmax <- max(mz, na.rm = TRUE) * 1.03
  plot(0, 0, type = "n", xlim = c(0, xmax), ylim = c(0, 1.15),
       xlab = "", ylab = "", yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i")
  title(ylab = "Intensity (%)", cex.lab = 0.75, line = 1.7)
  par(mgp = c(2.0, 0.8, 0))
  axis(2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 100, 20)),
       cex.axis = 0.83, las = 1)
  text(xmax * 0.995, 1.10, ann[[1]], cex = 0.83, adj = c(1, 0.5), xpd = NA)
  text(xmax * 0.995, 1.00, ann[[2]], cex = 0.83, adj = c(1, 0.5), xpd = NA)
  segments(mz, 0, mz, ints_n, lwd = 0.3, col = "grey60")
  if (!is.null(th)) {
    hit <- th[!is.na(th$peak), ]
    if (nrow(hit)) {
      oo <- order(ION_PRI[hit$ion], hit$num, hit$chg, na.last = TRUE)
      hit <- hit[oo, ]
      spl <- split(hit, hit$peak)
      pk_i <- as.integer(names(spl))
      labs <- vapply(spl, function(d) paste(d$lab, collapse = ""), "")
      cols <- vapply(spl, function(d) ION_COL[[d$ion[1]]], "")
      pmz <- mz[pk_i]; pin <- ints_n[pk_i]
      segments(pmz, 0, pmz, pin, lwd = 0.8, col = cols)
      # 标注直接画在各峰高 +0.03 处（不做防撞抬升）
      text(pmz, pin + 0.03, labs, cex = 8/12, srt = 90, col = cols, adj = 0)
    }
  }
  mtext(head, side = 3, line = 5.4, cex = 0.83, adj = 0)
  list(xmax = xmax, mz = mz, ints_n = ints_n)
}

# 倒置预测谱面板（isoPSV 绘图原语扩展：强度向下，标签 srt=90 自棒尖向下）
draw_pred_panel <- function(bmz, ymz, bint, yint, xmax, head, cos_line, geom) {
  par(mar = c(0.5, 2.2, 1.6, 0.3))
  plot(0, 0, type = "n", xlim = c(0, xmax), ylim = c(-1.15, 0.0),
       xlab = "", ylab = "", yaxt = "n", xaxt = "n", xaxs = "i", yaxs = "i")
  title(ylab = "Intensity (%)", cex.lab = 0.75, line = 1.7)
  par(mgp = c(2.0, 0.8, 0))
  axis(2, at = seq(0, -1, -0.2), labels = as.character(seq(0, 100, 20)),
       cex.axis = 0.83, las = 1)
  pmax <- max(max(bint, 0), max(yint, 0), 1e-9)
  pn <- c(bint, yint) / pmax
  pm <- c(bmz, ymz)
  pc <- c(rep(ION_COL[["b"]], length(bmz)), rep(ION_COL[["y"]], length(ymz)))
  plab <- c(paste0("b", seq_along(bmz)), paste0("y", seq_along(ymz)))
  segments(pm, 0, pm, -pn, lwd = 0.8, col = pc)
  sel <- pn > 0.10
  if (any(sel))
    text(pm[sel], -pn[sel] - 0.03, plab[sel], cex = 8/12, srt = 90,
         col = pc[sel], adj = 1)
  mtext(head, side = 3, line = 0.6, cex = 0.83, adj = 0, font = 2)
  if (nzchar(cos_line))
    mtext(cos_line, side = 3, line = 0.6, cex = 0.83, adj = 1)
}

# ppm 面板（isoPSV 原版）
draw_ppm_panel <- function(th, mz, xmax) {
  par(mar = c(2.6, 2.2, 0.5, 0.3), xpd = FALSE)
  plot(NA, xlim = c(0, xmax), ylim = c(-24, 24),
       xlab = "", ylab = "", type = "n", axes = FALSE, xaxs = "i", yaxs = "i")
  title(xlab = expression(italic("m/z")), cex.lab = 0.83, line = 1.8)
  title(ylab = "ppm", cex.lab = 0.75, line = 1.7)
  abline(h = 0, col = "#999999", lwd = 0.6)
  abline(h = c(-20, 20), col = "#BBBBBB", lty = 3, lwd = 0.6)
  if (!is.null(th)) {
    hm <- th[!is.na(th$peak), ]
    if (nrow(hm)) {
      px <- mz[hm$peak]
      py <- (px - hm$mz) / hm$mz * 1e6
      points(px, py, col = ION_COL[hm$ion], pch = 16, cex = 0.5)
    }
  }
  atk <- pretty(c(0, xmax), n = 10)
  atk <- atk[atk >= 0 & atk <= xmax]
  par(mgp = c(2.0, 0.8, 0))
  axis(1, at = atk, cex.axis = 0.83)
  axis(2, at = c(-20, 0, 20), las = 1, cex.axis = 0.83)
  box()
}

fmt_num <- function(v, dig) {
  x <- suppressWarnings(as.numeric(as.character(v)))
  if (is.na(x)) "-" else sprintf(paste0("%.", dig, "f"), x)
}

dir.create(dirname(OUT_PDF), showWarnings = FALSE, recursive = TRUE)
pdf(OUT_PDF, width = 8.27, height = 7.9, family = "Helvetica")
n_ok <- 0L
for (i in seq_len(nrow(pf))) {
  row <- pf[i, ]
  fn <- as.character(row$File_Name)
  spec_obj <- find_spec(fn)
  if (is.null(spec_obj)) { cat("skip (no spectrum):", fn, "\n"); next }
  mz0 <- as.numeric(spec_obj@spectrum[[1]]); ints0 <- as.numeric(spec_obj@spectrum[[2]])
  chg_v <- suppressWarnings(as.integer(spec_obj@charge))
  if (is.na(chg_v)) chg_v <- suppressWarnings(as.integer(row$Charge))
  if (is.na(chg_v) || chg_v < 1L) chg_v <- 0L
  base_pk <- max(ints0[is.finite(ints0) & ints0 > 0])
  bp_lab <- sub("E\\+0?", "E", sprintf("%.2E", base_pk))
  mz_pre <- suppressWarnings(as.numeric(sub("[^0-9.]", "", as.character(row$ExpMH))))
  mz_pre <- if (is.finite(mz_pre) && chg_v > 0) (mz_pre + (chg_v - 1) * H_M) / chg_v else NA
  ann <- c(sprintf("Base peak = %s", bp_lab),
           if (is.finite(mz_pre)) sprintf("Precursor m/z = %.4f (%d+)", mz_pre, chg_v) else "-")
  th <- all_ions(as.character(row$Sequence),
                 if (is.na(row$Modification) || row$Modification == "-") "" else as.character(row$Modification),
                 mz0[order(mz0, na.last = NA)], ints0[order(mz0, na.last = NA)])
  # 预测强度（b1..b_{n-1}, y1..y_{n-1} 1+）与 1+ m/z
  bip <- by_ions_position(as.character(row$Sequence),
                          if (is.na(row$Modification) || row$Modification == "-") "" else as.character(row$Modification),
                          "Unlabel", aa, aa_tab)
  nres <- nchar(as.character(row$Sequence))
  bmz <- ymz <- numeric(0); bint <- yint <- numeric(0)
  pv_ok <- !is.na(row$PredB) && nzchar(as.character(row$PredB))
  if (pv_ok && !is.null(bip$ions_table) && nrow(bip$ions_table) == nres - 1) {
    bmz <- bip$ions_table$b_1_charge
    ymz <- bip$ions_table$y_1_charge
    bint <- as.numeric(strsplit(as.character(row$PredB), ",", fixed = TRUE)[[1]])
    yint <- as.numeric(strsplit(as.character(row$PredY), ",", fixed = TRUE)[[1]])
  }
  cos_line <- ""
  if (!is.na(row$PredCos) && nzchar(as.character(row$PredCos)))
    cos_line <- sprintf("cosine(obs, pred) = %s (%s shared 1+ ions)",
                        as.character(row$PredCos),
                        if (is.na(row$PredN)) "-" else as.character(row$PredN))
  layout(matrix(c(1, 2, 3), 3, 1), heights = c(1.0, 0.85, 0.35))
  par(oma = c(2.0, 0.8, 3.2, 0.4)); par(cex = 1)
  g <- draw_obs_panel(spec_obj, th, chg_v, ann,
        sprintf("Observed spectrum - %s (scan %s)", as.character(row$Sample),
                as.character(row$Scan_No)))
  draw_ladder(as.character(row$Sequence), th,
              if (is.na(row$Modification)) "" else as.character(row$Modification), chg_v)
  if (pv_ok && length(bint) == length(bmz) && length(bint) > 0) {
    draw_pred_panel(bmz, ymz, bint, yint, g$xmax,
                    "AlphaPeptDeep prediction (inverted)", cos_line, g)
  } else {
    par(mar = c(0.5, 2.2, 1.6, 0.3))
    plot(NA, xlim = c(0, g$xmax), ylim = c(-1.15, 0), axes = FALSE, xlab = "", ylab = "")
    mtext("AlphaPeptDeep prediction unavailable for this precursor", side = 3,
          line = 0.6, cex = 0.83, adj = 0, font = 2)
  }
  draw_ppm_panel(th, g$mz, g$xmax)
  n_ok <- n_ok + 1L
  if (n_ok %% 20 == 0) cat("progress", n_ok, "/", nrow(pf), "\n")
}
dev.off()
cat(sprintf("done: %d pages -> %s\n", n_ok, OUT_PDF))
