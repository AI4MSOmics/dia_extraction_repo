# -*- coding: utf-8 -*-
"""
79_compose_a4.py —— 合成最终 A4：R 顶块（isoPSV 版式前五区块）+ 热图

每肽一页 A4（8.27x11.69in）：
  上 62%：78_s14_export.R 输出的顶块页（信息行/离子梯/原始谱/倒置预测/ppm）
  下 38%：热图（1 min 分箱 apex 确认离子靶向 XIC，白->红，无三角）matplotlib 渲染 PNG 后嵌入
封面与索引页由本脚本 matplotlib 生成。输出 S14_peptide_evidence_A4_isoPSV.pdf。
安全：路径常量；文件 open 经 _open() 白名单。
"""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pymupdf

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
TARGETS = ROOT + r"\work\tmp\s14_targets.parquet"
PROF_IDX = ROOT + r"\work\output\s14_evidence\profiles\index.parquet"
PROF_NPZ = ROOT + r"\work\output\s14_evidence\profiles\profiles.npz"
INTSV = ROOT + r"\work\risoPSV\in\s14_spectra.tsv"
TOPPDF = ROOT + r"\work\risoPSV\out\s14_topblock.pdf"
OUTPDF = ROOT + r"\work\output\s14_evidence\S14_peptide_evidence_A4_isoPSV.pdf"
OUTMANI = ROOT + r"\work\output\s14_evidence\pages_manifest_isoPSV.tsv"
OUTSUMM = ROOT + r"\work\output\s14_evidence\pages_summary_isoPSV.json"

plt.rcParams.update({"font.size": 7, "font.family": "Arial", "pdf.fonttype": 42})
HEAT = LinearSegmentedColormap.from_list("wr", ["#FFFFFF", "#FDD8B5", "#F68838", "#D7301f"])


def heatmap_png(seq, z, mz, gk, path, prof_idx=None, prof_npz=None):
    """热图渲染为 PNG（profiles 可指定主/补充两套；返回是否有内容）。"""
    pidx = pd.read_parquet(prof_idx or PROF_IDX)
    pz = np.load(prof_npz or PROF_NPZ)
    prof_map = {(r.seq, r.z, r.mz, r["sample"]): i for i, r in pidx.reset_index().iterrows()}
    samples = sorted(gk["sample"].unique(),
                     key=lambda s: (0 if s.startswith("HC") else 1, int(s[2:])))
    prof_rows = []
    for s in samples:
        i = prof_map.get((seq, int(z), float(mz), s))
        if i is None:
            continue
        r0, r1 = int(pz["off"][i]), int(pz["end"][i])
        prof_rows.append((s, pz["rt"][r0:r1], pz["val"][r0:r1]))
    if not prof_rows:
        return False
    fig = plt.figure(figsize=(8.0, 5.87), dpi=200)  # aspect 1.363 = 目标矩形 555x407pt 满铺
    ax = fig.add_axes([0.135, 0.108, 0.775, 0.816])  # 三带重分：标题31/色带332/底44pt，两侧余隙均衡
    b0 = math.floor(min(min(rt) for _, rt, _ in prof_rows))
    b1 = math.ceil(max(max(rt) for _, rt, _ in prof_rows)) + 1
    nb = b1 - b0
    M = np.zeros((len(prof_rows), nb))
    for ri, (s, rt, val) in enumerate(prof_rows):
        np.maximum.at(M[ri], np.clip((rt - b0).astype(int), 0, nb - 1), val)
    # 行级定量阈值（背景强度去除）：
    #   各行以该样本 S14/锚 RT ±2 min 为洗脱窗，窗内峰值 < 2×窗外中位基线
    #   → 该样本此肽无可定量峰（碎片埋于背景，相关剖面=背景纹理）→ 整行置白
    #   并从 apex RT 统计剔除。
    rt_ref = {r["sample"]: float(r.apex_rt) for _, r in gk.iterrows()}
    bin_c = b0 + np.arange(nb) + 0.5
    n_nq = 0
    base_row = np.zeros(M.shape[0])
    for ri, (s, _, _) in enumerate(prof_rows):
        r0 = rt_ref.get(s)
        if r0 is None:
            continue
        win = np.abs(bin_c - r0) <= 2.0
        if not win.any():
            continue
        peak = M[ri][win].max()
        base = np.median(M[ri][~win]) if (~win).any() else 0.0
        base_row[ri] = base
        if peak < 2.0 * base:
            M[ri, :] = 0.0
            n_nq += 1
    # 显示值 γ=2 锐化 + 双重卡阈：相对线（<20% 行最高不呈现）+ 绝对背景线
    #（原始值 < 1.5×行基线不算定量——动态范围仅 2–4 倍的弱峰行由此回归窄带）
    Mraw = M.copy()
    M = np.clip(M, 0.0, 1.0) ** 2
    row_max = M.max(axis=1, keepdims=True)
    M = np.where(M < 0.20 * row_max, 0.0, M)
    M = np.where(Mraw < 1.5 * base_row[:, None], 0.0, M)
    # 各样本最高定量峰的 RT（bin 中心）→ mean ± SD（仅计入有定量的行）
    apex_rts = []
    for ri in range(M.shape[0]):
        k = int(M[ri].argmax())
        if M[ri, k] > 0:
            apex_rts.append(b0 + k + 0.5)
    stat_txt = ""
    if len(apex_rts) >= 2:
        arr = np.array(apex_rts)
        stat_txt = "apex RT (max-quantity bin): %.2f ± %.2f min (mean ± SD, n=%d" % (
            float(arr.mean()), float(arr.std(ddof=1)), len(apex_rts))
        if n_nq:
            stat_txt += "; %d row(s) below quant threshold" % n_nq
        stat_txt += ")"
    elif len(apex_rts) == 1:
        stat_txt = "apex RT: %.2f min (n=1)" % apex_rts[0]
        if n_nq:
            stat_txt += "; %d row(s) below quant threshold" % n_nq
    vmax = M.max()
    if vmax > 0:
        M = M / vmax
    im = ax.imshow(M, aspect="auto", cmap=HEAT, vmin=0, vmax=1,
                   extent=[b0, b1, len(prof_rows) - 0.5, -0.5], interpolation="nearest")
    # 整行白（无定量）的行：行内居中标注该样本 Table S14 定量值（灰色）
    qty_ref = {}
    for _, r in gk.iterrows():
        try:
            q = float(r.qty)
            if q == q:  # not NaN
                qty_ref[r["sample"]] = q
        except (TypeError, ValueError):
            pass
    mid_rt = (b0 + b1) / 2.0
    for ri, (s, _, _) in enumerate(prof_rows):
        if M[ri].max() <= 0 and s in qty_ref:
            ax.text(mid_rt, ri, "%.2e" % qty_ref[s], ha="center", va="center",
                    fontsize=5.4, color="#666666")
    ax.set_yticks(range(len(prof_rows)))
    ax.set_yticklabels([x[0] for x in prof_rows], fontsize=6)
    ax.set_xlabel("retention time (min)", fontsize=8)
    ax.xaxis.set_label_coords(0.5, -0.062)
    ax.set_xlim(b0, b1)
    ax.tick_params(labelsize=6.5)
    ax.set_title("Fragment-ion XIC heatmap (1-min bins; white = 0, red = max; per-row cutoff = 20% of "
                 "row max; grey value in blank row = Table S14 quantity of that sample)\n" + stat_txt,
                 fontsize=7.2, loc="left", pad=4)
    cb = fig.colorbar(im, ax=ax, pad=0.01, fraction=0.04)
    cb.ax.tick_params(labelsize=6)
    cb.set_label("rel. XIC", fontsize=6.5)
    fig.savefig(path)
    plt.close(fig)
    return True


def fit_text(page, x, y, text, max_w, size, bold=False):
    """插入文字并自动缩字号至宽 max_w 内（pymupdf helv 量宽）。返回下一行 y。"""
    font = pymupdf.Font("hebo" if bold else "helv")
    s = size
    while s > 4 and font.text_length(text, fontsize=s) > max_w:
        s -= 0.25
    page.insert_text((x, y), text, fontsize=s, fontname="hebo" if bold else "helv",
                     color=(0, 0, 0))
    return y + s * 1.5


def main():
    only = ""
    for a in sys.argv[1:]:
        if a.startswith("--only="):
            only = a.split("=", 1)[1]
    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    ok = rows[rows.status == "ok"]
    tg = pd.read_parquet(TARGETS)
    seqs = list(ok.seq.unique())
    if only:
        seqs = [s for s in seqs if s == only]
        assert seqs, "not found"
    print("peptides:", len(seqs))
    toppdf = TOPPDF
    outpdf = OUTPDF
    if only:
        toppdf = ROOT + r"\work\risoPSV\out\sample_topblock.pdf"
        outpdf = ROOT + r"\work\output\s14_evidence\S14_evidence_sample_A4.pdf"

    tmpdir = ROOT + r"\work\output\s14_evidence\_heat_tmp"
    os.makedirs(tmpdir, exist_ok=True)
    info_map = {}
    if os.path.exists(INTSV):
        itab = pd.read_csv(INTSV, sep="\t", encoding="utf-8")
        info_map = {r.Sequence: (r.Info1, r.Info2, r.Info3)
                    for r in itab.drop_duplicates("Sequence").itertuples()}
    doc = pymupdf.open()
    topdoc = pymupdf.open(toppdf)
    pages = []
    for pi, seq in enumerate(seqs):
        g = ok[ok.seq == seq]
        b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
        gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
        m0 = tg[(tg.seq == seq) & (tg.z == b.z) & (tg.mz_s14 == b.mz_s14)].iloc[0]
        page = doc.new_page(width=595, height=842)  # A4 pt
        # 信息行（三行，自动缩字号防裁切）
        i1, i2, i3 = info_map.get(seq, ("", "", ""))
        y = fit_text(page, 14, 26, i1, 570, 11, bold=True)
        y = fit_text(page, 14, y + 1, i2, 570, 8)
        y = fit_text(page, 14, y + 1, i3, 570, 8)
        top_h = y + 2
        # 顶块（紧贴信息行；热图高度 ×1.25 → 顶块下缘相应上收，块间空隙压至 ~3pt）
        r = pymupdf.Rect(2, top_h, 593, 442)
        page.show_pdf_page(r, topdoc, pi if pi < topdoc.page_count else 0)
        # 热图（紧贴 ppm 轴题之下：顶 430pt，底 839pt）
        hp = os.path.join(tmpdir, "h%03d.png" % pi)
        has = heatmap_png(seq, b.z, b.mz_s14, gk, hp)
        if has:
            page.insert_image(pymupdf.Rect(20, 431, 575, 838), filename=hp)
        else:
            page.insert_text((30, 600), "no profiles", fontsize=10)
        pages.append(dict(page=len(pages) + 1, prot=m0.prot, gene=m0.gene, seq=seq,
                          z=int(b.z), mz=float(b.mz_s14),
                          flank="-", samples=len(gk),
                          rt_lo=float(gk.apex_rt.min()), rt_hi=float(gk.apex_rt.max()),
                          best=b["sample"], nmatched=int(b.nmatched), ms1=int(b.ms1_flag),
                          in_main=bool(m0.in_main)))
    # 封面插到最前
    cov = doc.new_page(width=595, height=842)
    topdoc2 = None
    # 用 matplotlib 生成封面/索引 PNG 再插入
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.5, 0.60, "Evidence for Mtb-derived peptides (Table S14)", ha="center",
             fontsize=16, weight="bold", family="Arial")
    fig.text(0.5, 0.555, "isoPSV-layout evidence pages: info lines, ion ladder, observed "
             "spectrum, inverted AlphaPeptDeep prediction, ppm panel, RT heatmap",
             ha="center", fontsize=9)
    for i, s in enumerate([
            "%d peptide sequences, %d samples (7 HC / 23 TB), hebei DIA cohort" %
            (len(seqs), rows["sample"].nunique()),
            "Top block rendered by the isoPSV R engine (single_spectrum_export.R drawing engine)",
            "Observed = RT-anchored targeted extraction, best sample by matched ions",
            "Prediction = AlphaPeptDeep 1.5.0 (inverted); heatmap = apex-confirmed-ion XIC, 1-min bins"]):
        fig.text(0.12, 0.49 - i * 0.028, s, fontsize=9)
    covpng = os.path.join(tmpdir, "_cover.png")
    fig.savefig(covpng, dpi=150); plt.close(fig)
    cov.insert_image(pymupdf.Rect(0, 0, 595, 842), filename=covpng)
    doc.move_page(len(doc) - 1, 0)
    # 索引页
    import math as _m
    chunk = 52
    npg = _m.ceil(len(pages) / chunk)
    for ci in range(npg):
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.965, "Index%s" % ("" if ci == 0 else " (cont'd %d/%d)" % (ci + 1, npg)),
                 fontsize=13, weight="bold", va="top", family="Arial")
        ax = fig.add_axes([0.04, 0.03, 0.94, 0.92]); ax.axis("off")
        hdr = ["Page", "Protein", "Gene", "Sequence", "z", "m/z", "N", "RT range", "Best PSM"]
        body = []
        for p in pages[ci * chunk:(ci + 1) * chunk]:
            body.append([p["page"], p["prot"], p["gene"], p["seq"], p["z"],
                         "%.4f" % p["mz"], p["samples"],
                         "%.2f-%.2f" % (p["rt_lo"], p["rt_hi"]),
                         "%s(%d)" % (p["best"], p["nmatched"])])
        body = [[(str(c)[:30] + "...") if len(str(c)) > 31 else c for c in row] for row in body]
        tb = ax.table(cellText=body, colLabels=hdr, loc="upper center", cellLoc="center",
                      bbox=[0, 0, 1, 0.94],
                      colWidths=[0.06, 0.10, 0.11, 0.31, 0.04, 0.10, 0.04, 0.12, 0.12])
        tb.auto_set_font_size(False); tb.set_fontsize(5)
        for (r_i, c_i), cell in tb.get_celld().items():
            cell.set_linewidth(0.2)
            if r_i == 0:
                cell.set_facecolor("#e8e8e8")
        ip = os.path.join(tmpdir, "_idx%d.png" % ci)
        fig.savefig(ip, dpi=150); plt.close(fig)
        pg = doc.new_page(width=595, height=842)
        pg.insert_image(pymupdf.Rect(0, 0, 595, 842), filename=ip)
    doc.save(outpdf, garbage=3, deflate=True)
    pd.DataFrame(pages).to_csv(OUTMANI, sep="\t", index=False, encoding="utf-8")
    summ = dict(peptide_pages=len(pages), total_pages=doc.page_count)
    with _open(OUTSUMM, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print(json.dumps(summ, ensure_ascii=False))


if __name__ == "__main__":
    main()
