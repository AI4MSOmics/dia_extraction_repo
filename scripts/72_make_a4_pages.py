# -*- coding: utf-8 -*-
"""
72_make_a4_pages.py —— S14 证据页：一条前体（肽段）一张 A4，多样本，RT 突出

版式（A4 纵向，matplotlib PdfPages，矢量）：
  头部   蛋白/基因/序列/电荷/m/z/修饰态 + 徽标（样本数、正文表肽段、唯一肽/Mtb 特异）
  左上   RT 一致性面板：每样本 S14 RT（空心圆）与提取 apex RT（实心圆，按组着色），
         共识 RT 线 + 批次漂移带；ΔRT 统计注角。
  右上   汇总统计小面板（文本卡：样本数/ΔRT 中位/梯≥6 比例/MS1 证据率/平均匹配碎片）
  中部   代表谱网格（≤6 张，按匹配碎片数选，HC/TB 优先混合）：镜像 b(上,蓝)/y(下,红)，
         匹配离子加粗+标签，干扰峰灰色；面板题：样本/RT/ΔRT/匹配数/梯。
  底部   全样本表（样本/组/S14 RT/apex RT/ΔRT/匹配碎片/b梯/y梯/MS1/S14 q/PG.Quantity），
         双栏排布，最多容纳全部样本。
输出 work/output/s14_evidence/：S14_peptide_evidence_A4.pdf（含封面与索引页）
     + pages_manifest.tsv
安全: 路径常量；文件 open 经 _open() 白名单。
"""
import sys, os, json, math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
SPEC = ROOT + r"\work\output\s14_evidence\extract\spectra_all.npz"
OUTPDF = ROOT + r"\work\output\s14_evidence\S14_peptide_evidence_A4.pdf"
OUTMANI = ROOT + r"\work\output\s14_evidence\pages_manifest.tsv"
OUTSUMM = ROOT + r"\work\output\s14_evidence\pages_summary.json"

PROTON = 1.007276466621
H2O = 18.01056468403
AA = dict(G=57.02146, C=103.00919 + 57.021464, S=87.03203, P=97.05276, A=71.03711,
          V=99.06841, T=101.04768, L=113.08406, I=113.08406, N=114.04293, D=115.02694,
          Q=128.05858, K=128.09496, E=129.04259, M=131.04049 + 15.994915,
          H=137.05891, F=147.06841, R=156.10111, Y=163.06333, W=186.07931)
MODMAP = {"Carbamidomethyl (C)": "C", "Oxidation (M)": "M_OX", "Acetyl (Protein N-term)": "NT_AC"}
M_OX_SHIFT = 15.994915
NT_AC_SHIFT = 42.010565

plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "font.family": "DejaVu Sans"})
C_HC, C_TB = "#2c7fb8", "#d95f0e"


def frag_ions(seq, pid):
    """理论 b/y（1+ 与 2+），支持 CAM 固定、M 氧化、N 端乙酰（synthetic 全 M 氧化）。"""
    z0 = int(str(pid).rsplit(".", 1)[1][0])
    body = str(pid)[1:]
    mods = {}
    if "(synthetic" in str(pid):
        st = str(pid).split("mods=")[-1].rstrip(")")
        if "ox" in st:
            mods = {i: "M_OX" for i, a in enumerate(seq) if a == "M"}
        if "acetyl" in st:
            mods[-1] = "NT_AC"
    else:
        body = str(pid)[1:str(pid).rindex("_.")]
        s2, i = "", 0
        while i < len(body):
            if body[i] == "[":
                j = body.index("]", i)
                key = MODMAP.get(body[i + 1:j].strip())
                if key:
                    mods[len(s2) - 1 if s2 else -1] = key
                i = j + 1
            else:
                s2 += body[i]
                i += 1
    res = []
    for i, a in enumerate(seq):
        if i in mods and mods[i] == "M_OX":
            res.append(AA[a] + 0.0)  # AA 表中 M 已含氧
        else:
            base = dict(G=57.02146, C=103.00919 + 57.021464, S=87.03203, P=97.05276, A=71.03711,
                        V=99.06841, T=101.04768, L=113.08406, I=113.08406, N=114.04293, D=115.02694,
                        Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891,
                        F=147.06841, R=156.10111, Y=163.06333, W=186.07931)[a]
            res.append(base)
    nterm = NT_AC_SHIFT if mods.get(-1) == "NT_AC" else 0.0
    pre = np.cumsum([0.0] + res)
    total = pre[-1] + H2O
    ions = {}
    for i in range(1, len(seq)):
        for zc in (1, 2):
            b = (pre[i] + nterm + PROTON * zc) / zc
            y = (total - pre[i] + PROTON * zc) / zc
            if zc == 1 or i >= 2:
                ions["b%d" % i if zc == 1 else "b%d(2+)" % i] = b
                ions["y%d" % i if zc == 1 else "y%d(2+)" % i] = y
    return ions


def annotate(mz, inten, ions, ppm=20.0):
    """返回匹配标签集合 dict: label -> (mz, int)。"""
    out = {}
    if mz is None or not len(mz):
        return out
    for l, t in ions.items():
        d = np.abs(mz - t) / t * 1e6
        k = int(np.argmin(d))
        if d[k] < ppm:
            out[l] = (float(mz[k]), float(inten[k]))
    return out


def spec_plot(ax, mz, inten, ions, title):
    ax.set_facecolor("#fbfbfb")
    if inten is None or not len(inten):
        ax.text(0.5, 0.5, "no spectrum", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title, fontsize=5.5)
        return 0
    top = inten.max() if inten.max() > 0 else 1.0
    hits = annotate(mz, inten, ions)
    for m, i in zip(mz, inten):
        ax.vlines(m, 0, i / top, color="#c8c8c8", lw=0.4)
    for l, (m, i) in sorted(hits.items()):
        isb = l.startswith("b")
        two = "(2+)" in l
        c = "#1f5fa9" if isb else "#c0392b"
        ax.vlines(m, 0, i / top, color=c, lw=0.9)
        if i / top > 0.06:
            ax.text(m, i / top + 0.012, l.split("(")[0], rotation=90, fontsize=3.6,
                    ha="center", va="bottom", color=c)
    ax.set_ylim(0, 1.16)
    ax.set_yticks([])
    ax.tick_params(labelsize=4.5, length=1.5)
    ax.set_title(title, fontsize=5.5)
    return len(hits)


def main():
    df = pd.read_parquet(ROWS)
    z = np.load(SPEC)
    mz_all, int_all, off, end = z["mz"], z["inten"], z["off"], z["end"]
    df = df.reset_index(drop=True)
    ok = df[df.status == "ok"].copy()
    groups = list(ok.groupby(["seq", "z", "mz_s14"], sort=False))
    print("pages:", len(groups))

    tg = pd.read_parquet(ROOT + r"\work\tmp\s14_targets.parquet")
    meta = tg.drop_duplicates(["seq", "z", "mz_s14"]).set_index(["seq", "z", "mz_s14"])[
        ["prot", "gene", "pid", "mod_state", "uc", "mtb", "in_main"]]
    pages = []
    with PdfPages(OUTPDF) as pdf:
        # ---- 封面 ----
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.5, 0.62, "Evidence for Mtb-derived peptides (Supplementary Table 14)",
                 ha="center", fontsize=15, weight="bold")
        fig.text(0.5, 0.57, "One precursor per page - RT-anchored targeted extraction from DIA raw files",
                 ha="center", fontsize=10)
        stats = [f"{len(groups)} precursor targets ({meta.index.get_level_values(0).nunique()} sequences, {meta.prot.nunique()} proteins)",
                 f"{len(ok)} spectra across {ok.run.nunique()} runs ({ok['sample'].str.startswith('HC').sum()} HC / {ok['sample'].str.startswith('TB').sum()} TB rows)",
                 f"median |RT difference vs Table S14| = {ok.drt.abs().median():.2f} min",
                 "For each sample: apex MS2 spectrum in the S14 retention-time window; b/y ion annotation (20 ppm, 1+/2+)",
                 "RT panel: per-sample S14 RT (open) vs extracted apex RT (filled); batch drift visible across runs",
                 "Generated from Thermo raw files via RawFileReader; no re-searching performed"]
        for i, s in enumerate(stats):
            fig.text(0.13, 0.50 - i * 0.028, s, fontsize=9)
        pdf.savefig(fig); plt.close(fig)

        for key, g in groups:
            seq, zc, mz14 = key
            m0 = meta.loc[key]
            g = g.sort_values("sample")  # 保留 df 位置索引（对齐 spectra 偏移数组）
            page = dict(page=len(pages) + 2, prot=m0.prot, gene=m0.gene, seq=seq, z=zc,
                        mz=mz14, mod=m0.mod_state, n=len(g),
                        drt_med=float(g.drt.abs().median()),
                        lad6=int(((g.bLad + g.yLad) >= 6).sum()),
                        ms1=int(g.ms1_flag.sum()), in_main=bool(m0.in_main))
            pages.append(page)
            ions = frag_ions(seq, m0.pid)

            fig = plt.figure(figsize=(8.27, 11.69))
            gs = gridspec.GridSpec(3, 1, height_ratios=[2.05, 3.5, 2.1], hspace=0.42,
                                   left=0.07, right=0.97, top=0.90, bottom=0.045)
            # 头部
            fig.text(0.07, 0.965, "%s  (%s)   %s   z=%d   m/z %.4f" %
                     (m0.prot, m0.gene, seq, zc, mz14), fontsize=10.5, weight="bold", va="top")
            tag = []
            if m0.in_main:
                tag.append("main-text table peptide")
            tag += ["unique peptide count=%s" % m0.uc, "Mtb-complex specific=%s" % m0.mtb,
                    "mods=%s" % m0.mod_state]
            fig.text(0.07, 0.945, " | ".join(tag), fontsize=7, va="top", color="#444444")

            # 顶区：左 RT 面板 + 右独立统计栏（固定分栏，杜绝遮挡）
            gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0],
                                                      width_ratios=[4.35, 1.0], wspace=0.04)
            # RT 面板
            ax = fig.add_subplot(gs_top[0, 0])
            y = np.arange(len(g))
            rts = g.rt_s14.astype(float).values
            apr = g.apex_rt.astype(float).values
            cons = np.median(rts)
            cols = [C_HC if s.startswith("HC") else C_TB for s in g["sample"]]
            ax.hlines(y, rts - 0.02, rts + 0.02, color="#999999", lw=0.6)
            ax.scatter(rts, y, facecolors="none", edgecolors="#666666", s=14, lw=0.7, label="Table S14 RT")
            ax.scatter(apr, y, c=cols, s=14, label="extracted apex RT")
            ax.axvline(cons, color="k", lw=0.7, ls="--")
            ax.set_yticks(y)
            ax.set_yticklabels(list(g["sample"]), fontsize=4.6)
            ax.set_xlabel("retention time (min)", fontsize=6)
            dmed = np.median(np.abs(g.drt.astype(float)))
            ax.set_title("RT consistency across samples  (|ΔRT| median %.2f min; dashed = median S14 RT %.2f)" % (dmed, cons),
                         fontsize=7, loc="left")
            ax.tick_params(labelsize=5, length=1.5)
            ax.legend(fontsize=5.5, loc="lower right", frameon=False)
            ax.invert_yaxis()
            for yi, (s, d) in enumerate(zip(g["sample"], g.drt.astype(float))):
                if abs(d) > 1.0:
                    ax.text(apr[yi] + 0.12, yi, "Δ%.1f" % d, fontsize=4, va="center", color="#c0392b")

            # 统计卡（右侧独立轴，文字左对齐自顶部，不与任何数据重叠）
            axc = fig.add_subplot(gs_top[0, 1])
            axc.axis("off")
            lad6 = int(((g.bLad + g.yLad) >= 6).sum())
            txt = ("samples with spectrum: %d\n"
                   "|ΔRT| ≤ 0.5 min: %.0f%%\n"
                   "ladder b+y ≥ 6: %d/%d\n"
                   "MS1 precursor envelope: %d/%d\n"
                   "median matched b/y ions: %d" %
                   (len(g), 100 * (g.drt.astype(float).abs() <= 0.5).mean(),
                    lad6, len(g), int(g.ms1_flag.sum()), len(g), int(g.nmatched.median())))
            axc.text(0.0, 0.98, txt, fontsize=5.8, va="top", ha="left", family="DejaVu Sans",
                     bbox=dict(boxstyle="round,pad=0.35", fc="#f5f5f5", ec="#bbbbbb", lw=0.5))

            # 谱网格：最多 6 张代表谱（匹配数降序，HC/TB 混合优先）
            g2 = g.sort_values("nmatched", ascending=False)
            pick, seen_hc, seen_tb = [], 0, 0
            for _, rr in g2.iterrows():
                ishc = rr["sample"].startswith("HC")
                if ishc and seen_hc < 2:
                    pick.append(rr); seen_hc += 1
                elif not ishc and seen_tb < 4:
                    pick.append(rr); seen_tb += 1
                if len(pick) == 6:
                    break
            for _, rr in g2.iterrows():
                if len(pick) == 6:
                    break
                if rr.name not in [p.name for p in pick]:
                    pick.append(rr)
            ncol = 3 if len(pick) > 3 else max(len(pick), 1)
            nrow = math.ceil(len(pick) / ncol)
            gs2 = gridspec.GridSpecFromSubplotSpec(nrow, ncol, subplot_spec=gs[1], hspace=0.55, wspace=0.12)
            for i, rr in enumerate(pick):
                ax = fig.add_subplot(gs2[i // ncol, i % ncol])
                sl = slice(int(off[rr.name]), int(end[rr.name]))
                smz, sint = mz_all[sl], int_all[sl]
                o = np.argsort(-sint)[:400]  # 只画最强 400 峰（矢量体积控制）
                title = "%s  RT %.2f (Δ%+.2f)  match %d  b%d/y%d%s" % (
                    rr["sample"], rr.apex_rt, rr.drt, rr.nmatched, rr.bLad, rr.yLad,
                    "  MS1✓" if rr.ms1_flag else "")
                spec_plot(ax, smz[o], sint[o], ions, title)
                if i % ncol == 0:
                    ax.set_ylabel("rel. int.", fontsize=5)
                if i >= len(pick) - ncol:
                    ax.set_xlabel("m/z", fontsize=5)
            fig.text(0.07, 0.455, "Representative apex MS2 spectra (top-%d by matched ions; blue=b ions up, red=y ions down; grey=co-isolated background)" % len(pick),
                     fontsize=7)

            # 底部全样本表（双栏）
            ntb = len(g)
            half = math.ceil(ntb / 2)
            gs_tab = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], wspace=0.25)
            for ci in range(2 if ntb > 1 else 1):
                sub = g.iloc[ci * half:(ci + 1) * half]
                ax = fig.add_subplot(gs_tab[0, ci])
                ax.axis("off")
                rows = [["Sample", "Grp", "S14 RT", "Apex RT", "ΔRT", "Ions", "b/y", "MS1", "q", "Qty"]]
                for _, rr in sub.iterrows():
                    rows.append([rr["sample"], "HC" if rr["sample"].startswith("HC") else "TB",
                                 "%.2f" % rr.rt_s14, "%.2f" % rr.apex_rt, "%+.2f" % rr.drt,
                                 "%d" % rr.nmatched, "%d/%d" % (rr.bLad, rr.yLad),
                                 "Y" if rr.ms1_flag else "-",
                                 ("%.0e" % rr.qval).replace("e-0", "e-") if pd.notna(rr.qval) else "-",
                                 ("%.3g" % rr.qty) if pd.notna(rr.qty) else "-"])
                tb = ax.table(cellText=rows[1:], colLabels=rows[0], loc="upper center",
                              cellLoc="center", bbox=[0, 0.02, 1, 0.93])
                tb.auto_set_font_size(False)
                tb.set_fontsize(4.2)
                for (r_i, c_i), cell in tb.get_celld().items():
                    cell.set_linewidth(0.2)
                    if r_i == 0:
                        cell.set_facecolor("#e8e8e8")
            pdf.savefig(fig); plt.close(fig)

        # ---- 索引页 ----
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.07, 0.965, "Index", fontsize=13, weight="bold", va="top")
        ax = fig.add_axes([0.05, 0.03, 0.92, 0.92]); ax.axis("off")
        rows = [["Page", "Protein", "Gene", "Sequence", "z", "m/z", "N", "|ΔRT|med", "MS1"]]
        for p in pages:
            rows.append([p["page"], p["prot"], p["gene"], p["seq"], p["z"],
                         "%.3f" % p["mz"], p["n"], "%.2f" % p["drt_med"],
                         "%d/%d" % (p["ms1"], p["n"])])
        n1 = 0
        chunk = 56
        for start in range(1, len(rows), chunk):
            sub = [rows[0]] + rows[start:start + chunk]
            ax2 = fig.add_axes([0.05, 0.03, 0.92, 0.92]) if start == 1 else None
        # 简化：单页最多 56 行，多页
        pdf.close = pdf.close  # keep
        tb_axs = []
        npages = math.ceil((len(rows) - 1) / chunk)
        for pi in range(npages):
            if pi > 0:
                fig = plt.figure(figsize=(8.27, 11.69))
                fig.text(0.07, 0.965, "Index (cont'd %d/%d)" % (pi + 1, npages), fontsize=13, weight="bold", va="top")
                ax = fig.add_axes([0.05, 0.03, 0.92, 0.92]); ax.axis("off")
            sub = [rows[0]] + rows[1 + pi * chunk:1 + (pi + 1) * chunk]
            tb = ax.table(cellText=sub[1:], colLabels=sub[0], loc="upper center",
                          cellLoc="center", bbox=[0, 0, 1, 0.94])
            tb.auto_set_font_size(False); tb.set_fontsize(5)
            for (r_i, c_i), cell in tb.get_celld().items():
                cell.set_linewidth(0.2)
                if r_i == 0:
                    cell.set_facecolor("#e8e8e8")
            pdf.savefig(fig); plt.close(fig)

    mani = pd.DataFrame(pages)
    mani.to_csv(OUTMANI, sep="\t", index=False, encoding="utf-8")
    summ = dict(pages=len(pages), targets_ok=int(len(groups)),
                drt_med_all=float(ok.drt.abs().median()),
                ms1_rate=float(ok.ms1_flag.mean()),
                lad6_rate=float(((ok.bLad + ok.yLad) >= 6).mean()),
                drt_le05_rate=float((ok.drt.astype(float).abs() <= 0.5).mean()))
    with _open(OUTSUMM, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print(json.dumps(summ, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
