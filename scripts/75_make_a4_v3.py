# -*- coding: utf-8 -*-
"""
75_make_a4_v3.py —— S14 证据页 v3：isoPSV 排版复刻（用户指定顺序）

顺序（A4 纵向，Arial，除热图外全部按 isoPSV single_spectrum_export.R /
spectrum_pdf_export_A.R 的 v7 版式）：
  1 信息行    （outer-title 风格：首行大字，后两行 0.83 号；内容=蛋白/侧翼肽段/价态/
              修饰/m/z/位点 与 RT 范围/最优样本）
  2 离子梯    （isoPSV draw_ladder：序列字母行、b 下绿标、y 上红标、切割位竖线、
              修饰残基橙、左端红 z+；独立窄轴、与谱共用 x 轴）
  3 原始谱图  （isoPSV 谱面板：x 从 0 起、灰60 背景 0.3 线宽、匹配彩棒 0.8 +
              旋转 90° 标签、右上 Base peak/前体注释、y 轴 0–100%）
  4 预测谱图（倒置）（同一画法镜像向下；面板题 AlphaPeptDeep；右上 cosine 注释）
  5 质量差分布（isoPSV ppm 面板：±20 虚线、匹配离子色点、pretty x 刻度、
              x 轴题 m/z 斜体）
  6 热图      （1 min 分箱碎片 XIC，白→红，无三角）
用法：python 75_make_a4_v3.py [--only=SEQ]；--only 时另存 300dpi PNG。
安全：路径常量；文件 open 经 _open() 白名单。
"""
import sys, os, json, math, pickle, importlib.util
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


TARGETS = ROOT + r"\work\tmp\s14_targets.parquet"
ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
SPECNPZ = ROOT + r"\work\output\s14_evidence\extract\spectra_all.npz"
PROF_IDX = ROOT + r"\work\output\s14_evidence\profiles\index.parquet"
PROF_NPZ = ROOT + r"\work\output\s14_evidence\profiles\profiles.npz"
PRED = ROOT + r"\work\output\summary\pred\AlphaPeptDeep.pkl"
FASTA = ROOT + r"\IPX0018588001_hebei\human_tuberculosis.fasta"
OUTPDF = ROOT + r"\work\output\s14_evidence\S14_peptide_evidence_A4_v3.pdf"
OUTMANI = ROOT + r"\work\output\s14_evidence\pages_manifest_v3.tsv"
OUTSUMM = ROOT + r"\work\output\s14_evidence\pages_summary_v3.json"

_spec = importlib.util.spec_from_file_location("ext71", ROOT + r"\work\scripts\71_extract_targeted.py")
m71 = importlib.util.module_from_spec(_spec)
sys.modules["ext71"] = m71
_spec.loader.exec_module(m71)

C_B, C_Y, C_MOD = "#008000", "#CA1D52", "#FF8C00"   # isoPSV ION_COL
GREY = "#999999"  # R grey60
plt.rcParams.update({
    "font.size": 7, "pdf.fonttype": 42, "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
HEAT = LinearSegmentedColormap.from_list("wr", ["#FFFFFF", "#FDD8B5", "#F68838", "#D7301f"])


def load_fasta_flanks():
    prots, name, buf = {}, None, []
    with _open(FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    prots[name] = "".join(buf)
                f = line[1:].split()[0].split("|")
                name = f[1] if len(f) >= 3 and f[0] == "Mtb" else f[-1]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        prots[name] = "".join(buf)
    return prots


def flank_of(prots, prot, seq):
    p = prots.get(prot)
    if not p:
        return "-", "-", None
    i = p.find(seq)
    if i < 0:
        return "-", "-", None
    return (p[i - 1] if i > 0 else "-",
            p[i + len(seq)] if i + len(seq) < len(p) else "-", i + 1)


def mod_positions(pid):
    if "(synthetic" in str(pid):
        st = str(pid).split("mods=")[-1].rstrip(")")
        pos = []
        if "ox" in st:
            seq = str(pid)[1:str(pid).index("_.")]
            pos += [i + 1 for i, a in enumerate(seq) if a == "M"]
        if "acetyl" in st:
            pos.append(0)
        return pos
    body = str(pid)[1:str(pid).rindex("_.")]
    s2, pos, i = "", [], 0
    while i < len(body):
        if body[i] == "[":
            j = body.index("]", i)
            pos.append(len(s2))
            i = j + 1
        else:
            s2 += body[i]
            i += 1
    return pos


def iso_ladder_ax(fig, spec_ax, seq, hits, mods_pos, z, xmax):
    """isoPSV draw_ladder 复刻：独立窄轴（谱上方），与谱共用 x 轴。"""
    ax = fig.add_subplot(spec_ax, sharex=None)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1)
    ax.axis("off")
    n = len(seq)
    x0, x1 = xmax * 0.015, xmax * 0.985
    xs = np.linspace(x0 + (x1 - x0) * 0.02, x1 - (x1 - x0) * 0.02, n)
    slot = xs[1] - xs[0] if n > 1 else 1.0
    ymid = 0.50
    b_sites = sorted({int(l[1:]) for l in hits if l.startswith("b") and l[1:].isdigit()})
    y_sites = sorted({n - int(l[1:]) for l in hits if l.startswith("y") and l[1:].isdigit()})
    for k in sorted(set(b_sites) | set(y_sites)):
        if k < n:
            ax.vlines(xs[k] + slot / 2, ymid - 0.16, ymid + 0.16, color="k", lw=1.2)
    for k in b_sites:
        if k < n:
            ax.hlines(ymid - 0.16, xs[k] + slot / 2 - slot * 0.45, xs[k] + slot / 2, color="k", lw=1.2)
            ax.text(xs[k], ymid - 0.22, r"$b_{%d}$" % k, fontsize=5.5, color=C_B,
                    ha="center", va="top")
    for k in y_sites:
        if k < n:
            ax.hlines(ymid + 0.16, xs[k] + slot / 2, xs[k] + slot / 2 + slot * 0.45, color="k", lw=1.2)
            ax.text(xs[k + 1] if k + 1 < n else xs[-1], ymid + 0.22,
                    r"$y_{%d}$" % (n - k), fontsize=5.5, color=C_Y, ha="center", va="bottom")
    ax.text(xs[0] - slot * 1.5, ymid, "%d+" % z, fontsize=7.5, color=C_Y,
            ha="right", va="center", weight="bold")
    mp = set(mods_pos)
    for i, (x, a) in enumerate(zip(xs, seq)):
        ax.text(x, ymid, a, fontsize=6.5, ha="center", va="center",
                color=C_MOD if (i + 1) in mp else "k",
                weight="bold" if (i + 1) in mp else "normal")


def iso_spectrum_ax(fig, spec_ax, mz, ints, hits, xmax, ann_lines, head,
                    inverted=False, ylab="Relative intensity (%)"):
    """isoPSV 谱面板复刻（可选倒置）。hits: label->(mz,int,ppm)。"""
    ax = fig.add_subplot(spec_ax)
    ax.set_xlim(0, xmax)
    top = 1.15
    ax.set_ylim(top, 0.0) if inverted else ax.set_ylim(0.0, top)
    ax.set_facecolor("white")
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=5.2)
    ax.set_ylabel(ylab, fontsize=6.2)
    ax.set_xticklabels([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2, pad=1.5)
    sgn = -1.0 if inverted else 1.0
    for m, i in zip(mz, ints):
        ax.vlines(m, 0, sgn * i, color=GREY, lw=0.3)
    # 标签防撞：相邻 2 Da 内交替抬升
    labs = sorted(hits.items(), key=lambda kv: kv[1][0])
    lift = 0.0
    prev_m = -9.0
    for l, (m, i, _) in labs:
        c = C_B if l.startswith("b") else C_Y
        h = max(i, 0.012)
        ax.vlines(m, 0, sgn * h, color=c, lw=0.8)
        if m - prev_m < 2.0:
            lift = 0.10 if lift == 0.0 else 0.0
        else:
            lift = 0.0
        ax.text(m, sgn * (h + 0.03 + lift), l.split("(")[0], rotation=90, fontsize=4.6,
                ha="center", va="bottom" if not inverted else "top", color=c)
        prev_m = m
    for k, line in enumerate(ann_lines):
        ax.text(xmax * 0.995, sgn * (1.10 - 0.10 * k), line, fontsize=5.8,
                ha="right", va="center", family="Arial",
                bbox=dict(fc="white", ec="none", pad=1.2))
    ax.text(0.008, 1.02 if not inverted else -1.02, head, transform=ax.transAxes,
            fontsize=6.2, ha="left", va="bottom" if not inverted else "top",
            family="Arial", weight="bold")
    return ax


def iso_ppm_ax(fig, spec_ax, hits, xmax):
    """isoPSV ppm 面板复刻。"""
    ax = fig.add_subplot(spec_ax)
    ax.set_xlim(0, xmax)
    ax.set_ylim(-24, 24)
    ax.axhline(0, color="#999999", lw=0.6)
    for hh in (-20, 20):
        ax.axhline(hh, color="#BBBBBB", ls=":", lw=0.6)
    for l, (m, i, ppm) in hits.items():
        ax.plot(m, ppm, "o", ms=1.6, color=C_B if l.startswith("b") else C_Y)
    ax.set_yticks([-20, 0, 20])
    ax.set_ylabel("ppm", fontsize=6.2)
    ax.set_xlabel(r"$m/z$", fontsize=6.5, style="italic")
    atk = np.linspace(0, xmax, 11)
    ax.set_xticks(atk)
    ax.tick_params(labelsize=5.2, length=2, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    return ax


def main():
    only = ""
    for a in sys.argv[1:]:
        if a.startswith("--only="):
            only = a.split("=", 1)[1]
    outpdf = OUTPDF
    if only:
        outpdf = ROOT + r"\work\output\s14_evidence\S14_evidence_sample_v3_%s.pdf" % only[:24]

    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    znp = np.load(SPECNPZ)
    mz_all, int_all, off, end = znp["mz"], znp["inten"], znp["off"], znp["end"]
    tg = pd.read_parquet(TARGETS)
    prots = load_fasta_flanks()
    pred = pickle.load(_open(PRED, "rb"))
    pidx = pd.read_parquet(PROF_IDX)
    pz = np.load(PROF_NPZ)
    prof_map = {(r.seq, r.z, r.mz, r["sample"]): i for i, r in pidx.reset_index().iterrows()}

    ok = rows[rows.status == "ok"]
    seqs = list(ok.seq.unique())
    if only:
        seqs = [s for s in seqs if s == only]
        assert seqs, "sequence not found: " + only
    print("sequences:", len(seqs))

    pages = []
    with PdfPages(outpdf) as pdf:
        if not only:
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.5, 0.60, "Evidence for Mtb-derived peptides (Table S14) — v3",
                     ha="center", fontsize=15, weight="bold", family="Arial")
            fig.text(0.5, 0.555, "isoPSV-style layout: info lines, ion ladder, observed spectrum, "
                     "inverted AlphaPeptDeep prediction, ppm panel, RT heatmap",
                     ha="center", fontsize=9)
            for i, s in enumerate([
                    f"{len(seqs)} peptide sequences ({tg.drop_duplicates(['prot']).prot.nunique()} proteins), "
                    f"{rows['sample'].nunique()} samples (7 HC / 23 TB)",
                    "Observed = RT-anchored targeted extraction, best sample by matched b/y ions",
                    "Prediction = AlphaPeptDeep 1.5.0 (inverted); cosine on shared 1+ b/y ions"]):
                fig.text(0.13, 0.49 - i * 0.028, s, fontsize=9)
            pdf.savefig(fig); plt.close(fig)

        for seq in seqs:
            g = ok[ok.seq == seq]
            b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
            gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
            m0 = tg[(tg.seq == seq) & (tg.z == b.z) & (tg.mz_s14 == b.mz_s14)].iloc[0]
            pid = m0.pid
            frags, _, _ = m71.parse_pid(pid)
            pre, post, pos = flank_of(prots, m0.prot, seq)
            rts_all = gk.apex_rt.astype(float)
            rt_lo, rt_hi = float(rts_all.min()), float(rts_all.max())
            mods_txt = (("Carbamidomethyl (C); " if seq.count("C") else "") +
                        ("Oxidation (M); " if "ox" in str(m0.mod_state) else "") +
                        ("Acetyl (N-term)" if "acetyl" in str(m0.mod_state) else "")).strip("; ") or "none"
            if str(m0.mod_state).startswith("unknown"):
                mods_txt += " [" + str(m0.mod_state) + "]"

            # 谱与匹配（匹配峰必进 + 背景 top-400）
            sl = slice(int(off[b.name]), int(end[b.name]))
            smz, sint = mz_all[sl], int_all[sl]
            hits = {}
            hit_idx = set()
            if len(smz):
                for l, t in frags.items():
                    d = np.abs(smz - t) / t * 1e6
                    k = int(np.argmin(d))
                    if d[k] < 20.0:
                        hits[l] = (float(smz[k]), float(sint[k]), (smz[k] - t) / t * 1e6)
                        hit_idx |= set(np.where(d < 20.0)[0].tolist())
            picked = list(hit_idx)
            for k in np.argsort(-sint):
                if len(picked) >= 400 + len(hit_idx):
                    break
                if int(k) not in hit_idx:
                    picked.append(int(k))
            sel = np.array(sorted(picked), dtype=int)
            omz, oint = smz[sel], sint[sel]
            top = oint.max() if len(oint) and oint.max() > 0 else 1.0
            oin = oint / top
            xmax = float(omz.max()) * 1.03 if len(omz) else 100.0

            fig = plt.figure(figsize=(8.27, 11.69))
            gs = gridspec.GridSpec(5, 1, height_ratios=[0.34, 1.05, 0.82, 0.34, 1.55],
                                   hspace=0.34, left=0.085, right=0.965,
                                   top=0.918, bottom=0.055)
            # 1 信息行（outer-title 风格）
            fig.text(0.085, 0.982, "| %s (%s) | (%s)%s(%s) | %d+ | %s |" %
                     (m0.prot, m0.gene, pre, seq, post, int(b.z), mods_txt),
                     fontsize=10.0, weight="bold", va="top", family="Arial")
            fig.text(0.085, 0.958, "Protein position %s   |   Precursor m/z %.4f   |   RT range "
                     "(%d samples): %.2f-%.2f min   |   main-text table peptide: %s" %
                     (("[%d]" % pos) if pos else "-", float(b.mz_s14), len(gk), rt_lo, rt_hi,
                      "yes" if m0.in_main else "no"),
                     fontsize=7.0, va="top", family="Arial", color="#222222")
            fig.text(0.085, 0.941, "Best spectrum: %s  |  RT %.2f min (S14 %.2f, Δ%+.2f)  |  matched b/y "
                     "%d (ladder b%d/y%d)  |  MS1 precursor envelope: %s  |  S14 q-value: %s" %
                     (b["sample"], float(b.apex_rt), float(b.rt_s14), float(b.drt), int(b.nmatched),
                      int(b.bLad), int(b.yLad), "yes" if b.ms1_flag else "not detected",
                      ("%.1e" % float(b.qval)) if pd.notna(b.qval) else "-"),
                     fontsize=7.0, va="top", family="Arial", color="#222222")

            # 2 离子梯
            iso_ladder_ax(fig, gs[0], seq, hits, mod_positions(pid), int(b.z), xmax)

            # 3 原始谱图（hits 强度须与背景同为归一化标度）
            hits_norm = {l: (v[0], v[1] / top if top > 0 else 0.0, v[2]) for l, v in hits.items()}
            iso_spectrum_ax(
                fig, gs[1], omz, oin, hits_norm, xmax,
                ann_lines=["Base peak = %.2E" % float(sint.max() if len(sint) else float("nan")),
                           "Precursor m/z = %.4f  (%d+)" % (float(b.mz_s14), int(b.z))],
                head="Observed spectrum - %s (scan %d)" % (b["sample"], int(b.apex_scan)))

            # 4 预测谱图（倒置）
            pv = pred.get(("hebei", pid))
            nres = len(seq)
            cos_txt = ""
            if pv is not None and len(np.asarray(pv["b"])) == nres - 1:
                xs_i = np.arange(1, nres)
                bm = np.array([frags["b%d" % i] for i in xs_i])
                ym = np.array([frags["y%d" % i] for i in xs_i])
                bi = np.asarray(pv["b"], dtype=float)
                yi = np.asarray(pv["y"], dtype=float)
                pmax = max(bi.max(), yi.max(), 1e-9)
                pmz = np.concatenate([bm, ym])
                pin = -np.concatenate([bi, yi]) / pmax
                phits = {}
                for m, v, lbl in zip(list(bm) + list(ym), list(pin),
                                     ["b%d" % i for i in xs_i] + ["y%d" % i for i in xs_i]):
                    phits[lbl] = (float(m), float(v), 0.0)
                obs1 = {l.split("(")[0]: v[1] * top for l, v in hits.items() if "(2+)" not in l}
                cp = ([(obs1["b%d" % i], bi[i - 1]) for i in xs_i if ("b%d" % i) in obs1] +
                      [(obs1["y%d" % i], yi[i - 1]) for i in xs_i if ("y%d" % i) in obs1])
                if len(cp) >= 4:
                    a1 = np.array([x[0] for x in cp]); a2 = np.array([x[1] for x in cp])
                    n1, n2 = float(np.linalg.norm(a1)), float(np.linalg.norm(a2))
                    cos_txt = ("cosine(obs, pred) = %.3f  (%d shared 1+ ions)" %
                               (float(np.dot(a1, a2) / (n1 * n2)), len(cp))) if n1 > 0 and n2 > 0 else ""
                # 倒置面板专用绘制（预测谱无背景灰峰）
                axp = fig.add_subplot(gs[2])
                axp.set_xlim(0, xmax)
                axp.set_ylim(-1.15, 0.0)
                axp.set_yticks([0.0, -0.5, -1.0])
                axp.set_yticklabels(["0", "50", "100"], fontsize=5.2)
                axp.set_ylabel("Relative intensity (%)", fontsize=6.2)
                axp.set_xticklabels([])
                axp.spines[["top", "right"]].set_visible(False)
                axp.tick_params(length=2, pad=1.5)
                for m, v, lbl in zip(pmz, pin, list(phits.keys())):
                    c = C_B if lbl.startswith("b") else C_Y
                    axp.vlines(m, 0, v, color=c, lw=0.8)
                    if -v > 0.10:
                        axp.text(m, v - 0.03, lbl, rotation=90, fontsize=4.6,
                                 ha="center", va="top", color=c)
                axp.text(0.008, 1.04, "AlphaPeptDeep prediction (inverted)", transform=axp.transAxes,
                         fontsize=6.2, ha="left", va="bottom", family="Arial", weight="bold")
                if cos_txt:
                    axp.text(0.992, 1.04, cos_txt, transform=axp.transAxes, fontsize=5.8,
                             ha="right", va="bottom", family="Arial")
            else:
                axp = fig.add_subplot(gs[2])
                axp.set_xlim(0, xmax); axp.set_ylim(-1.15, 0.0)
                axp.axis("off")
                axp.text(0.5, -0.5, "AlphaPeptDeep prediction unavailable for this precursor "
                         "(modified/synthetic variant not present in the prediction run)",
                         ha="center", fontsize=7, family="Arial")

            # 5 质量差分布
            axppm = iso_ppm_ax(fig, gs[3], hits, xmax)
            axppm.set_title("Mass error of matched ions (observed spectrum)", fontsize=6.0,
                            loc="left", pad=2)

            # 6 热图
            axh = fig.add_subplot(gs[4])
            samples = sorted(gk["sample"].unique(),
                             key=lambda s: (0 if s.startswith("HC") else 1, int(s[2:])))
            prof_rows = []
            for s in samples:
                i = prof_map.get((seq, int(b.z), float(b.mz_s14), s))
                if i is None:
                    continue
                r0, r1 = int(pz["off"][i]), int(pz["end"][i])
                prof_rows.append((s, pz["rt"][r0:r1], pz["val"][r0:r1]))
            if prof_rows:
                b0 = math.floor(min(min(rt) for _, rt, _ in prof_rows))
                b1 = math.ceil(max(max(rt) for _, rt, _ in prof_rows)) + 1
                nb = b1 - b0
                M = np.zeros((len(prof_rows), nb))
                for ri, (s, rt, val) in enumerate(prof_rows):
                    np.maximum.at(M[ri], np.clip((rt - b0).astype(int), 0, nb - 1), val)
                vmax = M.max()
                if vmax > 0:
                    M = M / vmax
                im = axh.imshow(M, aspect="auto", cmap=HEAT, vmin=0, vmax=1,
                                extent=[b0, b1, len(prof_rows) - 0.5, -0.5],
                                interpolation="nearest")
                axh.set_yticks(range(len(prof_rows)))
                axh.set_yticklabels([x[0] for x in prof_rows], fontsize=4.8)
                axh.set_xlabel("retention time (min)", fontsize=6.5)
                axh.set_xlim(b0, b1)
                axh.tick_params(labelsize=5.2, length=2, pad=1.5)
                axh.set_title("Fragment-ion XIC heatmap (1-min bins; white = 0, red = max)",
                              fontsize=6.5, loc="left", pad=2)
                cb = fig.colorbar(im, ax=axh, pad=0.008, fraction=0.035)
                cb.ax.tick_params(labelsize=4.5)
                cb.set_label("rel. XIC", fontsize=5)
            else:
                axh.text(0.5, 0.5, "no profiles", ha="center", va="center",
                         transform=axh.transAxes)

            pages.append(dict(page=len(pages) + 2, prot=m0.prot, gene=m0.gene, seq=seq,
                              z=int(b.z), mz=float(b.mz_s14), flank="%s.%s.%s" % (pre, seq, post),
                              samples=len(gk), rt_lo=rt_lo, rt_hi=rt_hi,
                              best=b["sample"], nmatched=int(b.nmatched),
                              ms1=int(b.ms1_flag), in_main=bool(m0.in_main),
                              pred_ok=bool(pv is not None and
                                           len(np.asarray(pv["b"])) == len(seq) - 1)))
            pdf.savefig(fig)
            if only:
                fig.savefig(ROOT + r"\work\output\s14_evidence\S14_evidence_sample_v3_%s.png" % only[:24],
                            dpi=300)
            plt.close(fig)

        if not only:
            rowsix = [["Page", "Protein", "Gene", "Peptide(flanks)", "z", "m/z", "N",
                       "RT range", "Best PSM"]]
            for p in pages:
                rowsix.append([p["page"], p["prot"], p["gene"], p["flank"], p["z"],
                               "%.4f" % p["mz"], p["samples"],
                               "%.2f-%.2f" % (p["rt_lo"], p["rt_hi"]),
                               "%s(%d)" % (p["best"], p["nmatched"])])
            chunk = 52
            npg = math.ceil((len(rowsix) - 1) / chunk)
            for pi in range(npg):
                fig = plt.figure(figsize=(8.27, 11.69))
                fig.text(0.06, 0.965, "Index%s" % ("" if pi == 0 else " (cont'd %d/%d)" % (pi + 1, npg)),
                         fontsize=13, weight="bold", va="top", family="Arial")
                ax = fig.add_axes([0.04, 0.03, 0.94, 0.92]); ax.axis("off")
                sub = [rowsix[0]] + rowsix[1 + pi * chunk:1 + (pi + 1) * chunk]
                sub = [[(c[:30] + "…") if isinstance(c, str) and len(c) > 31 else c
                        for c in row] for row in sub]
                tb = ax.table(cellText=sub[1:], colLabels=sub[0], loc="upper center",
                              cellLoc="center", bbox=[0, 0, 1, 0.94],
                              colWidths=[0.06, 0.10, 0.11, 0.33, 0.04, 0.09, 0.04, 0.11, 0.12])
                tb.auto_set_font_size(False); tb.set_fontsize(4.5)
                for (r_i, c_i), cell in tb.get_celld().items():
                    cell.set_linewidth(0.2)
                    if r_i == 0:
                        cell.set_facecolor("#e8e8e8")
                pdf.savefig(fig); plt.close(fig)

    pd.DataFrame(pages).to_csv(OUTMANI, sep="\t", index=False, encoding="utf-8")
    summ = dict(pages=len(pages), sequences=len(seqs),
                pred_available=int(sum(1 for p in pages if p["pred_ok"])))
    with _open(OUTSUMM, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print(json.dumps(summ, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
