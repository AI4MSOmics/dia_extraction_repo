# -*- coding: utf-8 -*-
"""
74_make_a4_v2.py —— S14 证据页 v2（按用户三条修改意见重排）

每条肽段序列一页（A4 纵向，Arial）：
  头部   一行或多行文字：蛋白(Rv+基因+位点) | 肽段(侧翼氨基酸 R.PEPTIDE.K，isoPSV 格式)
         | 价态 | 修饰 | RT 范围(该前体全部样本的实测 apex RT min–max)
  离子梯 isoPSV 风格：序列字母条，b 位点下方绿标 bₖ、y 位点上方红标 yₙ，切割位点竖线，
         修饰残基橙色加粗，左侧红色 z+
  实际谱 PSM 最优样本的 apex 谱：灰背景峰 + 匹配 b(绿)/y(红) 峰带旋转标签；
         右上注释 Base peak / Precursor m/z；下方 ±20 ppm 子面板
  理论谱 AlphaPeptDeep 预测（镜像向下，共用 m/z 轴），注释 obs-vs-pred 余弦（1+ 共有离子）
  热图   横轴 RT（1 min 分箱）、纵轴 S14 有定量的样本（HC 在上 TB 在下），
         色 = 该样本该 RT 箱的碎片匹配强度和（白=0 → 红=高）；每行右侧三角标 apex RT
输出：work/output/s14_evidence/S14_peptide_evidence_A4_v2.pdf + pages_manifest_v2.tsv
复用 71 的理论碎片函数（importlib）。安全：路径常量；open 经 _open() 白名单。
"""
import sys, os, json, math, pickle, re, importlib.util
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
OUTPDF = ROOT + r"\work\output\s14_evidence\S14_peptide_evidence_A4_v2.pdf"
OUTMANI = ROOT + r"\work\output\s14_evidence\pages_manifest_v2.tsv"
OUTSUMM = ROOT + r"\work\output\s14_evidence\pages_summary_v2.json"

PROTON = 1.007276466621
H2O = 18.01056468403
AA = dict(G=57.02146, C=103.00919, S=87.03203, P=97.05276, A=71.03711,
          V=99.06841, T=101.04768, L=113.08406, I=113.08406, N=114.04293, D=115.02694,
          Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
          R=156.10111, Y=163.06333, W=186.07931)

_spec = importlib.util.spec_from_file_location("ext71", ROOT + r"\work\scripts\71_extract_targeted.py")
m71 = importlib.util.module_from_spec(_spec)
sys.modules["ext71"] = m71
_spec.loader.exec_module(m71)

# isoPSV 配色
C_B, C_Y, C_MOD = "#008000", "#CA1D52", "#FF8C00"
plt.rcParams.update({
    "font.size": 7, "pdf.fonttype": 42, "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
})
HEAT = LinearSegmentedColormap.from_list("wr", ["#FFFFFF", "#FDD8B5", "#F68838", "#D7301F"])


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
    pre = p[i - 1] if i > 0 else "-"
    post = p[i + len(seq)] if i + len(seq) < len(p) else "-"
    return pre, post, i + 1


def mod_positions(pid):
    """pid 修饰位点（1-based，N 端乙酰=0）；synthetic 全 M 氧化。"""
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
            pos.append(len(s2))  # N 端修饰时 len=0
            i = j + 1
        else:
            s2 += body[i]
            i += 1
    return [p + (1 if p > 0 or True else 0) if False else p for p in pos]


def draw_ladder(ax, seq, hits, mods_pos, z, xpad):
    """isoPSV 风格离子梯（画在 ax 顶部 1.18–1.56 区域内，全部可见，x 轴与谱图共用）。"""
    n = len(seq)
    xs = np.linspace(xpad[0] + (xpad[1] - xpad[0]) * 0.06,
                     xpad[1] - (xpad[1] - xpad[0]) * 0.06, n)
    slot = xs[1] - xs[0] if n > 1 else 1
    ymid = 1.36
    b_sites = sorted({int(l[1:]) for l in hits if l.startswith("b") and l[1:].isdigit()})
    y_sites = sorted({n - int(l[1:]) for l in hits if l.startswith("y") and l[1:].isdigit()})
    for k in sorted(set(b_sites) | set(y_sites)):
        if k < len(xs):
            ax.vlines(xs[k] + slot / 2, ymid - 0.04, ymid + 0.04, color="k", lw=1.0)
    for k in b_sites:
        if k < len(xs):
            ax.hlines(ymid - 0.04, xs[k] + slot / 2 - slot * 0.45, xs[k] + slot / 2, color="k", lw=1.0)
            ax.text(xs[k], ymid - 0.07, r"$b_{%d}$" % k, fontsize=5.0, color=C_B,
                    ha="center", va="top")
    for k in y_sites:
        if k < len(xs):
            ax.hlines(ymid + 0.04, xs[k] + slot / 2, xs[k] + slot / 2 + slot * 0.45, color="k", lw=1.0)
            ax.text(xs[k + 1] if k + 1 < len(xs) else xs[-1], ymid + 0.07,
                    r"$y_{%d}$" % (n - k), fontsize=5.0, color=C_Y, ha="center", va="bottom")
    ax.text(xs[0] - slot * 1.4, ymid, "%d+" % z, fontsize=7, color=C_Y,
            ha="right", va="center", weight="bold")
    mp = set(mods_pos)
    for i, (x, a) in enumerate(zip(xs, seq)):
        col = C_MOD if (i + 1) in mp and (i + 1) != 0 else "k"
        ax.text(x, ymid, a, fontsize=5.6, ha="center", va="center",
                color=col, weight="bold" if (i + 1) in mp else "normal")


def main():
    only = ""
    for a in sys.argv[1:]:
        if a.startswith("--only="):
            only = a.split("=", 1)[1]
    outpdf = OUTPDF
    if only:
        outpdf = ROOT + r"\work\output\s14_evidence\S14_evidence_sample_%s.pdf" % only[:24]
    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    z = np.load(SPECNPZ)
    mz_all, int_all, off, end = z["mz"], z["inten"], z["off"], z["end"]
    tg = pd.read_parquet(TARGETS)
    prots = load_fasta_flanks()
    pred = pickle.load(_open(PRED, "rb"))

    pidx = pd.read_parquet(PROF_IDX)
    pz = np.load(PROF_NPZ)
    prof_map = {(r.seq, r.z, r.mz, r["sample"]): i for i, r in pidx.reset_index().iterrows()}

    ok = rows[rows.status == "ok"]
    pages, records = [], []
    seqs = list(ok.seq.unique())
    if only:
        seqs = [s for s in seqs if s == only]
        assert seqs, "sequence not found: " + only
    print("sequences:", len(seqs))

    with PdfPages(outpdf) as pdf:
        if not only:
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.5, 0.60, "Evidence for Mtb-derived peptides (Table S14) — v2",
                     ha="center", fontsize=15, weight="bold", family="Arial")
            fig.text(0.5, 0.555, "One peptide per page: sequence with flanking residues, ion ladder, "
                     "best-PSM observed spectrum vs AlphaPeptDeep prediction, RT heatmap across samples",
                     ha="center", fontsize=9)
            stats = [f"{len(seqs)} peptide sequences ({tg.drop_duplicates(['prot']).prot.nunique()} proteins), "
                     f"{rows['sample'].nunique()} samples (7 HC / 23 TB), hebei DIA cohort",
                     "Observed spectrum = RT-anchored targeted extraction, best sample by matched b/y ions",
                     "Prediction = AlphaPeptDeep 1.5.0 (QE, generic model); cosine on shared 1+ b/y ions",
                     "Heatmap = fragment-ion matched intensity per 1-min RT bin per sample (white = 0, red = max)"]
            for i, s in enumerate(stats):
                fig.text(0.13, 0.49 - i * 0.028, s, fontsize=9)
            pdf.savefig(fig); plt.close(fig)

        for seq in seqs:
            g = ok[ok.seq == seq]
            b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
            key3 = (seq, int(b.z), float(b.mz_s14))
            gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
            m0 = tg[(tg.seq == seq) & (tg.z == b.z) & (tg.mz_s14 == b.mz_s14)].iloc[0]
            pid = m0.pid
            frags, _, _ = m71.parse_pid(pid)
            pre, post, pos = flank_of(prots, m0.prot, seq)
            rts_all = gk.apex_rt.astype(float)
            rt_lo, rt_hi = float(rts_all.min()), float(rts_all.max())
            best_sl = b["sample"]

            fig = plt.figure(figsize=(8.27, 11.69))
            gs = gridspec.GridSpec(3, 1, height_ratios=[3.3, 0.95, 5.4], hspace=0.38,
                                   left=0.075, right=0.965, top=0.925, bottom=0.05)
            # 头部两行
            fig.text(0.075, 0.975, "%s (%s)   position %s   |   %s.%s.%s   |   %d+   |   mods: %s" %
                     (m0.prot, m0.gene, ("[%d]" % pos) if pos else "-",
                      pre, seq, post, int(b.z),
                      ("Carbamidomethyl (C)" if seq.count("C") else "none") +
                      ("; Oxidation (M)" if "ox" in str(m0.mod_state) else "") +
                      ("; Acetyl (N-term)" if "acetyl" in str(m0.mod_state) else "") +
                      ("; " + str(m0.mod_state) if str(m0.mod_state).startswith("unknown") else "")),
                     fontsize=9.0, weight="bold", va="top", family="Arial")
            fig.text(0.075, 0.951, "m/z %.4f   |   RT range (all %d samples): %.2f–%.2f min   |   "
                     "best spectrum: %s (matched b/y %d, ladder b%d/y%d, MS1 %s)   |   in main-text table: %s" %
                     (float(b.mz_s14), len(gk), rt_lo, rt_hi, best_sl, int(b.nmatched),
                      int(b.bLad), int(b.yLad), "yes" if b.ms1_flag else "n/a",
                      "yes" if m0.in_main else "no"),
                     fontsize=7.2, va="top", color="#333333", family="Arial")

            # ---- 实际谱（含离子梯 + ppm）----
            gsA = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0],
                                                   height_ratios=[3.2, 1.0], hspace=0.06)
            ax = fig.add_subplot(gsA[0])
            sl = slice(int(off[b.name]), int(end[b.name]))
            smz, sint = mz_all[sl], int_all[sl]
            # 先在全部峰上找匹配离子（弱匹配不能被 top-400 挤掉），再补 top-400 背景峰
            hits = {}
            hit_idx = set()
            if len(smz):
                for l, t in frags.items():
                    d = np.abs(smz - t) / t * 1e6
                    k = int(np.argmin(d))
                    if d[k] < 20.0:
                        hits[l] = (float(smz[k]), float(sint[k]), (smz[k] - t) / t * 1e6)
                        kk = np.where(d < 20.0)[0]
                        hit_idx |= set(kk.tolist())
            rest = np.argsort(-sint)
            picked = list(hit_idx)
            for k in rest:
                if len(picked) >= 400 + len(hit_idx):
                    break
                if int(k) not in hit_idx:
                    picked.append(int(k))
            sel = np.array(sorted(picked), dtype=int)
            omz, oint = smz[sel], sint[sel]
            top = oint.max() if len(oint) and oint.max() > 0 else 1.0
            ax.set_ylim(0, 1.62)
            for m, i in zip(omz, oint):
                ax.vlines(m, 0, i / top, color="#b8b8b8", lw=0.4)
            for l, (m, i, _) in sorted(hits.items()):
                c = C_B if l.startswith("b") else C_Y
                ax.vlines(m, 0, max(i / top, 0.012), color=c, lw=0.9)
                ax.text(m, max(i / top, 0.012) + 0.015, l.split("(")[0], rotation=90,
                        fontsize=3.8, ha="center", va="bottom", color=c)
            draw_ladder(ax, seq, hits, mod_positions(pid), int(b.z), (float(omz.min()), float(omz.max())))
            bp = float(sint.max()) if len(sint) else float("nan")
            ax.text(0.995, 1.56, "Base peak = %.2E" % bp, transform=ax.transAxes,
                    fontsize=6, ha="right", va="top", family="Arial")
            ax.text(0.995, 1.50, "Precursor m/z: %.4f   charge = %d+" % (float(b.mz_s14), int(b.z)),
                    transform=ax.transAxes, fontsize=6, ha="right", va="top", family="Arial")
            ax.set_xlim(min(float(omz.min()), 100) * 0.98 if len(omz) else 0,
                        float(omz.max()) * 1.02 if len(omz) else 100)
            ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=5)
            ax.set_ylabel("Relative intensity (%)", fontsize=6)
            ax.set_xticklabels([])
            ax.spines[["top", "right"]].set_visible(False)
            # ppm 面板
            axp = fig.add_subplot(gsA[1], sharex=ax)
            axp.axhline(0, color="#999999", lw=0.6)
            for hh in (-20, 20):
                axp.axhline(hh, color="#bbbbbb", ls=":", lw=0.5)
            for l, (m, i, ppm) in hits.items():
                axp.plot(m, ppm, "o", ms=1.4,
                         color=C_B if l.startswith("b") else C_Y)
            axp.set_ylim(-25, 25)
            axp.set_yticks([-20, 0, 20])
            axp.set_ylabel("ppm", fontsize=6)
            axp.set_xlabel("m/z", fontsize=6.5)
            axp.tick_params(labelsize=5)
            axp.spines[["top", "right"]].set_visible(False)
            ax_obs = ax

            # ---- AlphaPeptDeep 理论谱（镜像向下）----
            ax2 = fig.add_subplot(gs[1], sharex=ax)
            pv = pred.get(("hebei", pid))
            if pv is not None:
                nres = len(seq)
                xs = np.arange(1, nres)
                bm = np.array([frags["b%d" % i] for i in xs])
                ym = np.array([frags["y%d" % i] for i in xs])
                bi = np.asarray(pv["b"], dtype=float)
                yi = np.asarray(pv["y"], dtype=float)
                if len(bi) == nres - 1 and len(yi) == nres - 1:
                    bmax = max(bi.max(), yi.max(), 1e-9)
                    for m, v in zip(bm, bi):
                        ax2.vlines(m, 0, -v / bmax, color=C_B, lw=0.9)
                    for m, v in zip(ym, yi):
                        ax2.vlines(m, 0, -v / bmax, color=C_Y, lw=0.9)
                    for m, v, i in zip(bm, bi, xs):
                        if v / bmax > 0.10:
                            ax2.text(m, -v / bmax - 0.02, "b%d" % i, rotation=90, fontsize=3.6,
                                     ha="center", va="top", color=C_B)
                    for m, v, i in zip(ym, yi, xs):
                        if v / bmax > 0.10:
                            ax2.text(m, -v / bmax - 0.02, "y%d" % i, rotation=90,
                                     fontsize=3.6, ha="center", va="top", color=C_Y)
                    ax2.set_ylim(-1.12, 0.06)
                    # 余弦（1+ 共有离子）
                    obs1 = {l.split("(")[0]: v[1] for l, v in hits.items() if "(2+)" not in l}
                    cpairs = [(obs1["b%d" % i], bi[i - 1]) for i in xs if ("b%d" % i) in obs1] + \
                             [(obs1["y%d" % i], yi[i - 1]) for i in xs if ("y%d" % i) in obs1]
                    if len(cpairs) >= 4:
                        a1 = np.array([x[0] for x in cpairs]); a2v = np.array([x[1] for x in cpairs])
                        n1, n2 = float(np.linalg.norm(a1)), float(np.linalg.norm(a2v))
                        cosv = float(np.dot(a1, a2v) / (n1 * n2)) if n1 > 0 and n2 > 0 else float("nan")
                        ax2.text(0.995, 0.97, "AlphaPeptDeep prediction   |   cosine(obs, pred) = %.3f  (%d shared 1+ ions)" %
                                 (cosv, len(cpairs)), transform=ax2.transAxes, fontsize=6,
                                 ha="right", va="top", family="Arial",
                                 bbox=dict(fc="white", ec="none", pad=1.2))
                    else:
                        ax2.text(0.995, 0.97, "AlphaPeptDeep prediction", transform=ax2.transAxes,
                                 fontsize=6, ha="right", va="top", family="Arial",
                                 bbox=dict(fc="white", ec="none", pad=1.2))
                else:
                    ax2.text(0.5, -0.5, "prediction array length mismatch", ha="center",
                             fontsize=7, transform=ax2.transAxes)
            else:
                ax2.text(0.5, -0.5, "AlphaPeptDeep prediction unavailable for this precursor "
                         "(modified/synthetic variant not present in the prediction run)",
                         ha="center", fontsize=7, transform=ax2.transAxes, family="Arial")
            ax2.set_yticks([0, -0.5, -1.0])
            ax2.set_yticklabels(["0", "50", "100"], fontsize=5)
            ax2.set_ylabel("Predicted (%)", fontsize=6)
            ax2.tick_params(labelsize=5)
            ax2.spines[["top", "right"]].set_visible(False)

            # ---- RT 热图 ----
            axh = fig.add_subplot(gs[2])
            samples = sorted(gk["sample"].unique(), key=lambda s: (0 if s.startswith("HC") else 1, int(s[2:])))
            prof_rows = []
            for s in samples:
                i = prof_map.get((seq, int(b.z), float(b.mz_s14), s))
                if i is None:
                    continue
                r0, r1 = int(pz["off"][i]), int(pz["end"][i])
                rt = pz["rt"][r0:r1]; val = pz["val"][r0:r1]
                prof_rows.append((s, rt, val))
            if prof_rows:
                b0 = math.floor(min(min(rt) for _, rt, _ in prof_rows))
                b1 = math.ceil(max(max(rt) for _, rt, _ in prof_rows)) + 1
                nb = b1 - b0
                M = np.zeros((len(prof_rows), nb))
                for ri, (s, rt, val) in enumerate(prof_rows):
                    bi_ = np.clip(((rt - b0)).astype(int), 0, nb - 1)
                    np.maximum.at(M[ri], bi_, val)
                vmax = M.max()
                if vmax > 0:
                    M = M / vmax
                im = axh.imshow(M, aspect="auto", cmap=HEAT, vmin=0, vmax=1,
                                extent=[b0, b1, len(prof_rows) - 0.5, -0.5], interpolation="nearest")
                axh.set_yticks(range(len(prof_rows)))
                axh.set_yticklabels([x[0] for x in prof_rows], fontsize=4.6)
                axh.set_xlabel("retention time (min)", fontsize=6.5)
                axh.set_xlim(b0, b1)
                axh.tick_params(labelsize=5)
                axh.set_title("Fragment-ion XIC heatmap (1-min bins; white = 0, red = max)",
                              fontsize=6.5, loc="left")
                cb = fig.colorbar(im, ax=axh, pad=0.008, fraction=0.035)
                cb.ax.tick_params(labelsize=4.5)
                cb.set_label("rel. XIC", fontsize=5)
            else:
                axh.text(0.5, 0.5, "no profiles", ha="center", va="center", transform=axh.transAxes)

            pages.append(dict(page=len(pages) + 2, prot=m0.prot, gene=m0.gene, seq=seq,
                              z=int(b.z), mz=float(b.mz_s14), flank="%s.%s.%s" % (pre, seq, post),
                              samples=len(gk), rt_lo=rt_lo, rt_hi=rt_hi,
                              best=best_sl, nmatched=int(b.nmatched),
                              ms1=int(b.ms1_flag), in_main=bool(m0.in_main),
                              pred_ok=bool(pv is not None and len(np.asarray(pv["b"])) == len(seq) - 1)))
            pdf.savefig(fig)
            if only:
                fig.savefig(ROOT + r"\work\output\s14_evidence\S14_evidence_sample_%s.png" % only[:24],
                            dpi=300)
            plt.close(fig)

        # 索引
        rowsix = [["Page", "Protein", "Gene", "Peptide(flanks)", "z", "m/z", "N", "RT range", "Best PSM"]]
        for p in (pages if not only else []):
            rowsix.append([p["page"], p["prot"], p["gene"], p["flank"], p["z"],
                           "%.4f" % p["mz"], p["samples"], "%.2f-%.2f" % (p["rt_lo"], p["rt_hi"]),
                           "%s(%d)" % (p["best"], p["nmatched"])])
        chunk = 52
        np_ = math.ceil((len(rowsix) - 1) / chunk)
        for pi in range(np_):
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.06, 0.965, "Index%s" % ("" if pi == 0 else " (cont'd %d/%d)" % (pi + 1, np_)),
                     fontsize=13, weight="bold", va="top", family="Arial")
            ax = fig.add_axes([0.04, 0.03, 0.94, 0.92]); ax.axis("off")
            sub = [rowsix[0]] + rowsix[1 + pi * chunk:1 + (pi + 1) * chunk]
            sub = [[(c[:30] + "…") if isinstance(c, str) and len(c) > 31 else c for c in row] for row in sub]
            tb = ax.table(cellText=sub[1:], colLabels=sub[0], loc="upper center",
                          cellLoc="center", bbox=[0, 0, 1, 0.94],
                          colWidths=[0.06, 0.10, 0.11, 0.33, 0.04, 0.09, 0.04, 0.11, 0.12])
            tb.auto_set_font_size(False); tb.set_fontsize(4.5)
            for (r_i, c_i), cell in tb.get_celld().items():
                cell.set_linewidth(0.2)
                if r_i == 0:
                    cell.set_facecolor("#e8e8e8")
            pdf.savefig(fig); plt.close(fig)

    mani = pd.DataFrame(pages)
    mani.to_csv(OUTMANI, sep="\t", index=False, encoding="utf-8")
    summ = dict(pages=len(pages), sequences=len(seqs),
                pred_available=int(sum(1 for p in pages if p["pred_ok"])))
    with _open(OUTSUMM, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print(json.dumps(summ, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
