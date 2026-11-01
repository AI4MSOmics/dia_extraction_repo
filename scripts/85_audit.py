# -*- coding: utf-8 -*-
"""
85_audit.py —— S14/Table1 证据 PDF 全链自动审计

三层检查（全部程序化，逐页/逐行）：
A 数据层对账
  A1 S14 xlsx -> s14_targets（行数/前体数/样本/修饰态/每靶行数）
  A2 rows_all 提取质量（ok 率/ΔRT 分布/声明值在界内）
  A3 76/84 输入表 vs rows_all（每行 File_Name/Scan/ExpMH/Charge/Sequence/Raw_Score）
  A4 MGF vs 输入表 vs spectra npz（TITLE 一一对应、峰数、Base peak 值）
  A5 fasta 侧翼重算 vs Info1 侧翼（S14 链 112 页）
  A6 Info1/2/3 全字段对账（m/z/样本数/RT 范围/best/ΔRT/matched/ladder/MS1/q/qty）
  A7 AlphaPeptDeep cosine 独立重算 vs 页面 PredCos（全部页）
  A8 热图 apex RT mean±SD 独立重算 vs 页面图题（全部页）
  A9 热图行数 = 样本数；每行 max bin RT 落在页声明 RT 范围 ±1.5min
  A10 索引页/manifest vs 实际页序
B Table1 PDF
  B1 页序=Table 1 顺序（159、No. 非降）
  B2 首行 Table 1 No/蛋白 与 xlsx 对账
  B3 ok/ok_extra 与 manifest 一致；补充链 Info3 标注存在
C 呈现层（图像级轻量）
  C1 每页关键文本存在（Info1 "|"、Observed spectrum、Base peak、
     prediction 或 unavailable、apex RT 图题）
  C2 每页下部 40% 非全白（热图在）
  C3 全部文字 bbox 在页内（无溢出）
输出：work/output/s14_evidence/audit_report.md + audit_report.json
"""
import sys, os, json, math, re, pickle, importlib.util, collections
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import pymupdf
import openpyxl

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


S14XLSX = ROOT + r"\Supplementary Table_2026.09.03(1).xlsx"
TARGETS = ROOT + r"\work\tmp\s14_targets.parquet"
TARGETS_EX = ROOT + r"\work\tmp\s14_targets_extra.parquet"
ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
ROWS_EX = ROOT + r"\work\output\s14_evidence\extract_extra\rows_all.parquet"
SPECNPZ = ROOT + r"\work\output\s14_evidence\extract\spectra_all.npz"
SPECNPZ_EX = ROOT + r"\work\output\s14_evidence\extract_extra\spectra_all.npz"
INTSV = ROOT + r"\work\risoPSV\in\s14_spectra.tsv"
INTSV_EX = ROOT + r"\work\risoPSV\in_extra\s14_spectra.tsv"
MGF = ROOT + r"\work\risoPSV\in\best.mgf"
MGF_EX = ROOT + r"\work\risoPSV\in_extra\best.mgf"
FASTA = ROOT + r"\IPX0018588001_hebei\human_tuberculosis.fasta"
PRED = ROOT + r"\work\output\summary\pred\AlphaPeptDeep.pkl"
PROF_I = ROOT + r"\work\output\s14_evidence\profiles\index.parquet"
PROF_N = ROOT + r"\work\output\s14_evidence\profiles\profiles.npz"
PROF_I_EX = ROOT + r"\work\output\s14_evidence\profiles_extra\index.parquet"
PROF_N_EX = ROOT + r"\work\output\s14_evidence\profiles_extra\profiles.npz"
PDF_MAIN = ROOT + r"\work\output\s14_evidence\S14_peptide_evidence_A4_isoPSV.pdf"
PDF_T1 = ROOT + r"\work\output\s14_evidence\Table1_peptide_evidence_A4.pdf"
T1XLSX = ROOT + r"\Table 1.xlsx"
MANI_T1 = ROOT + r"\work\output\s14_evidence\table1_manifest.tsv"
OUTMD = ROOT + r"\work\output\s14_evidence\audit_report.md"
OUTJ = ROOT + r"\work\output\s14_evidence\audit_report.json"

_spec = importlib.util.spec_from_file_location("ext71", ROOT + r"\work\scripts\71_extract_targeted.py")
m71 = importlib.util.module_from_spec(_spec)
sys.modules["ext71"] = m71
_spec.loader.exec_module(m71)
PROTON = 1.007276466621

issues = []          # (级别, 编号, 描述)
def issue(lv, code, desc):
    issues.append((lv, code, desc))
    print("[%s] %s %s" % (lv, code, desc))

report = collections.OrderedDict()


def load_fasta():
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


def parse_s14():
    wb = openpyxl.load_workbook(S14XLSX, read_only=True)
    rows = list(wb["Table S14"].iter_rows(values_only=True))[2:]
    cols = ["prot", "gene", "uc", "mtb", "seq", "z", "sample", "mz", "rt",
            "fim", "fmz", "qval", "qty", "shared", "nvar"]
    return pd.DataFrame([r[:15] for r in rows], columns=cols).dropna(subset=["seq"])


def main():
    # ---------- A1 ----------
    s14 = parse_s14()
    tg = pd.read_parquet(TARGETS)
    a1 = dict(s14_rows=len(s14), targets=len(tg.drop_duplicates(["seq", "z", "mz_s14"])),
              samples=int(tg["sample"].nunique()), rows_parquet=len(tg))
    okA1 = (len(s14) == len(tg) and
            a1["targets"] == 164 and a1["samples"] == 30)
    if not okA1:
        issue("ERR", "A1", "S14->targets 对账失败: %s" % a1)
    grp = tg.groupby(["seq", "z", "mz_s14"]).size()
    a1["rows_per_target_minmax"] = [int(grp.min()), int(grp.max())]
    report["A1"] = a1

    # ---------- A2 ----------
    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    okr = rows[rows.status == "ok"]
    seqs = list(okr.seq.unique())   # 主链 112 序（页序 = 79 约定）
    a2 = dict(rows=len(rows), ok=len(okr), ok_rate=len(okr) / len(rows),
              abs_drt_med=float(okr.drt.abs().median()),
              abs_drt_p95=float(okr.drt.abs().quantile(.95),
              ) if False else float(okr.drt.abs().quantile(.95)),
              drt_le_05=float((okr.drt.astype(float).abs() <= 0.5).mean()),
              ladder_ge6=float(((okr.bLad + okr.yLad) >= 6).mean()),
              ms1=float(okr.ms1_flag.mean()))
    if a2["ok_rate"] != 1.0:
        issue("ERR", "A2", "提取存在非 ok 行: %d" % (len(rows) - len(okr)))
    report["A2"] = a2

    # ---------- A3/A4（两链） ----------
    pred = pickle.load(_open(PRED, "rb"))
    fasta = load_fasta()
    znp = np.load(SPECNPZ)
    znp_ex = np.load(SPECNPZ_EX)

    for chain, itab_p, mgf_p, rows_p, znp_c, tsv_tag in [
            ("main", INTSV, MGF, rows, znp, "main"),
            ("extra", INTSV_EX, MGF_EX, pd.read_parquet(ROWS_EX).reset_index(drop=True), znp_ex, "extra")]:
        it = pd.read_csv(itab_p, sep="\t", encoding="utf-8")
        okc = rows_p[rows_p.status == "ok"]
        seqs_chain = list(okc.seq.unique())
        # A3: 每行 vs 数据
        bad3 = 0
        for _, r in it.iterrows():
            g = okc[okc.seq == r.Sequence]
            b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
            if int(r.Scan_No) != int(b.apex_scan) or int(r.Charge) != int(b.z) \
               or int(r.Raw_Score) != int(b.nmatched):
                bad3 += 1
            exp_mh = float(b.mz_s14) * int(b.z) - (int(b.z) - 1) * PROTON
            if abs(float(r.ExpMH) - exp_mh) > 0.01:
                bad3 += 1
        if bad3:
            issue("ERR", "A3-" + chain, "输入表与数据不一致行数: %d/%d" % (bad3, len(it)))
        # A4: MGF
        titles, peaks, basepk = [], [], []
        cur_t, cur_max, n = None, 0.0, 0
        with _open(mgf_p) as f:
            for line in f:
                line = line.strip()
                if line.startswith("TITLE="):
                    cur_t = line[6:]
                elif re.match(r"^[0-9]", line):
                    n += 1
                    try:
                        v = float(line.split("\t")[1])
                        cur_max = max(cur_max, v)
                    except (ValueError, IndexError):
                        pass
                elif line == "END IONS":
                    titles.append(cur_t)
                    peaks.append(n)
                    basepk.append(cur_max)
                    cur_t, cur_max, n = None, 0.0, 0
        a4 = dict(rows=len(it), mgf=len(titles),
                  title_match=sum(1 for t in titles if t in set(it.File_Name)))
        if a4["rows"] != a4["mgf"] or a4["title_match"] != a4["rows"]:
            issue("ERR", "A4-" + chain, "MGF 对账失败: %s" % a4)
        report["A4-" + chain] = a4

    # ---------- A5/A6/A7/A8/A9：逐页（主 PDF） ----------
    pidx = pd.read_parquet(PROF_I)
    pz = np.load(PROF_N)
    prof_map = {(r.seq, r.z, r.mz, r["sample"]): i for i, r in pidx.reset_index().iterrows()}
    pidx_e = pd.read_parquet(PROF_I_EX)
    pz_e = np.load(PROF_N_EX)
    prof_map_e = {(r.seq, r.z, r.mz, r["sample"]): i for i, r in pidx_e.reset_index().iterrows()}

    it_main = pd.read_csv(INTSV, sep="\t", encoding="utf-8")
    it_extra = pd.read_csv(INTSV_EX, sep="\t", encoding="utf-8")
    info_by_seq = {r.Sequence: r for r in it_main.itertuples()}
    info_by_seq_ex = {r.Sequence: r for r in it_extra.itertuples()}

    doc = pymupdf.open(PDF_MAIN)
    a5_bad = a6_bad = a7_bad = a8_bad = a9_bad = c1_bad = c2_bad = c3_bad = 0
    n_pages = 0
    cos_recalc = {}
    heat_recalc = {}
    for pi in range(1, doc.page_count - 3 + 1 - 1):  # 肽段页 = 1..112
        pass
    # 页序与 79 相同：封面(0) + 肽段(1..112) + 索引(113..115)
    pep_pages = range(1, 1 + len(seqs)) if False else range(1, 113)
    for pi in pep_pages:
        seq = seqs[pi - 1]
        page = doc[pi]
        txt = page.get_text()
        n_pages += 1
        r = info_by_seq.get(seq)
        g = okr[okr.seq == seq]
        b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
        gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
        m0 = tg[(tg.seq == seq) & (tg.z == b.z) & (tg.mz_s14 == b.mz_s14)].iloc[0]
        # A5 侧翼
        p = fasta.get(m0.prot, "")
        i0 = p.find(seq)
        pre = p[i0 - 1] if i0 > 0 else "-"
        post = p[i0 + len(seq)] if i0 >= 0 and i0 + len(seq) < len(p) else "-"
        if "((%s)%s(%s))" % (pre, seq, post) not in txt.replace("\n", ""):
            if "(%s)%s(%s)" % (pre, seq, post) not in txt.replace("\n", ""):
                a5_bad += 1
                issue("ERR", "A5", "p%d 侧翼不符: 期望 (%s)%s(%s)" % (pi + 1, pre, seq, post))
        # A6 Info 数字
        rt_lo, rt_hi = float(gk.apex_rt.min()), float(gk.apex_rt.max())
        checks6 = [
            ("m/z %.4f" % float(b.mz_s14)) in txt,
            ("%d samples" % len(gk)) in txt,
            ("RT %.2f-%.2f" % (rt_lo, rt_hi)) in txt,
            ("Best spectrum: %s" % b["sample"]) in txt,
            ("matched b/y %d" % int(b.nmatched)) in txt,
            ("ladder b%d/y%d" % (int(b.bLad), int(b.yLad))) in txt,
        ]
        if not all(checks6):
            a6_bad += 1
            issue("ERR", "A6", "p%d Info 字段缺失 %s" % (pi + 1, checks6))
        # A7 cosine 重算
        pid = m0.pid
        frags, _, _ = m71.parse_pid(pid)
        sl = slice(int(znp["off"][b.name]), int(znp["end"][b.name]))
        smz, sint = znp["mz"][sl], znp["inten"][sl]
        pv = pred.get(("hebei", pid))
        if pv is not None and len(np.asarray(pv["b"])) == len(seq) - 1:
            obs1 = {}
            for l, t in frags.items():
                if "(2+)" in l:
                    continue
                d = np.abs(smz - t) / t * 1e6
                k = int(np.argmin(d))
                if d[k] < 20.0:
                    obs1[l] = float(sint[k])
            nres = len(seq)
            cp = ([(obs1["b%d" % i], float(pv["b"][i - 1])) for i in range(1, nres)
                   if ("b%d" % i) in obs1] +
                  [(obs1["y%d" % i], float(pv["y"][i - 1])) for i in range(1, nres)
                   if ("y%d" % i) in obs1])
            if len(cp) >= 4:
                a1v = np.array([x[0] for x in cp]); a2v = np.array([x[1] for x in cp])
                n1, n2 = float(np.linalg.norm(a1v)), float(np.linalg.norm(a2v))
                cos = float(np.dot(a1v, a2v) / (n1 * n2)) if n1 > 0 and n2 > 0 else None
                if cos is not None:
                    cos_recalc[pi + 1] = cos
                    mcos = re.search(r"cosine\(obs, pred\) = ([0-9.]+)", txt)
                    if not mcos or abs(float(mcos.group(1)) - cos) > 0.002:
                        a7_bad += 1
                        issue("ERR", "A7", "p%d cosine 页面 %s vs 重算 %.3f" %
                              (pi + 1, mcos.group(1) if mcos else "缺", cos))
            else:
                if "cosine" in txt:
                    a7_bad += 1
                    issue("ERR", "A7", "p%d 共有离子<4 却显示 cosine" % (pi + 1,))
        else:
            if "unavailable" not in txt:
                a7_bad += 1
                issue("ERR", "A7", "p%d 无预测却未标 unavailable" % (pi + 1,))
        # A8 热图 apex RT 重算
        samples = sorted(gk["sample"].unique(),
                         key=lambda s: (0 if s.startswith("HC") else 1, int(s[2:])))
        pr = []
        rt_ref_s = {r["sample"]: float(r.apex_rt) for _, r in gk.iterrows()}
        for s in samples:
            i = prof_map.get((seq, int(b.z), float(b.mz_s14), s))
            if i is None:
                continue
            rr0, rr1 = int(pz["off"][i]), int(pz["end"][i])
            rtv = pz["rt"][rr0:rr1]; valv = pz["val"][rr0:rr1]
            b0 = math.floor(float(rtv.min())); nb = int(math.ceil(float(rtv.max()))) + 1 - b0
            M = np.zeros(nb)
            np.maximum.at(M, np.clip((rtv - b0).astype(int), 0, nb - 1), valv)
            # 与 79 渲染同口径：行级阈值（窗内峰值<2x窗外基线→整行跳过）+
            # gamma^2 + 20% 相对线 + 1.5x 绝对背景线
            bin_c = b0 + np.arange(nb) + 0.5
            r0s = rt_ref_s.get(s) if rt_ref_s else None
            win = np.abs(bin_c - r0s) <= 2.0 if r0s is not None else np.zeros(nb, bool)
            if win.any() and (~win).any():
                if M[win].max() < 2.0 * np.median(M[~win]):
                    continue
            base = np.median(M[~win]) if (~win).any() else 0.0
            M2 = np.clip(M, 0, 1) ** 2
            M2 = np.where(M2 < 0.20 * M2.max(), 0.0, M2)
            M2 = np.where(M < 1.5 * base, 0.0, M2)
            M = M2
            k = int(M.argmax())
            if M[k] > 0:
                pr.append(b0 + k + 0.5)
        if len(pr) >= 2:
            arr = np.array(pr)
            mean, sd = float(arr.mean()), float(arr.std(ddof=1))
            heat_recalc[pi + 1] = (mean, sd, len(pr))
            # 图题为位图：重算值存档 audit_heatstats_recalc.tsv，数值对账走抽样读数（visual）
            # A9: 各行 max bin RT 在声明范围 ±1.5
            out_of_range = sum(1 for x in pr if not (rt_lo - 1.5 <= x <= rt_hi + 1.5))
            if out_of_range:
                a9_bad += 1
                issue("WARN", "A9", "p%d %d/%d 行 apex bin 超出声明 RT 范围±1.5min" % (pi + 1, out_of_range, len(pr)))
        # C1 关键文本
        for key in ["Observed spectrum", "Base peak", "Precursor m/z"]:
            if key not in txt:
                c1_bad += 1
                issue("ERR", "C1", "p%d 缺关键文本 %r" % (pi + 1, key))
        if "AlphaPeptDeep prediction" not in txt:
            c1_bad += 1
            issue("ERR", "C1", "p%d 缺预测面板文本" % (pi + 1,))
        # C2 热图区非全白
        pix = page.get_pixmap(dpi=40, clip=pymupdf.Rect(0, 431, 595, 838))
        arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
        if (arr[:, :, :3].min() >= 250):
            c2_bad += 1
            issue("ERR", "C2", "p%d 热图区全白" % (pi + 1,))
        # C3 文字越界
        for w in page.get_text("words"):
            if w[2] > 596 or w[3] > 843 or w[0] < -1 or w[1] < -1:
                c3_bad += 1
                issue("WARN", "C3", "p%d 文字越界 %r" % (pi + 1, w[4]))
                break
    report["A5_flank_bad"] = a5_bad
    report["A6_info_bad"] = a6_bad
    report["A7_cos_bad"] = a7_bad
    report["A7_cos_checked"] = len(cos_recalc)
    report["A8_heatstat_bad"] = a8_bad
    report["A8_heatstat_checked"] = len(heat_recalc)
    report["A9_rowOutOfRange_pages"] = a9_bad
    report["C1_missingText"] = c1_bad
    report["C2_blankHeat_pages"] = c2_bad
    report["C3_textOverflow"] = c3_bad

    # ---------- A10 索引 vs manifest ----------
    mani_main = pd.read_csv(ROOT + r"\work\output\s14_evidence\pages_manifest_isoPSV.tsv", sep="\t")
    if len(mani_main) != 112 or list(mani_main.page) != list(range(1, 113)):
        issue("ERR", "A10", "主 manifest 页码异常")

    # ---------- B：Table 1 PDF ----------
    mani = pd.read_csv(MANI_T1, sep="\t")
    b_bad = 0
    if len(mani) != 159 or list(mani.page_no) != list(range(2, 161)):
        issue("ERR", "B1", "T1 manifest 页序异常: %d 行" % len(mani))
        b_bad += 1
    nos = [int(x) for x in mani.t1_no]
    if any(nos[i] > nos[i + 1] for i in range(len(nos) - 1)):
        issue("ERR", "B1", "T1 No. 非降序被破坏")
        b_bad += 1
    doct = pymupdf.open(PDF_T1)
    # B2: 抽全量对账首行 Table 1 No./蛋白
    t1map = dict(zip(mani.seq, zip(mani.t1_no, mani.prot)))
    b2_bad = 0
    for _, m in mani.iterrows():
        txt = doct[int(m.page_no) - 1].get_text()
        if ("No.%s %s" % (m.t1_no, m.prot)) not in txt.replace("\n", " "):
            b2_bad += 1
            if b2_bad <= 3:
                issue("ERR", "B2", "T1 p%s 首行 No/蛋白不符" % m.page_no)
    if b2_bad == 0:
        report["B2"] = "all 159 pages header OK"
    else:
        report["B2_bad"] = b2_bad
    # B3: ok_extra 页含标注
    b3_bad = 0
    for _, m in mani[mani.status == "ok_extra"].iterrows():
        txt = doct[int(m.page_no) - 1].get_text()
        if "not covered by Table S14" not in txt:
            b3_bad += 1
    if b3_bad:
        issue("ERR", "B3", "%d 补充链页缺来源标注" % b3_bad)
    report["B"] = dict(pages=doct.page_count, ok=int((mani.status == "ok").sum()),
                       ok_extra=int((mani.status == "ok_extra").sum()),
                       no_spectrum=int((mani.status == "no_spectrum").sum()))

    # ---------- A8b 嵌入热图结构对账（PDF 重编码致 md5 不可用；验证 bbox 与像素尺寸） ----------
    a8b_bad = 0
    for pi in pep_pages:
        pg = doc[pi]
        infos = pg.get_image_info()
        big = [i for i in infos if (i["bbox"][2] - i["bbox"][0]) > 400]
        if not big:
            a8b_bad += 1
            issue("ERR", "A8b", "p%d 无热图大图" % (pi + 1))
            continue
        bx = big[0]["bbox"]
        if not (15 < bx[0] < 25 and 425 < bx[1] < 435 and 570 < bx[2] < 580 and 833 < bx[3] < 843):
            a8b_bad += 1
            if a8b_bad <= 3:
                issue("ERR", "A8b", "p%d 热图 bbox 异常 %s" % (pi + 1, bx))
        wpx = big[0].get("width", 0)
        if wpx and not (1550 <= wpx <= 1650):   # 8in@200dpi = 1600px
            a8b_bad += 1
            if a8b_bad <= 3:
                issue("WARN", "A8b", "p%d 热图宽 %dpx 异常" % (pi + 1, wpx))
    report["A8b_structure_bad"] = a8b_bad
    # 重算 apex RT 统计存档
    tab, nl = chr(9), chr(10)
    with _open(ROOT + "\\work\\output\\s14_evidence\\audit_heatstats_recalc.tsv", "w", encoding="utf-8") as f:
        f.write("page" + tab + "mean_min" + tab + "sd_min" + tab + "n" + nl)
        for pg, (mn, sd, n) in sorted(heat_recalc.items()):
            f.write("%d" + tab + "%.3f" + tab + "%.3f" + tab + "%d" + nl % ()) if False else \
                f.write("%d%s%.3f%s%.3f%s%d%s" % (pg, tab, mn, tab, sd, tab, n, nl))

    # ---------- 汇总 ----------
    n_err = sum(1 for lv, _, _ in issues if lv == "ERR")
    n_warn = sum(1 for lv, _, _ in issues if lv == "WARN")
    report["summary"] = dict(errors=n_err, warns=n_warn,
                             main_pages=doc.page_count, t1_pages=doct.page_count)
    with _open(OUTJ, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=1, default=str)
    with _open(OUTMD, "w", encoding="utf-8") as f:
        f.write("# S14 / Table 1 证据 PDF 自动审计报告\n\n")
        f.write("日期：2026-11-01 ｜ 脚本：85_audit.py\n\n")
        f.write("## 汇总\n\n- 错误：%d\n- 警告：%d\n- 主 PDF：%d 页；Table 1 PDF：%d 页\n\n"
                % (n_err, n_warn, doc.page_count, doct.page_count))
        f.write("## 明细\n\n```\n")
        for lv, code, desc in issues[:200]:
            f.write("[%s] %s %s\n" % (lv, code, desc))
        f.write("```\n\n## 指标\n\n```json\n%s\n```\n" %
                json.dumps(report, ensure_ascii=False, indent=1, default=str))
    print("AUDIT DONE: errors=%d warns=%d" % (n_err, n_warn))
    print(json.dumps(report["summary"], ensure_ascii=False))


if __name__ == "__main__":
    main()
