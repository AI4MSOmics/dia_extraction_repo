# -*- coding: utf-8 -*-
"""
84_extra_inputs.py —— 57 条补充肽的 isoPSV 输入构建（76 号的 extract_extra 版）

数据源 = extract_extra/rows_all.parquet（v9 主表 RT 锚定提取）+ s14_targets_extra.parquet。
输出：work/risoPSV/in_extra/{s14_spectra.tsv, best.mgf}（schema 与 76 相同），
供 78_s14_export.R 以 IN_TSV 环境变量渲染 topblock_extra.pdf。
Info1-3 合成：蛋白/基因来自 Table 1；RT 范围/最优样本/匹配信息来自 extra rows；
mods 文本从 pid 解析；预测强度/余弦从 AlphaPeptDeep pkl 取（键=主表真实 pid）。
安全：路径常量；open 经 _open() 白名单；pkl 为本项目自产。
"""
import sys, os, json, pickle, importlib.util
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


TARGETS = ROOT + r"\work\tmp\s14_targets_extra.parquet"
ROWS = ROOT + r"\work\output\s14_evidence\extract_extra\rows_all.parquet"
SPECNPZ = ROOT + r"\work\output\s14_evidence\extract_extra\spectra_all.npz"
MAPPING = ROOT + r"\work\data_cache\hebei_run_raw_mapping.tsv"
PRED = ROOT + r"\work\output\summary\pred\AlphaPeptDeep.pkl"
OUTDIR = ROOT + r"\work\risoPSV\in_extra"

_spec = importlib.util.spec_from_file_location("ext71", ROOT + r"\work\scripts\71_extract_targeted.py")
m71 = importlib.util.module_from_spec(_spec)
sys.modules["ext71"] = m71
_spec.loader.exec_module(m71)

PROTON = 1.007276466621


def mods_text(pid):
    parts = []
    body = str(pid)
    if "Carbamidomethyl (C)" in body:
        parts.append("Carbamidomethyl (C)")
    if "Oxidation (M)" in body:
        parts.append("Oxidation (M)")
    if "Acetyl (Protein N-term)" in body:
        parts.append("Acetyl (N-term)")
    return "; ".join(parts) if parts else "none"


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    znp = np.load(SPECNPZ)
    mz_all, int_all, off, end = znp["mz"], znp["inten"], znp["off"], znp["end"]
    tg = pd.read_parquet(TARGETS)
    mp = pd.read_csv(MAPPING, sep="\t")
    run2raw = dict(zip(mp.run_name, mp.raw_path))
    pred = pickle.load(_open(PRED, "rb"))

    ok = rows[rows.status == "ok"]
    recs, mgf = [], []
    n_pred = 0
    for seq in list(ok.seq.unique()):
        g = ok[ok.seq == seq]
        b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
        gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
        m0 = tg[(tg.seq == seq) & (tg.z == b.z) & (tg.mz_s14 == b.mz_s14)].iloc[0]
        pid = m0.pid
        raw_stem = os.path.basename(run2raw.get(b.run, "raw.raw"))[:-4]
        fn = "%s.%d.%d.%d.0.dta" % (raw_stem, int(b.apex_scan), int(b.apex_scan), int(b.z))
        mh = float(b.mz_s14) * int(b.z) - (int(b.z) - 1) * PROTON
        rt_lo, rt_hi = float(gk.apex_rt.min()), float(gk.apex_rt.max())
        # 修饰串（pFind 风格）
        modparts = []
        body = str(pid)[1:str(pid).rindex("_.")]
        s2, i = "", 0
        while i < len(body):
            if body[i] == "[":
                j = body.index("]", i)
                nm = body[i + 1:j].strip()
                if nm.startswith("Carbamidomethyl"):
                    modparts.append("%d,Carbamidomethyl[C]" % len(s2))
                elif nm.startswith("Oxidation"):
                    modparts.append("%d,Oxidation[M]" % len(s2))
                elif nm.startswith("Acetyl"):
                    modparts.append("0,Acetyl[AnyN-term]")
                i = j + 1
            else:
                s2 += body[i]
                i += 1
        mods = ";".join(modparts)
        # 预测 + 余弦
        frags, _, _ = m71.parse_pid(pid)
        pred_b = pred_y = ""
        cos_v = n_sh = ""
        pv = pred.get(("hebei", pid))
        if pv is not None and len(np.asarray(pv["b"])) == len(seq) - 1:
            n_pred += 1
            pred_b = ",".join("%.4f" % v for v in np.asarray(pv["b"], dtype=float))
            pred_y = ",".join("%.4f" % v for v in np.asarray(pv["y"], dtype=float))
            sl = slice(int(off[b.name]), int(end[b.name]))
            smz, sint = mz_all[sl], int_all[sl]
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
                a1 = np.array([x[0] for x in cp]); a2 = np.array([x[1] for x in cp])
                n1, n2 = float(np.linalg.norm(a1)), float(np.linalg.norm(a2))
                if n1 > 0 and n2 > 0:
                    cos_v = "%.3f" % float(np.dot(a1, a2) / (n1 * n2))
                    n_sh = str(len(cp))
        recs.append(dict(
            File_Name=fn, Scan_No=int(b.apex_scan), ExpMH="%.4f" % mh, Charge=int(b.z),
            Sequence=seq, Modification=mods if mods else "-",
            Raw_Score=int(b.nmatched), Final_Score=int(b.bLad + b.yLad),
            Proteins=m0.prot, Positions="-,-,-",
            Sample=b["sample"], RT_apex=float(b.apex_rt),
            Info1="| %s (%s) | %s | %d+ | mods: %s | m/z %.4f | %d samples, RT %.2f-%.2f min" %
                  (m0.prot, m0.gene, seq, int(b.z), mods_text(pid) if mods_text(pid) != "none" else "none",
                   float(b.mz_s14), len(gk), rt_lo, rt_hi),
            Info2="Best spectrum: %s  |  RT %.2f min  |  matched b/y %d (ladder b%d/y%d)  |  "
                  "MS1 precursor envelope: %s  |  RT anchor: v9 master (RTCorrected)" %
                  (b["sample"], float(b.apex_rt), int(b.nmatched), int(b.bLad), int(b.yLad),
                   "yes" if b.ms1_flag else "not detected"),
            Info3="Table 1 peptide not covered by Table S14 - evidence via v9-master "
                  "RT-anchored extraction  |  main-text table peptide: yes",
            PredB=pred_b, PredY=pred_y, PredCos=cos_v, PredN=n_sh,
            HeatSeq=seq, HeatZ=int(b.z), HeatMz=float(b.mz_s14),
        ))
        sl = slice(int(off[b.name]), int(end[b.name]))
        smz, sint = mz_all[sl], int_all[sl]
        mgf.append("BEGIN IONS")
        mgf.append("TITLE=%s" % fn)
        mgf.append("PEPMASS=%.4f\t%.1f" % (float(b.mz_s14), float(sint.max()) if len(sint) else 0.0))
        mgf.append("CHARGE=%d+" % int(b.z))
        mgf.append("RTINSECONDS=%.1f" % (float(b.apex_rt) * 60.0))
        mgf.append("SCANS=%d" % int(b.apex_scan))
        for m, i2 in zip(smz, sint):
            mgf.append("%.4f\t%.1f" % (m, i2))
        mgf.append("END IONS")
        mgf.append("")

    df = pd.DataFrame(recs)
    df.to_csv(_open(OUTDIR + r"\s14_spectra.tsv", "w", encoding="utf-8"), sep="\t", index=False)
    with _open(OUTDIR + r"\best.mgf", "w", encoding="utf-8") as f:
        f.write("\n".join(mgf))
    print(json.dumps(dict(rows=len(df), pred=n_pred), ensure_ascii=False))


if __name__ == "__main__":
    main()
