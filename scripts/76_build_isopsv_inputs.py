# -*- coding: utf-8 -*-
"""
76_build_isopsv_inputs.py —— 构造 isoPSV R 引擎的输入（每肽最优谱）

输出（work/risoPSV/in/）：
  s14_spectra.tsv   isoPSV/pFind 风格表：File_Name/Scan_No/Exp.MH+/Charge/Sequence/
                    Modification(pos,Name[AA];…, 1-based, 0=N-term)/Raw_Score(nmatched)/
                    Final_Score(b梯+y梯)/Proteins(Rv)/Positions(pos,pre,post 0-based)
                    + 附加信息列：Info1/Info2/Info3（信息行文字）、PredB/PredY（AlphaPeptDeep
                    强度逗号串，b1..b_{n-1}/y1..y_{n-1}）、PredCos/PredN（余弦与共有离子数）、
                    Sample/RT_apex/HeatIdx（热图行索引=行号）
  best.mgf          112 张最优谱（TITLE=File_Name，PEPMASS/CHARGE/RTINSECONDS/SCANS）
注：File_Name 与 MGF TITLE 一一对应（isoPSV find_spec 按 TITLE 定位）。
安全：路径常量；文件 open 经 _open() 白名单；pkl 为本项目自产预测（同 54 号口径）。
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


TARGETS = ROOT + r"\work\tmp\s14_targets.parquet"
ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
SPECNPZ = ROOT + r"\work\output\s14_evidence\extract\spectra_all.npz"
MAPPING = ROOT + r"\work\data_cache\hebei_run_raw_mapping.tsv"
PRED = ROOT + r"\work\output\summary\pred\AlphaPeptDeep.pkl"
FASTA = ROOT + r"\IPX0018588001_hebei\human_tuberculosis.fasta"
OUTDIR = ROOT + r"\work\risoPSV\in"
OUT_TSV = OUTDIR + r"\s14_spectra.tsv"
OUT_MGF = OUTDIR + r"\best.mgf"

_spec = importlib.util.spec_from_file_location("ext71", ROOT + r"\work\scripts\71_extract_targeted.py")
m71 = importlib.util.module_from_spec(_spec)
sys.modules["ext71"] = m71
_spec.loader.exec_module(m71)

PROTON = 1.007276466621


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


def pid_to_pfind_mods(pid, seq):
    """pid(Spectronaut) -> pFind 风格 'pos,Name[AA];…'（1-based；N 端=0）。synthetic: M 全氧化。"""
    parts = []
    if "(synthetic" in str(pid):
        st = str(pid).split("mods=")[-1].rstrip(")")
        if "ox" in st:
            parts += ["%d,Oxidation[M]" % (i + 1) for i, a in enumerate(seq) if a == "M"]
        if "acetyl" in st:
            parts.append("0,Acetyl[AnyN-term]")
        return ";".join(parts)
    body = str(pid)[1:str(pid).rindex("_.")]
    s2, i = "", 0
    while i < len(body):
        if body[i] == "[":
            j = body.index("]", i)
            nm = body[i + 1:j].strip()
            if nm.startswith("Carbamidomethyl"):
                parts.append("%d,Carbamidomethyl[C]" % (len(s2) if s2 else 0))
            elif nm.startswith("Oxidation"):
                parts.append("%d,Oxidation[M]" % len(s2))
            elif nm.startswith("Acetyl"):
                parts.append("0,Acetyl[AnyN-term]")
            i = j + 1
        else:
            s2 += body[i]
            i += 1
    return ";".join(parts)


def main():
    only = ""
    for a in sys.argv[1:]:
        if a.startswith("--only="):
            only = a.split("=", 1)[1]
    os.makedirs(OUTDIR, exist_ok=True)
    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    znp = np.load(SPECNPZ)
    mz_all, int_all, off, end = znp["mz"], znp["inten"], znp["off"], znp["end"]
    tg = pd.read_parquet(TARGETS)
    mp = pd.read_csv(MAPPING, sep="\t")
    run2raw = dict(zip(mp.run_name, mp.raw_path))
    prots = load_fasta_flanks()
    pred = pickle.load(_open(PRED, "rb"))

    ok = rows[rows.status == "ok"]
    seqs = list(ok.seq.unique())
    if only:
        seqs = [s for s in seqs if s == only]
    recs, mgf = [], []
    heat_idx = 0
    for seq in seqs:
        g = ok[ok.seq == seq]
        b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
        gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
        m0 = tg[(tg.seq == seq) & (tg.z == b.z) & (tg.mz_s14 == b.mz_s14)].iloc[0]
        pid = m0.pid
        raw_stem = os.path.basename(run2raw.get(b.run, "raw.raw"))[:-4]
        fn = "%s.%d.%d.%d.0.dta" % (raw_stem, int(b.apex_scan), int(b.apex_scan), int(b.z))
        # 侧翼（0-based pos）
        p_seq = prots.get(m0.prot, "")
        i0 = p_seq.find(seq)
        if i0 >= 0:
            posn = "%d,%s,%s" % (i0, p_seq[i0 - 1] if i0 > 0 else "-",
                                 p_seq[i0 + len(seq)] if i0 + len(seq) < len(p_seq) else "-")
            pre = p_seq[i0 - 1] if i0 > 0 else "-"
            post = p_seq[i0 + len(seq)] if i0 + len(seq) < len(p_seq) else "-"
        else:
            posn = "-,-,-"; pre = post = "-"
        mods = pid_to_pfind_mods(pid, seq)
        mh = float(b.mz_s14) * int(b.z) - (int(b.z) - 1) * PROTON
        rt_lo, rt_hi = float(gk.apex_rt.min()), float(gk.apex_rt.max())
        # 预测强度 + 余弦（1+ 共有离子，obs 用全谱匹配强度）
        frags, _, _ = m71.parse_pid(pid)
        pred_b = pred_y = ""
        cos_v = n_sh = ""
        pv = pred.get(("hebei", pid))
        if pv is not None and len(np.asarray(pv["b"])) == len(seq) - 1:
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
        rec = dict(
            File_Name=fn, Scan_No=int(b.apex_scan), ExpMH="%.4f" % mh, Charge=int(b.z),
            Sequence=seq, Modification=mods if mods else "-",
            Raw_Score=int(b.nmatched), Final_Score=int(b.bLad + b.yLad),
            Proteins=m0.prot, Positions=posn,
            Sample=b["sample"], RT_apex=float(b.apex_rt),
            Info1="| %s (%s) | (%s)%s(%s) | %d+ | mods: %s | m/z %.4f | %d samples, RT %.2f-%.2f min" %
                  (m0.prot, m0.gene, pre, seq, post, int(b.z),
                   mods.replace(";", "; ") if mods != "-" else "none",
                   float(b.mz_s14), len(gk), rt_lo, rt_hi),
            Info2="Best spectrum: %s  |  RT %.2f min (Table S14 %.2f, D%+.2f)  |  matched b/y %d "
                  "(ladder b%d/y%d)  |  MS1 precursor envelope: %s" %
                  (b["sample"], float(b.apex_rt), float(b.rt_s14), float(b.drt),
                   int(b.nmatched), int(b.bLad), int(b.yLad),
                   "yes" if b.ms1_flag else "not detected"),
            Info3="q-value (Table S14): %s  |  PG.Quantity: %s  |  main-text table peptide: %s  |  "
                  "unique peptide count %s / Mtb-complex specific %s" %
                  (("%.1e" % float(b.qval)) if pd.notna(b.qval) else "-",
                   ("%.3g" % float(b.qty)) if pd.notna(b.qty) else "-",
                   "yes" if m0.in_main else "no", m0.uc, m0.mtb),
            PredB=pred_b, PredY=pred_y, PredCos=cos_v, PredN=n_sh,
            HeatSeq=seq, HeatZ=int(b.z), HeatMz=float(b.mz_s14),
        )
        recs.append(rec)
        # MGF 谱
        sl = slice(int(off[b.name]), int(end[b.name]))
        smz, sint = mz_all[sl], int_all[sl]
        mgf.append("BEGIN IONS")
        mgf.append("TITLE=%s" % fn)
        mgf.append("PEPMASS=%.4f\t%.1f" % (float(b.mz_s14), float(sint.max()) if len(sint) else 0.0))
        mgf.append("CHARGE=%d+" % int(b.z))
        mgf.append("RTINSECONDS=%.1f" % (float(b.apex_rt) * 60.0))
        mgf.append("SCANS=%d" % int(b.apex_scan))
        for m, i in zip(smz, sint):
            mgf.append("%.4f\t%.1f" % (m, i))
        mgf.append("END IONS")
        mgf.append("")
        heat_idx += 1

    df = pd.DataFrame(recs)
    df.to_csv(_open(OUT_TSV, "w", encoding="utf-8"), sep="\t", index=False)
    with _open(OUT_MGF, "w", encoding="utf-8") as f:
        f.write("\n".join(mgf))
    summ = dict(rows=len(df), mgf_spectra=len(df),
                pred_rows=int((df.PredB != "").sum()))
    print(json.dumps(summ, ensure_ascii=False))


if __name__ == "__main__":
    main()
