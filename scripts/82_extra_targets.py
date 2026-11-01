# -*- coding: utf-8 -*-
"""
82_extra_targets.py —— Table 1 中缺 S14 页的 57 条肽：从 hebei v9 主表造补充靶点

对每条缺失肽段：
  1) 在 hebei 提取主表(rtqc 版)选最优前体（tier A>B、ApexSource ms1/consensus 优先、
     RTQC=0 优先、NMatchedFrags 大者优先）；
  2) 逐样本靶点 = 该前体在各 run 的主表行，rt_s14 = 该行 RTCorrected（与 S14 的
     逐样本 RT 同角色）；qty = hebei 报告缓存 PEP.Quantity（该 run 列）；
  3) pid 用主表 PrecursorId（真实修饰位点）；prot/gene 取 Table 1 行。
输出：work/tmp/s14_targets_extra.parquet（与 s14_targets 同 schema，供 71/73 以
     TGT 环境变量复用）+ extra_summary.json。
安全：路径常量；npz 为本项目自产缓存（路径经白名单根内校验）；open 经 _open()。
"""
import sys, os, json, collections
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import openpyxl

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


T1XLSX = ROOT + r"\Table 1.xlsx"
MASTER_TSV = ROOT + r"\work\output\summary\all_spectra_master.tsv"
ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
CACHE = ROOT + r"\work\data_cache\hebei_report_cache.npz"
MAPPING = ROOT + r"\work\data_cache\hebei_run_raw_mapping.tsv"
OUTP = ROOT + r"\work\tmp\s14_targets_extra.parquet"
OUTJ = ROOT + r"\work\tmp\extra_summary.json"


def main():
    # Table 1 肽 -> (prot, gene)
    wb = openpyxl.load_workbook(T1XLSX, read_only=True)
    t1rows = list(wb[wb.sheetnames[0]].iter_rows(values_only=True))[2:]
    t1pg, t1peps = {}, set()
    for r in t1rows:
        if not (r[0] and r[1]):
            continue
        if r[6]:
            for p in str(r[6]).replace("\n", ";").split(";"):
                p = p.strip()
                if p:
                    t1pg[p] = (str(r[1]).strip(), str(r[2]).strip() if r[2] else "")
                    t1peps.add(p)
    ok = pd.read_parquet(ROWS)
    have = set(ok[ok.status == "ok"].seq.unique())
    missing = sorted(t1peps - have)
    print("missing peptides:", len(missing))

    mcols = ["Dataset", "RunIdx", "PrecursorId", "StrippedSequence", "ModBody", "Charge",
             "PrecursorMz", "RTCorrected", "RTQC", "ApexSource", "QualityTier", "NMatchedFrags"]
    m = pd.read_csv(MASTER_TSV, sep="\t", usecols=mcols)
    m = m[m.Dataset == "hebei"]
    mp = pd.read_csv(MAPPING, sep="\t")
    mp["num"] = mp.run_name.str.extract(r"_ZhangYao_[HP]?_?(\d+)_ELU").astype(int)
    mp["sample"] = np.where(mp.group == "健康", "HC" + mp.num.astype(str), "TB" + mp.num.astype(str))
    ridx2sample = dict(zip(mp.run_idx, mp["sample"]))
    ridx2run = dict(zip(mp.run_idx, mp.run_name))
    c = np.load(CACHE, allow_pickle=True)
    pids = [str(p) for p in c["pid"]]
    pid_i = {p: i for i, p in enumerate(pids)}
    run_names = [str(x) for x in c["run_names"]]
    pep_qty = c["pep_qty"]  # (n_pid, n_run)

    recs = []
    stats = []
    for seq in missing:
        g = m[m.StrippedSequence == seq]
        if not len(g):
            stats.append((seq, "not_in_hebei_master"))
            continue
        g = g.copy()
        g["src_ok"] = g.ApexSource.isin(["ms1", "consensus"]).astype(int)
        g["rtqc0"] = (g.RTQC == 0).astype(int)
        g["tierA"] = (g.QualityTier == "A").astype(int)
        g = g.sort_values(["tierA", "rtqc0", "src_ok", "NMatchedFrags"], ascending=False)
        pid = g.PrecursorId.iloc[0]
        z = int(g.Charge.iloc[0])
        mz = float(g.PrecursorMz.iloc[0])
        gp = g[g.PrecursorId == pid]
        n_anchor = int(((gp.RTQC == 0) & gp.ApexSource.isin(["ms1", "consensus"])).sum())
        ci = pid_i.get(pid)
        for _, row in gp.iterrows():
            smp = ridx2sample.get(int(row.RunIdx))
            if not smp:
                continue
            qty = float(pep_qty[ci, int(row.RunIdx)]) if ci is not None else np.nan
            recs.append(dict(
                prot=t1pg[seq][0], gene=t1pg[seq][1], uc="", mtb="",
                seq=seq, z=z, mz_s14=mz, rt_s14=float(row.RTCorrected),
                fim="", fmz="", qval=np.nan, qty=qty, shared="No", nvar=0,
                mod_state="from_master", delta_mass=0.0, pid=pid,
                sample=smp, run_name=ridx2run[int(row.RunIdx)],
                group=None, in_main=True,
            ))
        stats.append((seq, "ok n_runs=%d anchor=%d" % (len(gp), n_anchor)))
    df = pd.DataFrame(recs)
    df.to_parquet(OUTP)
    summ = dict(missing=len(missing), targets=int(len(df)),
                seqs_ok=int(df.seq.nunique()),
                no_master=[s for s, r in stats if r == "not_in_hebei_master"],
                rows_per_seq_med=float(df.groupby("seq").size().median()) if len(df) else 0)
    with _open(OUTJ, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1, default=str)
    print(json.dumps(summ, ensure_ascii=False, indent=1, default=str))


if __name__ == "__main__":
    main()
