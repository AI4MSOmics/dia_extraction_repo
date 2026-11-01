# -*- coding: utf-8 -*-
"""
71_extract_targeted.py —— S14 RT 锚定靶向提取：一 (前体,样本) 一谱

对 S14 每行（肽段前体 × 样本 × S14 参考RT）：
  1) 由 PrecursorId（或合成）计算理论 b/y 碎片（电荷 1..min(2,z)，20 ppm）；
  2) raw 内定位含该 m/z 的 DIA 窗口，取 RT ∈ [rt_s14-0.6, +0.6] 的窗口 MS2 扫描
     （不足 2 张时放宽到 ±1.2 min）；
  3) 每候选扫描按 (bLad+yLad, 匹配数, 匹配强度) 择优 → apex；
  4) S14 'Fragment m/z' 列逐个在 apex 谱 ±20 ppm 复核（S14 碎片回收数）；
  5) apex 附近最近 MS1 扫描查 mono/M+1/M+2 包络（10 ppm）作前体证据。
输出（work/output/s14_evidence/extract/）：
  rows_all.parquet      每 (target,sample) 一行结果（含 drt/梯/匹配/包络/碎片回收）
  spectra_all.npz       apex 谱峰数组（concat + 偏移索引，无 pickle）
  summary.json          统计
安全: 路径全部硬编码常量；文件 open 一律经 _open()（root 白名单）。
"""
import sys, os, json, time
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

ROOT = r"F:\IPX0018588000"
BIN = r"D:\Bioinformatics\pFind3.1\bin"


def _open(path, mode="r", root=ROOT, **kw):
    """受 root 白名单限制的 open（同 xic_raw._open 惯用法）。"""
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


TARGETS = os.environ.get("TGT", ROOT + r"\work\tmp\s14_targets.parquet")
MAPPING = ROOT + r"\work\data_cache\hebei_run_raw_mapping.tsv"
OUTDIR = os.environ.get("OUTD", ROOT + r"\work\output\s14_evidence\extract")
OUT_ROWS = OUTDIR + r"\rows_all.parquet"
OUT_SPEC = OUTDIR + r"\spectra_all.npz"
OUT_SUMM = OUTDIR + r"\summary.json"

import clr
clr.AddReference("System")
from System.Globalization import CultureInfo
from System.Threading import Thread
Thread.CurrentThread.CurrentCulture = CultureInfo("en-US")
Thread.CurrentThread.CurrentUICulture = CultureInfo("en-US")
clr.AddReference(os.path.join(BIN, "ThermoFisher.CommonCore.Data.dll"))
clr.AddReference(os.path.join(BIN, "ThermoFisher.CommonCore.RawFileReader.dll"))
from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
from ThermoFisher.CommonCore.Data.Business import Device, Scan
from ThermoFisher.CommonCore.Data.Interfaces import IScanEventBase

PROTON = 1.007276466621
H2O = 18.01056468403
ISO = 1.0033548
AA = dict(G=57.02146, C=103.00919, S=87.03203, P=97.05276, A=71.03711,
          V=99.06841, T=101.04768, L=113.08406, I=113.08406, N=114.04293, D=115.02694,
          Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
          R=156.10111, Y=163.06333, W=186.07931, C_CAM=103.00919 + 57.021464,
          M_OX=131.04049 + 15.994915, NTerm_AC=42.010565)
MODMAP = {"Carbamidomethyl (C)": "C_CAM", "Oxidation (M)": "M_OX"}


def parse_pid(pid):
    """'_SEQ[Mod (X)]..._.z' -> (residue mass 列表 b 端起点, 序列, z, mods dict pos->res)
    返回 (prefix_masses, suffix_masses 前置表, n, z)。synthetic pid 按 mods 描述处理。"""
    z = int(pid.rsplit(".", 1)[1][0])
    body = pid[1:pid.rindex("_.")] if "_.\"" not in pid else pid[1:pid.index("_")]
    if "(synthetic" in pid:  # _SEQ_.z(synthetic,mods=ox1)
        seq = pid[1:pid.index("_.")]
        st = pid.split("mods=")[-1].rstrip(")")
        n_ox = 1 if "ox" in st else 0
        n_ac = 1 if "acetyl" in st else 0
        mods = {}
        if n_ox and seq.count("M") == 1:
            mods[seq.index("M")] = "M_OX"
        elif n_ox:
            mods = {i: "M_OX" for i, a in enumerate(seq) if a == "M"}  # 全置位（歧义注明）
        if n_ac:
            mods[-1] = "NTerm_AC"
        return build_fragments(seq, mods), seq, z
    seq = ""
    mods = {}
    i = 0
    while i < len(body):
        if body[i] == "[":
            j = body.index("]", i)
            name = body[i + 1:j].strip()
            key = MODMAP.get(name)
            if key is None:
                raise ValueError("unknown mod: " + name)
            pos = len(seq) - 1 if seq else -1
            mods[pos] = key
            i = j + 1
        else:
            seq += body[i]
            i += 1
    return build_fragments(seq, mods), seq, z


def build_fragments(seq, mods):
    """返回理论碎片 dict: label -> m/z（b1..bn, y1..yn, 电荷 1..2）。"""
    n = len(seq)
    residue = []
    for i, a in enumerate(seq):
        residue.append(AA[mods.get(i, a)])
    nterm = AA[mods.get(-1, "G")] - AA["G"]  # N 端修饰偏移
    pre = np.cumsum([0.0] + residue)
    total = pre[-1] + H2O + nterm
    frags = {}
    for i in range(1, n):  # i = N 端起残基数；b_i 与 y_{n-i}
        b = pre[i] + nterm + PROTON
        y = total - pre[i] + PROTON
        for zc in (1, 2):
            if zc == 2 and i < 2:
                continue
            frags["b%d" % i if zc == 1 else "b%d(2+)" % i] = (b + (zc - 1) * PROTON) / zc
            frags["y%d" % (n - i) if zc == 1 else "y%d(2+)" % (n - i)] = (y + (zc - 1) * PROTON) / zc
    return frags


def ladder_metrics(labels):
    """最长连续 b / y 梯（1+ 标签；'b12(2+)' 形式不参与）。"""
    bs = sorted(int(l[1:]) for l in labels if l.startswith("b") and l[1:].isdigit())
    ys = sorted(int(l[1:]) for l in labels if l.startswith("y") and l[1:].isdigit())

    def longest(seq):
        best = cur = 0
        prev = -9
        for v in seq:
            cur = cur + 1 if v == prev + 1 else 1
            best = max(best, cur)
            prev = v
        return best
    return longest(bs), longest(ys)


class RawIndex:
    """一次枚举：MS2 扫描 RT/窗口边界 + MS1 扫描 RT。"""

    def __init__(self, rf):
        self.rf = rf
        last = int(rf.RunHeaderEx.LastSpectrum)
        ms1_s, ms1_rt, ms2_s, ms2_rt, lo, hi = [], [], [], [], [], []
        for s in range(1, last + 1):
            ev = IScanEventBase(rf.GetScanEventForScanNumber(int(s)))
            if int(ev.MSOrder) == 1:  # MS1
                ms1_s.append(s)
                ms1_rt.append(rf.RetentionTimeFromScanNumber(s))
                continue
            try:
                l = float(ev.GetFirstPrecursorMass(0))
                h = float(ev.GetLastPrecursorMass(0))
            except Exception:
                continue
            if h <= 0 or h - l > 500:  # 防御：跳过异常窗口
                continue
            ms2_s.append(s)
            ms2_rt.append(rf.RetentionTimeFromScanNumber(s))
            lo.append(l)
            hi.append(h)
        self.ms1_s = np.array(ms1_s); self.ms1_rt = np.array(ms1_rt, dtype=float)
        self.ms2_s = np.array(ms2_s, dtype=np.int64); self.ms2_rt = np.array(ms2_rt, dtype=float)
        self.lo = np.array(lo); self.hi = np.array(hi)

    def centroid(self, scan):
        sc = Scan.FromFile(self.rf, int(scan))
        if not sc.HasCentroidStream:
            return None, None
        cs = self.rf.GetCentroidStream(int(scan), False)
        return (np.array(list(cs.Masses), dtype=float),
                np.array(list(cs.Intensities), dtype=float))

    def ms2_in(self, mz, rt, tol):
        m = (self.lo <= mz) & (self.hi >= mz) & (np.abs(self.ms2_rt - rt) <= tol)
        return self.ms2_s[m], self.ms2_rt[m]

    def ms1_near(self, rt, k=3):
        """最近的 k 张 MS1 扫描号（按 |RT-rt| 排序）。"""
        if len(self.ms1_rt) == 0:
            return []
        order = np.argsort(np.abs(self.ms1_rt - rt))[:k]
        return [int(self.ms1_s[i]) for i in order]


def match_ppm(mz, inten, targets, ppm=20.0):
    """返回 (labels, matched_mz, matched_int)；targets: label->mz。"""
    if mz is None or len(mz) == 0:
        return [], np.array([]), np.array([])
    lab, hit_m, hit_i = [], [], []
    for l, t in targets.items():
        d = np.abs(mz - t) / t * 1e6
        k = int(np.argmin(d)) if len(d) else 0
        if len(d) and d[k] < ppm:
            lab.append(l)
            hit_m.append(mz[k])
            hit_i.append(inten[k])
    return lab, np.array(hit_m), np.array(hit_i)


def open_raw(path):
    rf = RawFileReaderAdapter.FileFactory(path)
    for attempt in range(3):
        try:
            rf.SelectInstrument(Device.MS, 1)
            break
        except Exception:
            time.sleep(0.5)
    return rf


def main():
    limit_runs = 0
    for a in sys.argv[1:]:
        if a.startswith("--runs="):
            limit_runs = int(a.split("=", 1)[1])
    os.makedirs(OUTDIR, exist_ok=True)
    tg = pd.read_parquet(TARGETS)
    mp = pd.read_csv(MAPPING, sep="\t")
    tg = tg.merge(mp[["run_name", "raw_path"]], on="run_name", how="left")
    tg = tg[tg.raw_path.notna()]
    if limit_runs:
        keep = sorted(tg.run_name.unique())[:limit_runs]
        tg = tg[tg.run_name.isin(keep)]
    print("rows to extract:", len(tg), "| runs:", tg.run_name.nunique())

    # 预解析理论碎片（按 target 去重）
    frag_cache = {}
    for key, g in tg.groupby(["seq", "z", "mz_s14"]):
        pid = g.pid.iloc[0]
        try:
            frags, seq2, _ = parse_pid(pid)
        except ValueError as e:
            print("  pid parse fail:", pid, e)
            frags = None
        frag_cache[key] = frags

    rows_out, spec_mz, spec_int, spec_len = [], [], [], []
    t0 = time.time()
    for run_name, g in tg.groupby("run_name", sort=True):
        rawp = g.raw_path.iloc[0]
        if not os.path.exists(rawp):
            print("MISSING raw:", rawp)
            continue
        rf = open_raw(rawp)
        idx = RawIndex(rf)
        print("[%5.1fs] %s  MS2=%d MS1=%d  targets=%d" %
              (time.time() - t0, run_name[:40], len(idx.ms2_s), len(idx.ms1_s), len(g)))
        for _, r in g.iterrows():
            key = (r.seq, r.z, r.mz_s14)
            frags = frag_cache.get(key)
            rec = dict(prot=r.prot, gene=r.gene, seq=r.seq, z=r.z, mz_s14=r.mz_s14,
                       mod_state=r.mod_state, sample=r["sample"], group=r.group,
                       in_main=r.in_main, rt_s14=r.rt_s14, fim=r.fim, qval=r.qval,
                       qty=r.qty, run=run_name)
            fmz = []
            if isinstance(r.fmz, str) and r.fmz:
                try:
                    fmz = [float(x) for x in r.fmz.split(",") if x.strip()]
                except ValueError:
                    fmz = []
            if frags is None:
                rec.update(status="pid_parse_fail")
                rows_out.append(rec)
                spec_len.append(0)
                continue
            scans, rts = idx.ms2_in(r.mz_s14, r.rt_s14, 0.6)
            if len(scans) < 2:
                scans, rts = idx.ms2_in(r.mz_s14, r.rt_s14, 1.2)
            best = None
            for sc, rt in zip(scans, rts):
                mz, inten = idx.centroid(sc)
                if mz is None:
                    continue
                lab, hm, hi = match_ppm(mz, inten, frags)
                bL, yL = ladder_metrics(lab)
                score = (bL + yL, len(lab), float(hi.sum()))
                if best is None or score > best[0]:
                    best = (score, int(sc), float(rt), mz, inten, lab, bL, yL, len(lab), float(hi.sum()))
            if best is None:
                rec.update(status="no_scan")
                rows_out.append(rec)
                spec_len.append(0)
                continue
            _, sc, rt, mz, inten, lab, bL, yL, nm, msum = best
            tot = float(inten.sum())
            # S14 碎片回收
            s14hit = 0
            for t in fmz:
                if len(mz):
                    d = np.abs(mz - t) / t * 1e6
                    if d.min() < 20.0:
                        s14hit += 1
            # MS1 包络：最近 3 张 MS1 扫描里取 mono 最强的一张作证据
            mono = iso1 = iso2 = np.nan
            for ms1_scan in idx.ms1_near(rt, k=3):
                m1, i1 = idx.centroid(ms1_scan)
                if m1 is None or not len(m1):
                    continue
                d = np.abs(m1 - r.mz_s14) / r.mz_s14 * 1e6
                k = int(np.argmin(d))
                if d[k] < 10.0 and (np.isnan(mono) or i1[k] > mono):
                    mono = float(i1[k])
                    iso1 = iso2 = np.nan
                    for off, store in ((ISO / r.z, "iso1"), (2 * ISO / r.z, "iso2")):
                        t = r.mz_s14 + off
                        dd = np.abs(m1 - t) / t * 1e6
                        kk = int(np.argmin(dd))
                        if dd[kk] < 10.0:
                            if store == "iso1":
                                iso1 = float(i1[kk])
                            else:
                                iso2 = float(i1[kk])
                if np.isfinite(mono) and not np.isnan(i1[k]) and i1[k] == mono:
                    break
            spec_mz.append(mz)
            spec_int.append(inten)
            spec_len.append(len(mz))
            rec.update(status="ok", apex_scan=sc, apex_rt=rt, drt=rt - r.rt_s14,
                       bLad=bL, yLad=yL, nmatched=nm, matched_int_frac=(msum / tot if tot > 0 else np.nan),
                       s14_frag_hits=s14hit, s14_frag_n=len(fmz),
                       ms1_mono=mono, ms1_iso1=iso1, ms1_iso2=iso2,
                       ms1_flag=int(np.isfinite(mono)))
            rows_out.append(rec)
        rf.Dispose()

    df = pd.DataFrame(rows_out)
    df.to_parquet(OUT_ROWS)
    lengths = np.array(spec_len, dtype=np.int64)
    ends = np.cumsum(lengths)
    offs = ends - lengths
    np.savez(OUT_SPEC,
             mz=np.concatenate(spec_mz) if spec_mz else np.array([]),
             inten=np.concatenate(spec_int) if spec_int else np.array([]),
             off=offs, end=ends)
    ok = df[df.status == "ok"]
    summ = dict(
        rows=len(df), ok=len(ok),
        ok_rate=float(len(ok) / max(len(df), 1)),
        abs_drt_med=float(ok.drt.abs().median()) if len(ok) else None,
        abs_drt_p95=float(ok.drt.abs().quantile(0.95)) if len(ok) else None,
        drt_le_025=float((ok.drt.abs() <= 0.25).mean()) if len(ok) else None,
        ms1_flag_rate=float(ok.ms1_flag.mean()) if len(ok) else None,
        s14_frag_full_rate=float((ok.s14_frag_hits == ok.s14_frag_n).mean()) if len(ok) else None,
        ladder_ge6=float(((ok.bLad + ok.yLad) >= 6).mean()) if len(ok) else None,
        n_runs=int(df.run.nunique()),
    )
    with _open(OUT_SUMM, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print(json.dumps(summ, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
