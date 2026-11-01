# -*- coding: utf-8 -*-
"""
86_predict_missing.py —— 为 11 条无预测的肽补跑 AlphaPeptDeep 并并入主 pkl

原因：这 11 条肽不在 hebei/tianjin Spectronaut 报告内（S14 出自作者另一轮分析），
53 号全量预测只覆盖报告前体（37,814），故无预测条目。
流程：构造 apd 输入（precursor_id=我们页面的 pid；mods 用 AlphaPeptDeep 语法
Carbamidomethyl@C / Oxidation@M / Acetyl@Protein_N-term；mod_sites 1-based，N 端=0）
→ peptdeep library → 解析 HDF → 合并进 AlphaPeptDeep.pkl（键 (hebei, pid)）。
安全：路径常量；open 经 _open() 白名单。
"""
import sys, os, json, subprocess, pickle
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import h5py

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


INTSV = ROOT + r"\work\risoPSV\in\s14_spectra.tsv"
PRED = ROOT + r"\work\output\summary\pred\AlphaPeptDeep.pkl"
IN_TSV = ROOT + r"\work\tmp\apd_input_missing.tsv"
OUTDIR = ROOT + r"\work\tmp\apd_out_missing"
SETTINGS = ROOT + r"\work\tmp\apd_settings_missing.yaml"
PROTON = 1.007276466621
H2O = 18.01056468403
AA = dict(G=57.02146, C=103.00919, S=87.03203, P=97.05276, A=71.03711,
          V=99.06841, T=101.04768, L=113.08406, I=113.08406, N=114.04293, D=115.02694,
          Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
          R=156.10111, Y=163.06333, W=186.07931)


def parse_our_mods(mod_str, seq):
    """76 号输出的 pFind 风格 'pos,Name[AA];...' → AlphaPeptDeep 'Mod@AA;...' 与 mod_sites。"""
    if not mod_str or mod_str == "-" or mod_str == "nan":
        return "", ""
    mods, sites = [], []
    for part in str(mod_str).split(";"):
        part = part.strip()
        if not part:
            continue
        pos, name = part.split(",", 1)
        pos = int(pos)
        if name.startswith("Carbamidomethyl"):
            mods.append("Carbamidomethyl@C")
            sites.append(str(pos) if pos > 0 else "0")
        elif name.startswith("Oxidation"):
            mods.append("Oxidation@M")
            sites.append(str(pos))
        elif name.startswith("Acetyl"):
            mods.append("Acetyl@Protein_N-term")
            sites.append("0")
    return ";".join(mods), ";".join(sites)


def main():
    it = pd.read_csv(INTSV, sep="\t", encoding="utf-8")
    miss = it[it.PredB.isna() | (it.PredB == "")]
    print("missing predictions:", len(miss))
    if not len(miss):
        return
    recs = []
    for _, r in miss.iterrows():
        seq = r.Sequence
        z = int(r.Charge)
        # m/z 用理论值（CAM 固定 + 变基）
        mass = sum(AA[a] for a in seq) + H2O + 57.021464 * seq.count("C")
        mods_apd, sites = parse_our_mods(r.Modification, seq)
        for m in (mods_apd.split(";") if mods_apd else []):
            if m == "Oxidation@M":
                mass += 15.994915
            elif m == "Acetyl@Protein_N-term":
                mass += 42.010565
        mz = (mass + z * PROTON) / z
        recs.append(dict(precursor_id=r.File_Name.replace(".dta", "").rsplit(".", 1)[0]
                         if False else r.File_Name.split(".dta")[0],
                         sequence=seq, mods=mods_apd, mod_sites=sites,
                         charge=z, mz=round(mz, 5), rt=0.0,
                         _pid_target=None))
    df = pd.DataFrame(recs)
    # precursor_id 必须等于我们查 pkl 的键：Spectronaut 形 pid（76 号用 pid 列，此处以
    # Sequence+Modification 重构同形 pid）
    def pid_like(seq, charge, mods_str):
        body = seq
        if mods_str and mods_str not in ("-", "nan", ""):
            parts = []
            for part in str(mods_str).split(";"):
                p, n = part.split(",", 1)
                p = int(p)
                parts.append((p, n))
            out = ""
            prev = 0
            nterm = None
            for p, n in sorted(parts):
                if p == 0:
                    nterm = n
                    continue
                out += body[prev:p] + "[" + n + "]"
                prev = p
            out += body[prev:]
            if nterm:
                out = "[" + nterm + "]" + out
            return "_" + out + "_.%d" % int(charge)
        return "_%s_.%d" % (seq, int(charge))
    df["precursor_id"] = [pid_like(s, z, m) for s, z, m in zip(
        miss.Sequence, miss.Charge, miss.Modification)]
    df[["precursor_id", "sequence", "mods", "mod_sites", "charge", "mz", "rt"]] \
        .to_csv(_open(IN_TSV, "w", encoding="utf-8"), sep="\t", index=False)
    print(df[["precursor_id", "mods"]].to_string(index=False))
    # 配置（复制 full 版改输入输出）
    src = _open(ROOT + r"\work\tmp\apd_settings_full.yaml").read()
    src = src.replace("apd_input_all.tsv", "apd_input_missing.tsv")
    src = src.replace("apd_out_full", "apd_out_missing")
    with _open(SETTINGS, "w", encoding="utf-8") as f:
        f.write(src)
    # 运行 peptdeep
    exe = r"D:\Bioinformatics\peptdeep\peptdeep.exe"
    argv = [exe, "library", SETTINGS]
    ret = subprocess.run(argv, capture_output=True, text=True, shell=False,
                         cwd=r"D:\Bioinformatics\peptdeep")
    print("peptdeep rc=", ret.returncode)
    if ret.returncode != 0:
        print(ret.stdout[-2000:], ret.stderr[-2000:])
        return
    # 解析 HDF 并并入 pkl
    hdf = os.path.join(OUTDIR, "predict.speclib.hdf")
    pred = pickle.load(_open(PRED, "rb"))
    added = 0
    with h5py.File(hdf, "r") as f:
        p = f["library/precursor_df"]
        pid_arr = p["precursor_id"][:]
        frag = f["library/fragment_intensity_df"]
        b1 = np.asarray(frag["b_z1"], dtype=np.float32)
        y1 = np.asarray(frag["y_z1"], dtype=np.float32)
        starts = np.asarray(p["frag_start_idx"])
        stops = np.asarray(p["frag_stop_idx"])
        for i in range(len(pid_arr)):
            pid = pid_arr[i].decode() if isinstance(pid_arr[i], bytes) else str(pid_arr[i])
            s, e = int(starts[i]), int(stops[i])
            if e <= s:
                continue
            b = b1[s:e]           # 行 i=裂解位点 i → b_{i+1}
            yrow = y1[s:e]        # 行 i → y_{L-i}
            y = yrow[::-1]        # 反转成 y1..y_{L-1}（与 61 同口径）
            pred[("hebei", pid)] = {"b": b, "y": y}
            added += 1
    with _open(PRED + ".tmp", "wb") as f:
        pickle.dump(pred, f)
    os.replace(PRED + ".tmp", PRED)
    print(json.dumps(dict(added=added, total=len(pred))))


if __name__ == "__main__":
    main()
