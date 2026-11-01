# -*- coding: utf-8 -*-
"""
70_build_s14_targets.py —— 把 Supplementary Table S14 转成靶向提取靶点表

输出（work/tmp/s14_targets.parquet + s14_targets_summary.json）:
  每行 = 一个 (peptide target, sample) 待提取组合
  列:
    prot, gene, uc, mtb, seq, z, mz_s14, rt_s14, fim, fmz, qval, qty, shared, nvar,
    mod_state (unmod/oxN/acetyl/…), pid (hebei 报告匹配到的 PrecursorId 或合成 _SEQ_.z),
    sample, run_name, group(HC/TB), in_main (正文 Table 2 肽段标记)
  注: 不存 raw_path——下游脚本（71_extract_targeted.py）自行按 run_name 从
      hebei_run_raw_mapping.tsv 重查并经 _open() 白名单校验 raw 路径。
规则:
  1) S14 的 (seq, z, mz) 为前体单位（同 seq+z 可有氧化/未氧化两行）；
  2) 修饰态由 m/z 与理论质量的差值判定: 0 / +15.9949*nOx / +42.0106 (N-term Acetyl) 组合;
  3) 样本名 HCxx/TBxx -> hebei run（run 名 _xx_ + 健康/病例 组）;
  4) PrecursorId 优先取自提取主表（真实修饰位点）；否则合成 stripped pid（注明 synthetic）。
安全: 本脚本输入输出路径全部为硬编码常量（F 盘项目树内）；文件 open 一律经
     _open()（root 白名单 + PermissionError，校验与打开同函数，同 xic_raw 惯用法）；
     不使用 pickle 反序列化。
"""
import sys, os, json, collections, re
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import openpyxl

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    """受 root 白名单限制的 open（同 xic_raw._open 惯用法）：校验与打开同函数，越界显式拒绝。"""
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


S14XLSX = ROOT + r"\Supplementary Table_2026.09.03(1).xlsx"
DOCX = ROOT + r"\Manuscript_2026.09.04(1).docx"
MAPPING = ROOT + r"\work\data_cache\hebei_run_raw_mapping.tsv"
MASTER_SUBSET = ROOT + r"\work\tmp\hebei_master_subset.parquet"
MASTER_TSV = ROOT + r"\work\output\summary\all_spectra_master.tsv"
OUTP = ROOT + r"\work\tmp\s14_targets.parquet"
OUTJ = ROOT + r"\work\tmp\s14_targets_summary.json"

AA = dict(G=57.02146, C=103.00919 + 57.02146, S=87.03203, P=97.05276, A=71.03711,
          V=99.06841, T=101.04768, L=113.08406, I=113.08406, N=114.04293, D=115.02694,
          Q=128.05858, K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841,
          R=156.10111, Y=163.06333, W=186.07931)  # C 已含固定 Carbamidomethyl
H2O = 18.01056468403
PROTON = 1.007276466621
OX = 15.994915
ACETYL = 42.010565


def neutral_mass(seq: str) -> float:
    return sum(AA[a] for a in seq) + H2O


def resolve_mod_state(seq, z, mz):
    """返回 (mod_state 描述, delta_neutral)。C 已按固定 CAM 计入。"""
    m_theo = neutral_mass(seq)
    m_obs = mz * z - z * PROTON
    d = m_obs - m_theo
    cands = []
    for n_ox in range(seq.count("M") + 1):
        for n_ac in range(0, 2):
            dd = d - n_ox * OX - n_ac * ACETYL
            if abs(dd) < 0.02:
                cands.append((abs(dd), n_ox, n_ac))
    if not cands:
        return "unknown:%+.4f" % d, d
    cands.sort()
    _, n_ox, n_ac = cands[0]
    parts = []
    if n_ox:
        parts.append("ox%d" % n_ox)
    if n_ac:
        parts.append("acetyl")
    return "+".join(parts) if parts else "unmod", d


def load_precursor_index():
    """从提取主表(hebei 子集)构建 (seq,z)->mz->PrecursorId 索引。"""
    if os.path.exists(MASTER_SUBSET):
        h = pd.read_parquet(MASTER_SUBSET)
    else:
        h = pd.read_csv(MASTER_TSV, sep="\t",
                        usecols=["Dataset", "PrecursorId", "StrippedSequence", "Charge", "PrecursorMz"])
        h = h[h.Dataset == "hebei"]
    idx = collections.defaultdict(dict)
    for pid, seq, z, mz in zip(h.PrecursorId.astype(str), h.StrippedSequence.astype(str),
                               h.Charge.astype(int), h.PrecursorMz.astype(float)):
        idx[(seq, z)][round(mz, 4)] = pid
    return idx


def main():
    wb = openpyxl.load_workbook(S14XLSX, read_only=True)
    rows = list(wb["Table S14"].iter_rows(values_only=True))[2:]
    cols = ["prot", "gene", "uc", "mtb", "seq", "z", "sample", "mz", "rt",
            "fim", "fmz", "qval", "qty", "shared", "nvar"]
    s14 = pd.DataFrame([r[:15] for r in rows], columns=cols).dropna(subset=["seq"])
    s14["z"] = s14["z"].astype(int)
    s14["mz"] = s14["mz"].astype(float)
    s14["rt"] = s14["rt"].astype(float)

    # 正文 Table 2 肽段标记
    import docx
    doc = docx.Document(DOCX)
    main_peps = set()
    for row in doc.tables[0].rows[1:]:
        for p in row.cells[6].text.replace("\n", ";").split(";"):
            p = p.strip()
            if p:
                main_peps.add(p)
    print("main-text table peptides:", len(main_peps))

    # hebei run 映射: sample(HCxx/TBxx) -> run_name（不取 raw_path）
    mp = pd.read_csv(MAPPING, sep="\t")
    mp["num"] = mp.run_name.str.extract(r"_ZhangYao_[HP]?_?(\d+)_ELU").astype(int)
    mp["sample"] = np.where(mp.group == "健康", "HC" + mp.num.astype(str), "TB" + mp.num.astype(str))
    smap = mp.set_index("sample")[["run_name", "group"]].to_dict("index")
    print("hebei runs -> samples:", len(smap))

    pid_idx = load_precursor_index()

    recs = []
    for _, r in s14.iterrows():
        st, dlt = resolve_mod_state(r.seq, r.z, r.mz)
        pid = None
        for mz4, p in pid_idx.get((r.seq, r.z), {}).items():
            if abs(mz4 - r.mz) / r.mz * 1e6 < 15:
                pid = p
                break
        if pid is None:
            ox_n = 0
            m = re.match(r"ox(\d+)", st)
            if m:
                ox_n = int(m.group(1))
            if ox_n or "acetyl" in st or st.startswith("unknown"):
                pid = "_%s_.%d(synthetic,mods=%s)" % (r.seq, r.z, st)
            else:
                pid = "_%s_.%d" % (r.seq, r.z)
        s = smap.get(r["sample"])
        recs.append(dict(
            prot=r.prot, gene=r.gene, uc=r.uc, mtb=r.mtb, seq=r.seq, z=r.z,
            mz_s14=r.mz, rt_s14=r.rt, fim=r.fim, fmz=r.fmz, qval=r.qval,
            qty=r.qty, shared=r.shared, nvar=r.nvar,
            mod_state=st, delta_mass=dlt, pid=pid,
            sample=r["sample"],
            run_name=s["run_name"] if s else None,
            group=s["group"] if s else None,
            in_main=bool(r.seq in main_peps),
        ))
    df = pd.DataFrame(recs)
    df.to_parquet(OUTP)

    tg = df.drop_duplicates(["seq", "z", "mz_s14"])
    grp = df.groupby(["seq", "z", "mz_s14"]).size()
    summ = dict(
        s14_rows=len(df), targets=len(tg), sequences=int(tg.seq.nunique()),
        proteins=int(tg.prot.nunique()), samples=int(df["sample"].nunique()),
        mod_state=tg.mod_state.value_counts().to_dict(),
        pid_from_report=int((~tg.pid.str.contains("synthetic")).sum()),
        pid_synthetic=int(tg.pid.str.contains("synthetic").sum()),
        no_run=int(df.run_name.isna().sum()),
        in_main_targets=int(tg.in_main.sum()),
        per_target_rows_min=int(grp.min()),
        per_target_rows_med=float(grp.median()),
        per_target_rows_max=int(grp.max()),
    )
    with _open(OUTJ, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1, default=str)
    print(json.dumps(summ, ensure_ascii=False, indent=1, default=str))


if __name__ == "__main__":
    main()
