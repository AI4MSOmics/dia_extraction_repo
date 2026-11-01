# -*- coding: utf-8 -*-
"""
73_profiles.py —— 每个"序列最优前体"在其全部 S14 样本中的全 RT 碎片 XIC 剖面（热图数据）

对每条肽段序列取 PSM 最优前体（脚本内按 nmatched 选，见 select_best），对该前体的
全部 S14 样本 run：在含该 m/z 的 DIA 窗口的全部 MS2 扫描上算理论 b/y（20 ppm）
匹配强度和，输出逐扫描 (rt, matched_sum, nmatched) 剖面；窗口级缓存 centroid，
同窗多前体共享抓谱。74 号脚本把剖面分箱成 1 min 热图。
输出：work/output/s14_evidence/profiles/{index.parquet, profiles.npz}
  index: seq/z/mz/sample/run/best_flag 行；npz: 每行对应的 rt/val/nm 数组（concat+偏移）。
复用 71_extract_targeted 的 RawIndex/open_raw/parse_pid（importlib 加载）。
安全: 路径常量；文件 open 经 _open() 白名单；npz 仅存纯数组（无 pickle）。
"""
import sys, os, json, time, importlib.util
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


TARGETS = os.environ.get("TGT", ROOT + r"\work\tmp\s14_targets.parquet")
ROWS = os.environ.get("ROWSX", ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet")
MAPPING = ROOT + r"\work\data_cache\hebei_run_raw_mapping.tsv"
OUTDIR = os.environ.get("OUTD", ROOT + r"\work\output\s14_evidence\profiles")

_spec = importlib.util.spec_from_file_location(
    "ext71", ROOT + r"\work\scripts\71_extract_targeted.py")
m71 = importlib.util.module_from_spec(_spec)
sys.modules["ext71"] = m71
_spec.loader.exec_module(m71)  # 只用其函数/类，main() 有 __main__ 守卫


def select_best(rows):
    """每序列选最优前体与最优谱：key=(nmatched, bLad+yLad, matched_int_frac)。"""
    out = {}
    for seq, g in rows.groupby("seq"):
        b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
        key = (b.seq, int(b.z), float(b.mz_s14))
        out[seq] = key
    return out


def main():
    limit_runs = 0
    for a in sys.argv[1:]:
        if a.startswith("--runs="):
            limit_runs = int(a.split("=", 1)[1])
    os.makedirs(OUTDIR, exist_ok=True)
    rows = pd.read_parquet(ROWS)
    tg = pd.read_parquet(TARGETS)
    best = select_best(rows)
    best_keys = set(best.values())
    # 需要剖面的 (前体, 样本) 对：该前体的全部 S14 行
    need = tg[tg.apply(lambda r: (r.seq, int(r.z), float(r.mz_s14)) in best_keys, axis=1)]
    mp = pd.read_csv(MAPPING, sep="\t")
    need = need.merge(mp[["run_name", "raw_path"]], on="run_name", how="left")
    need = need[need.raw_path.notna()]
    if limit_runs:
        keep = sorted(need.run_name.unique())[:limit_runs]
        need = need[need.run_name.isin(keep)]
    print("profile pairs:", len(need), "| runs:", need.run_name.nunique(), flush=True)

    # 理论碎片缓存
    frag_cache = {}
    for key in best_keys:
        pid = tg[(tg.seq == key[0]) & (tg.z == key[1]) & (tg.mz_s14 == key[2])].pid.iloc[0]
        try:
            frags, _, _ = m71.parse_pid(pid)
            frag_cache[key] = frags
        except Exception as e:
            print("pid fail", key, e, flush=True)

    idx_rows, rt_cat, val_cat, nm_cat, lens = [], [], [], [], []
    t0 = time.time()
    for run_name, g in need.groupby("run_name", sort=True):
        rawp = g.raw_path.iloc[0]
        if not os.path.exists(rawp):
            print("MISSING raw:", rawp, flush=True)
            continue
        rf = m71.open_raw(rawp)
        idx = m71.RawIndex(rf)
        print("[%6.1fs] %s  windows=%d  pairs=%d" %
              (time.time() - t0, run_name[:38], len(set(zip(idx.lo, idx.hi))), len(g)), flush=True)
        # 窗口级缓存 centroid 块（同窗目标共享一次抓谱）
        win_cache = {}
        for _, r in g.iterrows():
            key = (r.seq, int(r.z), float(r.mz_s14))
            frags = frag_cache.get(key)
            if frags is None:
                continue
            m = (idx.lo <= key[2]) & (idx.hi >= key[2])
            if not m.any():
                continue
            scans = idx.ms2_s[m]
            rts = idx.ms2_rt[m]
            w = (round(float(idx.lo[m][0]), 3), round(float(idx.hi[m][0]), 3))
            blk = win_cache.get(w)
            if blk is None:
                rt_l, val_by_sc = [], []
                for sc, rt in zip(scans, rts):
                    mz, inten = idx.centroid(sc)
                    rt_l.append(float(rt))
                    if mz is None or not len(mz):
                        val_by_sc.append((np.array([]), np.array([])))
                        continue
                    o = np.argsort(mz)
                    val_by_sc.append((mz[o], inten[o]))
                blk = (np.array(rt_l), val_by_sc)
                win_cache[w] = blk
            targets = np.array(sorted(frags.values()))
            # 1) S14 RT±0.6 内按序列梯定 apex（同 71 口径），取确认匹配离子子集；
            # 2) 逐离子全 RT XIC（±20ppm 首峰强度），各自按 apex 值归一并 cap 1，
            #    取离子间均值 → 色谱剖面（0–1）。抗干扰：同窗强共隔离本底的随机
            #    ppm 命中被 cap 压平，真峰处全部确认离子同涨。
            rts0 = blk[0]
            cand = np.where(np.abs(rts0 - float(r.rt_s14)) <= 0.6)[0]
            best = None
            for si in cand:
                mzs, ins = blk[1][si]
                if not len(mzs):
                    continue
                hit_by_label = {}
                for k, t in frags.items():
                    d = np.abs(mzs - t) / t * 1e6
                    kk = int(np.argmin(d))
                    if d[kk] < 20.0:
                        hit_by_label[k] = kk
                bL, yL = m71.ladder_metrics(list(hit_by_label.keys()))
                score = (bL + yL, len(hit_by_label))
                if best is None or score > best[0]:
                    best = (score, si, list(hit_by_label.keys()))
            vals = np.zeros(len(blk[0])); nms = np.zeros(len(blk[0]), dtype=int)
            if best is not None and len(best[2]) >= 4:
                sub_vals = np.array(sorted(frags[l] for l in best[2]))
                series = np.zeros((len(sub_vals), len(blk[0])))
                for si in range(len(blk[0])):
                    mzs, ins = blk[1][si]
                    if not len(mzs):
                        continue
                    lo = np.searchsorted(mzs, sub_vals * (1 - 10e-6), "left")
                    hi = np.searchsorted(mzs, sub_vals * (1 + 10e-6), "right")
                    hit = hi > lo
                    if hit.any():
                        series[hit, si] = ins[lo[hit]]
                a = series[:, best[1]]
                oki = a > 0
                if oki.sum() >= 4:
                    # 背景强度去除门槛：离子 i 在扫描 j 计入 pres/cos 仅当
                    # 强度 >= 10% × 该离子 apex 强度（背景/噪声首峰与真实碎片
                    # 丰度无关，过不了自身门槛 → 偏峰处 pres→0，热图回归窄带）
                    S = series[oki]
                    S = np.where(S >= (0.10 * a[oki])[:, None], S, 0.0)
                    num = a[oki] @ S
                    den = np.linalg.norm(S, axis=0) * np.linalg.norm(a[oki])
                    cos = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
                    pres = (S > 0).mean(axis=0)
                    vals = np.clip(cos * pres, 0.0, 1.0)
                    nms = int(oki.sum()) * np.ones(len(blk[0]), dtype=int)
            if not vals.any() and len(blk[0]):
                # 兜底：无合格 apex 时退全离子和
                for si in range(len(blk[0])):
                    mzs, ins = blk[1][si]
                    if not len(mzs):
                        continue
                    lo = np.searchsorted(mzs, targets * (1 - 20e-6), "left")
                    hi = np.searchsorted(mzs, targets * (1 + 20e-6), "right")
                    hit = hi > lo
                    if hit.any():
                        vals[si] = float(ins[lo[hit]].sum())
                        nms[si] = int(hit.sum())
            idx_rows.append(dict(seq=key[0], z=key[1], mz=key[2], sample=r["sample"],
                                 run=run_name, n_scans=len(blk[0])))
            rt_cat.append(blk[0]); val_cat.append(vals); nm_cat.append(nms)
            lens.append(len(blk[0]))
        rf.Dispose()

    dfe = pd.DataFrame(idx_rows)
    dfe.to_parquet(OUTDIR + r"\index.parquet")
    ends = np.cumsum(lens); offs = ends - np.array(lens)
    np.savez(OUTDIR + r"\profiles.npz",
             rt=np.concatenate(rt_cat) if rt_cat else np.array([]),
             val=np.concatenate(val_cat) if val_cat else np.array([]),
             nm=np.concatenate(nm_cat) if nm_cat else np.array([]),
             off=offs, end=ends)
    print(json.dumps(dict(pairs=len(dfe), runs=int(dfe.run.nunique()),
                          elapsed_s=round(time.time() - t0, 1)), ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
