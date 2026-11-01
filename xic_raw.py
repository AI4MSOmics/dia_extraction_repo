# -*- coding: utf-8 -*-
"""rawdump --xic 封装：RAW 直读 XIC（无 .ms1 文件时替代 xic_stream）。

口径与 eval_ms_pipeline.xic_stream 对齐：±0.02 Da 窗内最强峰、半高线性插值、
FWHM 4–240 s 门、窗口 ±2.5 min；apex 顶点带 edge 标记（edge=1 弃用）。
文件读写一律经 _open（root 白名单校验，PermissionError 显式拒绝越界路径）。
"""
import os

import pandas as pd

RAWDUMP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'rawdump', 'bin', 'Debug', 'netcoreapp3.1', 'rawdump.exe')


def _open(path, mode, root, **kw):
    """受 root 白名单限制的 open（同 denovo_qa_common.open_text 惯用法）。"""
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError('path escapes allowed root: ' + str(path))
    return open(resolved, mode, **kw)


def xic(datadir, stem, targets):
    """targets: [(mz, rt_min, i, z)] → {i: dict(apex=(rt,int)[, fwhm=…])}。"""
    stem = os.path.basename(str(stem))          # 仅取文件名，禁 ../
    rawp = os.path.join(datadir, stem + '.raw')
    tf = os.path.join(datadir, stem + '_xic_targets.tsv')
    of = os.path.join(datadir, stem + '_xic_out.tsv')
    if not os.path.isfile(rawp):
        return {}
    with _open(tf, 'w', root=datadir, encoding='utf-8') as f:
        f.write('id\tmz\tz\trt_min\ttol\n')
        for mz, rt, i, z in targets:
            f.write('%d\t%.5f\t%d\t%.4f\t0.02\n' % (i, mz, int(z), rt))
    argv = [RAWDUMP, rawp, '--xic', tf, of]
    try:
        os.spawnv(os.P_WAIT, RAWDUMP, argv)
    except OSError:
        return {}
    out = {}
    try:
        with _open(of, 'r', root=datadir) as fh:
            x = pd.read_csv(fh, sep='\t')
        for _, r in x.iterrows():
            d = {}
            if pd.notna(r.get('apex_rt_min')) and pd.notna(r.get('max_int')):
                d['apex'] = (float(r['apex_rt_min']), float(r['max_int']))
            w, edge = r.get('fwhm_s'), r.get('edge')
            if (pd.notna(w) and edge == 0 and int(r.get('n_pts', 0)) >= 4
                    and 4.0 <= float(w) <= 240.0):
                d['fwhm'] = float(w)
            if d:
                out[int(r['id'])] = d
    except (PermissionError, OSError, ValueError):
        pass
    return out
