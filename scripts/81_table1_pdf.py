# -*- coding: utf-8 -*-
"""
81_table1_pdf.py —— 正文 Table 1 对应的逐肽证据 PDF

读取 F:\\IPX0018588000\\Table 1.xlsx（53 蛋白，Peptide sequence 列分号分隔多肽），
按 Table 1 顺序（No. → 蛋白内肽段顺序）为每条肽段出一页与主交付同版式的 A4
（信息行 + isoPSV R 顶块 + 热图）；Table 1 中无提取谱的肽段出占位页并在索引标注。
输出：work/output/s14_evidence/Table1_peptide_evidence_A4.pdf + table1_manifest.tsv
复用 79_compose_a4 的 heatmap_png / fit_text（importlib 加载）。
安全：路径常量；文件 open 经 _open() 白名单（校验与打开同函数）。
"""
import sys, os, json, math, importlib.util
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pymupdf
import openpyxl

ROOT = r"F:\IPX0018588000"


def _open(path, mode="r", root=ROOT, **kw):
    """受 root 白名单限制的 open（同 xic_raw._open 惯用法）：校验与打开同函数。"""
    resolved = os.path.normpath(os.path.abspath(path))
    root_r = os.path.normpath(os.path.abspath(root))
    if not (resolved == root_r or resolved.startswith(root_r + os.sep)):
        raise PermissionError("path escapes allowed root: " + str(path))
    return open(resolved, mode, **kw)


_spec = importlib.util.spec_from_file_location("e79", ROOT + r"\work\scripts\79_compose_a4.py")
m79 = importlib.util.module_from_spec(_spec)
sys.modules["e79"] = m79
_spec.loader.exec_module(m79)

T1XLSX = ROOT + r"\Table 1.xlsx"
INTSV = ROOT + r"\work\risoPSV\in\s14_spectra.tsv"
ROWS = ROOT + r"\work\output\s14_evidence\extract\rows_all.parquet"
TOPPDF = ROOT + r"\work\risoPSV\out\s14_topblock.pdf"
OUTPDF = ROOT + r"\work\output\s14_evidence\Table1_peptide_evidence_A4.pdf"
OUTMANI = ROOT + r"\work\output\s14_evidence\table1_manifest.tsv"
OUTSUMM = ROOT + r"\work\output\s14_evidence\table1_summary.json"


def parse_table1():
    wb = openpyxl.load_workbook(T1XLSX, read_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))[2:]
    out = []
    for r in rows:
        if r[0] is None or r[1] is None:
            continue
        peps = []
        if r[6]:
            for p in str(r[6]).replace("\n", ";").split(";"):
                p = p.strip()
                if p:
                    peps.append(p)
        out.append(dict(no=str(r[0]).strip(), prot=str(r[1]).strip(),
                        gene=str(r[2]).strip() if r[2] else "",
                        peptides=peps))
    return out


def main():
    t1 = parse_table1()
    print("Table1 proteins:", len(t1))
    rows = pd.read_parquet(ROWS).reset_index(drop=True)
    ok = rows[rows.status == "ok"]
    seqs = list(ok.seq.unique())          # 与 78/79 完全同序 → topblock 页索引
    seq2page = {s: i for i, s in enumerate(seqs)}
    itab = pd.read_csv(INTSV, sep="\t", encoding="utf-8")
    info_map = {r.Sequence: (r.Info1, r.Info2, r.Info3)
                for r in itab.drop_duplicates("Sequence").itertuples()}
    topdoc = pymupdf.open(TOPPDF)
    # 补充链（Table 1 有而 S14 无的肽）：extract_extra + topblock_extra + profiles_extra
    EX_ROWS = ROOT + r"\work\output\s14_evidence\extract_extra\rows_all.parquet"
    EX_TOP = ROOT + r"\work\risoPSV\out\topblock_extra.pdf"
    EX_INTSV = ROOT + r"\work\risoPSV\in_extra\s14_spectra.tsv"
    EX_PROF_I = ROOT + r"\work\output\s14_evidence\profiles_extra\index.parquet"
    EX_PROF_N = ROOT + r"\work\output\s14_evidence\profiles_extra\profiles.npz"
    ex_ok = None
    ex_seq2page = {}
    ex_topdoc = None
    ex_info = {}
    if os.path.exists(EX_ROWS) and os.path.exists(EX_TOP):
        exr = pd.read_parquet(EX_ROWS).reset_index(drop=True)
        ex_ok = exr[exr.status == "ok"]
        ex_seqs = list(ex_ok.seq.unique())   # 与 84→78(in_extra) 同序
        ex_seq2page = {s: i for i, s in enumerate(ex_seqs)}
        ex_topdoc = pymupdf.open(EX_TOP)
        eit = pd.read_csv(EX_INTSV, sep="\t", encoding="utf-8")
        ex_info = {r.Sequence: (r.Info1, r.Info2, r.Info3)
                   for r in eit.drop_duplicates("Sequence").itertuples()}
        print("extra chain ready:", len(ex_seqs), "seqs")
    tmpdir = ROOT + r"\work\output\s14_evidence\_heat_tmp"
    os.makedirs(tmpdir, exist_ok=True)

    def compose(page, rec, pep, src_ok, pi, topd, info3, gk_row, gk, prof_i, prof_n, idx_pi):
        i1, i2, i3 = info3
        y = m79.fit_text(page, 14, 26,
                         "Table 1 No.%s %s (%s)  |  " % (rec["no"], rec["prot"], rec["gene"]) + i1,
                         570, 10.5, bold=True)
        y = m79.fit_text(page, 14, y + 1, i2, 570, 8)
        y = m79.fit_text(page, 14, y + 1, i3, 570, 8)
        r = pymupdf.Rect(2, y + 2, 593, 442)
        page.show_pdf_page(r, topd, pi)
        hp = os.path.join(tmpdir, "t1_%03d.png" % idx_pi)
        b = gk_row
        has = m79.heatmap_png(pep, b.z, b.mz_s14, gk, hp, prof_i, prof_n)
        if has:
            page.insert_image(pymupdf.Rect(20, 431, 575, 838), filename=hp)
        else:
            page.insert_text((30, 600), "no profiles", fontsize=10)

    doc = pymupdf.open()
    mani = []
    idx_pi = 0
    for rec in t1:
        for pep in rec["peptides"]:
            idx_pi += 1
            page = doc.new_page(width=595, height=842)
            if pep in seq2page:
                pi = seq2page[pep]
                g = ok[ok.seq == pep]
                b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
                gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
                compose(page, rec, pep, ok, pi, topdoc, info_map.get(pep, ("", "", "")),
                        b, gk, None, None, idx_pi)
                mani.append(dict(page_no=idx_pi + 1, t1_no=rec["no"], prot=rec["prot"],
                                  gene=rec["gene"], seq=pep, status="ok"))
            elif pep in ex_seq2page:
                g = ex_ok[ex_ok.seq == pep]
                b = g.sort_values(["nmatched", "bLad", "yLad"], ascending=False).iloc[0]
                gk = g[(g.z == b.z) & (g.mz_s14 == b.mz_s14)]
                compose(page, rec, pep, ex_ok, ex_seq2page[pep], ex_topdoc,
                        ex_info.get(pep, ("", "", "")), b, gk, EX_PROF_I, EX_PROF_N, idx_pi)
                mani.append(dict(page_no=idx_pi + 1, t1_no=rec["no"], prot=rec["prot"],
                                  gene=rec["gene"], seq=pep, status="ok_extra"))
            else:
                page.insert_text((30, 60), "Table 1 No.%s %s (%s)" % (rec["no"], rec["prot"], rec["gene"]),
                                 fontsize=12, fontname="hebo")
                page.insert_text((30, 90), "Peptide: %s" % pep, fontsize=11, fontname="helv")
                page.insert_text((30, 120), "No extracted spectrum: this peptide is not present in "
                                 "Table S14 of the revision (no RT-anchored extraction was possible).",
                                 fontsize=9, fontname="helv")
                mani.append(dict(page_no=idx_pi + 1, t1_no=rec["no"], prot=rec["prot"],
                                  gene=rec["gene"], seq=pep, status="no_spectrum"))

    # 封面
    cov = doc.new_page(width=595, height=842)
    fig = plt.figure(figsize=(8.27, 11.69))
    n_ok = sum(1 for m in mani if m["status"] == "ok")
    n_ex = sum(1 for m in mani if m["status"] == "ok_extra")
    fig.text(0.5, 0.60, "Evidence pages for peptides of Table 1 (Mtb-annotated proteins)",
             ha="center", fontsize=14, weight="bold", family="Arial")
    fig.text(0.5, 0.555, "Same isoPSV-layout page per peptide as the main S14 evidence PDF, "
             "ordered by Table 1 (protein No., peptide order)", ha="center", fontsize=9)
    for i, s in enumerate([
            "%d proteins / %d peptides from Table 1" % (len(t1), len(mani)),
            "%d peptides with Table S14 evidence; %d additional via v9-master RT-anchored "
            "extraction (not covered by Table S14); %d without any spectrum" %
            (n_ok, n_ex, len(mani) - n_ok - n_ex),
            "Page content: info lines, ion ladder, observed spectrum, inverted AlphaPeptDeep "
            "prediction, ppm panel, RT heatmap"]):
        fig.text(0.12, 0.49 - i * 0.028, s, fontsize=9)
    covpng = os.path.join(tmpdir, "_t1_cover.png")
    fig.savefig(covpng, dpi=150); plt.close(fig)
    cov.insert_image(pymupdf.Rect(0, 0, 595, 842), filename=covpng)
    doc.move_page(len(doc) - 1, 0)
    # 索引
    chunk = 52
    npg = math.ceil(len(mani) / chunk)
    for ci in range(npg):
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.965, "Index%s" % ("" if ci == 0 else " (cont'd %d/%d)" % (ci + 1, npg)),
                 fontsize=13, weight="bold", va="top", family="Arial")
        ax = fig.add_axes([0.04, 0.03, 0.94, 0.92]); ax.axis("off")
        sub = mani[ci * chunk:(ci + 1) * chunk]
        body = [[m["page_no"], m["t1_no"], m["prot"], m["gene"], m["seq"][:30],
                 m["status"]] for m in sub]
        tb = ax.table(cellText=body,
                      colLabels=["Page", "No.", "Protein", "Gene", "Peptide", "Status"],
                      loc="upper center", cellLoc="center", bbox=[0, 0, 1, 0.94],
                      colWidths=[0.07, 0.06, 0.13, 0.12, 0.44, 0.18])
        tb.auto_set_font_size(False); tb.set_fontsize(5)
        for (r_i, c_i), cell in tb.get_celld().items():
            cell.set_linewidth(0.2)
            if r_i == 0:
                cell.set_facecolor("#e8e8e8")
        ip = os.path.join(tmpdir, "_t1_idx%d.png" % ci)
        fig.savefig(ip, dpi=150); plt.close(fig)
        pg = doc.new_page(width=595, height=842)
        pg.insert_image(pymupdf.Rect(0, 0, 595, 842), filename=ip)

    doc.save(OUTPDF, garbage=3, deflate=True)
    pd.DataFrame(mani).to_csv(OUTMANI, sep="\t", index=False, encoding="utf-8")
    summ = dict(proteins=len(t1), peptides=len(mani), ok=n_ok, ok_extra=n_ex,
                no_spectrum=len(mani) - n_ok - n_ex, total_pages=doc.page_count)
    with _open(OUTSUMM, "w", encoding="utf-8") as f:
        json.dump(summ, f, ensure_ascii=False, indent=1)
    print(json.dumps(summ, ensure_ascii=False))


if __name__ == "__main__":
    main()
