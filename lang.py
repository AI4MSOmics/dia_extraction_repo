# -*- coding: utf-8 -*-
"""
lang.py — 语言常量与切换（默认英文 en，支持中文 zh）

用法：
  import os, lang
  L = lang.make(os.environ.get("DIA_LANG", "en"))
  ... L("Observed spectrum") ...            # zh 模式返回 "实测谱图"

切换：环境变量 DIA_LANG=zh（或 zh_CN）；默认/未设/en → 英文。
所有面向用户的文字（PDF 信息行/热图图题/封面/索引/审计报告）一律走本表；
代码注释与脚本内日志保持英文（不影响交付物语言）。
"""
import os

LANG = {
    # ---- PDF info lines (76/84 → 78 R topblock / 79 compose) ----
    "Observed spectrum": "实测谱图",
    "Best spectrum": "最优谱图",
    "matched b/y": "匹配 b/y",
    "ladder": "序列梯",
    "MS1 precursor envelope": "MS1 前体包络",
    "yes": "有",
    "not detected": "未检测到",
    "samples": "样本",
    "RT range": "RT 范围",
    "Table 1 peptide not covered by Table S14 - evidence via v9-master RT-anchored extraction":
        "Table 1 肽段未收录于 Table S14 —— 证据经 v9 主表 RT 锚定提取",
    "main-text table peptide": "正文表肽段",
    "RT anchor: v9 master (RTCorrected)": "RT 锚: v9 主表 (RTCorrected)",
    "mods": "修饰",
    "none": "无",
    "Protein position": "蛋白位点",
    "Precursor m/z": "前体 m/z",
    "charge": "电荷",
    "AlphaPeptDeep prediction (inverted)": "AlphaPeptDeep 理论谱（倒置）",
    "AlphaPeptDeep prediction unavailable for this precursor":
        "该前体无 AlphaPeptDeep 预测",
    "cosine(obs, pred)": "余弦(实测, 预测)",
    "shared 1+ ions": "共有 1+ 离子",
    "Base peak": "基峰",
    "Mass error of matched ions (observed spectrum)": "匹配离子质量差（实测谱）",

    # ---- Heatmap (79) ----
    "Fragment-ion XIC heatmap (1-min bins; white = 0, red = max; per-row cutoff = 20% of "
    "row max; grey value in blank row = Table S14 quantity of that sample)":
        "碎片离子 XIC 热图（1 min 分箱; 白=0, 红=最高; 逐行截止=行最高的 20%; "
        "白行中灰色数值 = 该样本 Table S14 定量值）",
    "apex RT (max-quantity bin)": "顶点 RT（最高定量分箱）",
    "mean ± SD": "均值 ± 标准差",
    "row(s) below quant threshold": "行低于定量阈值",
    "retention time (min)": "保留时间 (min)",
    "rel. XIC": "相对 XIC",

    # ---- Cover & index (79/81) ----
    "Evidence for Mtb-derived peptides (Table S14)": "Mtb 来源肽段证据 (Table S14)",
    "Evidence pages for peptides of Table 1 (Mtb-annotated proteins)":
        "Table 1 肽段证据页（Mtb 注释蛋白）",
    "Same isoPSV-layout page per peptide as the main S14 evidence PDF, "
    "ordered by Table 1 (protein No., peptide order)":
        "每肽段一页 isoPSV 版式（与主 S14 证据 PDF 相同），按 Table 1 顺序（蛋白编号、肽段序）",
    "peptide sequences": "肽段序列",
    "proteins": "蛋白",
    "samples": "样本",
    "hebei DIA cohort": "hebei DIA 队列",
    "peptides with extracted evidence pages": "有提取证据页的肽段",
    "peptides with Table S14 evidence": "有 Table S14 证据的肽段",
    "additional via v9-master RT-anchored extraction (not covered by Table S14)":
        "经 v9 主表 RT 锚定提取补充（未收录于 Table S14）",
    "without any spectrum": "无任何谱图",
    "Page content: info lines, ion ladder, observed spectrum, inverted AlphaPeptDeep "
    "prediction, ppm panel, RT heatmap":
        "页面内容：信息行、离子梯、实测谱图、倒置 AlphaPeptDeep 理论谱、ppm 面板、RT 热图",
    "Observed = RT-anchored targeted extraction, best sample by matched b/y ions":
        "实测谱 = RT 锚定靶向提取，按匹配 b/y 离子数取最优样本",
    "Prediction = AlphaPeptDeep 1.5.0 (inverted); cosine on shared 1+ b/y ions":
        "理论谱 = AlphaPeptDeep 1.5.0（倒置）; 共有 1+ b/y 离子余弦",
    "Top block rendered by the isoPSV R engine":
        "顶块由 isoPSV R 引擎渲染",
    "Index": "索引",
    "cont'd": "续",
    "Page": "页",
    "No." : "编号",
    "Protein": "蛋白",
    "Gene": "基因",
    "Peptide": "肽段",
    "Sequence": "序列",
    "Peptide(flanks)": "肽段(侧翼)",
    "Best PSM": "最优 PSM",
    "Status": "状态",
    "ok": "正常",
    "ok_extra": "补充链",
    "no_spectrum": "无谱图",
    "No extracted spectrum: this peptide is not present in "
    "Table S14 of the revision (no RT-anchored extraction was possible).":
        "无提取谱图：该肽段未收录于修回 Table S14（无法进行 RT 锚定提取）。",
    "Table 1 peptide not covered by Table S14":
        "Table 1 肽段未收录于 Table S14",

    # ---- Audit (85) ----
    "S14 / Table 1 evidence PDF automated audit report":
        "S14 / Table 1 证据 PDF 自动审计报告",
    "Errors": "错误",
    "Warnings": "警告",
    "Details": "明细",
    "Metrics": "指标",
    "Summary": "汇总",

    # ---- ppm / misc ----
    "ppm": "ppm",
    "m/z": "m/z",
    "Intensity (%)": "强度 (%)",
    "Predicted (%)": "预测 (%)",
}


def make(lang_code=None):
    """返回查找函数 L(text) → 翻译文本。lang_code 从 DIA_LANG 环境变量或参数取。"""
    code = (lang_code or os.environ.get("DIA_LANG", "en") or "en").lower()
    if code.startswith("zh"):
        def L(text):
            return LANG.get(text, text)
    else:
        def L(text):
            return text
    return L


def current():
    """当前生效语言代码。"""
    code = (os.environ.get("DIA_LANG", "en") or "en").lower()
    return "zh" if code.startswith("zh") else "en"
