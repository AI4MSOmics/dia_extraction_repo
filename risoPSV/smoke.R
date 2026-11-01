# isoPSV 引擎冒烟：理论 b/y 质量与 Python 侧对账
source("compat.R", encoding = "UTF-8")
source("ReadMGF_script.R", encoding = "UTF-8")
setwd("F:/IPX0018588000/work/risoPSV")
aa <- read.delim("AA_mass_dataframe.txt", sep = "\t", header = TRUE,
                 check.names = FALSE, stringsAsFactors = FALSE)
b <- by_ions_position("FTPSGAAVANFTVASTPR", "", "Unlabel", aa, amino_acid_table)
cat("cols:", paste(names(b), collapse = ","), "\n")
cat("nrow:", nrow(b), "\n")
print(head(b, 6))
# 带修饰：SQTVFMGDFPMMTEKGTFIINGTER 无修饰版本 + 一个 CAM 例
b2 <- by_ions_position("ALCTVATLR", "2,Carbamidomethyl[C]", "Unlabel", aa, amino_acid_table)
print(head(b2, 4))
cat("SMOKE_OK\n")
