#!/usr/bin/env python3
"""
compute_tau.py — 用 GTEx v8 median TPM 计算每基因组织特异性 tau 分数（学生跟跑版扩展）

tau = Σ_i (1 - x_i / max(x)) / (n - 1)，Yanai et al. 2004
  tau→1: 表达高度组织特异；tau→0: 广泛表达

输出 data/features/gene_tau.csv: gene_id, symbol, tau, max_tpm_tissue, max_tpm,
      cortex_tpm, cortex_specific (TPM > 5× 全组织均值)
"""
import pandas as pd
import numpy as np

GCT = "./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
OUT = "./data/features/gene_tau_cortex.csv"      #文件输出

df = pd.read_csv(GCT, sep="\t", skiprows=2)
tissues = [c for c in df.columns if c not in ("Name", "Description")]  
# 查看54个组织名称
print(tissues)
print(f"genes: {len(df)}, tissues: {len(tissues)}")

df["gene_id"] = df["Name"].str.split(".").str[0]      # 去版本号，统一主键
X = df[tissues].clip(lower=0)  
mx = X.max(axis=1)
tau = (1 - X.div(mx.replace(0, np.nan), axis=0)).sum(axis=1) / (len(tissues) - 1)   #Tau公式

out = pd.DataFrame({
    "gene_id": df["gene_id"], "symbol": df["Description"], "tau": tau.round(4),
    "max_tpm_tissue": X.idxmax(axis=1).where(mx > 0),
    "max_tpm": mx.round(2),
    "cortex_tpm": X["Brain - Cortex"].round(2),   #组织替换
})
out["cortex_specific"] = out["cortex_tpm"] > 5 * X.mean(axis=1)
# 去版本号后可能出现重复（如 PAR 区基因副本）：保留 max_tpm 最大的一行
out = out.sort_values("max_tpm", ascending=False).drop_duplicates("gene_id")  #排序
out = out.sort_values("gene_id")
assert out["gene_id"].is_unique   #确保没有重复
out.to_csv(OUT, index=False)

print("\ntau 分布:")
print(out["tau"].describe().round(3).to_string())
print("\n高组织特异基因（tau>0.95）数量:", (out["tau"] > 0.95).sum())
print("脑特异基因数量:", out["cortex_specific"].sum())
cb = out[(out["cortex_specific"])].nlargest(10, "cortex_tpm")
print("\nTPM 最高的脑特异基因:\n", cb[["symbol", "tau", "cortex_tpm", "max_tpm_tissue"]].to_string(index=False))
print(f"\n→ {OUT}")
