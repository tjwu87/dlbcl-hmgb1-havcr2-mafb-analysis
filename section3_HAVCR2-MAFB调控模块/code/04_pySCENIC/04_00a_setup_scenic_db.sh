#!/usr/bin/env bash
# =============================================================================
#  pySCENIC 参考数据库下载（论文 section3 上游）
#
#  来源：Biomni 平台分析记录（19_PAGA与注释重做）
#  三个文件均来自 Aerts Lab 官方资源站，约 400 MB，下载后放到 $SCENIC_DB。
#
#  用法：
#      export SCENIC_DB=/path/to/scenic_db      # 默认 <数据根>/external/scenic_db
#      bash 04_00a_setup_scenic_db.sh
#
#  下载完成后，可直接用官方 CLI 跑三步（等价于 04_00c_pyscenic_run.py）：
#
#      pyscenic grn  $SCENIC_DB/mono_ifntam_latam.loom $SCENIC_DB/allTFs_hg38.txt \
#          -o $SCENIC_DB/adjacencies.csv -m grnboost2 --num_workers 4 --seed 42
#
#      pyscenic ctx  $SCENIC_DB/adjacencies.csv \
#          $SCENIC_DB/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
#          --annotations_fname $SCENIC_DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
#          --expression_mtx_fname $SCENIC_DB/mono_ifntam_latam.loom \
#          --mode dask_multiprocessing --output $SCENIC_DB/regulons.csv \
#          --num_workers 4 --mask_dropouts
#
#      pyscenic aucell $SCENIC_DB/mono_ifntam_latam.loom $SCENIC_DB/regulons.csv \
#          --output $SCENIC_DB/auc_matrix.csv --num_workers 4
#
#  注：本论文实际采用 Python API 形式（04_00c），原因是 CLI 的 dask 后端在
#      容器里不稳定（ctx 长时间 0 字节输出）；API 版用 custom_multiprocessing
#      绕开 dask，结果一致。两条路径的阈值参数完全相同。
# =============================================================================
set -euo pipefail

SCENIC_DB="${SCENIC_DB:-./scenic_db}"
mkdir -p "$SCENIC_DB"

CISTARGET="https://resources.aertslab.org/cistarget"

# ① 人类 TF 列表（1 892 个 TF）
curl -sSL "$CISTARGET/tf_lists/allTFs_hg38.txt" \
     -o "$SCENIC_DB/allTFs_hg38.txt"

# ② hg38 motif 排名数据库（~300 MB）
curl -L "$CISTARGET/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/\
hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
     -o "$SCENIC_DB/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
     --progress-bar

# ③ motif → TF 注释表（~95 MB）
curl -L "$CISTARGET/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
     -o "$SCENIC_DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
     --progress-bar

echo "--- 校验 ---"
wc -l "$SCENIC_DB/allTFs_hg38.txt"
ls -lh "$SCENIC_DB"/*.feather "$SCENIC_DB"/*.tbl
echo "[OK] 数据库就绪 -> $SCENIC_DB"
