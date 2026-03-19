
output_dir= '/home/zhangfeng/projects/QRLu/data/Xenium_14/'


##########

import pandas as pd
import anndata as ad  # 若使用AnnData格式

# 读取文件
cells = pd.read_parquet(output_dir+'/'+"cells.parquet")  # 细胞元信息（cell_id、centroid_x/y等）
transcripts = pd.read_parquet(output_dir+'/'+"transcripts.parquet")  # 转录本（cell_id、gene、x/y等）

# 查看关键列（根据实际数据调整列名）
print("cells列名:", cells.columns)
print("transcripts列名:", transcripts.columns)


transcripts_filtered = transcripts.dropna(subset=["cell_id"])

# 统计每个细胞-基因的转录本数量（构建count矩阵）
count_matrix = (
    transcripts_filtered
    .groupby(["cell_id", "feature_name"])  # 按细胞和基因分组
    .size()  # 统计数量
    .unstack(fill_value=0)  # 转换为矩阵（行：cell_id，列：gene）
)

# 确保count_matrix的行是cell_id，列是基因名
count_matrix.head()



cell_metadata = cells.set_index("cell_id")[["x_centroid", "y_centroid"]]  # 空间坐标
# 可选：添加其他元数据（如细胞面积、类型等）
# cell_metadata["area"] = cells.set_index("cell_id")["area"]

# 确保坐标与count_matrix的细胞ID一致（按count_matrix的行排序）
cell_metadata = cell_metadata.reindex(count_matrix.index)



count_matrix.to_csv(output_dir+'/'+"xenium_count_matrix.csv")

# 导出细胞坐标（行：cell_id，列：x/y）
cell_metadata.to_csv(output_dir+'/'+"xenium_cell_metadata.csv")



#adata = ad.AnnData(
#    X=count_matrix.values,  # 表达矩阵（数值部分）
#    obs=cell_metadata,  # 细胞元数据（包含空间坐标）
#    var=pd.DataFrame(index=count_matrix.columns)  # 基因信息（此处仅保留基因名）
#)

# 添加空间坐标到obsm（Seurat会识别"spatial"键）
#adata.obsm["spatial"] = cell_metadata[["centroid_x", "centroid_y"]].values

# 保存为AnnData格式（.h5ad）
#adata.write("xenium_data.h5ad")






