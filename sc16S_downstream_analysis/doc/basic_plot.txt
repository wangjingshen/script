结果说明：
featureplot_total_genus：genus 的总 UMI 的 UMAP 图；
featureplot_total_genus_group：genus 的总 UMI 按组拆分的 UMAP 图；
featureplot_genus_detect：检测到 genus 的 UMAP 图，其中，检测到 genus 的总 UMI 数 > 2 定义为 genus+；反之，定义为 genus-；
featureplot_genus_detect_group：检测到 genus 按组拆分的 UMAP 图；
barplot_total_genus_umi: 各细胞类型各组检测到 genus 的总UMI数的 条形图；
barplot_total_genus_umi_group: 各组检测到 genus 的总UMI数的 条形图；
featureplot_top10_genus_*: *数据中检出细胞数前 10 的 genus 的 UMAP 图（*为 Main：整体数据 or groupX：X组数据）；
top10_genus_split：*数据中检出细胞数前 10 的 genus 的 UMAP 图拆分版，其中数字后缀为检出细胞数的降序排序号；
barplot_genus_count_sample_cluster：横坐标为样本+细胞类型，纵坐标为检测到的细菌属个数；
barplot_genus_umi_sample_cluster：横坐标为样本+细胞类型，纵坐标为检测到的细菌属 UMI 总数；
barplot_plot_data_cluster：作图用到的数据；