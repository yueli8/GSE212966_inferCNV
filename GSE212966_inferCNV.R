# https://zhuanlan.zhihu.com/p/625589597
#原理大意就是用染色体上相邻的101个基因的平均表达值表示中间位置基因的CNV值

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")

#​​​​​
#还需要安装JAGS4.3.1（R.Version()>4.2.0）
#https://sourceforge.net/projects/mcmc-jags/

library(Seurat)
library(infercnv)
library(ggplot2)
library(ggpubr)
library(tidyverse)
dir.create("./data/output/inferCNV")

#0 seurat_object####
##0.1 seurat_refer####
hms_cluster_id <- readRDS("./data/temp/hms_cluster_id_test1.3.rds")
p1 <- DimPlot(hms_cluster_id, reduction = "umap", label = TRUE)
p1
levels(hms_cluster_id)
#labels存放的自动注释的分类名，这里替换掉，细胞分群名，其他分群也这样修改⭐
hms_cluster_id$labels = as.character(Idents(hms_cluster_id))
table(hms_cluster_id@meta.data$labels) 

#取T cells, macrophages, endothelial cells, and stellate cells为参考细胞
idents_id <- c("Macrophage", "T_cell","Endothelial_cell","Stellate_cell")
seurat_refer <- subset(hms_cluster_id, idents = idents_id)

#⭐下采样200，试下运行时间
seurat_refer <- subset(seurat_refer, downsample = 200)#这里下采样减少了很多数据,各亚群取200个细胞，存疑

p1 <- DimPlot(seurat_refer, reduction = "umap", label = TRUE)
p1

##0.2 seurat_target####
#Duct_epithelial_cell inferCNV对象
###1. 采用分群好的作为对象####
seurat_target <- readRDS('./data/temp/Duct_epithelial_cell_cluster_id_test.rds')
levels(seurat_target)
#labels存放的自动注释的分类名，这里替换掉，细胞分群名，其他分群也这样修改⭐
seurat_target$labels = as.character(Idents(seurat_target))

###2. 采用未分群的作为对象####

idents_id <- c("Duct_epithelial_cell")
seurat_target <- subset(hms_cluster_id, idents = idents_id)
table(seurat_target@meta.data$labels) 
p1 <- DimPlot(seurat_target, reduction = "umap", label = TRUE)
p1

#⭐筛选分析需要的细胞类型（备用）
# idents_id <- c("Group_1", "Group_2","Group_3","Group_4","Group_5","Group_6","Group_7","Group_8","Group_9")
# seurat_target <- subset(seurat_target, idents = idents_id)

#⭐下采样200，试下运行时间
seurat_target <- subset(seurat_target, downsample = 200)#这里下采样减少了很多数据,各亚群取200个细胞，存疑
# 

##0.3 seurat_object####
seurat_object <- merge(seurat_target,seurat_refer)
levels(seurat_object)
seurat_object$labels = as.character(Idents(seurat_object))
saveRDS(seurat_object, file = "./data/temp/seurat_object_refer+target-分群下采样200.rds")
# seurat_object<- readRDS("./data/temp/seurat_object_refer+target-分群下采样200.rds")

#这里没办法umap图
# p1 <- DimPlot(seurat_object, reduction = "umap", label = TRUE)
# p1


#1 文件准备####
seurat_object <- readRDS('./data/temp/seurat_object_refer+target.rds')
## 1.1 expr_matrix表达矩阵文件####
expr_matrix  <-as.matrix(seurat_object@assays$RNA@counts)
head(expr_matrix)

## 1.2 cell_annotation细胞类型注释文件####
cell_annotation  <-data.frame(seurat_object$labels)
head(cell_annotation)

## 1.3 gene_position基因定位文件####
gene_position <-read.csv("./data/input/hg38_gencode_v27.txt",header=F,row.names=1,sep = "\t")
head(gene_position)
# CreateInfercnvObject报错：Error in FUN(left, right) : 二进列运算符中有非数值参数，header=F,row.names=1没有设置原因
#https://www.jianshu.com/p/268463268851

# 2.创建inferCNV对象####
infercnv_obj = CreateInfercnvObject(delim = "\t",
                                    #min_max_counts_per_cell = c(100, +Inf),
                                    raw_counts_matrix = expr_matrix,
                                    annotations_file = cell_annotation,
                                    gene_order_file = gene_position,
                                    ref_group_names = c("Macrophage", "T_cell","Endothelial_cell","Stellate_cell"))#参考细胞，优先选择免疫细胞作为参考
#问题3在这里就要思考一下了，参考细胞必须提供吗？不是。不提供并不会影响整个程序的运行，程序会计算所有基因的表达均值作为参考细胞的表达矩阵，但是分析的合理性就要打个问号了。
saveRDS(infercnv_obj, file = "./data/temp/seurat_obj_refer+target-分群下采样200.rds")
# infercnv_obj<- readRDS("./data/temp/seurat_obj_refer+target-分群下采样200.rds")

#3.运行inferCNV####

##inferCNV1####
infercnv_obj <- readRDS('./data/temp/seurat_obj_refer+target.rds')
dir.create("./data/output/inferCNV1")
infercnv_obj = infercnv::run( infercnv_obj,
                              cutoff = 0.1,#cutoff：默认0.1，过滤低表达基因。Smart-seq2数据选“1”, 10x Genomics数据选“0.1”
                              #min_cells_per_gene：默认3，过滤表达比例较低基因
                              #window_length：默认101，滑窗中包含的基因数目
                              out_dir = "./data/output/inferCNV1",#输出结果位置
                              
                              cluster_by_groups = T,#选择"TRUE"时，先对分析细胞类型各自聚类，区分细胞来源，再做层次聚类。
                                                    #选择"FALSE"时，对分析细胞类型整体聚类，并按照的数值进行聚类。
                                                    #默认"FALSE"，对肿瘤细胞进行聚类，展示为不同的亚群
                              #k_obs_groups = 8,#自定义分析细胞类型分群数，与为“FALSE"时联合使用。
                              #k_obs_groups：默认1，肿瘤细胞聚类数目
                              
                              # write_exper_matrix=T,
                              
                              HMM = T,#是否使用隐马尔可夫模型(Hidden Markov Model,,HMM)进行预测，默认使用。
                              #analysis_mode=samples#使用HMM进行预测的分析模式，包括samples,subclusters和cells.
                              #samples:默认，对每个样本进行分析，不考虑样本内的细胞差异。运行速度快，但是可能忽略一些细胞的拷贝数变异。
                              #subclusters:推荐，先对每个样本内的细胞进行聚类，再对每个聚类进行分析，考虑细胞之间的差异。能够更好地检测细胞的拷贝数变异，但是运行速度较慢。
                              #cells:对每个细胞进行分析，不进行聚类。能够得到最详细的拷贝数变异结果，但是运行速度最慢。
                              denoise = TRUE,#默认FALSE(更快)，对CNV矩阵进行降噪
                              num_threads = 10)#线程数

saveRDS(infercnv_obj,file="./data/temp/inferCNV_infercnv_obj_Duct_epithelial_cell1.rds")



##inferCNV1 全采样####
dir.create("./data/output/inferCNV1")
infercnv_obj <- readRDS('./data/temp/seurat_obj_refer+target.rds')
infercnv_obj = infercnv::run( infercnv_obj,
                              cutoff = 0.1,                           #cutoff：默认0.1，过滤低表达基因。Smart-seq2数据选“1”, 10x Genomics数据选“0.1”
                              #min_cells_per_gene=,                   #默认3，过滤表达比例较低基因
                              #window_length=,                        #默认101，滑窗中包含的基因数目
                              out_dir = "./data/output/inferCNV1",    #输出结果位置
                              cluster_by_groups = T,                  #选择"TRUE"时，先对分析细胞类型各自聚类，区分细胞来源，再做层次聚类。
                              #选择"FALSE"时，对分析细胞类型整体聚类，并按照的数值进行聚类。
                              #默认"FALSE"，对肿瘤细胞进行聚类，展示为不同的亚群
                              #k_obs_groups = 8,                      #自定义分析细胞类型分群数，与为“FALSE"时联合使用。#k_obs_groups：默认1，肿瘤细胞聚类数目
                              HMM = T,                                #是否使用隐马尔可夫模型(Hidden Markov Model,,HMM)进行预测，默认使用。
                              analysis_mode="subclusters",            #使用HMM进行预测的分析模式，包括samples,subclusters和cells.
                              #samples:默认，对每个样本进行分析，不考虑样本内的细胞差异。运行速度快，但是可能忽略一些细胞的拷贝数变异。
                              #subclusters:推荐，先对每个样本内的细胞进行聚类，再对每个聚类进行分析，考虑细胞之间的差异。能够更好地检测细胞的拷贝数变异，但是运行速度较慢。
                              #cells:对每个细胞进行分析，不进行聚类。能够得到最详细的拷贝数变异结果，但是运行速度最慢。
                              denoise = T,                            #默认FALSE(更快)，对CNV矩阵进行降噪
                              num_threads = 10,                       #线程数
                              write_expr_matrix = T,                  #输出矩阵
                              output_format = "pdf",
                              save_final_rds = T
)


saveRDS(infercnv_obj,file="./data/temp/inferCNV_infercnv_obj_Duct_epithelial_cell2.rds")

##inferCNV2下采样200####
dir.create("./data/output/inferCNV2")
infercnv_obj <- readRDS('./data/temp/seurat_object_refer+target-分群下采样200')
infercnv_obj = infercnv::run( infercnv_obj,
                              cutoff = 0.1,                           #cutoff：默认0.1，过滤低表达基因。Smart-seq2数据选“1”, 10x Genomics数据选“0.1”
                              #min_cells_per_gene=,                   #默认3，过滤表达比例较低基因
                              #window_length=,                        #默认101，滑窗中包含的基因数目
                              out_dir = "./data/output/inferCNV2",    #输出结果位置
                              cluster_by_groups = T,                  #选择"TRUE"时，先对分析细胞类型各自聚类，区分细胞来源，再做层次聚类。
                                                                      #选择"FALSE"时，对分析细胞类型整体聚类，并按照的数值进行聚类。
                                                                      #默认"FALSE"，对肿瘤细胞进行聚类，展示为不同的亚群
                              #k_obs_groups = 8,                      #自定义分析细胞类型分群数，与为“FALSE"时联合使用。#k_obs_groups：默认1，肿瘤细胞聚类数目
                              HMM = T,                                #是否使用隐马尔可夫模型(Hidden Markov Model,,HMM)进行预测，默认使用。
                              analysis_mode="subclusters",            #使用HMM进行预测的分析模式，包括samples,subclusters和cells.
                                                                      #samples:默认，对每个样本进行分析，不考虑样本内的细胞差异。运行速度快，但是可能忽略一些细胞的拷贝数变异。
                                                                      #subclusters:推荐，先对每个样本内的细胞进行聚类，再对每个聚类进行分析，考虑细胞之间的差异。能够更好地检测细胞的拷贝数变异，但是运行速度较慢。
                                                                      #cells:对每个细胞进行分析，不进行聚类。能够得到最详细的拷贝数变异结果，但是运行速度最慢。
                              denoise = T,                            #默认FALSE(更快)，对CNV矩阵进行降噪
                              num_threads = 10,                       #线程数
                              write_expr_matrix = T,                  #输出矩阵
                              output_format = "pdf",
                              save_final_rds = T
                              )


saveRDS(infercnv_obj,file="./data/temp/inferCNV_infercnv_obj_Duct_epithelial_cell2.rds")


##inferCNV3-未分群####
dir.create("./data/output/inferCNV3-未分群")
infercnv_obj = infercnv::run( infercnv_obj,
                              cutoff = 0.1,#cutoff：默认0.1，过滤低表达基因。Smart-seq2数据选“1”, 10x Genomics数据选“0.1”
                              #min_cells_per_gene：默认3，过滤表达比例较低基因
                              #window_length：默认101，滑窗中包含的基因数目
                              out_dir = "./data/output/inferCNV3-未分群",#输出结果位置
                              
                              cluster_by_groups = F,#选择"TRUE"时，先对分析细胞类型各自聚类，区分细胞来源，再做层次聚类。
                              #选择"FALSE"时，对分析细胞类型整体聚类，并按照的数值进行聚类。
                              #默认FALSE，对肿瘤细胞进行聚类，展示为不同的亚群
                              k_obs_groups = 3,#自定义分析细胞类型分群数，与为“FALSE"时联合使用。
                              #k_obs_groups：默认1，肿瘤细胞聚类数目
                              
                              HMM = FALSE,#是否使用隐马尔可夫模型(Hidden Markov Model,,HMM)进行预测，默认使用。
                              #analysis_mode=samples#使用HMM进行预测的分析模式，包括samples,subclusters和cells.
                              #samples:默认，对每个样本进行分析，不考虑样本内的细胞差异。运行速度快，但是可能忽略一些细胞的拷贝数变异。
                              #subclusters:推荐，先对每个样本内的细胞进行聚类，再对每个聚类进行分析，考虑细胞之间的差异。能够更好地检测细胞的拷贝数变异，但是运行速度较慢。
                              #cells:对每个细胞进行分析，不进行聚类。能够得到最详细的拷贝数变异结果，但是运行速度最慢。
                              denoise = TRUE,#默认FALSE(更快)，对CNV矩阵进行降噪
                              num_threads = 10)#线程数

saveRDS(infercnv_obj,file="./data/temp/inferCNV_infercnv_obj_Duct_epithelial_cell-未分群.rds")


##inferCNV4-未分群####
infercnv_obj <- readRDS('./data/temp/seurat_obj_refer+target-未分群.rds')
dir.create("./data/output/inferCNV4-未分群")
infercnv_obj = infercnv::run( infercnv_obj,
                              cutoff = 0.1,#cutoff：默认0.1，过滤低表达基因。Smart-seq2数据选“1”, 10x Genomics数据选“0.1”
                              #min_cells_per_gene：默认3，过滤表达比例较低基因
                              #window_length：默认101，滑窗中包含的基因数目
                              out_dir = "./data/output/inferCNV4-未分群",#输出结果位置
                              cluster_by_groups = T,#选择"TRUE"时，先对分析细胞类型各自聚类，区分细胞来源，再做层次聚类。
                              #选择"FALSE"时，对分析细胞类型整体聚类，并按照的数值进行聚类。
                              #默认"FALSE"，对肿瘤细胞进行聚类，展示为不同的亚群
                              #k_obs_groups = 8,#自定义分析细胞类型分群数，与为“FALSE"时联合使用。
                              #k_obs_groups：默认1，肿瘤细胞聚类数目
                              HMM = T,#是否使用隐马尔可夫模型(Hidden Markov Model,,HMM)进行预测，默认使用。
                              #analysis_mode=samples#使用HMM进行预测的分析模式，包括samples,subclusters和cells.
                              #samples:默认，对每个样本进行分析，不考虑样本内的细胞差异。运行速度快，但是可能忽略一些细胞的拷贝数变异。
                              analysis_mode = "subclusters",#推荐，先对每个样本内的细胞进行聚类，再对每个聚类进行分析，考虑细胞之间的差异。能够更好地检测细胞的拷贝数变异，但是运行速度较慢。
                              #cells:对每个细胞进行分析，不进行聚类。能够得到最详细的拷贝数变异结果，但是运行速度最慢。
                              denoise = T,#默认FALSE(更快)，对CNV矩阵进行降噪
                              num_threads = 10,
                              write_expr_matrix = T,#输出矩阵
                              output_format = "pdf",
                              save_final_rds = T
                              )#线程数

saveRDS(infercnv_obj,file="./data/temp/inferCNV_infercnv_obj_Duct_epithelial_cell1_inferCNV4-未分群.rds")
#参数和输出文件解读####
# https://www.jianshu.com/p/6b44e511f641
# infercnv.preliminary.png : 初步的inferCNV展示结果（未经去噪或HMM预测）
# 
# infercnv.png : 最终inferCNV产生的去噪后的热图.
# 
# infercnv.references.txt : 正常细胞矩阵.
# 
# infercnv.observations.txt : 肿瘤细胞矩阵.
# 
# infercnv.observation_groupings.txt : 肿瘤细胞聚类后的分组关系.
# 
# infercnv.observations_dendrogram.txt : NEWICK格式，展示细胞间的层次关系.
##plot####
library(RColorBrewer)
#每一步都会输出_obj文件，可以话不同步骤下的cnv图
infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "better_plot",
                   output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色

infercnv_obj <- readRDS("./data/temp/inferCNV_infercnv_obj_Duct_epithelial_cell.rds")



library(scales)
cnvScore < function(data) {
  data <- data %>% as.matrix %>%
    t() %>%
    scale() %>%
    scales::rescale(to = c(-1, 1)) %>%
    t()
  cnv_score < as.data.frame(colSums(data,data))
  return(cnv_score)
}
data <- read.table("./data/inferCNV/infercnv.observations.txt", header = T)
expr = data %>% as.matrix()
expr.scale = scale(t(expr))
tmpl = sweep(expr.scale, 2, apply(expr.scale, 2, min),'-')
tmp2 = apply(expr.scale, 2, max) - apply(expr.scale, 2, min)
expr_1 = t(2 * sweep(tmp1, 2, tmp2, "/") - 1)
cnv_score = as.data.frame(colSums(expr_1 * expr_1))
colnames(cnv_score)="cnv_score"
cnv_score=rownames_to_column(cnv_score,var='cell')
gene_data2$cell=rownames (gene_data2)
colnames(gene_data2)[1]="cluster"
test=merge(cnv_score,gene_data2,by="cell",all=F)
ggplot2:ggplot(test,aes (x=cluster,y=cnv_score))+
  geom_vio1in(aes(fi11=cluster),cex=1.2)+#根据Ancestry的不同因子使用不同颜色，其实用R
  scale_fill_manual(values=c("#FB5554", "#868B31","#42F203","#579ABB","#B978AE","#FFA500","#b5aa82","#de9d3d","#347852","#ca8399","#296097","#564b84"))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size=20)+
  theme(axis.text, element_text(color="black",legend.position='none'))

#4.出图#####

library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")
seurat_object<- readRDS("./data/temp/seurat_object_refer+target-分群下采样200.rds")
# 目标
obs <- read.table("./data/output/inferCNV2/infercnv.observations.txt",header = T,check.names = F) 
# 参考
ref <- read.table("./data/output/inferCNV2/infercnv.references.txt",header = T,check.names = F) 

infercnv_obj = readRDS("./data/output/inferCNV2/run.final.infercnv_obj")

#提取obs和ref的inferCNV的结果   
# min(obs)
# max(obs)
# if(F){              #可以通过定义obs的元素数值来让差异变大   在后面画图就能够更大的差异   也可以不运行，更真实体现拷贝数（但是可能没有啥差异）
#   max(obs)          #根据最大最小值来定义
#   min(obs)     
#   obs[obs > 0.6 & obs < 0.7] <- -2   #把0.6-0.7定义为数值2  后面依此类推
#   obs[obs >= 0.7 & obs < 0.8] <- -1
#   obs[obs >= 0.8 & obs < 1.0] <- 0
#   obs[obs >= 1.0 & obs <= 1.1] <- 1
#   obs[obs > 1.1 & obs <= 1.3] <- 2
#   obs[obs > 1.3] <- 2
# }

if(T){  
  obs=obs-1
  obs=obs ^ 2
  ref=ref-1
  ref=ref ^ 2
}

score_obs=as.data.frame(colSums(obs))   #把obs的每一个基因拷贝量加起来，就是这个细胞的总拷贝数obs
score_ref=as.data.frame(colSums(ref))   #把obs的每一个基因拷贝量加起来，就是这个细胞的总拷贝数obs
colnames(score_obs)=c("CNV_score")
colnames(score_ref)=c("CNV_score")
score_all <- rbind(score_obs,score_ref)

#提取meta信息  celltype,seurat_clusters,orig.ident
meta <- subset(seurat_object@meta.data,select = c("labels","celltype","tech","seurat_clusters"))  #提取该细胞的其他的meta信息

#将meta信息添加给score
meta <- rownames_to_column(meta)
score <- rownames_to_column(score_all)
meta_score <- merge(score,meta,by.x = "rowname",by.y = "rowname")    #这里会可能损失一些细胞   为什么呢，因为在前面infer时，有一些细胞无法推断，会被删掉，但是总体上问题不大
#我们现在可以用ggplot画图了     但是直接这样出图很丑   因为他会根据你的Y轴的大小  所以建议定义y轴范围
meta_score$labels <- factor(meta_score$labels, levels = c("Group_1", "Group_2", "Group_3", "Group_4", "Group_5", "Group_6", "Group_7", "Group_8", "Group_9", "Endothelial_cell", "Macrophage", "Stellate_cell", "T_cell"))
meta_score$tech <- factor(meta_score$tech, levels = c("Tumor", "Normal"))


##4.1 CNV_all####
p1<- ggplot(meta_score, aes(x = labels, y = CNV_score, fill = labels, color=tech)) +
  # geom_boxplot(geom = "errorbar",width=0.2) +
  stat_boxplot(geom = "errorbar",width=0.2, position = position_dodge(0.8),aes(color=tech))+
  geom_boxplot(width=0.7, position = position_dodge(0.8))+
  # facet_wrap(~tech, scales = "free_y") +
  # theme(strip.text.x = element_text(angle = 90))+
  scale_fill_manual(values=c("#EB746A", "#CB8F0D", "#8FA526","#2CAD3F","#29B396","#23AED4","#6792CD","#AA79B3","#DA68A3",
                             "#DCDDDD","#DCDDDD","#DCDDDD","#DCDDDD"))+
  theme_bw()+                                        #去除背景
  theme(panel.grid=element_blank())+                 #去除网格
  coord_cartesian (ylim = c(0, 40.5))#xlim =c(5, 20), #坐标范围
p1
pdf("./data/output/Duct_epithelial_CNV_all.pdf",width = 12,height = 6)
p1
dev.off()
ttest_results <- t.test(CNV_score ~ tech, data = meta_score[meta_score$labels=="Endothelial_cell",])
ttest_results <- t.test(CNV_score ~ tech, data = meta_score[meta_score$labels=="G",])


groups <- split(meta_score$CNV_score, meta_score$labels)

# 初始化一个空的数据框，用于存储t检验结果
t_test_results <- data.frame(group = character(), t_statistic = numeric(), p_value = numeric())

# 对每个组进行t检验，并将结果存储到t_test_results数据框中
for (i in 1:length(groups)) {
  for (j in (i + 1):length(groups)) {
    t_test_result <- t.test(groups[[i]], groups[[j]])
    t_test_results <- rbind(t_test_results, data.frame(group = paste(c(names(groups)[i], names(groups)[j]), collapse = " vs "),
                                                       t_statistic = t_test_result$statistic,
                                                       p_value = t_test_result$p.value))
  }
}





##4.2 CNV_tumor_normal####
p2 <- ggboxplot(meta_score,"tech","CNV_score",fill="tech")
labely = max(meta_score$CNV_score)# 用的max值做的T检验
compare_means(CNV_score ~ tech,  data = meta_score)
my_comparisons <- list( c("Tumor", "Normal") )# 修改
p3 = p2+ 
  stat_compare_means(label = "p.format",label.x=2.12)+#size = 3,label.y = 0.6,label.x = 1.7
  stat_compare_means(comparisons = my_comparisons,
                     # label = "p.format",
                     aes(label = paste0("p = ", after_stat(p.format))),
                     # size = 3,
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))
p3
pdf("./data/output/Duct_epithelial_CNV_normal_tumor.pdf",width = 6,height = 6)
p3
dev.off()
#小提琴图（未使用）
meta_score%>%ggplot(aes(labels,score_sum))+geom_violin(aes(fill=labels),color="NA")
