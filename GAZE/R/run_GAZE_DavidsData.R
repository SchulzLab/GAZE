set.seed(0)
source("/projects/triangulate/work/GAZE/GAZE/R/GAZE_parallel_complete.R")

## Instantiate GAZE object
file_name <- "/projects/triangulate/work/data/seurat_object_hca_as_harmonized_AS_SP_nuc_refined_cells_NAsRemoved.rds"
gaze.opt <- GAZE.options(ATAC_exists= F, MetacellaR= T, gtf_file= "/projects/abcpp/work/base_data/gencode.v38lift37.basic.annotation.gtf")

## Set MetaCellaR's options and create the GAZE object
mc.opts <- MetaCellaR.options(output_file= "/projects/triangulate/work/GAZE/GAZE/metacell_David_e100_d20",
															RNA_count_slot= "assays$RNA@data",
															celltype_info= "meta.data$cell_type",
															umap_dim= 20,
															expected_cells= 100,
															threshold= 300,
															umap_flag= T,
															metadata= c("cell_type", "orig.ident", "dataset", "condition"),
															normalization_flag= F
)
gaze.obj <- CreateGAZEobject(gaze.options= gaze.opt, seurat_file= file_name, metacell.options= mc.opts)



## Set the STARE options and run
x <- read.table("/projects/triangulate/work/GAZE/GAZE/data/Davids_genes.bed", header= F)
dups_res <- duplicated(sapply(seq(nrow(x)), function(i)paste(x[i,1:3], collapse= "_")))
xx <- x[!dups_res,]
common_genes <- intersect(xx[, 5], rownames(gaze.obj@expr_mat))
subset_x <- xx[which(xx[, 5] %in% common_genes), seq(3)]
write.table(subset_x, "/projects/triangulate/work/GAZE/GAZE/data/Davids_genes_subset.bed", quote= F, row.names= F, col.names= F, sep= "\t")


## create the file for the -u flag in STARE for limiting its computation for a subset of genes only
writeLines(text= rownames(gaze.obj@expr_mat), "/projects/triangulate/work/GAZE/GAZE/data/Davids_genes_subset.txt")
stare_opt <- new("STARE.options",
								 g= "/projects/abcpp/work/base_data/GRCh37.p13.genome.fa",
#								 g= "/projects/abcpp/work/base_data/hg38.fa",
#b= paste0("/projects/triangulate/work/GAZE/GAZE/data/Davids_genes_subset.bed"), ## have to get this one from TSS_Fethcer
o= "Davids_data_STARE",
u= "/projects/triangulate/work/GAZE/GAZE/data/Davids_genes_subset.txt",
s= "/projects/triangulate/work/STARE/PWMs/2.2/Jaspar_Hocomoco_Kellis_human_transfac.txt",
a= "/projects/abcpp/work/base_data/gencode.v38lift37.basic.annotation.gtf",
#a= "/projects/abcpp/work/base_data/gencode.v38.annotation.gtf",
w= 5000,
e= "False", c= 10)

gaze.obj@STARE_options <- stare_opt

stare_res <- STARE.run(stare_opt, mode= "static")

gtf <- as.data.frame(rtracklayer::import(gaze.opt@gtf_file))
gtf_genes <- subset(gtf, type == "gene")
ID_mapping <- gtf_genes[, c("gene_id", "gene_name")]

gaze.obj@ID_mapping <- ID_mapping

GAZE.build(gaze.obj, "static")





high_cors <- which(sapply(seq(length(all_cors)), function(i) sum(all_cors[[i]])) > 30)
high_gene_hits <- which(gene_names %in% names(all_cors)[high_cors])
pdf("/projects/triangulate/work/GAZE/GAZE/data/GAZE_outputs/mojitoo_metacell_30_mm10_STARE_UMAP.pdf");
for(i in high_gene_hits[seq(20)]){
	umap_res <- uwot::umap(as.matrix(feature_mat_list[[i]], scale = T));
	cell_names_df <- sapply(cell_names, function(j){strsplit(j, "_")[[1]][2]})
	df <- data.frame(umap_res, cell= cell_names_df);
	print(ggplot2::ggplot(df, aes(x= X1, y= X2, shape= cell)) + geom_point(aes(colour= cell)) + scale_shape_manual(values= seq(length(unique(cell_names_df)))) + ggtitle(paste0("gene: ", i, " or ", gene_names[i])))
};
dev.off()


###  Add the metadata slots for metacell and singlecell
metadata_metacell <- data.frame(celltype= sapply(cell_names_df$group_names, function(i)strsplit(i, "_")[[1]][1]), donor= sapply(cell_names_df$group_names, function(i)strsplit(i, "_")[[1]][2]), condition= sapply(cell_names_df$group_names, function(i)strsplit(i, "_")[[1]][3]))
gaze.obj@metadata_metacell <- metadata_metacell

metadata_singlecell <- data.frame(celltype= gaze.obj@celltype, condition= gaze.obj@condition)
gaze.obj@metadata_singlecell <- metadata_singlecell
gaze.lite <- gaze.obj
gaze.lite@shap_values <- list()
saveRDS(gaze.lite, "/projects/triangulate/work/GAZE/GAZE/metacell_David_e100_d20_peakFeatures_FALSE_GAZE_object_lite_noSHAPs.rds")
