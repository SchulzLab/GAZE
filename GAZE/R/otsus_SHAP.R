library(image.Otsu)
library(pheatmap)
load("/Users/fatemehbehjati/Documents/GAZE/SHAP/Davids_shap_summary.RData")

my_scale <- function(x, a, b){
	return((x - min(x))/(max(x) - min(x)) * (b - a) + a)
}

bin_res <- list();
if(F){
for(i in seq(length(gene_summary))){
	mat <- abs(gene_summary[[i]]);
	mat_scaled <- my_scale(mat, 0, 255);
	bin_res[[i]] <- image_otsu(round(mat_scaled))$x
}
}

for(i in seq(length(bin_res))){bin_res[[i]] <- apply(bin_res[[i]], 2, as.logical)}
names(bin_res) <- names(gene_summary)

TF_SHAP_names <- list()
for(i in seq(length(bin_res))){
	if(!is.null(dim(gaze.obj@shap_values[[i]]))){
		TF_SHAP_names[[i]] <- rownames(gene_summary[[i]])
		rownames(gaze.obj@shap_values[[i]]) <- rownames(gene_summary[[i]])
	}
	else{
		TF_SHAP_names[[i]] <- names(gene_summary[[i]])
		tmp <- t(as.matrix(gaze.obj@shap_values[[i]]))
		rownames(tmp) <- rownames(gene_summary[[i]])
		gaze.obj@shap_values[[i]] <- tmp
		#names(gaze.obj@shap_values[[i]]) <- rownames(gene_summary[[i]])
	}
}
names(TF_SHAP_names) <- names(bin_res)
save(bin_res, file= "bin_res_bool.RData")
pdf("SHAP_otsu.pdf");
for(i in seq(5)){
	bin_mat <- image_otsu(round(my_scale(abs(gene_summary[[i]]), 0, 255)));
	pheatmap(bin_mat$x, cluster_rows= F, cluster_cols= F, fontsize_row= 3, main= paste0("Otsu's thresholding\nfor ", names(gene_summary)[i]));
	pheatmap(gene_summary[[i]], cluster_rows= F, cluster_cols= F, fontsize_row= 3, main= paste0("SHAP average\n for ", names(gene_summary)[i]))
};
dev.off();
