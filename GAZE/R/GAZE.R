#library(MetaCellaR)
library(triangulate)
## triangulate options
triang_opts <- list(file_name= "",
										RNA_count= "",
										ATAC_count= "",
										celltype= "",
										assay= "",
										umap_flag= T,
										expected_cells= 30,
										pca_components= 30,
										umap_components= 20,
										conditions= ""
										)
stare_opts <- list(g= "",
									 b= "",
									 o= "",
									 p= "",
									 a= "",
									 w= 5000000,
									 e= F,
									 n= "4+",
									 j= "",
									 f= "",
									 k= "",
									 t= .02,
									 z= T,
									 c= 10,
									 y= .41,
									 x= "",
									 q= T,
									 m= 5*10^6,
									 d= T,
									 r= ""
									 )
setClass("STARE.options", slots= list(g= "character",
                                     b= "character",
                                     o= "character",
                                     p= "character",
																		 a= "character",
																		 w= "numeric",
																		 e= "character",
																		 n= "character",
																		 j= "character",
																		 f= "character",
																		 k= "numeric",
																		 t= "numeric",
																		 z= "character",
																		 c= "numeric",
																		 y= "numeric",
																		 x= "character",
																		 q= "character", ### In fact a logical
																		 m= "numeric",
																		 d= "character", ### In fact a logical
																		 r= "character"
									 								)
)
setClass("MetaCellaR.options", slots= list(umap_flag= "logical",
																					 RNA_slot= "character",
																					 ATAC_slot= "character",
																					 celltype_slot= "character",
																					 assay_slot= "character",
																					 expected_cells= "numeric",
																					 pca_components= "numeric",
																					 umap_components= "numeric",
																					 output_file= "character",
																					 condition= "character")
)

MetaCellaR.options <- function(umap_flag= T, expected_cells= 30, pca_components= 20, umap_components= 2, ...){
	return(new("MetaCellaR.options", umap_flag= umap_flag, expected_cells= expected_cells, pca_components= pca_components, umap_components= umap_components, ...))
}

MetaCellaR.run <- function(gaze.obj, obj){
	mc_command <- paste0("time Rscript /projects/triangulate/work/MetaCellaR/MetaCellaR/MetaCellaR.R -file ", gaze.obj@seurat_file, " -RNA ", obj@RNA_slot, " -celltype ", obj@celltype_slot, " -e ", obj@expected_cells, " -umap ", obj@umap_flag, " -output ", obj@output_file)
	print("Running MetaCellaR started...")
	sys_res <- system(mc_command)
	print("MetaCellaR finished!")
	if(!sys_res){
		mc_res <- read.csv(paste0(obj@output_file, "/cellSummarized_kmed_means_sum.csv"), row.names= 1)
		gaze.obj@metacell_mat <- mc_res

		mc_map <- read.csv(paste0(obj@output_file, "/RNA_cell2metacell_info_kmed_means.csv"))
		gaze.obj@RNA_MC_map <- mc_map
		print("metacell results successfully added to the metacell_mat slot of the GAZE object!")
		print(paste("There are", ncol(gaze.obj@metacell_mat), "metacells created!"))

	}
	return(gaze.obj)
}

setGeneric("STARE.run", function(obj) standardGeneric("STARE.run"))

#setMethod("STARE.run", "STARE", function(g, b, o, p, a, s= 5000000, e= F, n= "4+", j, f, k, t= .02, z= T, c= 10){
setMethod("STARE.run", "STARE.options", function(obj){
						exceptional_args <- ""
						if(obj@j){
							exceptional_args <- paste0(exceptional_args, " -j")
						}
						if(obj@z){
              exceptional_args <- paste0(exceptional_args, " -z")
            }
						command <- paste0("time bash /projects/triangulate/work/GAZE/GAZE/R/Code/STARE.sh -g ", obj@g, " -b ", obj@b, " -p ", obj@p, " -o ", obj@o, " -a ", obj@a, " -w ", obj@w, " -e ", obj@e, " -n ", obj@n, " -f ", obj@f, " -k ", obj@k, " -t ", obj@t, " -c ", obj@c, " -z ", obj@z, " -j ", obj@j)
						system(command)
					}
)

## Only when sc-RNA data is available
setClass("Gene_Learning.options", slots= list(inputX= "character",
																			inputY= "character",
																			trainingSet_percent= "numeric",
																			cvFolds= "numeric",
																			maxIter= "numeric",
																			output= "character")
)

## Gene_Learning.options constructor ##
Gene_Learning.options <- function(inputX, inputY, trainingSet_percent= .6, cvFolds= 5, maxIter= 200, output){
	return(new("Gene_Learning.options", inputX= inputX, inputY= inputY, trainingSet_percent= trainingSet_percent,
						 cvFolds= cvFolds, output= output))
}
#######################################

## Only when sc-ATAC & sc-RNA data is available
setClass("Cell_Learning.options", slots= list(STARE_options= "STARE.options")
)

Cell_Learning.options <- function(STARE_options, nfolds= 5, train_split= .8)

setClass("Learning", slots= list(mode= "character",
																 Gene_Learning= "Gene_Learning"),

)

setClass("GAZE.options", slots= list(TPM= "logical",
																		 MetacellaR= "logical",
																		 genome_version= "character", ## MUST
																		 organism= "character", ## MUST
																		 gtf_file= "character", ## MUST
																		 ATAC_exists= "logical",
																		 filtering= "logical",
																		 rna_file= "character",
																		 atac_file= "character")
)

GAZE.options <- function(TPM= F, MetacellaR= F, ATAC_exists= F, filtering= F, ...){
	return(new("GAZE.options", TPM= TPM, MetacellaR= MetacellaR, ATAC_exists= ATAC_exists, filtering= filtering, ...))
}

library(Matrix) # have to run this for the "dgCMatrix" to be recognized during the compliation of a gaze object
setClass("GAZE", slots= list(STARE_options= "STARE.options",
														 triang_options= "list",
														 gaze.options= "GAZE.options",
														 seurat_file= "character",
														 expr_mat= "dgCMatrix",
														 metacell_mat= "data.frame",
														 RNA_MC_map= "data.frame",
														 condition= "character",
														 celltype= "character",
														 ATAC_mat= "matrix",
														 ATAC_peaks= "character",
														 lost_genes= "data.frame",
														 lost_TFs= "data.frame",
														 models= "list",
														 shap_values= "list",
														 ID_mapping= "data.frame",
														 TF_cell_mat= "matrix"
														 )
)

CreateGAZEobject <- function(gaze.options= gaze.opt, seurat_file= seurat_file, ...){
	obj <- new("GAZE", gaze.options= gaze.opt, seurat_file= seurat_file, ...)

	if(is.null(obj@seurat_file)){
		tryCatch({
			obj@expr_mat <- read.csv(obj@rna_file)
			if(gaze.opt@ATAC_exists){
				obj@ATAC_mat <- read.csv(atac_file)
			}
		},
		error= function(e){
			message("GAZE initialize failed on reading the `rna_file` due to the error below:")
			print(e)
		})
	}else{## The Seurat object exists:
		seu_obj <- tryCatch({
			seu_obj <- readRDS(obj@seurat_file)
			return(seu_obj)
		},
		error= function(e){
			print("Trying to use the load function instead of readRDS to load the Seurat object!")
			seu_obj <- load(obj@seurat_file)
			if(class(seu_obj) == "character"){
				seu_obj <- get(seu_obj) ## get the actual Seurat object from the name of the stored object in obj@seurat_file
			}
			return(seu_obj)
		}
		)
		seu_obj_names <- names(attributes(seu_obj))
		if(length(which(seu_obj_names %in% c("meta.data", "assays"))) < 2){
			stop("The given Seurat object doesn't contain either `meta.data` or `assays`! These entries must contain `meta.data` and `assays`.")
		}
		obj@expr_mat <- seu_obj[["RNA"]]@counts
		obj@celltype <- as.character(seu_obj@meta.data$celltype)
		if(!is.null(seu_obj@meta.data$condition)){
			obj@condition <- seu_obj@meta.data$condition
		}
		if(gaze.opt@ATAC_exists){
			obj@ATAC_mat <- seu_obj[["ATAC"]]@counts
			obj@ATAC_peaks <- seu_obj[["ATAC"]]@peaks
		}

	}
	return(obj)
}
if(F){
setGeneric("GAZE.initialize", function(obj) standardGeneric("GAZE.initialize"))
setGeneric("GAZE.options.initialize", function(obj) standardGeneric("GAZE.options.initialize"))

setMethod("GAZE.options.initialize", "GAZE.options", function(obj){
						if(is.null(obj@seurat_file)){
							tryCatch({
								expr_mat <- read.csv(obj@rna_file)
								if(obj@ATAC_exists){
									atac_mat <- read.csv(atac_file)
								}
							},
							error= function(e){
								message("GAZE.options.initialize failed on reading the `rna_file` due to the error below:")
								print(e)
							})
						}
						else{## The Seurat object exists:
							seu_obj <- readRDS(obj@seurat_file)
							seu_obj_names <- names(seu_obj)
							if(length(which(seu_obj_names %in% c("meta.data", "RNA")))){
								error("The given Seurat object doesn't contain either `meta.data` or `RNA`! These entries must contain `meta.data` and `RNA`.")
							}
							expr_mat <- seu_obj[["RNA"]]@counts
							celltype <- seu_obj[["meta.data"]]$celltype
							if(!is.null(seu_obj[["meta.data"]]$condition)){
								condition <- seu_obj[["meta.data"]]$condition
							}
							if(obj@ATAC_exists){
								atac_mat <- seu_obj[["ATAC"]]@counts
							}

						}

						

})


#setGeneric("get.GAZE.options", function(obj)standardGeneric("get.GAZE.options"))
setMethod("get.GAZE.options", "GAZE.options",
					function(obj){
						cat("TPM:", obj@TPM, "\n")
						cat("MetacellaR:", obj@MetacellaR, "\n")
						cat("ATAC_exists:", obj@ATAC_exists, "\n")
						cat("filtering:", obj@filtering, "\n")
					}
)
}
GAZE_constructor <- function(){
	GAZE_obj <- list()
	GAZE_obj[["options"]][["TPM"]] <- F
	GAZE_obj[["options"]][["MetacellaR"]] <- F
	GAZE_obj[["options"]][["ATAC_exists"]] <- F
	GAZE_obj[["options"]][["filtering"]] <- T

	GAZE_obj[["expr_mat"]] <- NULL
}



if(F){
### Test STARE ###

	stare_opt <- new("STARE.options", g= "/projects/abcpp/work/base_data/GRCh37.p13.genome.fa", b= "/projects/abcpp/work/base_data/K562_10xRandomActivity.bed", o= "Test_K562_10xRandom_GAZE3", p= "/projects/triangulate/work/STARE/PWMs/2.0/human_jaspar_hoc_kellis.PSEM", a= "/projects/abcpp/work/base_data/gencode.v38lift37.basic.annotation.gtf", w= 5000000, e= "False", n= "4+", j= "True", f= "/projects/abcpp/work/base_data/K562_ChrWise_Contacts/", k= 5000, t= .2, z= "True", c= 10)
	stare_res <- STARE.run(stare_opt)

	if(stare_res){
		stop("STARE module failed to finish successfully. Please review the corresponding arguments more carefully based on the error messages from the STARE run!")
	}

## read the reshape matrix and load it into the GAZE object
	reshaped <- read.table(paste0(stare_opt@o, "/Stacked_Matrix/", stare_opt@o, "_StackedAffinityMatrix.txt.gz"))
}

test_mode <- F

if(test_mode){

## Test MetaCellaR ##
	file_name <- "../data/3008_.Rds"
	mc.opts <- MetaCellaR.options(assay_slot= "sc-ATAC", condition= "Age", output_file= "test_out", RNA_slot= "assays\\$RNA@counts", celltype_slot= "meta.data\\$celltype")

	gaze.opt <- GAZE.options(ATAC_exists= F, MetacellaR= T, gtf_file= "/projects/abcpp/work/base_data/gencode.v38lift37.basic.annotation.gtf")
#	gaze.obj <- new("GAZE", gaze.options= gaze.opt, seurat_file= file_name)
	gaze.obj <- CreateGAZEobject(gaze.options= gaze.opt, seurat_file= file_name)

	gaze.obj <- MetaCellaR.run(gaze.obj, mc.opts)

	## Add STARE's output to the GAZE object
	#input.x <- read.csv("../data/scMTL_heartecCellSumKmedmeans3008_notImputed_epigenetic_feature_doubleReduced.txt")
	input.x <- readr::read_tsv("../../Test_K562_10xRandom_GAZE4/Gene_TF_matrices/Test_K562_10xRandom_GAZE4_TF_Gene_Affinities_Apple.txt.gz")
	gtf <- as.data.frame(rtracklayer::import(gaze.opt@gtf_file))
	gtf_genes <- subset(gtf, type == "gene")


	## Merge X and Y variables based on their gene names and keep track of the lost genes
	ID_mapping <- gtf_genes[, c("gene_id", "gene_name")]
	merged_input <- merge(ID_mapping, input.x, by.x= "gene_id", by.y= "geneID")
	x_y_merged <- merge(gaze.obj@metacell_mat, merged_input, by.x= 0, by.y= "gene_name")
	ENS_IDs <- x_y_merged[, ncol(gaze.obj@metacell_mat) + 2]
	x_y_merged <- x_y_merged[, -(ncol(gaze.obj@metacell_mat) + 2)]
	x_y_hits <- intersect(rownames(gaze.obj@metacell_mat), merged_input$gene_name)
	## Separate X and Y from the merged matrix
	Y <- x_y_merged[, seq(2, ncol(gaze.obj@metacell_mat) + 1)]
	rownames(Y) <- ENS_IDs
	Y <- as.matrix(Y)
	
	X <- x_y_merged[, seq(ncol(gaze.obj@metacell_mat) + 2, ncol(x_y_merged))]
	rownames(X) <- ENS_IDs
	X <- as.matrix(X)

	## Update the GAZE object


	lost_genes <- NULL;
	lost_genes <- data.frame(genes= setdiff(merged_input$gene_name, x_y_hits), source= "Gene-TF input matrix")
	lost_genes <- data.frame(genes= setdiff(rownames(gaze.obj@metacell_mat), x_y_hits), source= "Gene expression input matrix (gaze.object)")
	gaze.obj@lost_genes <- lost_genes
	gaze.obj@TF_cell_mat <- X

	tri_res <- triangulate_main(input.x= X, input.y= Y, percent= .8, output= "triangulate_out_GAZE")
	gaze.obj@TF_cell_mat <- tri_res$TGL.model$B

	sec_res <- second_level_learnings(TGL.model_m= tri_res$TGL.model, partition_m= tri_res$partition, ID_mapping)

	gaze.obj@models <- sec_res[["models"]]
	gaze.obj@shap_values <- sec_res[["shap"]]
	gaze.obj@ID_mapping <- ID_mapping

##########
	#gaze.opt <- new("GAZE.options", ATAC_exists= ATAC_exists, MetacellaR= MetacellaR,filtering= filtering)

	#gaze <- new("GAZE", gaze.opt, expr_mat= expr_mat, metacell_mat= as.matrix(metacell_mat), condition= condition, celltype= as.character(celltype), models= models$models, shap_values= models$shap)

}
