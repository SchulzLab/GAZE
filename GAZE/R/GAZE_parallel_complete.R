library(MetaCellaR)
library(triangulate)
source("/projects/triangulate/work/triangulate_pkg/triangulate/R/one_level_learning_parallel_pkg.R")
options(scipen=999)
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
setClass("STARE.options", slots= list(g= "character",
                                     b= "character",
                                     o= "character",
                                     p= "character",
																		 s= "character",
																		 a= "character",
																		 w= "numeric",
																		 e= "character",
																		 n= "character",
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
																		 r= "character",
																		 u= "character"
									 								)
)
setClass("MetaCellaR.options", slots= list(umap_flag= "logical", file_name= "character",
																					 RNA_count_slot= "character",
																					 ATAC_count_slot= "character",
																					 celltype_info= "character",
																					 assay_slot= "character",
																					 expected_cells= "numeric",
																					 umap_dim= "numeric",
																					 reduction= "character",
																					 threshold= "numeric",
																					 k= "numeric",
																					 output_file= "character",
																					 condition= "character",
																					 normalization_flag= "logical",
																					 merge_flag= "logical",
																					 csv_flag= "logical",
																					 gtf_file= "character",
																					 metadata= "character")
)

MetaCellaR.options <- function(normalization_flag = T, csv_flag = F, umap_flag= T, expected_cells= 30, umap_dim= 20, merge_flag= T,reduction= "umap", threshold = 90, output_file = getwd(), ...){
	return(new("MetaCellaR.options", normalization_flag= normalization_flag, csv_flag= csv_flag, umap_flag= umap_flag, expected_cells= expected_cells,
	umap_dim= umap_dim, merge_flag= merge_flag, reduction= reduction, threshold= threshold, output_file=output_file, ...))
}

MetaCellaR.run <- function(gaze.obj, obj){
	#nonemptySlotNames <- unlist(sapply(slotNames(obj), function(sn)if(length(slot(obj, sn)))return(sn)))
	#nonemptySlots <- sapply(nonemptySlotNames, function(sn) slot(obj, sn))

	#flags_dt <- data.table::data.table(keys= nonemptySlotNames, values= nonemptySlots)

	#mc_command <- paste0("time Rscript /projects/triangulate/work/MetaCellaR/MetaCellaR/MetaCellaR.R -file ", gaze.obj@seurat_file, paste0(" -", flags_dt$keys, " ", flags_dt$values, collapse= ""))
	#print(mc_command)

	print("Running MetaCellaR started...")
	#sys_res <- system(mc_command)
	sys_res <- 0
	metacellar.run(file_name= gaze.obj@seurat_file,
merge_flag = obj@merge_flag, normalization_flag = obj@normalization_flag, reduction = obj@reduction, gtf_file= obj@gtf_file,
umap_flag = obj@umap_flag, RNA_count_slot= obj@RNA_count_slot, ATAC_count_slot= obj@ATAC_count_slot, celltype_info= obj@celltype_info , assay_slot= obj@assay_slot, output_file= obj@output_file,
expected_cells = obj@expected_cells, threshold = obj@threshold, umap_dim = obj@umap_dim, k= obj@k, metadata= obj@metadata)
	if(!sys_res){
		print("MetaCellaR finished!")
		print(obj@output_file)
		if(obj@normalization_flag){
		mc_res <- read.csv(paste0(obj@output_file, "/results/RNA_cellSummarized_normalized_TPM.csv"), row.names= 1, check.names= F) ## first DESeq2 and then exonic length normalization
		}else{
			mc_res <- read.csv(paste0(obj@output_file, "/results/RNA_cellSummarized_normalized.csv"), row.names=1, check.names= F) ## the data was already normalized, so we just take it as it is!
		}
		gaze.obj@metacell_mat <- mc_res

		mc_map <- read.csv(paste0(obj@output_file, "/results/RNA_cell2metacell_info.csv"))
		gaze.obj@RNA_MC_map <- mc_map

		mc_umap <- read.csv(paste0(obj@output_file, "/results/RNA_metacell_umap.csv"), row.names= 1)
		gaze.obj@RNA_MC_UMAP <- mc_umap

		print("metacell results successfully added to the metacell_mat slot of the GAZE object!")
		print(paste("There are", ncol(gaze.obj@metacell_mat), "metacells created!"))
		if(gaze.obj@gaze.options@ATAC_exists){
			atac_res <- read.csv(paste0(obj@output_file, "/results/cellSummarized_ATAC.csv"), row.names= 1)
			tokenzied_names <- sapply(rownames(atac_res), function(rn) strsplit(rn, "[: -]")[[1]])
			tokenzied_names <- t(tokenzied_names)
			colnames(tokenzied_names) <- c("#chr", "start", "end")

			new_atac_res <- cbind(tokenzied_names, atac_res)
			write.table(new_atac_res, paste0(obj@output_file, "/results/cellSummarized_ATAC_modifiedByGAZE.txt"), row.names= F, quote= F, sep= "\t")
			gaze.obj@ATAC_MC_map <- read.csv(paste0(obj@output_file, "/results/ATAC_cell2metacell_info.csv"))

			## Not all RNA metacells are present in ATAC metacells. Therefore, one needs to amend such discrepency through taking the intersect of surviving metacells between both assays.
			common_metacells <- intersect(gaze.obj@ATAC_MC_map$metacell, gaze.obj@RNA_MC_map$metacell)
			## fix the mismatch in actual given cell types in the data and those placed in the column names of dataframes by R
			common_metacells <- gsub(c(" |-"), ".", common_metacells)
			gaze.obj@metacell_mat <- dplyr::select(gaze.obj@metacell_mat, common_metacells[seq(common_metacells)]) #this "[seq(common_metacells)]" is a workaround I found to resolve this mysterious bug in the dplyr package!
		}

	}
	return(gaze.obj)
}

setGeneric("STARE.run", function(obj, mode= NULL) standardGeneric("STARE.run"))

setMethod("STARE.run", "STARE.options", function(obj, mode= NULL){
						nonemptySlotNames <- unlist(sapply(slotNames(obj), function(sn)if(length(slot(obj, sn)))return(sn)))
						nonemptySlots <- sapply(nonemptySlotNames, function(sn) slot(obj, sn))

						flags_dt <- data.table::data.table(keys= nonemptySlotNames, values= nonemptySlots)
						if(mode == "static"){
						command <- paste0("time bash /projects/triangulate/work/STARE/Code/STARE.sh", paste0(" -", flags_dt$keys, " ", flags_dt$values, collapse= ""))
						}else{
							command <- paste0("time bash /projects/triangulate/work/STARE/Code/STARE.sh", paste0(" -", flags_dt$keys, " ", flags_dt$values, collapse= ""))
						}
						res <- system(command)
						if(res){
							stop("STARE failed! Please refer to the error messages printed on the prompt and make changes accordingly!")
						}
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
						 cvFolds= cvFolds, output= output)) #TODO
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
																		 output_file= "character",
																		 atac_file= "character")
)

GAZE.options <- function(TPM= F, MetacellaR= F, ATAC_exists= F, filtering= F, output_file= paste0("GAZE_output_", gsub(" ", "_", date())), ...){
	return(new("GAZE.options", TPM= TPM, MetacellaR= MetacellaR, ATAC_exists= ATAC_exists, filtering= filtering, output_file= output_file, ...))
}

library(Matrix) # have to run this for the "dgCMatrix" to be recognized during the compliation of a gaze object
setClassUnion("data.frameORmatrix", c("data.frame", "matrix", "dgCMatrix"))
setClassUnion("MetaCellaR.optionsORnull", c("MetaCellaR.options", "NULL"))
setClass("GAZE", slots= list(STARE_options= "STARE.options",
														 triang_options= "list",
														 gaze.options= "GAZE.options",
														 metacell.options= "MetaCellaR.optionsORnull",
														 seurat_file= "character",
														 rna_file= "character",
														 expr_mat= "data.frameORmatrix",
														 metacell_mat= "data.frame",
														 RNA_MC_UMAP= "data.frame",
														 RNA_raw_UMAP= "data.frame",
														 RNA_MC_map= "data.frame",
														 ATAC_MC_map= "data.frame",
														 condition= "character",
														 celltype= "character",
														 ATAC_mat= "dgCMatrix",
														 ATAC_peaks= "character",
														 lost_genes= "data.frame",
														 lost_TFs= "data.frame",
														 models= "list",
														 shap_values= "list",
														 ID_mapping= "data.frame",
														 TF_cell_mat= "matrix",
														 coefficient_map= "matrix",
														 model_accuracy= "data.frame",
														 stat_test_TF= "list",
														 stat_test_gene= "list",
														 shiny_store= "data.frame",
														 metadata_metacell= "data.frame",
														 metadata_singlecell= "data.frame"
														 )
)

CreateGAZEobject <- function(gaze.options, seurat_file, metacell.options, ...){
#CreateGAZEobject <- function(gaze.options, seurat_file= "", ...){
	#obj <- new("GAZE", gaze.options= gaze.options, seurat_file= seurat_file, metacell.options= metacell.options, ...)
	obj <- new("GAZE", gaze.options= gaze.options, seurat_file=  seurat_file, ...)
	if(!is.null(metacell.options)){
		obj@metacell.options <- metacell.options
	}
	gtf <- as.data.frame(rtracklayer::import(gaze.options@gtf_file))
	gtf_genes <- subset(gtf, type == "gene")
	## Merge X and Y variables based on their gene names and keep track of the lost genes
	ID_mapping <- gtf_genes[, c("gene_id", "gene_name")]
	obj@ID_mapping <- ID_mapping

	if(obj@seurat_file == ""){ ## This part doesn't  work, and I dunno why... changed the type of expr_mat to dataframOrmatrix but it still throws an error 
		seu_obj <- tryCatch({
			file_extension <- substr(gaze.options@rna_file, nchar(gaze.opt@rna_file) - 2, nchar(gaze.options@rna_file) )
			if(file_extension == "csv"){
				obj@expr_mat <- read.csv(obj@rna_file)
			}else{
				print("inside !csv check")
				print(gaze.options@rna_file)
				obj@metacell_mat <- read.table(gaze.options@rna_file, header= T, row.names = 1)## I have to assign it to metacell_mat because I later access only this slot for training models
				if(length(strsplit(obj@ID_mapping$gene_id[1], "\\.")[[1]]) > length(strsplit(rownames(obj@metacell_mat)[1], "\\.")[[1]])){
					obj@ID_mapping$gene_id <- sapply(obj@ID_mapping$gene_id, function(g) strsplit(g, "\\.")[[1]][1])
				}else if (length(strsplit(obj@ID_mapping$gene_id[1], "\\.")[[1]]) < length(strsplit(rownames(obj@metacell_mat)[1], "\\.")[[1]])){
					rownames(obj@metacell_mat) <- sapply(rownames(obj@metacell_mat), function(g) strsplit(g, "\\.")[[1]][1])
				}

			}
			print("read RNA file successfully!")
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
#			return(seu_obj)
		},
		error = function(e){ # nolint
			print("Trying to use the load function instead of readRDS to load the Seurat object!")
			seu_obj <- load(obj@seurat_file)
			if(class(seu_obj) == "character"){
				seu_obj <- get(seu_obj) ## get the actual Seurat object from the name of the stored object in obj@seurat_file
			}
#			return(seu_obj)
		}
		)
		seu_obj_names <- names(attributes(seu_obj))
		if(length(which(seu_obj_names %in% c("meta.data", "assays"))) < 2){
			stop("The given Seurat object doesn't contain either `meta.data` or `assays`! These entries must contain `meta.data` and `assays`.")
		}
		if(length(obj@metacell.options@assay_slot)){
			assays <- eval(parse(text= paste0("seu_obj", "@",obj@metacell.options@assay_slot)))
			print(paste("class(obj):", class(obj)))
				rna_hits <- which(tolower(assays) == "scrna-seq" | tolower(assays) == "scrna" | tolower(assays) == "scrna" | tolower(assays) == "rna")
			obj@expr_mat <- seu_obj[["RNA"]]@counts[, rna_hits]
		}else{
			obj@expr_mat <- eval(parse(text= paste0("seu_obj", "@",obj@metacell.options@RNA_count_slot)))
		}
		RNA_raw_UMAP <- uwot::umap(t(as.matrix(obj@expr_mat)), pca= 30, pca_center= T, n_components= 2, scale= T)
		rownames(RNA_raw_UMAP) <- colnames(obj@expr_mat)
		colnames(RNA_raw_UMAP) <- c("UMAP1", "UMAP2")
		obj@RNA_raw_UMAP <- as.data.frame(RNA_raw_UMAP)

		obj@celltype <- as.character(seu_obj@meta.data$celltype)
		if(!is.null(seu_obj@meta.data$condition)){
			obj@condition <- seu_obj@meta.data$condition
		}
		if(gaze.opt@ATAC_exists){
			if(length(which(names(seu_obj) == "ATAC"))){
				obj@ATAC_mat <- (seu_obj[["ATAC"]]@counts)
				obj@ATAC_peaks <- seu_obj[["ATAC"]]@peaks
			}else{
				obj@ATAC_mat <- (seu_obj[["peaks"]]@counts)
				obj@ATAC_peaks <- rownames(seu_obj[["peaks"]]@counts)
			}
		}

	}
	if(gaze.opt@MetacellaR){
		print("Running MetacellaR from GAZE...")
		obj <- MetaCellaR.run(obj, obj@metacell.options)
		print("Done running MetacellaR from GAZE!")
	}
	return(obj)
}



check_endianness <- function(parent_path){
	end_bin_path <- paste0(parent_path, "_Endianness_check.bin")
	end_txt_path <- paste0(parent_path, "_Endianness_check.txt") 
	endianness <- "little"
	end_bin <- readBin(file(end_bin_path, "rb"), n= 10, what= "double", size= 4, endian= "little")
	end_txt <- as.numeric(readLines(end_txt_path, n= 10)[-1])

	if(!all(end_bin == end_txt)){## The check failed, the endianness must be chagned
		end_bin <- readBin(file(end_bin_path, "rb"), n= 10, what= "double", size= 4, endian= "big")
		if(!all(end_bin == end_txt)){
			stop("The binary file could not be read properly. This error has occurred due to the inconsistency between the content of the binary file and its corresponding text file.")
		}
		endianness <- "big" ## Setting the right endianness for later when the real data is read
	}
	print(paste("Endianness check was successful! The endianness is set to", endianness))
	return(endianness)
}
setGeneric("GAZE.build", function(gaze.obj, mode) standardGeneric("GAZE.build"))
setMethod("GAZE.build", "GAZE",
	function(gaze.obj, mode){
		stare_opt <- gaze.obj@STARE_options
		if(tolower(mode) == "static"){
### Triangulate ###
			peak_flag <- F

			input.x <- readr::read_tsv(paste0(stare_opt@o, "/Gene_TF_matrices/", stare_opt@o, "_TF_Gene_Affinities.txt.gz"))
			### Remove the peak-related features: NumPeaks, AvgPeakDistance, AvgPeakSize
			input.x <- input.x[, -c((ncol(input.x)-2):ncol(input.x))]


			## Merge X and Y variables based on their gene names and keep track of the lost genes
			merged_input <- merge(gaze.obj@ID_mapping, input.x, by.x= "gene_id", by.y= "geneID")
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
			genediff <- setdiff(merged_input$gene_name, x_y_hits)
			if(length(genediff)){
				lost_genes <- data.frame(genes= setdiff(merged_input$gene_name, x_y_hits), source= "Gene-TF input matrix")
			}
			genediff <- setdiff(rownames(gaze.obj@metacell_mat), x_y_hits)
			if(length(genediff)){
				lost_genes <- data.frame(genes= setdiff(rownames(gaze.obj@metacell_mat), x_y_hits), source= "Gene expression input matrix (gaze.object)")
			}
			gaze.obj@lost_genes <- lost_genes
			gaze.obj@TF_cell_mat <- X

			print(paste("dim(X)", dim(X)))
			print(paste("dim(Y)", dim(Y)))
			tri_res <- triangulate_main(input.x= X, input.y= Y, percent= .8, output= "triangulate_out_GAZE")
			reduced_ID_mapping <- merge(gaze.obj@ID_mapping, gaze.obj@metacell_mat, by.x= "gene_name", by.y= 0)[, 1:2]
			gene_res <- list()

			run_one_level_learning <- function(TGL.model_m, gene_names){
				print(class(gaze.obj))
				gene_hit <- which(gene_names == gaze.obj@ID_mapping$gene_name)
				if(!length(gene_hit)){
					gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "STARE output (before per-gene learning starting)"))
					return(NULL);
				}
				gaze_gene_hit <- which(reduced_ID_mapping$gene_name %in% gene_names)
				if(length(gaze_gene_hit)){
					if(sum(gaze.obj@metacell_mat[which(rownames(gaze.obj@metacell_mat) == reduced_ID_mapping$gene_name[gaze_gene_hit]), ]) & sum(TGL.model_m)){
						gene_res[[gene_names]] <- one_level_learnings(TGL.model_m= TGL.model_m, partition_m=  log2(1 + t(gaze.obj@metacell_mat)), gene_name= gene_names, id_mapping= gaze.obj@ID_mapping)
							return(gene_res)
					}else if(sum(TGL.model_m)){
						gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "Zero expression across all cells"))
					}else{
						gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "Zero affinity across all cells and TFs"))
					}
				}else{
					gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "Gene did not exist in the expression data!"))
				}
			}
			#res <- one_level_learnings(t(tri_res$TGL.model$B), partition_m= t(rbind(tri_res$partition$train$y, tri_res$partition$test$y)), gene_name= "A1BG", gaze.obj@ID_mapping)
			gene_names <- rownames(gaze.obj@metacell_mat)

			library(future.apply)
			options(future.globals.maxSize= 5.5 *1024^3)
			#plan(list(tweak(multicore, workers = 10), tweak(multicore, workers = 5)))
			memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern= T)) ## in kB
			gaze_size <- object.size(gaze.obj)/10^3
			workers <- floor(memfree/(gaze_size * 10))

			print(paste("workers=", workers))
			print(paste("assigning workers:", min((workers -1 ), 20)))
			#plan(tweak(multicore, workers = ifelse(min((workers -1 ), 20) < 1, 1, min((workers -1 ), 20)) ))# 40))
			plan(tweak(multicore, workers = ifelse(min((workers -1 ), 10) < 1, 1, min((workers -1 ), 10)) ))# 40))
			#models <- future.apply::future_lapply(seq(length(gene_names)), function(j){run_one_level_learning(t(tri_res$TGL.model$B[!is.na(rowSums(tri_res$TGL.model$B)), ]), gene_names[j])}, future.seed = 0xBEEF)
			chunk_size <- 10
			total_chunks <- floor(length(gene_names) / chunk_size)
			chunk_rem <- length(gene_names) %% chunk_size
			models <- list()
			models <- future.apply::future_lapply(gene_names, function(j){
												run_one_level_learning(t(tri_res$TGL.model$B[!is.na(rowSums(tri_res$TGL.model$B)), ]), j)
						}, future.seed = 0xBEEF
			)#, future.globals= c("tri_res", "j"))
			names(models) <- gene_names

			if(F){
				for(i in seq(1, total_chunks)){
					models_chunk <- future.apply::future_lapply(gene_names[seq((i - 1) * chunk_size + 1, i * chunk_size)], function(j){
														run_one_level_learning(t(tri_res$TGL.model$B[!is.na(rowSums(tri_res$TGL.model$B)), ]), j)
								}, future.seed = 0xBEEF)#, future.globals= c("tri_res", "j"))
					names(models_chunk) <- gene_names[seq((i - 1) * chunk_size + 1, i * chunk_size)]
					print(paste0("========================== ", i, " ==========================="))
					models <- c(models, models_chunk)
				}
				for(i in seq(chunk_rem)){
					models_chunk <- future.apply::future_lapply(gene_names[seq((i - 1) * chunk_size + 1, i * chunk_size)], function(j){
														run_one_level_learning(t(tri_res$TGL.model$B[!is.na(rowSums(tri_res$TGL.model$B)), ]), i)
								}, future.seed = 0xBEEF)
					names(models_chunk) <- gene_names[seq((i - 1) * chunk_size + 1, i * chunk_size)]
					models <- c(models, models_chunk)
				}
			}
		}else{
### STARE ###


			## Set if the peak-related features should be included:
			peak_flag <- T
			## check Endianness
			STARE_tokens <- unlist(strsplit(stare_opt@o, "/"))
			endianness <- check_endianness(paste0(stare_opt@o, "/Gene_TF_matrices/", STARE_tokens[length(STARE_tokens)]))

			## get the TF, cell, and gene names from STARE metadata file
			print(paste("Reading the Reshape metadata file:", paste0(stare_opt@o, "/Gene_TF_matrices/", STARE_tokens[length(STARE_tokens)], "_Reshape_Meta.txt")))
			meta <- readLines(paste0(stare_opt@o, "/Gene_TF_matrices/", STARE_tokens[length(STARE_tokens)], "_Reshape_Meta.txt"))
			print(paste("Done reading the Reshape metadata file:", paste0(stare_opt@o, "/Gene_TF_matrices/", STARE_tokens[length(STARE_tokens)], "_Reshape_Meta.txt")))

			number_of_entries <- as.numeric(strsplit(meta[1], "\t")[[1]])
			TF_names <- strsplit(meta[2], "\t")[[1]]
			cell_names <- strsplit(meta[3], "\t")[[1]]
			gene_names <- strsplit(meta[4], "\t")[[1]]
			
			if(length(strsplit(gene_names[1], "\\.")[[1]]) > length(strsplit(rownames(gaze.obj@metacell_mat)[1], "\\.")[[1]])){
				gene_names <- sapply(gene_names, function(g) strsplit(g, "\\.")[[1]][1])
			}else{
				rownames(gaze.obj@metacell_mat) <- sapply(rownames(gaze.obj@metacell_mat), function(g) strsplit(g, "\\.")[[1]][1])
			}
			
			## reshaped_bin <- file(paste0(stare_opt@o, "_stare2reshape.bin", "rb")) ## This is the file for my wrapper version. But now that Dennis incorporated the reshape into STARE, I can get the correct shape directly out of STARE (without requiring the intermediate wrapper script).
			all_ps <- list()
			all_qs <- list()
			feature_mat_list <- list()
			gene_res <- list()
			var_thresh <- quantile(apply(gaze.obj@metacell_mat, 1, FUN= var), probs= seq(0, 1, .1))[5] ## 40% quantile
			gaze.obj@ID_mapping$gene_id <- sapply(gaze.obj@ID_mapping$gene_id, function(g) strsplit(g, "\\.")[[1]][1]) ## remove the subversion as it's causing a lot of trouble for theIHEC bulk data
			if(substr(rownames(gaze.obj@metacell_mat)[1], 1, 6) == "ENSG00"){
				reduced_ID_mapping <- merge(gaze.obj@ID_mapping, gaze.obj@metacell_mat, by.x= "gene_id", by.y= 0)[, 1:2]
			}else{
				reduced_ID_mapping <- merge(gaze.obj@ID_mapping, gaze.obj@metacell_mat, by.x= "gene_name", by.y= 0)[, 1:2]
			}


			run_one_level_learning <- function(TGL.model_m, gene_names){
				print(class(gaze.obj))
				print(gene_names)
				gene_hit <- which(gene_names == gaze.obj@ID_mapping$gene_id)
				if(!length(gene_hit)){
					gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "STARE output (before per-gene learning starting)"))
					return(NULL);
				}
				gaze_gene_hit <- which(reduced_ID_mapping$gene_id %in% gene_names)
				if(length(gaze_gene_hit)){
					#if(sum(gaze.obj@metacell_mat[which(rownames(gaze.obj@metacell_mat) == reduced_ID_mapping$gene_name[gaze_gene_hit]), ]) & sum(feature_mat)){
					if(sum(gaze.obj@metacell_mat[which(rownames(gaze.obj@metacell_mat) == reduced_ID_mapping$gene_name[gaze_gene_hit]), ]) & sum(TGL.model_m)){
						gene_res[[gene_names]] <- one_level_learnings(TGL.model_m= log2(1+ TGL.model_m), partition_m=  log2(1 + t(gaze.obj@metacell_mat)), gene_name= gene_names, id_mapping= gaze.obj@ID_mapping)
							return(gene_res)
						print(range(TGL.model_m))
					}else if(sum(TGL.model_m)){
						gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "Zero expression across all cells"))
					}else{
						gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "Zero affinity across all cells and TFs"))
					}
				}else{
					gaze.obj@lost_genes <- rbind(gaze.obj@lost_genes, data.frame(genes= gene_names, source= "Gene did not exist in the expression data!"))
				}
			}
			### perform parallelization of the binary file
			chunk_size <- 10 ## number of genes to be read from file
			rem_genes <- length(gene_names) - (length(gene_names) %/%  chunk_size) * chunk_size
			library(future.apply)
			#options(future.globals.maxSize= 2 *1024^3)
			## Has to be much larger for the MPAL dataset.. Increasing from 2 *1024^3 to 4.5 *1024^3
			options(future.globals.maxSize= 4.5 *1024^3)
			#plan(list(tweak(multicore, workers = 10), tweak(multicore, workers = 5)))
			memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern= T)) ## in kB
			gaze_size <- object.size(gaze.obj)/10^3
			workers <- floor(memfree/(gaze_size * 10))

			print(paste("workers=", workers))
			print(paste("assigning workers:", min((workers -1 ), 20)))
			#plan(tweak(multicore, workers = ifelse(min((workers -1 ), 20) < 1, 1, min((workers -1 ), 20)) ))# 40))
			plan(tweak(multicore, workers = ifelse(min((workers -1 ), 10) < 1, 1, min((workers -1 ), 10)) ))# 40))

			models <- vector(mode= "list", length= ceiling(length(gene_names) / chunk_size))
			print(paste("Reading the binary file:", paste0(stare_opt@o, "/Gene_TF_matrices/", STARE_tokens[length(STARE_tokens)], "_Reshape_Binary.bin.gz")))
			reshaped_bin <- gzfile(paste0(stare_opt@o, "/Gene_TF_matrices/", STARE_tokens[length(STARE_tokens)], "_Reshape_Binary.bin.gz"), "rb")
			i <- 0
			if(chunk_size < length(gene_names)){
				for(i in seq(length(gene_names) %/%  chunk_size)){
					print(paste0("Current chunk: ", i))
					feature_mat_list <- list()
					for(chnk in seq(chunk_size)){
						feature_mat <- matrix(readBin(reshaped_bin, n= length(TF_names) * length(cell_names), size= 4, what= "double", endian= endianness), byrow= T, nrow= length(cell_names), ncol= length(TF_names))
						rownames(feature_mat) <- cell_names
						colnames(feature_mat) <- TF_names
						if(!peak_flag){
							### remove the last three features that related to peak info
							feature_mat <- feature_mat[, -seq(ncol(feature_mat) - 2, ncol(feature_mat))]
						}
						feature_mat_list[[chnk]] <- feature_mat


					}
					gene_chunk <- gene_names[seq((i - 1) * chunk_size + 1, i*chunk_size)]

					models_chunk <- future.apply::future_lapply(seq(length(feature_mat_list)), function(j){
												run_one_level_learning(feature_mat_list[[j]], gene_chunk[j])
						}, future.seed = 0xBEEF)
					print(length(gene_chunk))
					print(length(models_chunk))
					#gaze.obj@models <- models_chunk
					#gaze.obj@TF_cell_mat <- feature_mat_list[[1]]
					#return(gaze.obj)
					names(models_chunk) <- gene_chunk

					models[[i]] <- models_chunk
					if(i == 5000){
						saveRDS(models, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_5000_models.rds"))
					}
					if(i == 60000){
						saveRDS(models, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_60000_models.rds"))
					}
				}
			}
			saveRDS(gaze.obj, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_GAZE_object.rds"))
			gene_chunk <- gene_names[seq(i * chunk_size + 1, length(gene_names))]
			feature_mat_list <- list()
			for(chnk in seq(rem_genes)){
				feature_mat <- matrix(readBin(reshaped_bin, n= length(TF_names) * length(cell_names), size= 4, what= "double", endian= endianness), byrow= T, nrow= length(cell_names), ncol= length(TF_names))
				rownames(feature_mat) <- cell_names
				colnames(feature_mat) <- TF_names
				if(!peak_flag){
					### remove the last three features that related to peak info
					feature_mat <- feature_mat[, -seq(ncol(feature_mat) - 2, ncol(feature_mat))]
				}
				feature_mat_list[[chnk]] <- feature_mat
			}

			models_chunk <- future.apply::future_lapply(seq(length(feature_mat_list)), function(j){
										run_one_level_learning(feature_mat_list[[j]], gene_chunk[j])
				}, future.seed = 0xBEEF)
			names(models_chunk) <- gene_chunk

			models[[i + 1]] <- models_chunk

			close(reshaped_bin)
		}
	#####################################################################
		model_hits <- which(unlist(sapply(seq(length(models)), function(m) sapply(seq(length(models[[m]])), function(mm)return(class(models[[m]][[mm]]))))) == "list")

		gene_res <- NULL
		ctr <- 1;
		if(tolower(mode) == "static"){
			for(m in seq(length(models))){
				if(class(models[[m]]) == "list"){
					gene_res <- c(gene_res, models[[m]])
					ctr <- ctr + 1
				}
			}
		}else{
			for(m in seq(length(models))){
				for(mm in seq(length(models[[m]]))){
					if(class(models[[m]][[mm]]) == "list"){
						gene_res <- c(gene_res, models[[m]][[mm]])
						ctr <- ctr + 1
					}
				}
			}
		}
		gaze.obj@models <- gene_res
		#saveRDS(gaze.obj, paste0(gaze.obj@STARE_options@o, "_peakFeatures_", as.character(peak_flag), "_GAZE_object.rds"))
		saveRDS(gaze.obj, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_GAZE_object.rds"))
	#####################################################################
		null_idx <- sapply(gene_res, function(g)is.null(g$models))
		test_cor <- list()
		train_cor <- list()
		coef_tmp <- glmnet::coef.glmnet(gene_res[!null_idx][[1]]$models, s= gene_res[!null_idx][[1]]$alphas$lambda.min)

		coef_map <- matrix(NA, nrow= length(gene_res[!null_idx]), ncol= length(coef_tmp))
		rownames(coef_map) <- names(gene_res[!null_idx])
		colnames(coef_map) <- rownames(coef_tmp)

		cors_all <- lapply(gene_res[!null_idx], function(gg)cor(gg$d$test$y, predict(gg$models, gg$d$test$x, s= gg$alphas$lambda.min), method= "p"))
		cors_all_s <- lapply(gene_res[!null_idx], function(gg)cor(gg$d$test$y, predict(gg$models, gg$d$test$x, s= gg$alphas$lambda.min), method= "s"))

		cors_train_all <- lapply(gene_res[!null_idx], function(gg)cor(gg$d$train$y, predict(gg$models, gg$d$train$x, s= gg$alphas$lambda.min), method= "p"))
		cors_train_all_s <- lapply(gene_res[!null_idx], function(gg)cor(gg$d$train$y, predict(gg$models, gg$d$train$x, s= gg$alphas$lambda.min), method= "s"))

		cors_df <- data.frame(gene= names(cors_all), cor= unlist(cors_all), method= "Pearson", partition= "test"); cors_df <- rbind(cors_df, data.frame(gene= names(cors_all), cor= unlist(cors_all_s), method= "Spearman", partition= "test"));
		cors_df <- rbind(cors_df, data.frame(gene= names(cors_train_all), cor= unlist(cors_train_all), method= "Pearson", partition= "train")); cors_df <- rbind(cors_df, data.frame(gene= names(cors_train_all), cor= unlist(cors_train_all_s), method= "Spearman", partition= "train"));
		#write.table(cors_df, paste0(gaze.obj@STARE_options@o, "_test_correlation.txt"), quote= F, row.names= F)
		#write.table(cors_df, paste0(gaze.obj@STARE_options@o, "_peakFeatures_", as.character(peak_flag), "_correlation.txt"), quote= F, row.names= F)
		write.table(cors_df, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_correlation.txt"), quote= F, row.names= F)
		gaze.obj@model_accuracy <- cors_df

		preds_all <- NULL
		responses_all <- NULL

		for(j in seq(length(gene_res[!null_idx]))){
			preds <- predict(gene_res[!null_idx][[j]]$models, gene_res[!null_idx][[j]]$d$test$x, s= gene_res[!null_idx][[j]]$alphas$lambda.min)
			train_preds <- predict(gene_res[!null_idx][[j]]$models, gene_res[!null_idx][[j]]$d$train$x, s= gene_res[!null_idx][[j]]$alphas$lambda.min)
			preds_all <- c(preds_all, preds)
			responses_all <- c(responses_all, gene_res[!null_idx][[j]]$d$test$y)
			test_cor[[names(gene_res[!null_idx][j])]] <- cor(preds, gene_res[!null_idx][[j]]$d$test$y)
			train_cor[[names(gene_res[!null_idx][j])]] <- cor(train_preds, gene_res[!null_idx][[j]]$d$train$y)

			coef_map[j, ] <- glmnet::coef.glmnet(gene_res[!null_idx][[j]]$models, s= gene_res[!null_idx][[j]]$alphas$lambda.min)[, 1]
		}
		### Fill the GAZE object and save it ###
		gaze.obj@coefficient_map <- coef_map
		gaze.obj@shap_values <- lapply(gene_res[!null_idx], function(gr)gr$shap)
		gaze.obj@models <- list()
		#saveRDS(gaze.obj, paste0(gaze.obj@STARE_options@o, "_peakFeatures_", as.character(peak_flag), "_GAZE_object_noModels.rds"))
		saveRDS(gaze.obj, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_GAZE_object_noModels.rds"))
		######################################

		informative_genes <- which(rowSums(abs(coef_map)) > ncol(coef_map))
		informative_genes <- which(test_cor > .7)
		rownames(gaze.obj@ID_mapping) <- gaze.obj@ID_mapping$gene_id
		rownames(coef_map) <- gaze.obj@ID_mapping[rownames(coef_map), ]$gene_name
		#coef_map_filtered <- coef_map[informative_genes, -seq(ncol(coef_map), 1)[seq(3)]] ## remove peak-related features
		coef_map_filtered <- coef_map[informative_genes, -1]

		paletteLength <- 50
		myBreaks <- c(seq(min(coef_map_filtered), 0, length.out=ceiling(paletteLength/2) + 1),
								seq(max(coef_map_filtered)/paletteLength, max(coef_map_filtered), length.out=floor(paletteLength/2)))
		myColor <- colorRampPalette(c("black", "white", "green"))(paletteLength)


		pdf(paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag),"_prediction_results.pdf"))
		hist(unlist(test_cor), main= "prediction correlation on test", xlab= "Pearson correlation")
		par(mfrow= c(5,2))
		for(i in seq(nrow(coef_map_filtered))[seq(30)]){
			plot(coef_map_filtered[i, ], main= rownames(coef_map_filtered)[i], pch= 20, xlab= "", ylab= "")
		}
		pheatmap::pheatmap(coef_map_filtered[order(unlist(test_cor[informative_genes]), decreasing= T)[seq(20)], ], cluster_rows= T, cluster_cols= T, color=myColor, breaks=myBreaks, main= "regression coefficients for top genes")
		dev.off()
		####### scale coefficients ##########
		#coef_filtered <- coef_map[informative_genes, -c(1, seq(ncol(coef_map), 1)[seq(3)])] ## remove intercept and the peak-related features
		coef_filtered <- coef_map[informative_genes, -c(1)] ## remove intercept only, because at some point we decided to run the models without the peak-related features

		useless_TFs <- which(colSums(abs(coef_filtered)) == 0)
		useless_TFs2 <- which(colSums(abs(coef_filtered[order(unlist(test_cor[informative_genes]), decreasing= T)[seq(20)],])) == 0)

		coef_scaled_filtered <- t(apply(abs(coef_filtered[, -useless_TFs2]), 1, FUN= function(i)return(if(min(i) == max(i)){return(i)}else{(i-min(i))/(max(i) - min(i))})))
		coef_scaled_filtered <- coef_scaled_filtered * sign(coef_filtered[, -useless_TFs2])

		myBreaks <- c(seq(min(coef_scaled_filtered, na.rm= T), 0, length.out=ceiling(paletteLength/2) + 1),
								seq(max(coef_scaled_filtered, na.rm= T)/paletteLength, max(coef_scaled_filtered, na.rm= T), length.out=floor(paletteLength/2)))

		low_coef_cutoff <- .2 #.3
		#zero_idx <- which(colSums(abs(coef_scaled_filtered[seq(20), ])) == 0)
		zero_idx <- which(colSums(abs(coef_scaled_filtered[seq(20), ])) < low_coef_cutoff)

		coef_scaled_filtered_final <- coef_scaled_filtered[seq(min(20, nrow(coef_scaled_filtered))), -zero_idx]
		pdf(paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag),"_scaled_coefficients_uselessTFsRemoved_", low_coef_cutoff, ".pdf"))
		pheatmap::pheatmap(coef_scaled_filtered_final[!(rowSums(abs(coef_scaled_filtered_final)) == 0), ], fontsize_col= 4, cluster_rows= T, cluster_cols= T, color=myColor, breaks=myBreaks, main= "scaled regression coefficients for top genes")
		dev.off()
		##########################################
	#####################################################################
		####### Run the statistical tests #######
		table_cols <- lapply(gaze.obj@metadata_metacell, function(l)unique(l))
		summarize_shaps <- function(i, md, table_cols){
			md_vals <- table_cols[[md]];
			res <- NULL;
			for(j in md_vals){
				hit <- which(metadata_metacell_reordered[, md] == j);
				res <- cbind(res, colMeans(gaze.obj@shap_values[[i]][hit, ]))
			};
			return(res)
		}
		gene_summary <- list();
		for(g in names(gaze.obj@shap_values)){
					gene_summary_i <- NULL;
					for(md in seq(length(table_cols))){
						gene_summary_i <- cbind(gene_summary_i, summarize_shaps(g, md, table_cols))
					}
					gene_summary[[g]] <- gene_summary_i
		}
		print("Started computing the differential tests!")
		source("compute_stat_test_indivTFs.R")
		test_res <- compute_TF_view_tests(gaze.obj)
		gaze.obj@stat_test_gene <- test_res
		crossing_res <- do.call(tidyr::crossing, gaze.obj@metadata_metacell)
		saveRDS(gaze.obj, paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_GAZE_object_noModels.rds"))
		#### Prepare the test table for Dennis' heat map and other plots
		stat_trimmed <- lapply(gaze.obj@stat_test_gene, function(st){if(!nrow(st)) return(NULL); return(st[, c(1, 4, 7, 8, 9, 10)])})
		stat_trimmed_unlist <- do.call(rbind, stat_trimmed)
		stat_trimmed_unlist_uuniq <- stat_trimmed_unlist[!duplicated(stat_trimmed_unlist) , ]
		write.table(stat_trimmed_unlist_uuniq, paste0(gaze.obj@gaze.options@outout_file, "_peakFeatures_", as.character(peak_flag), "_GAZE_stat_table_trimmed.txt"), quote= F, sep= "\t")
		##########################################
		#pdf("GAZE_debug2.pdf")
		pdf(paste0(gaze.obj@gaze.options@output_file, "_peakFeatures_", as.character(peak_flag), "_GAZE_debug.pdf"))
		for(i in seq(length(all_ps), 2)){
			print(ggplot2::ggplot(all_ps[[i]], aes(x= X1, y= X2, shape= cell)) + geom_point(aes(colour= cell)) + scale_shape_manual(values= seq(length(unique(all_ps[[i]]$cell)))) + ggtitle(paste0("gene: ", i, " or ", gene_names[i])))
			print(ggplot(data.frame(correlation= ), aes(x= correlation)) + geom_histogram() + ggtitle(paste0("gene: ", i, " or ", gene_names[i])))

		}
		dev.off()
		if(stare_res){
			stop("STARE module failed to finish successfully. Please review the corresponding arguments more carefully based on the error messages from the STARE run!")
		}

	## read the reshape matrix and load it into the GAZE object
	#	reshaped <- read.table(paste0(stare_opt@o, "/Stacked_Matrix/", stare_opt@o, "_StackedAffinityMatrix.txt.gz"))
	}
)
test_mode <- F

if(test_mode){

## Test MetaCellaR ##
	file_name <- "../data/3008_.Rds"
	file_name <- "/projects/triangulate/work/mojitoo/rna_atac_merged_coembedded_mojitoo.rds"
	# mc.opts <- MetaCellaR.options(assay_slot= "sc-ATAC", condition= "Age", output_file= "test_out", RNA_slot= "assays\\$RNA@counts", celltype_slot= "meta.data\\$celltype") #old implementation
	mc.opts <- MetaCellaR.options(assay= "meta.data\\$datasets2", ATAC= "assays\\$peaks@counts", output= "test_out", RNA= "assays\\$RNA@counts", celltype= "meta.data\\$celltype", reduction= "MOJITOO_UMAP")

	gaze.opt <- GAZE.options(ATAC_exists= T, MetacellaR= T, gtf_file= "/projects/abcpp/work/base_data/gencode.v38lift37.basic.annotation.gtf")
#	gaze.obj <- new("GAZE", gaze.options= gaze.opt, seurat_file= file_name)
	gaze.obj <- CreateGAZEobject(gaze.options= gaze.opt, seurat_file= file_name, metacell.options= mc.opts)

	#gaze.obj <- MetaCellaR.run(gaze.obj, mc.opts)

	## Add STARE's output to the GAZE object
	#input.x <- read.csv("../data/scMTL_heartecCellSumKmedmeans3008_notImputed_epigenetic_feature_doubleReduced.txt")
	## Now with the binary and reshaped STARE, I can directly get the output and run the pipeline, so the read_tsv function below is not required anymore (commenting it out)!
	#input.x <- readr::read_tsv("../../Test_K562_10xRandom_GAZE4/Gene_TF_matrices/Test_K562_10xRandom_GAZE4_TF_Gene_Affinities_Apple.txt.gz")
	gtf <- as.data.frame(rtracklayer::import(gaze.opt@gtf_file))
	gtf_genes <- subset(gtf, type == "gene")


	## Merge X and Y variables based on their gene names and keep track of the lost genes
	ID_mapping <- gtf_genes[, c("gene_id", "gene_name")]
	gene_id_noversion <- sapply(ID_mapping$gene_id, function(i) unlist(strsplit(i, "\\."))[1])
	ID_mapping$gene_id <- gene_id_noversion
	merged_input <- merge(ID_mapping, input.x, by.x= "gene_id", by.y= "geneID")
	x_y_merged <- merge(gaze.obj@metacell_mat, merged_input, by.x= 0, by.y= "gene_name")
	ENS_IDs <- x_y_merged[, ncol(gaze.obj@metacell_mat) + 2]
	## spot the duplicated row names
	dup.hit <- which(duplicated(ENS_IDs))
	if(length(dup.hit)){
		x_y_merged <- x_y_merged[-dup.hit, ]
		ENS_IDs <- ENS_IDs[-dup.hit]
		merged_input <- merged_input[-dup.hit, ]
	}
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
	stop()
	tri_res <- triangulate_main(input.x= X, input.y= Y, percent= .8, output= "triangulate_out_GAZE")
	gaze.obj@TF_cell_mat <- tri_res$TGL.model$B

	stop()

	sec_res <- second_level_learnings(TGL.model_m= tri_res$TGL.model, partition_m= tri_res$partition, ID_mapping)

	gaze.obj@models <- sec_res[["models"]]
	gaze.obj@shap_values <- sec_res[["shap"]]
	gaze.obj@ID_mapping <- ID_mapping

##########
	#gaze.opt <- new("GAZE.options", ATAC_exists= ATAC_exists, MetacellaR= MetacellaR,filtering= filtering)

	#gaze <- new("GAZE", gaze.opt, expr_mat= expr_mat, metacell_mat= as.matrix(metacell_mat), condition= condition, celltype= as.character(celltype), models= models$models, shap_values= models$shap)

}

