
#' Project PCs
#'
#' @param bfile a plink1 binary fileset containing genotypes for individuals whose ancestral superpopulation you want to estimate.
#' @param out a path to the output file
#' @param path_to_plink_exe the path to your plink2 executable
#' @param n_threads
#'
#' @return a '.sscore' file with the projected PC values for individuals in your dataset in 1kg PC space.
#' @export
#'
#' @examples

project_pcs = function(bfile = "ukb_merged_chroms", out = "projected_pcs", n_threads = 10, path_to_plink_exe = "module unload plink; module load plink/2.0-20200328; plink2 "){
  # arg checking
  if(!(file.exists(paste0(bfile,".bed")) & file.exists(paste0(bfile,".bim")) & file.exists(paste0(bfile,".fam")))){
    stop(bfile, "is not a valid PLINK1 .bed .bim .fam fileset")
  }



  # load references
  ref_snp_freqs = data_list$ref_freqs
  pc_loadings = data_list$ref_loadings

  # write to files for plink
  readr::write_tsv(ref_snp_freqs,"ref_snp_freqs.tsv")
  readr::write_tsv(pc_loadings,"pc_loadings.tsv")

  # run plink
  cmd = paste0(path_to_plink_exe,"--bfile ",bfile," --read-freq ref_snp_freqs.tsv --score pc_loadings.tsv 2 5 header-read no-mean-imputation variance-standardize --threads ",n_threads," --score-col-nums 6-55 --out ",out)
  message("Trying to run plink")
  system(cmd)
}

#' Transform projected PCs onto 1kg scale
#'
#' @param projected_pcs
#'
#' @return transformed PC scores on the same co-ordinates as the 1kg reference samples
#' @export
#'
#' @examples
transform_pcs = function(projected_pcs = "projected_pcs.sscore"){

  # references
  eigenvals = data_list$ref_eigenval

  # read in files
  message("Reading in projected scores")
  new_scores = readr::read_tsv(projected_pcs)
  message("Read in projected scores")
  # calculate scale factors from PCs
  scale_factors = 1/(sqrt(eigenvals)/-2)

  # transform scores
  just_pc_scores = new_scores[,grepl("PC",colnames(new_scores))]
  transformed_pcs = lapply(c(1:ncol(just_pc_scores)),function(i){
    just_pc_scores[,i] * scale_factors[i,]
  })
  transformed_pcs = do.call("cbind",transformed_pcs) %>% data.frame()
  transformed_pcs$IID = new_scores$IID
  transformed_pcs = transformed_pcs %>% select(IID,1:ncol(transformed_pcs))
  return(transformed_pcs)
}



#' Title Plot PCs with 1kg reference
#'
#' @param projected_pcs projected PC scores - a '.sscore' file produced by project_pcs & transformed with transform_pcs
#'
#' @return a plot with the new samples projected on 1kg PC space
#' @export
#'
#' @examples
plot_with_1kg = function(projected_pcs = newpcs){

  # references pheno
  pheno = data_list$ref_pheno
  ref_pcs = data_list$ref_eigenvec
  # join
  ref_pcs_pheno = ref_pcs %>% dplyr::rename(sample = IID) %>% dplyr::left_join(pheno,by="sample") %>% rename(IID = sample) %>% select(IID,contains("PC"),super_pop)

  # strip _AVG from newpcs name
  colnames(projected_pcs) = stringr::str_remove(colnames(projected_pcs),"_AVG")

  # covert IID to char
  projected_pcs$IID = as.character(projected_pcs$IID)

  # superpop
  projected_pcs$super_pop = "Unknown"

  # join
  joint_df = projected_pcs %>% dplyr::bind_rows(ref_pcs_pheno %>% select(colnames(projected_pcs)))

  ggplot2::ggplot(data=joint_df,aes(PC1,PC2,col=super_pop))+geom_point()
  +ggtitle("Projected PCs with 1000 Genomes reference")

}


#' Title
#'
#' @param projected_pcs
#'
#' @return
#' @export
#'
#' @examples
predict_ancestry = function(projected_pcs = newpcs){

  # strip _AVG from newpcs name
  colnames(projected_pcs) = stringr::str_remove(colnames(projected_pcs),"_AVG")

  # predict ancestry
  predictions = predict(data_list$rf,newdata = projected_pcs)

  # append to IIDs
  pred_output = data.frame("IID" = projected_pcs$IID,"ancestry_projection" = predictions)

  return(pred_output)
}

#' Title
#'
#' @param projected_pcs
#' @param projected_ancestry
#'
#' @return
#' @export
#'
#' @examples
plot_predicted_ancestry = function(projected_pcs = newpcs, projected_ancestry = pred_output){
  # strip _AVG from newpcs name
  colnames(projected_pcs) = stringr::str_remove(colnames(projected_pcs),"_AVG")

  # join
  joint_df = projected_pcs %>% left_join(pred_output,by="IID")

  ggplot2::ggplot(joint_df,aes(PC1,PC2,col=ancestry_projection))+geom_point()+
    ggtitle("Projected PCs with ancestry estimates")

}

#' Title
#'
#' @return
#' @export
#'
#' @examples
compare_1kg_with_predicted_ancestry = function(){
  p1=plot_with_1kg()
  p2=plot_predicted_ancestry()
  gridExtra::grid.arrange(p1,p2)
}
