
#' Project UKB individual-level genotypes onto 1kg principal component space
#'
#' @param bfile a plink1 binary fileset containing genotypes for individuals whose ancestral superpopulation you want to estimate.
#' @param out a path to the output file
#' @param path_to_plink_exe the path to your plink2 executable
#' @param n_threads number of computing threads to use, defaults to 1
#'
#' @return a '.sscore' file with the projected PC values for individuals in your dataset in 1kg PC space.
#' @export
#'
#' @examples
#' bfile = "ukb_merged_chroms"
#' path_to_plink = "plink.exe"
#' project_pcs()

project_pcs = function(
  bfile = "ukb_merged_chroms",
  out = "projected_pcs",
  n_threads = 1,
  path_to_plink_exe = "module unload plink; module load plink/2.0-20200328; plink2 "){

  # arg checking
  if(!(file.exists(paste0(bfile,".bed")) & file.exists(paste0(bfile,".bim")) & file.exists(paste0(bfile,".fam")))){
    stop(bfile, "is not a valid PLINK1 .bed .bim .fam fileset")
  }

  if(!(is.character(bfile) & is.character(out) & is.character(path_to_plink_exe) & is.numeric(n_threads))){
      stop("Argument ", bfile," must be a valid path to a PLINK1 binary fileset \n",
           "Argument ", out," must be a valid output file path \n",
           "Argument ", path_to_plink_exe," must be a valid path to a PLINK1 executable \n",
           "Argument ", threads," must be a valid integer with the number of computing threads")
  }

  # remove whitespace from args
  path_to_plink_exe = stringr::str_trim(path_to_plink_exe,side = "right")
  bfile = stringr::str_trim(bfile)

  # load references
  message("Loading reference PCs and SNP loadings.")
  ref_snp_freqs = data_list$ref_freqs
  pc_loadings = data_list$ref_loadings

  # write to files for plink
  readr::write_tsv(ref_snp_freqs,"ref_snp_freqs.tsv")
  readr::write_tsv(pc_loadings,"pc_loadings.tsv")

  # run plink
  cmd = paste0(path_to_plink_exe," --bfile ",bfile," --read-freq ref_snp_freqs.tsv --score pc_loadings.tsv 2 5 header-read no-mean-imputation variance-standardize --threads ",n_threads," --score-col-nums 6-55 --out ",out)
  message("Trying to run plink")
  system(cmd)

  # output message
  outfile = paste0(out,".sscore")
  if(file.exists(outfile)){
    message("Finished successfully. Projected PC scores written to ",outfile)
  } else {
    stop("Something went wrong.")
  }
}

#' Transform projected PCs onto 1kg scale
#'
#' @param projected_pcs transformed PC scores on the same co-ordinates as the 1kg reference samples. Default output from project_pcs().
#'
#' @return A data frame with raw PC scores on the same co-ordinates as the 1kg reference samples. Default output from project_pcs().
#' @export
#'
#' @examples
#'
#' # First - project PCs:
#' project_pcs(out = "projected_pcs")
#' # Then run:
#' new_pcs = transform_pcs("projected_pcs.sscore")

transform_pcs = function(projected_pcs = "projected_pcs.sscore"){

  # references
  eigenvals = data_list$ref_eigenval
  if(!is.data.frame(eigenvals)){
      stop("Issue loading internal data. Aborting. ")
  }

  # read in files
  message("Reading in projected scores")
  new_scores = readr::read_tsv(projected_pcs)
  message("Read in projected scores")

  # calculate scale factors from PCs
  message("Calculating scale factors")
  scale_factors = 1/(sqrt(eigenvals)/-2)

  # transform scores
  message("Transforming scores by scale factors")
  just_pc_scores = new_scores[,grepl("PC",colnames(new_scores))]
  transformed_pcs = lapply(c(1:ncol(just_pc_scores)),function(i){
    just_pc_scores[,i] * scale_factors[i,]
  })
  transformed_pcs = do.call("cbind",transformed_pcs) %>% data.frame()
  transformed_pcs$IID = new_scores$IID
  transformed_pcs = transformed_pcs %>% select(IID,1:ncol(transformed_pcs))
  message("Finished")
  return(transformed_pcs)
}



#' Title Plot PCs with 1kg reference
#'
#' @param projected_pcs a data frame with transformed PC scores. Default output from transform_pcs().
#' @param save_plot TRUE/FALSE, whether to save the plot to file
#' @param outfile the output filename for the plot
#'
#' @return a plot with the new samples projected on 1kg PC space
#' @export
#'
#' @examples
#' # First - project PCs:
#' project_pcs(out = "projected_pcs")
#' # Then run:
#' new_pcs = transform_pcs("projected_pcs.sscore")
#' # Then plot with:
#' plot_with_1kg(projected_pcs = new_pcs, save_plot = TRUE, outfile = "niceplot")

plot_with_1kg = function(projected_pcs = new_pcs, save_plot = TRUE, outfile = "projected_pcs"){

  # check args
  if(!is.data.frame(projected_pcs)){
    stop("Aborting. Argument ",projected_pcs," must be a data frame with individual PC projections.")
  }

  if(!is.logical(save_plot)){
    stop("Aborting. Argument ",save_plot," must be TRUE or FALSE.")
  }

  if(!is.character(outfile)){
    stop("Aborting. Argument ",outfile," must be a valid path for the output plot.")
  }

  # references pheno
  pheno = data_list$ref_pheno
  ref_pcs = data_list$ref_eigenvec
  # join
  ref_pcs_pheno = ref_pcs %>% dplyr::rename(sample = IID) %>% dplyr::left_join(pheno,by="sample") %>% rename(IID = sample) %>% select(IID,contains("PC"),super_pop)

  # strip _AVG from newpcs name
  colnames(projected_pcs) = stringr::str_remove(colnames(projected_pcs),"_AVG")

  # covert IID to char
  projected_pcs$IID = as.character(projected_pcs$IID)

  # define target superpop as unknown
  projected_pcs$super_pop = "Unknown"

  # join
  joint_df = projected_pcs %>%
    dplyr::bind_rows(ref_pcs_pheno %>%
                       select(colnames(projected_pcs))
                     )

  p1=ggplot2::ggplot(data=joint_df,aes(PC1,PC2,col=super_pop))+geom_point()+ggtitle("Projected PCs with 1000 Genomes reference")+theme_classic()

  # save to file
  if(save_plot == TRUE){
    message("Saving PC plot to file ",outfile,".png")
    png(filename = paste0(outfile,".png"),height=8,width=8,res=300,units="in")
    print(p1)
    dev.off()
  } else {
    return(p1)
  }
}


#' Predict ancestral superpopulation from individual-level PC projections using a pre-trained classifier from 1kg
#'
#' @param projected_pcs a data frame with transformed PC scores. Default output from transform_pcs().
#' @param outfile a path for the output file
#'
#' @return a ".tsv" file with best-guess ancestry estimates for each individual in the dataset
#' @export
#'
#' @examples
#'
#' # First - project PCs:
#' project_pcs(out = "projected_pcs")
#' # Then run:
#' new_pcs = transform_pcs("projected_pcs.sscore")
#' # Then plot with:
#' plot_with_1kg(projected_pcs = new_pcs, save_plot = TRUE, outfile = "niceplot")
#' # Predict ancestry with:
#' predict_ancestry(new_pcs)

predict_ancestry = function(projected_pcs = new_pcs, outfile = "ancestry_estimates"){

  # check args
  if(!is.data.frame(projected_pcs)){
    stop("Aborting. Argument ",projected_pcs," must be a data frame with individual PC projections.")
  }

  # strip _AVG from newpcs name
  colnames(projected_pcs) = stringr::str_remove(colnames(projected_pcs),"_AVG")

  # predict ancestry
  message("Predicting ancestral superpopulation for each individual in the dataset.")
  predictions = predict(data_list$rf,newdata = projected_pcs)

  # append to IIDs
  pred_output = data.frame("IID" = projected_pcs$IID,"ancestry_projection" = predictions)
  outpath = paste0(outfile,".tsv")
  readr::write_tsv(pred_output,path=outpath)
  message("Finished - written ancestry estimates to ",outpath)
}

#' Plot projected PCs by predicted ancestral superpopulation, and compare with 1kg reference
#'
#' @param projected_pcs a data frame with transformed PC scores. Default output from transform_pcs().
#' @param outfile a path for the output plot
#' @param projected_ancestry a path to a .tsv file with ancestry estimates per person, as generated by predict_ancestry().
#'
#' @return a plot which will display the projected PC values for each individual in your sample, grouped by predicted ancestral superpopulation of origin, alongside reference data from 1kg in the same PC co-ordinate space
#' @export
#'
#' @examples
#' #' # First - project PCs:
#' project_pcs(out = "projected_pcs")
#' # Then run:
#' new_pcs = transform_pcs("projected_pcs.sscore")
#' # Then plot with:
#' plot_with_1kg(projected_pcs = new_pcs, save_plot = TRUE, outfile = "niceplot")
#' # Predict ancestry with:
#' predict_ancestry(new_pcs, outfile = "ancestry_estimates")
#' # Now plot predictions alongside 1kg superpopulations for reference
#' plot_predicted_ancestry(projected_pcs = new_pcs, projected_ancestry = "ancestry_estimates.tsv")

plot_predicted_ancestry = function(projected_pcs = new_pcs,
                                   projected_ancestry = "ancestry_estimates.tsv",
                                   outfile = "ancestry_estimates_plot"){

  # check args
  if(!is.data.frame(projected_pcs)){
    stop("Aborting. Argument ",projected_pcs," must be a data frame with individual PC projections.")
  }
  if(!file.exists(projected_ancestry)){
    stop("Argument ",projected_ancestry," must be a full path to a .tsv file with ancestry estimates, as generated by predict_ancestry()")
  }


  # strip _AVG from newpcs name
  colnames(projected_pcs) = stringr::str_remove(colnames(projected_pcs),"_AVG")

  # read in ancestry predictions
  ancestry_preds = readr::read_tsv(projected_ancestry)

  # join
  joint_df = projected_pcs %>% left_join(ancestry_preds,by="IID")


  # do 1kg plot
  # references pheno
  pheno = data_list$ref_pheno
  ref_pcs = data_list$ref_eigenvec
  # join
  ref_pcs_pheno = ref_pcs %>% dplyr::rename(sample = IID) %>% dplyr::left_join(pheno,by="sample") %>% rename(IID = sample) %>% select(IID,contains("PC"),super_pop)
  p2=ggplot2::ggplot(data=ref_pcs_pheno,aes(PC1,PC2,col=super_pop))+geom_point()+ggtitle("Projected PCs for 1000 Genomes reference")+theme_classic()

  # plot
  p1=ggplot2::ggplot(joint_df,aes(PC1,PC2,col=ancestry_projection))+geom_point()+
    ggtitle("Projected PCs with ancestry estimates")+
    theme_classic()
  message("Saving PC plot to file ",outfile,".png")
  png(filename = paste0(outfile,".png"),height=8,width=16,res=300,units="in")
  print(gridExtra::grid.arrange(p1,p2,nrow=1))
  dev.off()
}




