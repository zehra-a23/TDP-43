library(tidyverse)
library(ggplot2)
my_clean_reader <- function(file){
  df = data.table::as.data.table(janitor::clean_names(data.table::fread(file)))
  return(df)
}
fix_column_names = function(df){
  ori_colnames = df |>  colnames()
  ori_colnames[12] = "baseline_PSI"
  ori_colnames[13] = "contrast_PSI"
  colnames(df) = ori_colnames
  
  return(df)
  
}

read_in_and_combine_data_sets = function(folder_path, suffix = "_annotated_junctions.csv"){
  
  
  
  estimate_files = list.files(folder_path,
                              pattern = suffix,
                              full.names = TRUE)
  print(estimate_files)
  sf = purrr::map(estimate_files,my_clean_reader)
  samp_ids = base::tolower(purrr::simplify(purrr::map(estimate_files, basename)))
  samp_ids = gsub(suffix,"",samp_ids)
  ##add on the name as an additional column
  sf = purrr::map2(sf, samp_ids, ~cbind(.x, comparison = .y))
  sf = purrr::map(sf,fix_column_names)
  sf = data.table::rbindlist(sf) 
  
  return(sf)
  
  
  
}

hallegger_kds = read_in_and_combine_data_sets("/Users/recklessreverie/Desktop/UCL internship/majiq/delta_psi_voila_tsv", 
                                              suffix = "_annotated_junctions.csv")

hallegger_kds = hallegger_kds |>  
  dplyr::select(seqnames:gene_id,comparison,mean_dpsi_per_lsv_junction:lsv_type,paste_into_igv_junction,junc_in_ref,junc_cat,exon_rank_start,exon_rank_end) |> 
  dplyr::select(-lsv_type)
