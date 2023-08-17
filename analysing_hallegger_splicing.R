    library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)
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
# Find all the significantly spliced junctions, use probability_changing > 0.9 and the absolute value mean_dpsi_per_lsv_junction >= 0.2 
significantly_spliced_junctions <- hallegger_kds |>
  filter(probability_changing > 0.9, mean_dpsi_per_lsv_junction >= 0.2 )
significantly_spliced_junctions
# Which comparison has the most number of significant junctions, use the n_distinct function and the paste_into_igv_junction column 
number_of_sig_junctions <- significantly_spliced_junctions |>
  group_by(comparison) |>
  summarize(n_junctions = n_distinct(paste_into_igv_junction))   
# Make a plot of the number of significant junctions by the comparison
number_of_sig_junctions |>
  mutate(comparison = forcats::fct_reorder(comparison,n_junctions)) |>
  ggplot(aes(x=comparison, y=n_junctions)) + 
  geom_col() + 
  coord_flip()
# What splicing events are found in every mutation?
number_of_comparisons <- significantly_spliced_junctions |>
  group_by(paste_into_igv_junction) |>
  summarize(n_comparison = n_distinct(comparison))   
number_of_comparisons
# What splicing event occurs most frequently?
 max_splicing_events <- number_of_comparisons %>%
  group_by(n_comparison) %>%
 summarize(max_junction= max(paste_into_igv_junction)) 
max_splicing_events
install.packages("clipr")
library(clipr)
max_splicing_events |> colnames()  |> clipr::write_clip()
# What splicing event occurs most frequently? (using slice_max)
maximum_splicing_events <- number_of_comparisons %>%
  group_by(n_comparison) %>%
slice_max(paste_into_igv_junction, n = 1)
maximum_splicing_events

# Venn diagram
install.packages("ggvenn")
install.packages("ggVennDiagram")
library(ggvenn)
library(ggVennDiagram)

significantly_spliced_junctions |> colnames()  |> clipr::write_clip()



first_comparison <- significantly_spliced_junctions |>
filter(comparison == "notdpkdrescuenotinduceda326p-sitdp43rescuenotinduceda326p") %>%
pull(paste_into_igv_junction) 
 

second_comparison <- significantly_spliced_junctions |>
  filter(comparison == "sitdp43rescueinduceda326p-sitdp43rescuenotinduceda326p") %>%
pull(paste_into_igv_junction) 

junction_list <- list(junctions_normal_tdpkd = first_comparison, junctions_mutant_tdp= second_comparison)
venn_diagram_1 <- ggVennDiagram(junction_list)
venn_diagram_1

#The  delta psi of the junction chr9:120404113-120408347 across all the different comparisons?

delta_psi_chr_9 <- significantly_spliced_junctions %>%
  group_by(mean_dpsi_per_lsv_junction) %>%
  filter(paste_into_igv_junction == "chr9:120404113-120408347")
delta_psi_chr_9
