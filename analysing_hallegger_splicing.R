library(tidyverse)
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
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
head(significantly_spliced_junctions)
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

shared_junctions <- intersect(first_comparison, second_comparison)

junction_list <- list(junctions_normal_tdpkd = first_comparison, junctions_mutant_tdp= second_comparison)
venn_diagram_1 <- ggVennDiagram(junction_list)
venn_diagram_1

first_comparison_df <- as.data.frame(first_comparison)
second_comparison_df <- as.data.frame(second_comparison)


delta_psi_shared_events <- average_table %>%


#The  delta psi of the junction chr9:120404113-120408347 across all the different comparisons?

average_table <- significantly_spliced_junctions |>
  group_by(paste_into_igv_junction,comparison) |>
  summarize(mean_contrast_PSI = mean(contrast_PSI),
            mean_baseline_PSI = mean(baseline_PSI),mean_delta_psi = mean(mean_dpsi_per_lsv_junction),gene_name,gene_id,junc_cat) |> unique()

delta_psi_chr_9 <- average_table %>%
  group_by(mean_delta_psi) %>%
  filter(paste_into_igv_junction == "chr9:120404113-120408347")
delta_psi_chr_9

delta_psi_chr_9 %>%
  ggplot(aes(x = comparison, y = mean_delta_psi)) + 
  geom_col() + 
  coord_flip()


install.packages(data.table)
library(data.table)
library(data.table)


library(tidyr)
dpsi_shared_events <- subset(significantly_spliced_junctions, paste_into_igv_junction %in% c("chr7:92242137-92244902",  "chr1:33010833-33013207" ,  "chr6:15248964-15374117"  ,  "chr10:3099352-3099557" ,  
                "chr10:3099819-3101365",  "chrX:53625258-53625786" ,  "chr22:42599793-42602770"  , "chr20:33680416-33686004"  ,
                 "chr7:32576072-32580802" ,   "chr10:122427031-122428276" ,"chr17:50750694-50752200" ,  "chr12:27676599-27677064"  ,
               "chr12:27677096-27679489"  , "chr5:647916-648100"   ,     "chr5:75399387-75400205"  ,  "chr1:51756358-51760690"  , 
                "chr16:50224994-50225239"  , "chr7:44765579-44766164"  ,  "chr19:14446561-14450088"  , "chr5:7878322-7885701"  ,                   "chr7:22946163-22960256"  ,  "chr20:43458509-43459771" ,  "chr3:15096710-15099035"  ,  "chr7:87194102-87195031"  , 
                "chr13:106557131-106559432" , "chr15:72266822-72267475" ,  "chr9:120404113-120408347",  "chr15:22838921-22845146" , 
                 "chr17:18858261-18865802" ,  "chr18:12429437-12431149"  , "chr12:57757084-57758108"  , "chr13:113461136-113474124",
              "chr10:14878919-14896846"  , "chr1:222649791-222650292" , "chr5:115269000-115271213" , "chr5:134738792-134750759" ,
              "chr13:52661570-52662649" ,  "chr3:9464672-9468519",      "chr3:121790215-121807344",  "chr12:94402337-94402778"  ,
              "chr12:94403262-94410621"   ,"chr1:28506101-28508127"  ,  "chr12:131992230-132005077" ,"chr16:16333884-16340306" ,                  "chr12:102182627-102195951" ,"chr19:36572392-36573043" ,  "chr1:1044440-1044887"   ,  "chr1:1044987-1045161"   ,  
               "chr1:1045081-1045160"    ,  "chr20:38421098-38426419"  , "chr14:100293377-100298842", "chr19:37465473-37466739"  ,
                "chr1:244842258-244843041" , "chr20:64257051-64259942"  , "chr6:30884710-30885625" ,   "chr6:30885680-30888688"   ,
           "chr1:16576330-16582608"  ,  "chr6:30199153-30204656"  ,  "chr12:57615017-57615512"  , "chr19:36773526-36775386"  ,
               "chr19:36773538-36775386" ,  "chr7:45735670-45736787"  ,  "chr7:45728313-45735607"   , "chr15:92883186-92885515" , 
              "chr6:57949068-57959028"))


pivoted_junctions <- significantly_spliced_junctions |>
  filter(paste_into_igv_junction %in% shared_junctions) |> 
  select(comparison, mean_dpsi_per_lsv_junction, paste_into_igv_junction) %>%
pivot_wider(names_from = "comparison",
            values_from = "comparison") 
pivoted_junctions

comparison_1 <- pivoted_junctions %>%
  filter(`notdpkdrescuenotinduceda326p-sitdp43rescuenotinduceda326p` == "notdpkdrescuenotinduceda326p-sitdp43rescuenotinduceda326p", paste_into_igv_junction %in% shared_junctions ) %>%
pull(mean_dpsi_per_lsv_junction) 
comparison_1 

comparison_2 <- pivoted_junctions %>%
  filter(`sitdp43rescueinduceda326p-sitdp43rescuenotinduceda326p` == "sitdp43rescueinduceda326p-sitdp43rescuenotinduceda326p") %>%
  pull(mean_dpsi_per_lsv_junction) 
comparison_2

# The two vectors
comparison_one <- c(0.4194, 0.2151, 0.3242, 0.2745, 0.3094, 0.3763, 0.5860, 0.4212, 0.7300, 0.8042, 0.4996, 0.2317, 0.3432,
                    0.4788, 0.2804, 0.3508, 0.7509, 0.7269, 0.7830, 0.3876, 0.2538, 0.3569, 0.5105, 0.2464, 0.3851, 0.4080,
                    0.3426, 0.2993, 0.2126, 0.2767, 0.4772, 0.2443, 0.4118, 0.2935, 0.2199, 0.2337, 0.5652, 0.5127, 0.3200,
                    0.3489, 0.3429, 0.3066, 0.3726, 0.2974, 0.2784, 0.2939, 0.3021, 0.3205, 0.2683, 0.3233, 0.3107, 0.3969,
                    0.5774, 0.5214, 0.3166, 0.4737, 0.2590, 0.3067, 0.3691, 0.3110, 0.4260, 0.4420, 0.2596, 0.2523, 0.2189,
                    0.3140, 0.2444, 0.2907, 0.3091, 0.7911, 0.3166, 0.4068, 0.3525, 0.4557, 0.3595, 0.2848, 0.2889, 0.3263,
                    0.3336, 0.2032, 0.3132, 0.2619, 0.4136, 0.3357, 0.3250, 0.3469, 0.3617, 0.3912, 0.3990)

comparison_two <- c(0.2804, 0.2827, 0.2401, 0.3586, 0.5795, 0.2835, 0.3339, 0.6682, 0.7005, 0.2289, 0.2827, 0.2022, 0.3537,
                    0.2247, 0.5068, 0.4493, 0.7055, 0.2186, 0.2354, 0.2095, 0.3834, 0.4104, 0.2035, 0.2311, 0.2591, 0.3857,
                    0.3097, 0.3482, 0.4696, 0.4344, 0.2171, 0.2788, 0.2877, 0.3559, 0.2374, 0.3397, 0.2826, 0.4128, 0.3057,
                    0.2627, 0.2861, 0.2683, 0.2404, 0.3018, 0.3299, 0.3918, 0.6330, 0.2994, 0.2750, 0.3913, 0.4090, 0.2422,
                    0.2518, 0.2003, 0.2266, 0.3991, 0.3182, 0.4764, 0.2354, 0.4502, 0.2392, 0.2958, 0.2615, 0.2100, 0.3674,
                    0.2095, 0.2982, 0.2475, 0.4556, 0.2925, 0.2428, 0.2429, 0.3939, 0.3370, 0.2150, 0.2492, 0.3512, 0.3664)

# Ensuring both vectors have the same length (take the minimum length)
min_length <- min(length(comparison_one), length(comparison_two))

# Creating the scatterplot
plot(comparison_one[1:min_length], comparison_two[1:min_length], 
     xlab = "Δ PSI of shared junctions in TDP-43 knockdown", 
     ylab = "Δ PSI of shared junctions in  mutant TDP-43",
     pch = 19) 
