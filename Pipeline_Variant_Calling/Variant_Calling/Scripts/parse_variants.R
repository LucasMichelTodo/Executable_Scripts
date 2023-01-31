library(tidyverse)

wd <- '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
setwd(wd)

variants <- read_tsv('parsed_variants_new.tsv') %>%
  mutate(Var_id = paste0('Variant_', `...1`)) %>%
  select(Var_id, everything(), `...1`)

parse_variants_bystrain <- function(strain, depth_filter, refratio_filter, impact_filter){
  depthcol <- paste0('depth_', strain)
  ratiocol <- paste0('RefRatio_', strain)
  outname <- paste0('./Parsed_by_Strain/',
                    strain,
                    '_variants_depth_', depth_filter,
                    '_refratio_', refratio_filter,
                    '_impactfilter_', impact_filter,
                    '.tsv'
                    )

  if (impact_filter){
    variants <- variants %>%
      filter(IMPACT == 'HIGH')
  }

  variants %>%
    filter(get(depthcol) >= depth_filter &
           get(ratiocol) <= refratio_filter) %>%
    select(Var_id, contains(strain), Gene, Annot, Consequence,
           Chrom, Pos, Ref, Alt, everything()) %>%
    mutate(Annot = gsub('\"', '', Annot)) %>%
    write_tsv(outname)
}

#### Parse variants per strain

depth_filter <- 20
refratio_filter <- 0.5
impact_filter <- F

strains <- c('12B', '10G', 'A7', 'E5', 'B11')
for (strain in strains){
  parse_variants_bystrain(
    strain, depth_filter, refratio_filter, impact_filter
  )
}

#### Parse variants, all strains together

## depth_filter <- 20
## refratio_filter <- 0.5
## impact_filter <- F

## outname <- paste0(
##   'allstrains_variants_depth_', depth_filter,
##   '_refratio_', refratio_filter,
##   '_impactfilter_', impact_filter,
##   '.tsv'
## )

## filtered_vars <- variants
## if (impact_filter){
##   filtered_vars <- variants %>%
##     filter(IMPACT == 'HIGH')
## }

## filtered_vars <- variants %>%
##   filter(
##   (RefRatio_12B <= refratio_filter & depth_12B >= depth_filter) |
##   (RefRatio_10G <= refratio_filter & depth_10G >= depth_filter) |
##   (RefRatio_A7 <= refratio_filter & depth_A7 >= depth_filter) |
##   (RefRatio_E5 <= refratio_filter & depth_E5 >= depth_filter) |
##   (RefRatio_B11 <= refratio_filter & depth_B11 >= depth_filter)
##   ) %>%
##   select(Chrom, Pos, Ref, Alt, everything()) %>%
##   mutate(Annot = gsub('\"', '', Annot)) %>%
##   write_tsv(outname)
