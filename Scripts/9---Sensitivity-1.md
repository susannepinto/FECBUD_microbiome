Sensitivity 1
================
S. Pinto
December 2022

## Load required packages

``` r
library( dplyr )
library( magrittr )
library( knitr )
library( tidyverse )
library( reshape2 )
```

``` r
library( phyloseq )
library( microbiome )
library( devtools )
library( stringr )
```

``` r
library( aod )
library( nlme )
library( lme4 )
library( splines )
library( jtools )
```

## Load the data

Here we use the relative abundance dataset.

## Select ‘average’ donor mOTUs

``` r
donor_A_phyloseq <- subset_samples( physeq_mOTU.new, subject_id == "Donor A" )
donor_A_phyloseq
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1552 taxa and 13 samples ]
    ## sample_data() Sample Data:       [ 13 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1552 taxa by 7 taxonomic ranks ]

``` r
## Now determine the actual 'core' or persistent species per donor
donor_A_core_species <- core_members(
  x = donor_A_phyloseq,
  detection = 0.001, # detection limit is the same as the general limit (0.1%)
  prevalence = 6 / ( length( sample_names( donor_A_phyloseq ))),
  include.lowest = TRUE
  # select only bacteria that occur at at least 2 different timepoints
)

donor_B_phyloseq <- subset_samples( physeq_mOTU.new, subject_id == "Donor B" )
donor_B_phyloseq
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1552 taxa and 14 samples ]
    ## sample_data() Sample Data:       [ 14 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1552 taxa by 7 taxonomic ranks ]

``` r
donor_B_core_species <- core_members(
  x =  donor_B_phyloseq,
  detection = 0.001, # detection limit is the same as the general limit (0.1%)
  prevalence = 6 / ( length( sample_names( donor_B_phyloseq ))),
  include.lowest = TRUE
  # select only bacteria that occur at at least 2 different timepoints
)
```

## Select patient mOTUs per timepoint

``` r
# Collect data for patients per timepoint
# and convert to dataframe for easier processing
baseline_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Baseline" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
baseline_bacteria.select <- select( baseline_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
baseline_per_recipient <- unstack(baseline_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Pre_FMT_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Pre-FMT" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Pre_FMT_bacteria.select <- select( Pre_FMT_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Pre_FMT_per_recipient <- unstack( Pre_FMT_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Post_1_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Post-1" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Post_1_bacteria.select <- select( Post_1_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Post_1_per_recipient <- unstack( Post_1_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Post_2_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Post-2" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Post_2_bacteria.select <- select( Post_2_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Post_2_per_recipient <- unstack( Post_2_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Post_3_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Post-3" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Post_3_bacteria.select <- select( Post_3_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Post_3_per_recipient <- unstack( Post_3_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Post_4_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Post-4" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Post_4_bacteria.select <- select( Post_4_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Post_4_per_recipient <- unstack( Post_4_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Week8_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Week8" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Week8_bacteria.select <- select( Week8_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Week8_per_recipient <- unstack( Week8_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Week10_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Week10" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Week10_bacteria.select <- select( Week10_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Week10_per_recipient <- unstack( Week10_bacteria.select[ ,1:2 ])

# Collect data for patients per timepoint
# and convert to dataframe for easier processing
Week14_bacteria <- psmelt( physeq_mOTU.new ) %>%
  filter( timepoint.new == "Week14" ) %>%
  filter( Abundance > 0 ) %>%
  arrange( subject_id )

# Now list all species present in each patient before FMT
Week14_bacteria.select <- select( Week14_bacteria, OTU, subject_id, treated_with_donor )
# Make a list per recipient
Week14_per_recipient <- unstack( Week14_bacteria.select[ ,1:2 ])
```

## Obtain results per patient per sample

``` r
# Remove donors
physeq_mOTU.new = subset_samples( physeq_mOTU.new, subject_id != "Donor A" )
physeq_mOTU.new = subset_samples( physeq_mOTU.new, subject_id != "Donor B" )
```

``` r
sub.id <- physeq_mOTU.new@sam_data[[ "subject_id" ]] %>% unique()
mOTUs.all <- physeq_mOTU.new@otu_table %>% rownames()
times <- physeq_mOTU.new@sam_data[[ "timepoint.new" ]] %>% unique()
```

``` r
data.recipient <- matrix( nrow = length( mOTUs.all ), ncol = 11 )  
colnames( data.recipient ) <- c( "mOTUs", "Donor", "Baseline" ,"Pre_FMT", "Post_1", "Post_2", "Post_3",
                          "Post_4", "Week8", "Week10", "Week14" )
data.recipient[ ,1 ] <- mOTUs.all 
head( data.recipient )
```

    ##      mOTUs               Donor Baseline Pre_FMT Post_1 Post_2 Post_3 Post_4
    ## [1,] "ref_mOTU_v3_00077" NA    NA       NA      NA     NA     NA     NA    
    ## [2,] "ref_mOTU_v3_00084" NA    NA       NA      NA     NA     NA     NA    
    ## [3,] "ref_mOTU_v3_00085" NA    NA       NA      NA     NA     NA     NA    
    ## [4,] "ref_mOTU_v3_00086" NA    NA       NA      NA     NA     NA     NA    
    ## [5,] "ref_mOTU_v3_00087" NA    NA       NA      NA     NA     NA     NA    
    ## [6,] "ref_mOTU_v3_00095" NA    NA       NA      NA     NA     NA     NA    
    ##      Week8 Week10 Week14
    ## [1,] NA    NA     NA    
    ## [2,] NA    NA     NA    
    ## [3,] NA    NA     NA    
    ## [4,] NA    NA     NA    
    ## [5,] NA    NA     NA    
    ## [6,] NA    NA     NA

``` r
custom_min <- function(x) {if (length(x)>0) min(x) else 0}
custom_max <- function(x) {if (length(x)>0) max(x) else 0}

# The list to store the results
listofresults <- list()

for( i in 1:length( physeq_mOTU.new@sam_data[[ "subject_id" ]] %>% unique())){
  patient = sub.id[[ i ]]
  patient.characteristics <- physeq_mOTU.new@sam_data[ paste( patient, "-FMT1", sep = "" )] %>% as.matrix()
  
  # select donor bacteria
  if( patient.characteristics[ ,"treated_with_donor" ] == "Donor A" ){
  Donor = donor_A_core_species
  } else if( patient.characteristics[ ,"treated_with_donor" ] == "Donor B" ) {
  Donor = donor_B_core_species }
  
  # select patient characteristics of interest
  patient.characteristics <- paste( patient.characteristics[ ,"subject_id" ] %>% as.character(), "-" , patient.characteristics[ ,"clinical_outcome_wk14" ] %>%
                                     as.character(), "-", patient.characteristics[ ,"treated_with_donor" ] %>% as.character(), "-",
                                   patient.characteristics[ ,"pretreatment" ] %>% as.character(), "-", patient.characteristics[ ,"age" ] %>% as.character(), "-",
                                   patient.characteristics[ ,"sex" ] %>% as.character())

  # mOTUs present/absent data
  for( j in 1:nrow( data.recipient )) {
    data.recipient[ j, 2 ] <- ifelse( data.recipient[ j, 1 ] %in% Donor, 1, 0)
    
    data.recipient[ j, 3 ] <- ifelse( data.recipient[ j, 1 ] %in% baseline_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 4 ] <- ifelse( data.recipient[ j, 1 ] %in% Pre_FMT_per_recipient[[ patient ]], 1, 0 )
  
    data.recipient[ j, 5 ] <- ifelse( data.recipient[ j, 1 ] %in% Post_1_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 6 ] <- ifelse( data.recipient[ j, 1 ] %in% Post_2_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 7 ] <- ifelse( data.recipient[ j, 1 ] %in% Post_3_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 8 ] <- ifelse( data.recipient[ j, 1 ] %in% Post_4_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 9 ] <- ifelse( data.recipient[ j, 1 ] %in% Week8_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 10 ] <- ifelse( data.recipient[ j, 1 ] %in% Week10_per_recipient[[ patient ]], 1, 0 )
    data.recipient[ j, 11 ] <- ifelse( data.recipient[ j, 1 ] %in% Week14_per_recipient[[ patient ]], 1, 0 )
    }
  
  data.recipient.sel <- data.recipient[apply(data.recipient[ ,-1 ], 1, function(x) !all( x == 0 )),] %>% as.data.frame()
  data.recipient.sel %>% head()
  data.recipient.sel[ ,2:11 ] <- sapply(data.recipient.sel[ ,2:11 ],as.numeric )
  data.recipient.sel[ ,2:11 ] %>% colSums()
  # some patients have missing data (colsums = NA)
  data.recipient.sel[ ,2:11 ][ data.recipient.sel[ , 2:11 ] %>% colSums() == 0 ] <- NA
  data.recipient.sel %>% head()
  
  # Possible combinations
  #expand.grid( rep( list( 0:1 ), 3 ))
  #expand.grid( rep( list( 0:1 ), 4 ))
  #expand.grid( rep( list( 0:1 ), 7 ))
  #expand.grid( rep( list( 0:1 ), 8 ))
  
  data.recipient.results <- matrix( ncol = ncol( data.recipient ), nrow = nrow( data.recipient.sel ))
  colnames( data.recipient.results ) <- colnames( data.recipient )
  data.recipient.results <- data.recipient.results %>% as.data.frame()
  data.recipient.results[ ,"mOTUs" ] <- data.recipient.sel[ ,"mOTUs" ]
  data.recipient.results %>% head()
  
  for( j in 1:nrow( data.recipient.results )) {
    for (k in 5:ncol( data.recipient.results )) {
    data.recipient.results[ j, k ] <- 
      # If the pre FMT samples or the current sample is missing, name 'NA'
      if( ( is.na( data.recipient.sel[ j, 3 ]) & is.na( data.recipient.sel[ j, 4 ]) ) | is.na( data.recipient.sel[ j, k ]) ) {
        NA
        # If present at Pre-FMT, present at the current sample and been present till now (1x 0 possible), name Resident
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) != 0 &
                 data.recipient.sel[ j, k ] == 1 & 
                 sum( data.recipient.sel[ j, 5:k ], na.rm = TRUE ) >= length(data.recipient.sel[ j, 5:k ])  ) {
        "Resident"
        # If present at Pre-FMT, present at the current sample and have not always been present till now, name host transient
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) != 0 &
                 data.recipient.sel[ j, k ] == 1 & 
                 sum( data.recipient.sel[ j, 5:k ], na.rm = TRUE ) < length(data.recipient.sel[ j, 5:k ])  ) {
        "Host Transient"
        # If present at donor (not pre-fmt), present at the current sample and have been present since first 1 (1x 0 possible), name colonisation
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) == 0 &
                 data.recipient.sel[ j, 2 ] == 1 & 
                 data.recipient.sel[ j, k ] == 1 & 
                 (sum( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))], na.rm = TRUE ) >= length( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))])  )) {
        "Colonisation" 
        # If present at donor (not pre-fmt), present at the current sample and have not been present since first 1 (1x 0 possible), name donor transient
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) == 0 &
                 data.recipient.sel[ j, 2 ] == 1 & 
                 data.recipient.sel[ j, k ] == 1 & 
                 (sum( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))], na.rm = TRUE ) < length( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))])) ) {
       "Donor Transient"
      # If not present at pre-fmt or donor, but current present and have been present since first 1 (1x 0 possible), name Novel
        } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) == 0 &
                 data.recipient.sel[ j, 2 ] == 0 & 
                 data.recipient.sel[ j, k ] == 1 & 
                 (sum( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))], na.rm = TRUE ) >= length( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))]) )) {
       "Novel"
          # If not present at pre-fmt or donor,  and have not been present since first 1 (1x 0 possible), name Novel transient
        } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) == 0 &
                 data.recipient.sel[ j, 2 ] == 0 & 
                 data.recipient.sel[ j, k ] == 1 &
                 (sum( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))], na.rm = TRUE ) < length( data.recipient.sel[ j, 5:k ][ custom_min( which( data.recipient.sel[ j, 5:k ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:k ] != 0 ))]) ) ) {
       "Novel Transient"
# ************************************************************          
        # # the first 0 after the start of 1 is NA
        #   } else if( 
        #          data.recipient.sel[ j, k ] == 0 & 
        #          
        #          ( data.recipient.sel[ j, ( k - 1 ) ] == 1 &&
        #          !is.na( data.recipient.sel[ j, ( k - 1 ) ] )) &
        #          
        #          ( (data.recipient.sel[ j, ( k + 1 ) ] == 1 &&
        #            !is.na( data.recipient.sel[ j, ( k + 1 ) ] )) ||
        #             k == 11 ) &
        #          
        #          (sum( data.recipient.sel[ j, 5:(k-1) ][ custom_min( which( data.recipient.sel[ j, 5:(k-1) ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:(k-1) ] != 0 ))], na.rm = TRUE ) >= length( data.recipient.sel[ j, 5:(k-1) ][ custom_min( which( data.recipient.sel[ j, 5:(k-1) ] != 0 )) : custom_max( which( data.recipient.sel[ j, 5:(k-1) ] != 0 ))]) )) {
        # NA
        # 
            # till presence and absent in donor an pre FMT, name absent
          } else if( 
            sum( data.recipient.sel[ j, 2:k ], na.rm = TRUE ) == 0 &
                       data.recipient.sel[ j, k ] == 0 ) {
        "Absent" 
            # if present in pre FMT samples, but not present in the current sample (except the first one), name Species loss
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) != 0 &
                 data.recipient.sel[ j, k ] == 0 ) {
        "Species loss" 
    # if present in donor sample, but not present in the current sample (except the first one), name Rejection
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) == 0 &
                 data.recipient.sel[ j, 2 ] == 1 & 
                 data.recipient.sel[ j, k ] == 0 ) {
        "Rejection" 
      # if not present in pre-FMT or donor samples, not present currently, but have been present before as novel, novel loss
      } else if( sum( data.recipient.sel[ j, 3:4 ], na.rm = TRUE ) == 0 &
                 data.recipient.sel[ j, 2 ] == 0 & 
                 data.recipient.sel[ j, k ] == 0 ) {
       "Novel loss"
        # if not been classified according to the rules above
      } else {
        "Unclassified" 
      }
    }
  }
  data.recipient.results.meta <- cbind( data.recipient.results, patient.characteristics )
  listofresults[[ i ]]  <- data.recipient.results.meta

  # print( i )
  }
  
head( listofresults[[ i ]])
```

    ##               mOTUs Donor Baseline Pre_FMT       Post_1       Post_2
    ## 1 ref_mOTU_v3_00077    NA       NA      NA Species loss Species loss
    ## 2 ref_mOTU_v3_00086    NA       NA      NA Species loss Species loss
    ## 3 ref_mOTU_v3_00095    NA       NA      NA     Resident     Resident
    ## 4 ref_mOTU_v3_00312    NA       NA      NA Species loss Species loss
    ## 5 ref_mOTU_v3_00321    NA       NA      NA Species loss Species loss
    ## 6 ref_mOTU_v3_00322    NA       NA      NA Species loss Species loss
    ##         Post_3       Post_4        Week8       Week10         Week14
    ## 1 Species loss Species loss Species loss Species loss Host Transient
    ## 2 Species loss Species loss Species loss Species loss   Species loss
    ## 3     Resident     Resident     Resident     Resident   Species loss
    ## 4 Species loss Species loss Species loss Species loss   Species loss
    ## 5 Species loss Species loss Species loss Species loss Host Transient
    ## 6 Species loss Species loss Species loss Species loss Host Transient
    ##                      patient.characteristics
    ## 1 125 - Good - Donor B - budesonide - 32 - M
    ## 2 125 - Good - Donor B - budesonide - 32 - M
    ## 3 125 - Good - Donor B - budesonide - 32 - M
    ## 4 125 - Good - Donor B - budesonide - 32 - M
    ## 5 125 - Good - Donor B - budesonide - 32 - M
    ## 6 125 - Good - Donor B - budesonide - 32 - M

``` r
# count data
resultscounts <- list()
for(i in 1:length( listofresults )){
  resultscounts[[ i ]] <- listofresults[[ i ]] %>% 
    pivot_longer( cols = c( Donor, Baseline, Pre_FMT, Post_1, Post_2, Post_3,
                          Post_4, Week8, Week10, Week14 ), names_to = "timepoint", values_to = "category" ) %>% 
    group_by( timepoint ) %>% 
    mutate( category = factor( category, levels = c( "Resident", "Host Transient", "Species loss", 
                                                     "Colonisation", "Donor Transient", "Rejection",
                                                     "Novel", "Novel Transient", "Novel loss",
                                                     "Absent" ))) %>% 
    filter( category != "Absent" ) %>%
    summarise( n = c(( table( category )))) %>% 
    mutate( category = c( "Resident", "Host Transient", "Species loss", 
                                                     "Colonisation", "Donor Transient", "Rejection",
                                                     "Novel", "Novel Transient", "Novel loss",
                                                     "Absent" )) 
  
  resultscounts[[ i ]]$timepoint = factor( resultscounts[[ i ]]$timepoint, levels = c(  "Donor", "Baseline", "Pre_FMT","Post_1", "Post_2", "Post_3",
                                                        "Post_4", "Week8", "Week10", "Week14" ))
  resultscounts[[ i ]]$patient.characteristics <- listofresults[[ i ]]$patient.characteristics %>% unique()
  resultscounts[[ i ]][c( 'Patient_ID', 'State', 'Donor', 'Pretreatment', 'Age', 'Sex' )] <- str_split_fixed(resultscounts[[ i ]]$patient.characteristics, ' - ', 6)
  resultscounts[[ i ]] <- resultscounts[[ i ]][ c( 'Patient_ID', 'State', 'Donor', 'Pretreatment', 'Age', 'Sex' ,'timepoint', 'category', 'n' )]
}
```

``` r
results.all.count <- bind_rows( resultscounts, .id = "column_label" )
results.all.count <- results.all.count[ c( 2:10 )]
```

## Save results

``` r
save( results.all.count, file = "Scenario2_sens.Rda" )
```

## Investigate results

``` r
dim( results.all.count )
```

    ## [1] 1330    9

``` r
head( results.all.count )
```

    ## # A tibble: 6 × 9
    ## # Groups:   timepoint [1]
    ##   Patient_ID State Donor   Pretreatment Age   Sex   timepoint category         n
    ##   <chr>      <chr> <chr>   <chr>        <chr> <chr> <fct>     <chr>        <int>
    ## 1 101        Good  Donor A placebo      58    F     Post_1    Resident       199
    ## 2 101        Good  Donor A placebo      58    F     Post_1    Host Transi…     0
    ## 3 101        Good  Donor A placebo      58    F     Post_1    Species loss    72
    ## 4 101        Good  Donor A placebo      58    F     Post_1    Colonisation    12
    ## 5 101        Good  Donor A placebo      58    F     Post_1    Donor Trans…     0
    ## 6 101        Good  Donor A placebo      58    F     Post_1    Rejection       27

``` r
typeof( results.all.count )
```

    ## [1] "list"

``` r
hist( results.all.count$n, n = 50 )
```

![](9---Sensitivity-1_files/figure-gfm/Investigate%20results-1.png)<!-- -->

``` r
results.all.count$timepoint %>% unique()
```

    ## [1] Post_1 Post_2 Post_3 Post_4 Week10 Week14 Week8 
    ## 10 Levels: Donor Baseline Pre_FMT Post_1 Post_2 Post_3 Post_4 Week8 ... Week14

``` r
results.all.count$category %>% unique()
```

    ##  [1] "Resident"        "Host Transient"  "Species loss"    "Colonisation"   
    ##  [5] "Donor Transient" "Rejection"       "Novel"           "Novel Transient"
    ##  [9] "Novel loss"      "Absent"
