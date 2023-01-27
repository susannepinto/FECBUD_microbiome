Average Microbiota Compostition
================
December 2022

## Load required packages

``` r
library( dplyr )
library( magrittr )
library( knitr )
library( tidyverse )
library( lattice )
library( latticeExtra)
```

``` r
library( phyloseq )
library( microbiome )
library( microViz )
library( ggplot2 )
library( microbiomeutilities )
library( RColorBrewer )
library( patchwork )
library( ggpubr )
library( gridExtra )
library( grid )
```

## Load the data

Here we use the relative abundance dataset.

## Preprocessing of the data

``` r
# Make compositional
physeq_mOTU.rel = microbiome::transform(physeq_mOTU, "compositional")

# Aggregate to family level
physeq_mOTU.rel = aggregate_taxa( physeq_mOTU.rel, "family" )
physeq_mOTU.rel
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 95 taxa and 207 samples ]
    ## sample_data() Sample Data:       [ 207 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 95 taxa by 6 taxonomic ranks ]

``` r
get_taxa_unique(physeq_mOTU.rel, "family")
```

    ##  [1] "Actinobacteriafam.[Sanguibacteraceae/Actinomycetaceae/Promicromonosporaceae]"
    ##  [2] "Actinobacteriafam.[Actinomycetaceae/Streptomycetaceae]"                      
    ##  [3] "Actinomycetaceae"                                                            
    ##  [4] "Bifidobacteriaceae"                                                          
    ##  [5] "Corynebacteriaceae"                                                          
    ##  [6] "Micrococcaceae"                                                              
    ##  [7] "Atopobiaceae"                                                                
    ##  [8] "Coriobacteriaceae"                                                           
    ##  [9] "Coriobacterialesfam.incertaesedis"                                           
    ## [10] "Coriobacteriiafam.incertaesedis"                                             
    ## [11] "Eggerthellaceae"                                                             
    ## [12] "Eggerthellalesfam.incertaesedis"                                             
    ## [13] "Bacteriafam.[Lactobacillaceae/Rhodospirillaceae]"                            
    ## [14] "Bacteriafam.incertaesedis"                                                   
    ## [15] "Bacteroidaceae"                                                              
    ## [16] "Bacteroidalesfam.[Porphyromonadaceae/Barnesiellaceae]"                       
    ## [17] "Bacteroidalesfam.incertaesedis"                                              
    ## [18] "Barnesiellaceae"                                                             
    ## [19] "Muribaculaceae"                                                              
    ## [20] "Odoribacteraceae"                                                            
    ## [21] "Porphyromonadaceae"                                                          
    ## [22] "Prevotellaceae"                                                              
    ## [23] "Rikenellaceae"                                                               
    ## [24] "Tannerellaceae"                                                              
    ## [25] "Flavobacteriaceae"                                                           
    ## [26] "Flavobacteriiafam.incertaesedis"                                             
    ## [27] "Dehalococcoidalesfam.incertaesedis"                                          
    ## [28] "Elusimicrobiafam.incertaesedis"                                              
    ## [29] "Bacillalesfam.incertaesedis"                                                 
    ## [30] "Staphylococcaceae"                                                           
    ## [31] "Carnobacteriaceae"                                                           
    ## [32] "Enterococcaceae"                                                             
    ## [33] "Lactobacillaceae"                                                            
    ## [34] "Leuconostocaceae"                                                            
    ## [35] "Streptococcaceae"                                                            
    ## [36] "Catabacteriaceae"                                                            
    ## [37] "Christensenellaceae"                                                         
    ## [38] "Clostridiaceae"                                                              
    ## [39] "Clostridialesfam.[Clostridiaceae/Christensenellaceae]"                       
    ## [40] "Clostridialesfam.[Clostridiaceae/Heliobacteriaceae]"                         
    ## [41] "Clostridialesfam.[Clostridiaceae/Ruminococcaceae]"                           
    ## [42] "Clostridialesfam.[Eubacteriaceae/Clostridiaceae]"                            
    ## [43] "Clostridialesfam.[Eubacteriaceae/Ruminococcaceae]"                           
    ## [44] "Clostridialesfam.[Lachnospiraceae/Clostridiaceae]"                           
    ## [45] "Clostridialesfam.[Lachnospiraceae/Clostridiaceae/Ruminococcaceae]"           
    ## [46] "Clostridialesfam.[Lachnospiraceae/Ruminococcaceae]"                          
    ## [47] "Clostridialesfam.incertaesedis"                                              
    ## [48] "ClostridialesFamilyXIII.IncertaeSedis"                                       
    ## [49] "Eubacteriaceae"                                                              
    ## [50] "Lachnospiraceae"                                                             
    ## [51] "Oscillospiraceae"                                                            
    ## [52] "Peptostreptococcaceae"                                                       
    ## [53] "Ruminococcaceae"                                                             
    ## [54] "Clostridiafam.incertaesedis"                                                 
    ## [55] "Erysipelotrichaceae"                                                         
    ## [56] "Firmicutesfam.[Lachnospiraceae/Lactobacillaceae]"                            
    ## [57] "Firmicutesfam.[Erysipelotrichaceae/Clostridiaceae]"                          
    ## [58] "Firmicutesfam.[Erysipelotrichaceae/Eubacteriaceae/Clostridiaceae]"           
    ## [59] "Firmicutesfam.[Erysipelotrichaceae/Lachnospiraceae]"                         
    ## [60] "Firmicutesfam.[Erysipelotrichaceae/Peptostreptococcaceae]"                   
    ## [61] "Firmicutesfam.[Veillonellaceae/Ruminococcaceae]"                             
    ## [62] "Firmicutesfam.[Carnobacteriaceae/Acidaminococcaceae]"                        
    ## [63] "Firmicutesfam.incertaesedis"                                                 
    ## [64] "Acidaminococcaceae"                                                          
    ## [65] "Selenomonadaceae"                                                            
    ## [66] "Veillonellaceae"                                                             
    ## [67] "Peptoniphilaceae"                                                            
    ## [68] "Fusobacteriaceae"                                                            
    ## [69] "Leptotrichiaceae"                                                            
    ## [70] "Lentisphaeraefam.incertaesedis"                                              
    ## [71] "Alphaproteobacteriafam.incertaesedis"                                        
    ## [72] "Acetobacteraceae"                                                            
    ## [73] "Rhodospirillaceae"                                                           
    ## [74] "Burkholderialesfam.incertaesedis"                                            
    ## [75] "Comamonadaceae"                                                              
    ## [76] "Oxalobacteraceae"                                                            
    ## [77] "Sutterellaceae"                                                              
    ## [78] "Neisseriaceae"                                                               
    ## [79] "Desulfovibrionaceae"                                                         
    ## [80] "Myxococcaceae"                                                               
    ## [81] "Campylobacteraceae"                                                          
    ## [82] "Enterobacteriaceae"                                                          
    ## [83] "Hafniaceae"                                                                  
    ## [84] "Morganellaceae"                                                              
    ## [85] "Pasteurellaceae"                                                             
    ## [86] "Moraxellaceae"                                                               
    ## [87] "Pseudomonadaceae"                                                            
    ## [88] "Proteobacteriafam.incertaesedis"                                             
    ## [89] "Brachyspiraceae"                                                             
    ## [90] "Synergistaceae"                                                              
    ## [91] "Mollicutesfam.incertaesedis"                                                 
    ## [92] "Mycoplasmataceae"                                                            
    ## [93] "Verrucomicrobiafam.incertaesedis"                                            
    ## [94] "Akkermansiaceae"

## PCOA plot of microbiome composition

``` r
# Change the data for the PCOA plot
physeq_PCOA <- physeq_mOTU.rel

cl_outcome <- sample_data( physeq_PCOA ) %>%
    data.frame() %>%
    select( "subject_id", "clinical_outcome_wk14" ) %>%
    mutate_if( is.factor, as.character )

cl_outcome <- cl_outcome %>% 
    mutate( clinical_outcome_wk14 = coalesce( clinical_outcome_wk14, subject_id ))

physeq_PCOA@sam_data[["clinical_outcome_wk14"]] <- cl_outcome$clinical_outcome_wk14
```

``` r
physeq_PCOA %>% tax_transform( "identity", rank = "unique" ) %>% # best graph
  ord_calc( method = "PCA" ) %>%
  dist_calc( "aitchison" ) %>%
  ord_plot(
    axes = c( 1, 2 ), 
    colour = "clinical_outcome_wk14", fill = "clinical_outcome_wk14", plot_taxa = 1:3
  ) +
  scale_color_manual( values = c("#85C1E9", "#BD86EB", "#82E0AA", "#EC7063" )) +
  ggplot2::stat_ellipse(
    ggplot2::aes( colour = clinical_outcome_wk14 )
  )
```

![](Average_Microbiota_Composition_files/figure-gfm/PCOA%20plot-1.png)<!-- -->

## Subset data and select core patient micorbiota

``` r
temp <- aggregate_top_taxa2( physeq_mOTU.rel, 15, "family" )
temp@otu_table %>% rownames()
```

    ##  [1] "Bacteroidaceae"                 "Bacteroidalesfam.incertaesedis"
    ##  [3] "Bifidobacteriaceae"             "Clostridiaceae"                
    ##  [5] "Clostridialesfam.incertaesedis" "Coriobacteriaceae"             
    ##  [7] "Eubacteriaceae"                 "Firmicutesfam.incertaesedis"   
    ##  [9] "Lachnospiraceae"                "Oscillospiraceae"              
    ## [11] "Other"                          "Prevotellaceae"                
    ## [13] "Rikenellaceae"                  "Ruminococcaceae"               
    ## [15] "Sutterellaceae"                 "Veillonellaceae"

``` r
temp@otu_table <- temp@otu_table[ c( 1:10, 12:16, 11 ), ]

temp@tax_table %>% rownames()
```

    ##  [1] "Bacteroidaceae"                 "Bacteroidalesfam.incertaesedis"
    ##  [3] "Bifidobacteriaceae"             "Clostridiaceae"                
    ##  [5] "Clostridialesfam.incertaesedis" "Coriobacteriaceae"             
    ##  [7] "Eubacteriaceae"                 "Firmicutesfam.incertaesedis"   
    ##  [9] "Lachnospiraceae"                "Oscillospiraceae"              
    ## [11] "Other"                          "Prevotellaceae"                
    ## [13] "Rikenellaceae"                  "Ruminococcaceae"               
    ## [15] "Sutterellaceae"                 "Veillonellaceae"

``` r
temp@tax_table <- temp@tax_table[ c( 1:10, 12:16, 11 ), ]

temp
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16 taxa and 207 samples ]
    ## sample_data() Sample Data:       [ 207 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16 taxa by 2 taxonomic ranks ]

``` r
temp.donorA = temp %>%
  ps_filter( subject_id == "Donor A", .keep_all_taxa = TRUE )

temp.donorA@sam_data[[ "timepoint.new" ]] <- temp.donorA@sam_data[[ "timepoint.new" ]] %>% as.factor()
temp.donorA@sam_data[[ "timepoint.new" ]] <- factor( temp.donorA@sam_data[[ "timepoint.new" ]], levels = c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"  ))

temp.donorB = temp %>%
  ps_filter( subject_id == "Donor B", .keep_all_taxa = TRUE )

temp.donorB@sam_data[[ "timepoint.new" ]] <- temp.donorB@sam_data[[ "timepoint.new" ]] %>% as.factor()
temp.donorB@sam_data[[ "timepoint.new" ]] <- factor( temp.donorB@sam_data[[ "timepoint.new" ]], levels = c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"  ))

temp.base = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Baseline", .keep_all_taxa = TRUE )

temp.fmt1 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Pre-FMT", .keep_all_taxa = TRUE ) 

temp.fmt2 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Post-1", .keep_all_taxa = TRUE ) 

temp.fmt3 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Post-2", .keep_all_taxa = TRUE )

temp.fmt4 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Post-3", .keep_all_taxa = TRUE )

temp.wk7 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Post-4", .keep_all_taxa = TRUE )

temp.wk8 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Week8", .keep_all_taxa = TRUE )

temp.wk10 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Week10", .keep_all_taxa = TRUE )

temp.wk14 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( timepoint.new == "Week14", .keep_all_taxa = TRUE )
```

## Average plots per timepoint

``` r
# Choose colors
mycolors <- colorRampPalette( brewer.pal( 8, "Set1" ))( 16 )
```

``` r
pdonorA = plot_composition( temp.donorA, average_by = "treated_with_donor" ) +
  theme( text = element_text(size = 20),
         legend.position = "none", axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        # axis.text = element_text( size = 8 ),
        # axis.title.x = element_text( size = 8 ),
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Average", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "Donor A" )) 

pdonorB = plot_composition( temp.donorB, average_by = "treated_with_donor" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_line( "white" ),
        axis.text.y = element_blank(),
        # axis.title.x = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Average", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "Donor B" ))

p1 = plot_composition( temp.base, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_line( "white" ),
        axis.text.y = element_blank(),
        # axis.title.x = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Baseline", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 

p2 = plot_composition( temp.fmt1, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_line( "white" ),
        axis.text.y = element_blank(),
        # axis.title.x = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Pre-FMT" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
# + theme( aspect.ratio = 2/1 )

p3 = plot_composition( temp.fmt2, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        # axis.title.x = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Post-1" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
# + theme(aspect.ratio = 2/1)

p4 = plot_composition( temp.fmt3, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        # axis.title.x = element_text( size = 8 ),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Post-2" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete(labels=c( "NR", "R" )) 
# + theme( aspect.ratio = 2/1 ) 

p5 = plot_composition( temp.fmt4, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        # axis.title.x = element_text( size = 8 ),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Post-3" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
# + theme( aspect.ratio = 2/1 )

p6 = plot_composition( temp.wk7, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        # axis.title.x = element_text( size = 8 ),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Post-4" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
# + theme( aspect.ratio = 2/1 )
 
p7 = plot_composition( temp.wk8, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x=element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        # axis.title.x = element_text( size = 8 ),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Week 8" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
# + theme( aspect.ratio = 2/1 )

p8 = plot_composition( temp.wk10, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        # axis.title.x = element_text( size = 8 ),
        legend.position = "none",
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Week 10" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" ))  
# + theme( aspect.ratio = 2/1 ) 

p9 = plot_composition( temp.wk14, average_by = "clinical_outcome_wk14" ) +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        #legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.x = element_text( size = 12 ),
        legend.title = element_blank(),  
        # legend.key.height = unit( 0.5, "cm" ), 
        # legend.key.width = unit( 0.5, "cm" ),
         panel.grid.major = element_blank(),
          legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Week 14" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
          # + theme( aspect.ratio = 2/1 )
```

``` r
plot.compositions = (pdonorA + pdonorB + p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9) + plot_layout( nrow = 1 ) 
plot.compositions
```

![](Average_Microbiota_Composition_files/figure-gfm/Final%20plot-1.png)<!-- -->

### Compostion per timepoint per patient

``` r
temp@otu_table %>% rownames()
```

    ##  [1] "Bacteroidaceae"                 "Bacteroidalesfam.incertaesedis"
    ##  [3] "Bifidobacteriaceae"             "Clostridiaceae"                
    ##  [5] "Clostridialesfam.incertaesedis" "Coriobacteriaceae"             
    ##  [7] "Eubacteriaceae"                 "Firmicutesfam.incertaesedis"   
    ##  [9] "Lachnospiraceae"                "Oscillospiraceae"              
    ## [11] "Prevotellaceae"                 "Rikenellaceae"                 
    ## [13] "Ruminococcaceae"                "Sutterellaceae"                
    ## [15] "Veillonellaceae"                "Other"

``` r
# Change names so the other category ends up in the bottom of the graph
rownames( temp@otu_table )[rownames( temp@otu_table ) == "Other"] <- "z_Other"
rownames( temp@tax_table )[rownames( temp@tax_table ) == "Other"] <- "z_Other"
temp = microbiome::merge_taxa2( temp, pattern = "Other", name = "z_Other" )

# temp@otu_table %>% rownames()
# temp@tax_table %>% rownames()
# temp@tax_table

temp
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16 taxa and 207 samples ]
    ## sample_data() Sample Data:       [ 207 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16 taxa by 2 taxonomic ranks ]

``` r
temp.donorA = temp %>%
  ps_filter( subject_id == "Donor A", .keep_all_taxa = TRUE )

temp.donorA@sam_data[[ "timepoint.new" ]] <- temp.donorA@sam_data[[ "timepoint.new" ]] %>% as.factor()
temp.donorA@sam_data[[ "timepoint.new" ]] <- factor( temp.donorA@sam_data[[ "timepoint.new" ]], levels = c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"  ))

temp.donorB = temp %>%
  ps_filter( subject_id == "Donor B", .keep_all_taxa = TRUE )

temp.donorB@sam_data[[ "timepoint.new" ]] <- temp.donorB@sam_data[[ "timepoint.new" ]] %>% as.factor()
temp.donorB@sam_data[[ "timepoint.new" ]] <- factor( temp.donorB@sam_data[[ "timepoint.new" ]], levels = c( "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"  ))

temp.101 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "101", .keep_all_taxa = TRUE )

temp.101@sam_data[["timepoint.new"]] <- factor( temp.101@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.102 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "102", .keep_all_taxa = TRUE ) 

temp.102@sam_data[["timepoint.new"]] <- factor( temp.102@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.103 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "103", .keep_all_taxa = TRUE ) 

temp.103@sam_data[["timepoint.new"]] <- factor( temp.103@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.104 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "104", .keep_all_taxa = TRUE )

temp.104@sam_data[["timepoint.new"]] <- factor( temp.104@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.105 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "105", .keep_all_taxa = TRUE )

temp.105@sam_data[["timepoint.new"]] <- factor( temp.105@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.106 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "106", .keep_all_taxa = TRUE )

temp.106@sam_data[["timepoint.new"]] <- factor( temp.106@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.107 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "107", .keep_all_taxa = TRUE )

temp.107@sam_data[["timepoint.new"]] <- factor( temp.107@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.108 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "108", .keep_all_taxa = TRUE )

temp.108@sam_data[["timepoint.new"]] <- factor( temp.108@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.109 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "109", .keep_all_taxa = TRUE ) 

temp.109@sam_data[["timepoint.new"]] <- factor( temp.109@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.110 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "110", .keep_all_taxa = TRUE ) 

temp.110@sam_data[["timepoint.new"]] <- factor( temp.110@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.111 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "111", .keep_all_taxa = TRUE )

temp.111@sam_data[["timepoint.new"]] <- factor( temp.111@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.112 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "112", .keep_all_taxa = TRUE )

temp.112@sam_data[["timepoint.new"]] <- factor( temp.112@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.113 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "113", .keep_all_taxa = TRUE )

temp.113@sam_data[["timepoint.new"]] <- factor( temp.113@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.114 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "114", .keep_all_taxa = TRUE )

temp.114@sam_data[["timepoint.new"]] <- factor( temp.114@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.115 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "115", .keep_all_taxa = TRUE )

temp.115@sam_data[["timepoint.new"]] <- factor( temp.115@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.117 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "117", .keep_all_taxa = TRUE )

temp.117@sam_data[["timepoint.new"]] <- factor( temp.117@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.118 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "118", .keep_all_taxa = TRUE )

temp.118@sam_data[["timepoint.new"]] <- factor( temp.118@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.119 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "119", .keep_all_taxa = TRUE ) 

temp.119@sam_data[["timepoint.new"]] <- factor( temp.119@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.120 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "120", .keep_all_taxa = TRUE ) 

temp.120@sam_data[["timepoint.new"]] <- factor( temp.120@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.121 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "121", .keep_all_taxa = TRUE )

temp.121@sam_data[["timepoint.new"]] <- factor( temp.121@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.122 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "122", .keep_all_taxa = TRUE )

temp.122@sam_data[["timepoint.new"]] <- factor( temp.122@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.123 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "123", .keep_all_taxa = TRUE )

temp.123@sam_data[["timepoint.new"]] <- factor( temp.123@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.124 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "124", .keep_all_taxa = TRUE )

temp.124@sam_data[["timepoint.new"]] <- factor( temp.124@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))

temp.125 = temp %>%
  ps_select( subject_id, clinical_outcome_wk14, timepoint.new ) %>%
  ps_filter( subject_id == "125", .keep_all_taxa = TRUE )

temp.125@sam_data[["timepoint.new"]] <- factor( temp.125@sam_data[["timepoint.new"]],
                                               levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))
```

``` r
# here we use plot_bar because we want to show samples that are not available
pdonorA = plot_bar( temp.donorA, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none", axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        panel.grid.major = element_blank(),
        axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Donor A" ) 

pdonorB = plot_bar( temp.donorB, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Donor B" ) 

p101 = plot_bar( temp.101, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 101" ) 

p102 = plot_bar( temp.102, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 102" ) 

p103 = plot_bar( temp.103, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE )+ 
  ggtitle( "Patient 103" ) 

p104 = plot_bar( temp.104, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 104" )  

p105 = plot_bar( temp.105, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 105" ) 

p106 = plot_bar( temp.106, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE )+ 
  ggtitle( "Patient 106" ) 
 
p107 = plot_bar( temp.107, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 107" )  

p108 = plot_bar( temp.108, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 108" ) 

p109 = plot_bar( temp.109, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 109" ) 

p110 = plot_bar( temp.110, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.title = element_blank(),  
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 110" ) 

p111 = plot_bar( temp.111, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 111" ) 

p112 = plot_bar( temp.112, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 112" ) 

p113 = plot_bar( temp.113, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 113" ) 

p114 = plot_bar( temp.114, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.title = element_blank(),  
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 114" ) 

p115 = plot_bar( temp.115, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 115" ) 

p117 = plot_bar( temp.117, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 117" ) 

p118 = plot_bar( temp.118, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.title = element_blank(),  
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 118" ) 

p119 = plot_bar( temp.119, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 119" ) 

p120 = plot_bar( temp.120, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 120" ) 

p121 = plot_bar( temp.121, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 45, hjust = 1 ),
        legend.title = element_blank(),  
        legend.position = "none",
        axis.title.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE )+ 
  ggtitle( "Patient 121" ) 

p122 = plot_bar( temp.122, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 45, hjust = 1 ),
        legend.title = element_blank(),  
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 122" ) 

p123 = plot_bar( temp.123, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank() ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 123" ) 

p124 = plot_bar( temp.124, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 45, hjust = 1 ),
        legend.title = element_blank(),  
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()  ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 124" ) 

p125 = plot_bar( temp.125, "timepoint.new", fill = "family" ) +
  theme_bw() +
  theme( text = element_text(size = 20),
         axis.text.x = element_text( angle = 45, hjust = 1 ),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()  ) + 
  labs( x = "Timepoint", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) +  
  scale_x_discrete( drop = FALSE ) + 
  ggtitle( "Patient 125" ) 
```

``` r
# plot.compositions = (( pdonorA | pdonorB )  / ( p101 | p102 | p103 | p104 ) / ( p105 | p106 | p107 | p108 ) / ( p109 | p110 | p111 | p112 ) / ( p113 | p114 | p115 | p117 ) / ( p118 | p119 | p120 | p121 ) / ( p122 | p123 | p124 | p125 )) + plot_layout( guides = 'collect' )

# plot.compositions = (( pdonorA + labs( tag = 'A' ) | pdonorB + labs( tag = 'B' ))  / ( p101 + labs( tag = 'C' ) | p105 + labs( tag = 'D' ) | p102 + labs( tag = 'E' ) | p103 + labs( tag = 'F' ) ) / ( p108 + labs( tag = 'G' ) | p109 + labs( tag = 'H' ) | p104 + labs( tag = 'I' ) | p106 + labs( tag = 'J' )) / ( p110 + labs( tag = 'K' ) | p111 + labs( tag = 'L' ) | p107 + labs( tag = 'M' ) | p112 + labs( tag = 'N' ) ) / ( p114 + labs( tag = 'O' ) | p117 + labs( tag = 'P' ) | p113 + labs( tag = 'Q' ) | p115 + labs( tag = 'R' ) ) / ( p118 + labs( tag = 'S' ) | p119 + labs( tag = 'T' ) | p120 + labs( tag = 'U' ) | p123 + labs( tag = 'V' ) ) / ( p121 + labs( tag = 'W' ) | p122 + labs( tag = 'X' ) | p124 + labs( tag = 'Y' ) | p125 + labs( tag = 'Z' ) )) + plot_layout( guides = 'collect' )

plot.compositions = (( pdonorA | pdonorB )  / ( p101 | p105 | p102 | p103 ) / ( p108 | p109 | p104 | p106 ) / ( p110 | p111 | p107 | p112 ) / ( p114 | p117 | p113 | p115 ) / ( p118 | p119 | p120 | p123 ) / ( p121 | p122 | p124 | p125 )) + plot_layout( guides = 'collect' )

grid.arrange( patchworkGrob( plot.compositions ), left = textGrob( "Relative abundances", gp = gpar( fontsize = 20 ), rot = 90 ) )
```

![](Average_Microbiota_Composition_files/figure-gfm/Final%20plot%20comp-1.png)<!-- -->

## Relative abundances of top 15 species over time

``` r
# Sample data
sample.data = as.data.frame( as.matrix( temp@sam_data ))
sample.data$timepoint.new = as.factor( sample.data$timepoint.new )
sample.data$timepoint.new = factor( sample.data$timepoint.new, levels = c( "Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14" ))
sample.data$timepoint.new.num = as.numeric( sample.data$timepoint.new )
sample.data$clinical_outcome_wk14 = as.factor( sample.data$clinical_outcome_wk14 )
sample.data$clinical_outcome_wk14 = factor( sample.data$clinical_outcome_wk14, levels = c( "None", "Good" ))

# Abundance data
abund = as.data.frame( as.matrix( t( temp@otu_table )))

# Combine
df = data.frame( abund, sample.data )
colnames( df )
```

    ##  [1] "Bacteroidaceae"                 "Bacteroidalesfam.incertaesedis"
    ##  [3] "Bifidobacteriaceae"             "Clostridiaceae"                
    ##  [5] "Clostridialesfam.incertaesedis" "Coriobacteriaceae"             
    ##  [7] "Eubacteriaceae"                 "Firmicutesfam.incertaesedis"   
    ##  [9] "Lachnospiraceae"                "Oscillospiraceae"              
    ## [11] "Prevotellaceae"                 "Rikenellaceae"                 
    ## [13] "Ruminococcaceae"                "Sutterellaceae"                
    ## [15] "Veillonellaceae"                "z_Other"                       
    ## [17] "subject_id"                     "filename"                      
    ## [19] "file_id"                        "timepoint"                     
    ## [21] "timepoint.new"                  "sample_id.new"                 
    ## [23] "days_offset"                    "treated_with_donor"            
    ## [25] "age"                            "sex"                           
    ## [27] "pretreatment"                   "clinical_outcome_wk10"         
    ## [29] "clinical_outcome_wk14"          "raw_reads"                     
    ## [31] "human_reads"                    "human_percentage"              
    ## [33] "high_quality_reads"             "timepoint.new.num"

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Bacteroidaceae

``` r
a = xyplot( Bacteroidaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Bacteroidaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Bacteroidaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Bacteroidaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Bacteroidaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Bacteroidaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Bacteroidaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Bacteroidaceae-1.png)<!-- -->

## Plot Bacteroidalesfam.incertaesedis

``` r
a = xyplot( Bacteroidalesfam.incertaesedis ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Bacteroidalesfam.incertaesedis", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Bacteroidalesfam.incertaesedis ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Bacteroidalesfam.incertaesedis ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Bacteroidalesfam.incertaesedis ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Bacteroidalesfam.incertaesedis.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Bacteroidalesfam.incertaesedis.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Bacteroidalesfam.incertaesedis-1.png)<!-- -->

## Plot Bifidobacteriaceae

``` r
a = xyplot( Bifidobacteriaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Bifidobacteriaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Bifidobacteriaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Bifidobacteriaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Bifidobacteriaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Bifidobacteriaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Bifidobacteriaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Bifidobacteriaceae-1.png)<!-- -->

## Plot Clostridiaceae

``` r
a = xyplot( Clostridiaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Clostridiaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Clostridiaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Clostridiaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Clostridiaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Clostridiaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Clostridiaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Clostridiaceae-1.png)<!-- -->

## Plot Clostridialesfam.incertaesedis

``` r
a = xyplot( Clostridialesfam.incertaesedis ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Clostridialesfam.incertaesedis", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Clostridialesfam.incertaesedis ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Clostridialesfam.incertaesedis ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Clostridialesfam.incertaesedis ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Clostridialesfam.incertaesedis.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Clostridialesfam.incertaesedis.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Clostridialesfam.incertaesedis-1.png)<!-- -->

## Plot Coriobacteriaceae

``` r
a = xyplot( Coriobacteriaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Coriobacteriaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Coriobacteriaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Coriobacteriaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Coriobacteriaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Coriobacteriaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Coriobacteriaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Coriobacteriaceae-1.png)<!-- -->

## Plot Eubacteriaceae

``` r
a = xyplot( Eubacteriaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Eubacteriaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Eubacteriaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Eubacteriaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Eubacteriaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Eubacteriaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Eubacteriaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Eubacteriaceae-1.png)<!-- -->

## Plot Firmicutesfam.incertaesedis

``` r
a = xyplot( Firmicutesfam.incertaesedis ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Firmicutesfam.incertaesedis", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Firmicutesfam.incertaesedis ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Firmicutesfam.incertaesedis ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Firmicutesfam.incertaesedis ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Firmicutesfam.incertaesedis.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Firmicutesfam.incertaesedis.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Firmicutesfam.incertaesedis-1.png)<!-- -->

## Plot Lachnospiraceae

``` r
a = xyplot( Lachnospiraceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Lachnospiraceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Lachnospiraceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Lachnospiraceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Lachnospiraceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Lachnospiraceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Lachnospiraceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Lachnospiraceae-1.png)<!-- -->

## Plot Oscillospiraceae

``` r
a = xyplot( Oscillospiraceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Oscillospiraceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Oscillospiraceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Oscillospiraceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Oscillospiraceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Oscillospiraceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Oscillospiraceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Oscillospiraceae-1.png)<!-- -->

## Plot Prevotellaceae

``` r
a = xyplot( Prevotellaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Prevotellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Prevotellaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Prevotellaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Prevotellaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Prevotellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Prevotellaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Prevotellaceae-1.png)<!-- -->

## Plot Rikenellaceae

``` r
a = xyplot( Rikenellaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Rikenellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Rikenellaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Rikenellaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Rikenellaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Rikenellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Rikenellaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Rikenellaceae-1.png)<!-- -->

## Plot Ruminococcaceae

``` r
a = xyplot( Ruminococcaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Ruminococcaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Ruminococcaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Ruminococcaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Ruminococcaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Ruminococcaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Ruminococcaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Ruminococcaceae-1.png)<!-- -->

## Plot Sutterellaceae

``` r
a = xyplot( Sutterellaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Sutterellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Sutterellaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Sutterellaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Sutterellaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Sutterellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Sutterellaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Sutterellaceae-1.png)<!-- -->

## Plot Veillonellaceae

``` r
a = xyplot( Veillonellaceae ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Veillonellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Veillonellaceae ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Veillonellaceae ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Veillonellaceae ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

Veillonellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Veillonellaceae.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20Veillonellaceae-1.png)<!-- -->

## Plot z_Other

``` r
a = xyplot( z_Other ~ timepoint.new, data = df.none, xlab = "Timepoint", ylab = "Relative abundance", main="Other", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "none", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( z_Other ~ timepoint.new, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( z_Other ~ timepoint.new, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( z_Other ~ timepoint.new, data = df.good, type = "p", col = "#82E0AA" )

z_Other.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
z_Other.plot
```

![](Average_Microbiota_Composition_files/figure-gfm/Plot%20z_Other-1.png)<!-- -->

``` r
grid.arrange( Bacteroidaceae.plot, Bacteroidalesfam.incertaesedis.plot, Bifidobacteriaceae.plot, Clostridiaceae.plot, Clostridialesfam.incertaesedis.plot, Coriobacteriaceae.plot, Eubacteriaceae.plot , Firmicutesfam.incertaesedis.plot, Lachnospiraceae.plot, Oscillospiraceae.plot, Prevotellaceae.plot, Rikenellaceae.plot, Ruminococcaceae.plot, Sutterellaceae.plot, Veillonellaceae.plot, z_Other.plot, ncol = 4) 
```

![](Average_Microbiota_Composition_files/figure-gfm/Final%20plot%20abuns-1.png)<!-- -->
