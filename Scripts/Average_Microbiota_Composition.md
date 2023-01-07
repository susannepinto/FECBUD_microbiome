Average Microbiota Compostition
================
December 2022

## Load required packages

``` r
library( dplyr )
library( magrittr )
library( knitr )
library( tidyverse )
```

``` r
library( phyloseq )
library( microbiome )
library( microViz )
library( ggplot2 )
library( microbiomeutilities )
library( RColorBrewer )
library( patchwork )
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
# Remove donors
physeq_mOTU.rel = subset_samples( physeq_mOTU.rel, subject_id != "Donor A" )
physeq_mOTU.rel = subset_samples( physeq_mOTU.rel, subject_id != "Donor B" )
```

``` r
temp <- aggregate_top_taxa2( physeq_mOTU.rel, 15, "family" )
# temp = microbiome::merge_taxa2( temp, taxa = c( "Unknown", "Other" ), name = "Other" )
temp@otu_table %>% rownames()
```

    ##  [1] "Bacteroidaceae"                 "Bacteroidalesfam.incertaesedis"
    ##  [3] "Bifidobacteriaceae"             "Clostridiaceae"                
    ##  [5] "Clostridialesfam.incertaesedis" "Coriobacteriaceae"             
    ##  [7] "Eubacteriaceae"                 "Firmicutesfam.incertaesedis"   
    ##  [9] "Lachnospiraceae"                "Oscillospiraceae"              
    ## [11] "Other"                          "Prevotellaceae"                
    ## [13] "Rikenellaceae"                  "Ruminococcaceae"               
    ## [15] "Streptococcaceae"               "Sutterellaceae"

``` r
temp@otu_table <- temp@otu_table[ c( 1:10, 12:16, 11 ), ]
temp
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16 taxa and 180 samples ]
    ## sample_data() Sample Data:       [ 180 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16 taxa by 2 taxonomic ranks ]

``` r
temp.base = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "Baseline" )

temp.fmt1 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "FMT1" ) 

temp.fmt2 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "FMT2" ) 

temp.fmt3 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "FMT3" )

temp.fmt4 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "FMT4" )

temp.wk7 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "Week7" )

temp.wk8 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "Week8" )

temp.wk10 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "Week10" )

temp.wk14 = temp %>%
  ps_select( clinical_outcome_wk14, timepoint ) %>%
  ps_filter( timepoint == "Week14" )
```

``` r
# Choose colors
mycolors <- colorRampPalette( brewer.pal( 8, "Set1" ))( 16 )
```

``` r
p1 = plot_composition( temp.base, average_by = "clinical_outcome_wk14" ) +
  theme( legend.position = "none", axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        # axis.text = element_text( size = 8 ),
        # axis.title.x = element_text( size = 8 ),
        panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs(x = "Baseline", y = "Relative Abundance" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) + 
  ylim( 0, 1.0 )

p2 = plot_composition( temp.fmt1, average_by = "clinical_outcome_wk14" ) +
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x=element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
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
  theme( axis.text.x = element_text( angle = 0, hjust = 0.6 ),
        legend.text = element_text( face = "italic" ),
        axis.title.y = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),
        # axis.title.x = element_text( size = 12 ),
        legend.title = element_blank(),  
        # legend.key.height = unit( 0.5, "cm" ), 
        # legend.key.width = unit( 0.5, "cm" ),
         panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) + 
  labs( x = "Week 14" ) + 
  scale_fill_manual( values = mycolors ) + 
  scale_x_discrete( labels = c( "NR", "R" )) 
          # + theme( aspect.ratio = 2/1 )
```

``` r
plot.compositions = p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout( nrow = 1 )
plot.compositions
```

![](Average_Microbiota_Composition_files/figure-gfm/Final%20plot-1.png)<!-- -->
