Top 10 families Mixed models
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
library( lattice )
library( latticeExtra )
library( nlme )
library( glmmTMB )
library( splines )
library( lme4 )
library( MuMIn ) 
library( aod )
library( DHARMa )
library( jtools )
library( dotwhisker )
```

``` r
library( broom )
library( broom.mixed )
library( dotwhisker )
library( tidyverse )
library( tibble )
library( RColorBrewer )
```

## Load the data

Here we use the compositional count abundance dataset, else problems
with convergence occur.

## Preprocessing of the data

``` r
# Make compositional
physeq_mOTU = microbiome::transform(physeq_mOTU, "compositional")

# Aggregate to family level
physeq_mOTU = aggregate_taxa( physeq_mOTU, "family" )
physeq_mOTU
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 94 taxa and 207 samples ]
    ## sample_data() Sample Data:       [ 207 samples by 17 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 94 taxa by 6 taxonomic ranks ]

``` r
# Remove donors from data
physeq_mOTU = subset_samples( physeq_mOTU, subject_id != "Donor A" )
physeq_mOTU = subset_samples( physeq_mOTU, subject_id != "Donor B" )
```

``` r
# Sample data
sample.data = as.data.frame(as.matrix(physeq_mOTU@sam_data))
sample.data$timepoint.new = as.factor(sample.data$timepoint.new)
sample.data$timepoint.new = factor(sample.data$timepoint.new, levels = c("Baseline", "Pre-FMT" , "Post-1", "Post-2", "Post-3", "Post-4", "Week8", "Week10", "Week14"))
sample.data$timepoint.new.num = as.numeric(sample.data$timepoint.new)
sample.data$clinical_outcome_wk14 = as.factor(sample.data$clinical_outcome_wk14)
sample.data$clinical_outcome_wk14 = factor(sample.data$clinical_outcome_wk14, levels = c("None", "Good"))
sample.data$sex = as.factor(sample.data$sex)
sample.data$age = as.numeric(sample.data$age)
sample.data$subject_id = as.factor(sample.data$subject_id)
sample.data$pretreatment = as.factor(sample.data$pretreatment)
sample.data$treated_with_donor = as.factor(sample.data$treated_with_donor)

# Abundance data
abund = as.data.frame(as.matrix(t(physeq_mOTU@otu_table)))

# Combine
df = data.frame(abund, sample.data)
```

# Prevotellaceae

## Transform the data

``` r
hist( df$Prevotellaceae )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Prevotellaceae-1.png)<!-- -->

``` r
df$Prevotellaceae.new = asin( sqrt( df$Prevotellaceae ))
hist( df$Prevotellaceae.new )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Prevotellaceae-2.png)<!-- -->

``` r
sum( df$Prevotellaceae.new == 0 ) / length( df$Prevotellaceae.new ) # proportion of zeros
```

    ## [1] 0.3611111

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Prevotellaceae

``` r
a = xyplot( Prevotellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Prevotellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Prevotellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Prevotellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Prevotellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Prevotellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Prevotellaceae.plot
```

![](Significant-mixed-models_files/figure-gfm/Plot%20Prevotellaceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Prevotellaceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Prevotellaceae = lmer( Prevotellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), data = df )

# lme.slope.Prevotellaceae = lmer( Prevotellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), data = df )

glmm.int.Prevotellaceae = glmmTMB( Prevotellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               family = gaussian, 
               data = df, 
               zi = ~ 1,
               REML = TRUE )

# glmm.slope.Prevotellaceae = glmmTMB( Prevotellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), 
               #family = gaussian, 
               #data = df, 
               #zi = ~ 1,
               #REML = TRUE )

AICc( lme.int.Prevotellaceae ) # lower
```

    ## [1] 20.47816

``` r
AICc( glmm.int.Prevotellaceae )
```

    ## [1] 22.80267

# Diagnostics Prevotellaceae

``` r
lme.Prevotellaceae.diag = DHARMa::simulateResiduals( lme.int.Prevotellaceae )
plot( lme.Prevotellaceae.diag )

plotResiduals( lme.Prevotellaceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( lme.Prevotellaceae.diag, form = df$timepoint.new )
plotResiduals( lme.Prevotellaceae.diag, form = df$treated_with_donor )
plotResiduals( lme.Prevotellaceae.diag, form = df$sex ) 
plotResiduals( lme.Prevotellaceae.diag, form = df$age )
plotResiduals( lme.Prevotellaceae.diag, form = df$pretreatment )

testZeroInflation( lme.Prevotellaceae.diag )
testDispersion( lme.Prevotellaceae.diag )

output1 = recalculateResiduals( lme.Prevotellaceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

As there were no major differences in diagnostics between LMM and ZIGMM,
the LMM model was chosen due to its lower AICc value.

``` r
plot_summs( lme.int.Prevotellaceae )
```

![](Significant-mixed-models_files/figure-gfm/Results%20Prevotellaceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( lme.int.Prevotellaceae ), b = fixef( lme.int.Prevotellaceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##        chi2          df           P 
    ## 13.39249451  3.00000000  0.00386031

``` r
wald.test( Sigma = vcov( lme.int.Prevotellaceae ), b = fixef( lme.int.Prevotellaceae ), Terms = 2 )$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.5899848 1.0000000 0.2073283

``` r
wald.test( Sigma = vcov( lme.int.Prevotellaceae ), b = fixef( lme.int.Prevotellaceae ), Terms = 3 )$result # age
```

    ## $chi2
    ##     chi2       df        P 
    ## 1.089428 1.000000 0.296598

``` r
wald.test( Sigma = vcov( lme.int.Prevotellaceae ), b = fixef( lme.int.Prevotellaceae ), Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##        chi2          df           P 
    ## 0.008037087 1.000000000 0.928565509

``` r
wald.test( Sigma = vcov( lme.int.Prevotellaceae ), b = fixef( lme.int.Prevotellaceae ), Terms = 5 )$result # donor
```

    ## $chi2
    ##     chi2       df        P 
    ## 0.305256 1.000000 0.580606

# Lachnospiraceae

## Transform the data

``` r
hist( df$Lachnospiraceae )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Lachnospiraceae-1.png)<!-- -->

``` r
df$Lachnospiraceae.new = asin( sqrt( df$Lachnospiraceae ))
hist( df$Lachnospiraceae.new )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Lachnospiraceae-2.png)<!-- -->

``` r
sum( df$Lachnospiraceae.new == 0 ) / length( df$Lachnospiraceae.new ) # proportion of zeros
```

    ## [1] 0.01111111

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Lachnospiraceae

``` r
a = xyplot( Lachnospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Lachnospiraceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Lachnospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Lachnospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Lachnospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Lachnospiraceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Lachnospiraceae.plot
```

![](Significant-mixed-models_files/figure-gfm/Plot%20Lachnospiraceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Lachnospiraceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Lachnospiraceae = lmer( Lachnospiraceae.new~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), data = df )

# lme.slope.Lachnospiraceaee = lmer( Lachnospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ),data = df )

glmm.int.Lachnospiraceae = glmmTMB( Lachnospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               family = gaussian, 
               data = df, 
               zi = ~ 1,
               REML = TRUE )

# glmm.slope.Lachnospiraceae = glmmTMB( Lachnospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), 
               #family = gaussian, 
               #data = df, 
               #zi = ~ 1,
               #REML = TRUE )

AICc( lme.int.Lachnospiraceae ) # lower
```

    ## [1] -248.4765

``` r
AICc( glmm.int.Lachnospiraceae )
```

    ## [1] -246.1519

# Diagnostics Lachnospiraceae

``` r
lme.Lachnospiraceae.diag = DHARMa::simulateResiduals( lme.int.Lachnospiraceae )
plot( lme.Lachnospiraceae.diag )

plotResiduals( lme.Lachnospiraceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( lme.Lachnospiraceae.diag, form = df$timepoint.new )
plotResiduals( lme.Lachnospiraceae.diag, form = df$treated_with_donor )
plotResiduals( lme.Lachnospiraceae.diag, form = df$sex ) 
plotResiduals( lme.Lachnospiraceae.diag, form = df$age )
plotResiduals( lme.Lachnospiraceae.diag, form = df$pretreatment )

testZeroInflation( lme.Lachnospiraceae.diag )
testDispersion( lme.Lachnospiraceae.diag )

output1 = recalculateResiduals( lme.Lachnospiraceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

LMM was chosen due to its lower AICc value compared to ZIGMM as well as
its adequate diagnostics.

``` r
plot_summs( lme.int.Lachnospiraceae )
```

![](Significant-mixed-models_files/figure-gfm/Results%20Lachnospiraceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( lme.int.Lachnospiraceae ), b = fixef( lme.int.Lachnospiraceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##        chi2          df           P 
    ## 10.66533649  3.00000000  0.01368026

``` r
wald.test( Sigma = vcov( lme.int.Lachnospiraceae ), b = fixef( lme.int.Lachnospiraceae ), Terms = 2 )$result # sex
```

    ## $chi2
    ##       chi2         df          P 
    ## 3.55444847 1.00000000 0.05938604

``` r
wald.test( Sigma = vcov( lme.int.Lachnospiraceae ), b = fixef( lme.int.Lachnospiraceae ), Terms = 3 )$result # age
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.01468133 1.00000000 0.90355912

``` r
wald.test( Sigma = vcov( lme.int.Lachnospiraceae ), b = fixef( lme.int.Lachnospiraceae ), Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.2188204 1.0000000 0.6399402

``` r
wald.test( Sigma = vcov( lme.int.Lachnospiraceae ), b = fixef( lme.int.Lachnospiraceae ), Terms = 5 )$result # donor
```

    ## $chi2
    ##     chi2       df        P 
    ## 0.115300 1.000000 0.734189

# Ruminococcaceae

## Transform the data

``` r
hist( df$Ruminococcaceae )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Ruminococcaceae-1.png)<!-- -->

``` r
# No transformation needed
```

``` r
sum( df$Ruminococcaceae == 0 ) / length( df$Ruminococcaceae ) # proportion of zeros
```

    ## [1] 0

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Ruminococcaceae

``` r
a = xyplot( Ruminococcaceae ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Ruminococcaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Ruminococcaceae ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Ruminococcaceae ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Ruminococcaceae ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Ruminococcaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Ruminococcaceae.plot
```

![](Significant-mixed-models_files/figure-gfm/Plot%20Ruminococcaceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Ruminococcaceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Ruminococcaceae = lme( Ruminococcaceae ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
               data = df, 
               random = ~ 1 | subject_id )

# lme.slope.Ruminococcaceae = lme( Ruminococcaceae ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
               #data = df, 
               #random = ~ ns( timepoint.new.num, knots = 7 ) | subject_id )

summary( lme.int.Ruminococcaceae )
```

    ## Linear mixed-effects model fit by REML
    ##   Data: df 
    ##         AIC       BIC   logLik
    ##   -344.7528 -307.1232 184.3764
    ## 
    ## Random effects:
    ##  Formula: ~1 | subject_id
    ##         (Intercept)  Residual
    ## StdDev:  0.05669035 0.0668476
    ## 
    ## Fixed effects:  Ruminococcaceae ~ sex + age + pretreatment + treated_with_donor +      clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) 
    ##                                                                   Value
    ## (Intercept)                                                  0.23843953
    ## sexM                                                        -0.00129266
    ## age                                                         -0.00031894
    ## pretreatmentplacebo                                         -0.00374556
    ## treated_with_donorDonor B                                    0.02594808
    ## clinical_outcome_wk14Good                                   -0.05993587
    ## ns(timepoint.new.num, knots = 7)1                           -0.10750425
    ## ns(timepoint.new.num, knots = 7)2                           -0.02545539
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1  0.18376382
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 -0.01081229
    ##                                                              Std.Error  DF
    ## (Intercept)                                                 0.04545174 152
    ## sexM                                                        0.02789324  18
    ## age                                                         0.00085208  18
    ## pretreatmentplacebo                                         0.02743016  18
    ## treated_with_donorDonor B                                   0.02960734  18
    ## clinical_outcome_wk14Good                                   0.03706142  18
    ## ns(timepoint.new.num, knots = 7)1                           0.03710392 152
    ## ns(timepoint.new.num, knots = 7)2                           0.02836876 152
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1 0.05530239 152
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 0.03806099 152
    ##                                                               t-value p-value
    ## (Intercept)                                                  5.245994  0.0000
    ## sexM                                                        -0.046343  0.9635
    ## age                                                         -0.374307  0.7125
    ## pretreatmentplacebo                                         -0.136549  0.8929
    ## treated_with_donorDonor B                                    0.876407  0.3924
    ## clinical_outcome_wk14Good                                   -1.617204  0.1232
    ## ns(timepoint.new.num, knots = 7)1                           -2.897383  0.0043
    ## ns(timepoint.new.num, knots = 7)2                           -0.897303  0.3710
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1  3.322891  0.0011
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 -0.284078  0.7767
    ##  Correlation: 
    ##                                                             (Intr) sexM  
    ## sexM                                                        -0.310       
    ## age                                                         -0.761 -0.009
    ## pretreatmentplacebo                                         -0.008  0.043
    ## treated_with_donorDonor B                                   -0.168 -0.302
    ## clinical_outcome_wk14Good                                   -0.209  0.313
    ## ns(timepoint.new.num, knots = 7)1                           -0.253 -0.009
    ## ns(timepoint.new.num, knots = 7)2                            0.042  0.035
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1  0.170  0.006
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 -0.031 -0.026
    ##                                                             age    prtrtm
    ## sexM                                                                     
    ## age                                                                      
    ## pretreatmentplacebo                                         -0.341       
    ## treated_with_donorDonor B                                    0.067  0.028
    ## clinical_outcome_wk14Good                                   -0.087  0.026
    ## ns(timepoint.new.num, knots = 7)1                           -0.015  0.039
    ## ns(timepoint.new.num, knots = 7)2                           -0.017  0.050
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1  0.010 -0.026
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2  0.012 -0.037
    ##                                                             tr__DB c__14G
    ## sexM                                                                     
    ## age                                                                      
    ## pretreatmentplacebo                                                      
    ## treated_with_donorDonor B                                                
    ## clinical_outcome_wk14Good                                   -0.422       
    ## ns(timepoint.new.num, knots = 7)1                           -0.042  0.342
    ## ns(timepoint.new.num, knots = 7)2                           -0.033 -0.038
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1  0.028 -0.538
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2  0.025  0.054
    ##                                                             n(..,k=7)1
    ## sexM                                                                  
    ## age                                                                   
    ## pretreatmentplacebo                                                   
    ## treated_with_donorDonor B                                             
    ## clinical_outcome_wk14Good                                             
    ## ns(timepoint.new.num, knots = 7)1                                     
    ## ns(timepoint.new.num, knots = 7)2                            0.018    
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1 -0.671    
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 -0.014    
    ##                                                             n(..,k=7)2
    ## sexM                                                                  
    ## age                                                                   
    ## pretreatmentplacebo                                                   
    ## treated_with_donorDonor B                                             
    ## clinical_outcome_wk14Good                                             
    ## ns(timepoint.new.num, knots = 7)1                                     
    ## ns(timepoint.new.num, knots = 7)2                                     
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1 -0.012    
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 -0.745    
    ##                                                             c__14G:(..,k=7)1
    ## sexM                                                                        
    ## age                                                                         
    ## pretreatmentplacebo                                                         
    ## treated_with_donorDonor B                                                   
    ## clinical_outcome_wk14Good                                                   
    ## ns(timepoint.new.num, knots = 7)1                                           
    ## ns(timepoint.new.num, knots = 7)2                                           
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)1                 
    ## clinical_outcome_wk14Good:ns(timepoint.new.num, knots = 7)2 -0.029          
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -2.3145998 -0.5266265 -0.0244018  0.4820287  4.2885897 
    ## 
    ## Number of Observations: 180
    ## Number of Groups: 24

# Diagnostics Ruminococcaceae

``` r
lmer.int.Ruminococcaceae = lmer( Ruminococcaceae ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ),
               data = df )

lmer.Ruminococcaceae.diag = DHARMa::simulateResiduals( lmer.int.Ruminococcaceae )
plot( lmer.Ruminococcaceae.diag )

plotResiduals( lmer.Ruminococcaceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( lmer.Ruminococcaceae.diag, form = df$timepoint.new )
plotResiduals( lmer.Ruminococcaceae.diag, form = df$treated_with_donor )
plotResiduals( lmer.Ruminococcaceae.diag, form = df$sex ) 
plotResiduals( lmer.Ruminococcaceae.diag, form = df$age )
plotResiduals( lmer.Ruminococcaceae.diag, form = df$pretreatment )

testZeroInflation( lmer.Ruminococcaceae.diag )
testDispersion( lmer.Ruminococcaceae.diag )

output1 = recalculateResiduals( lmer.Ruminococcaceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

LMM was chosen due to its lower AICc value and adequate diagnostics.

``` r
plot_summs( lme.int.Ruminococcaceae )
```

![](Significant-mixed-models_files/figure-gfm/Results%20Ruminococcaceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( lme.int.Ruminococcaceae ), b = fixef( lme.int.Ruminococcaceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##        chi2          df           P 
    ## 11.12211779  3.00000000  0.01108354

``` r
wald.test( Sigma = vcov( lme.int.Ruminococcaceae ), b = fixef( lme.int.Ruminococcaceae ), Terms = 2 )$result # sex
```

    ## $chi2
    ##        chi2          df           P 
    ## 0.002147698 1.000000000 0.963036659

``` r
wald.test( Sigma = vcov( lme.int.Ruminococcaceae ), b = fixef( lme.int.Ruminococcaceae ), Terms = 3 )$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.1401060 1.0000000 0.7081757

``` r
wald.test( Sigma = vcov( lme.int.Ruminococcaceae ), b = fixef( lme.int.Ruminococcaceae ), Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.01864563 1.00000000 0.89138730

``` r
wald.test( Sigma = vcov( lme.int.Ruminococcaceae ), b = fixef( lme.int.Ruminococcaceae ), Terms = 5 )$result # donor
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.7680893 1.0000000 0.3808088

# Oscillospiraceae

## Transform the data

``` r
hist( df$Oscillospiraceae )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Oscillospiraceae-1.png)<!-- -->

``` r
df$Oscillospiraceae.new = asin( sqrt( df$Oscillospiraceae ))
hist( df$Oscillospiraceae.new )
```

![](Significant-mixed-models_files/figure-gfm/Histogramm%20Oscillospiraceae-2.png)<!-- -->

``` r
sum( df$Oscillospiraceae.new == 0 ) / length( df$Oscillospiraceae.new ) # proportion of zeros
```

    ## [1] 0.1

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Oscillospiraceae

``` r
a = xyplot( Oscillospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Oscillospiraceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Oscillospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Oscillospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Oscillospiraceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Oscillospiraceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Oscillospiraceae.plot
```

![](Significant-mixed-models_files/figure-gfm/Plot%20Oscillospiraceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Oscillospiraceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Oscillospiraceae = lme( Oscillospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7), 
               data = df, 
               random = ~ 1 | subject_id )

# lme.slope.Oscillospiraceae = lme( Oscillospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7), 
               #data = df, 
               #random = ~ ns( timepoint.new.num, knots = 7 ) | subject_id )

glmm.int.Oscillospiraceae = glmmTMB( Oscillospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               family = gaussian, 
               data = df, 
               zi= ~ 1 )

#glmm.slope.Oscillospiraceae. = glmmTMB( Oscillospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ 1 )

AICc( lme.int.Oscillospiraceae )
```

    ## [1] -399.2128

``` r
AICc( glmm.int.Oscillospiraceae ) # lower
```

    ## [1] -459.8495

# Diagnostics Oscillospiraceae

``` r
lmer.int.Oscillospiraceae = lmer(Oscillospiraceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

lme.Oscillospiraceae.diag = DHARMa::simulateResiduals( lmer.int.Oscillospiraceae )
plot( lme.Oscillospiraceae.diag )

plotResiduals( lme.Oscillospiraceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( lme.Oscillospiraceae.diag, form = df$timepoint.new )
plotResiduals( lme.Oscillospiraceae.diag, form = df$treated_with_donor )
plotResiduals( lme.Oscillospiraceae.diag, form = df$sex ) 
plotResiduals( lme.Oscillospiraceae.diag, form = df$age )
plotResiduals( lme.Oscillospiraceae.diag, form = df$pretreatment )

testZeroInflation( lme.Oscillospiraceae.diag )
testDispersion( lme.Oscillospiraceae.diag )

output1 = recalculateResiduals( lme.Oscillospiraceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

LMM was chosen due to better diagnostics than ZIGMM despite ZIGMM having
a lower AICc value.

``` r
plot_summs( lme.int.Oscillospiraceae )
```

![](Significant-mixed-models_files/figure-gfm/Results%20Oscillospiraceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( lme.int.Oscillospiraceae ), b = fixef( lme.int.Oscillospiraceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##       chi2         df          P 
    ## 9.89128300 3.00000000 0.01951324

``` r
wald.test( Sigma = vcov( lme.int.Oscillospiraceae ), b = fixef( lme.int.Oscillospiraceae ), Terms = 2 )$result # sex
```

    ## $chi2
    ##     chi2       df        P 
    ## 0.549351 1.000000 0.458583

``` r
wald.test( Sigma = vcov( lme.int.Oscillospiraceae ), b = fixef( lme.int.Oscillospiraceae ), Terms = 3 )$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 2.2284132 1.0000000 0.1354929

``` r
wald.test( Sigma = vcov( lme.int.Oscillospiraceae ), b = fixef( lme.int.Oscillospiraceae ), Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.3570830 1.0000000 0.5501307

``` r
wald.test( Sigma = vcov( lme.int.Oscillospiraceae ), b = fixef( lme.int.Oscillospiraceae ), Terms = 5 )$result # donor
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.4209052 1.0000000 0.2332545

# Coefficients Plot - Comparison

``` r
# Prevotellaceae
prevo = broom.mixed::tidy( lme.int.Prevotellaceae, conf.int = TRUE )
prevo = prevo[ -c( 11:12 ), -c( 7:8 )]
model = "Prevotellaceae"
prevo = data.frame( model = model, prevo )
prevo = as_tibble( prevo )

# Lachnospiraceae
lachno = broom.mixed::tidy( lme.int.Lachnospiraceae, conf.int = TRUE )
lachno = lachno[ -c( 11:12 ), -c( 7:8 )]
model = "Lachnospiraceae"
lachno = data.frame( model = model, lachno )
lachno = as_tibble( lachno )

# Ruminococcaceae
rumino = broom.mixed::tidy( lme.int.Ruminococcaceae, conf.int = TRUE )
rumino = rumino[ -c( 11:12 ), -c( 6,8 )]
model = "Ruminococcaceae"
rumino = data.frame( model = model, rumino )
rumino = as_tibble( rumino )

# Oscillospiraceae
oscillo = broom.mixed::tidy( lme.int.Oscillospiraceae, conf.int = TRUE )
oscillo = oscillo[ -c( 11:12 ), -c( 6,8 )]
model = "Oscillospiraceae"
oscillo = data.frame( model = model, oscillo )
oscillo = as_tibble( oscillo )
```

``` r
# combine
common <- intersect( colnames( prevo ), colnames( rumino ))
all_bacteria = rbind( prevo[ common ], lachno[ common ], rumino[ common ], oscillo[ common ])

# pdf( "coefficientsplot.pdf", width = 15, height = 11 )
# dwplot( all_bacteria, dodge_size = 0.8, dot_args = list( size = 3 ),
#   whisker_args = list( size = 2 )) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw( base_size = 4 ) + theme( legend.title = element_blank(), text = element_text( size = 20 ))

mycolors <- colorRampPalette( brewer.pal( 8, "Set1" ))( 16 )

dwplot( all_bacteria, dodge_size = 0.6, dot_args = list(size = 4) ) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw() + theme( text = element_text( size = 20 ), legend.title = element_blank()) + scale_color_manual(values=c("#FFD422", "#E3712B", "#C9992C", "#FF980A"))
```

![](Significant-mixed-models_files/figure-gfm/Coefficient%20Plot-1.png)<!-- -->

``` r
# dev.off()
```
