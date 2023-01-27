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

# make new class - Coriobacteriaceae
physeq_mOTU = merge_taxa2(physeq_mOTU, pattern = "_Coriobacteriaceae", name = "Coriobacteriaceae")
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

# Bacteroidaceae

## Transform the data

``` r
hist( df$Bacteroidaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Bacteroidaceae-1.png)<!-- -->

``` r
df$Bacteroidaceae.new = asin( sqrt( df$Bacteroidaceae ))
hist( df$Bacteroidaceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Bacteroidaceae-2.png)<!-- -->

``` r
sum( df$Bacteroidaceae.new == 0 ) / length( df$Bacteroidaceae.new ) # proportion of zeros
```

    ## [1] 0.02222222

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Bacteroidaceae

``` r
a = xyplot( Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Bacteroidaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Bacteroidaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Bacteroidaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Bacteroidaceae.plot
```

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Bacteroidaceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Bacteroidaceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Bacteroidaceae = lme4::lmer(Bacteroidaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

#lme.int.Bacteroidaceae = lme4::lmer(Bacteroidaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

glmm.int.Bacteroidaceae = glmmTMB(Bacteroidaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi= ~ 1) 

#glmm.slope.Bacteroidaceae = glmmTMB(Bacteroidaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ 1) 

AICc(lme.int.Bacteroidaceae)
```

    ## [1] -274.0412

``` r
AICc(glmm.int.Bacteroidaceae) # lower
```

    ## [1] -328.5049

# Diagnostics Bacteroidaceae

``` r
lme.Bacteroidaceae.diag = DHARMa::simulateResiduals( lme.int.Bacteroidaceae )
plot( lme.Bacteroidaceae.diag )

plotResiduals( lme.Bacteroidaceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( lme.Bacteroidaceae.diag, form = df$timepoint.new )
plotResiduals( lme.Bacteroidaceae.diag, form = df$treated_with_donor )
plotResiduals( lme.Bacteroidaceae.diag, form = df$sex ) 
plotResiduals( lme.Bacteroidaceae.diag, form = df$age )
plotResiduals( lme.Bacteroidaceae.diag, form = df$pretreatment )

testZeroInflation( lme.Bacteroidaceae.diag )
testDispersion( lme.Bacteroidaceae.diag )

output1 = recalculateResiduals( lme.Bacteroidaceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

LMM model was chosen due to slightly better diagnostics.

``` r
plot_summs( lme.int.Bacteroidaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Bacteroidaceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( lme.int.Bacteroidaceae ), b = fixef( lme.int.Bacteroidaceae ), Terms = 2 )$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.2508746 1.0000000 0.2633855

``` r
wald.test( Sigma = vcov( lme.int.Bacteroidaceae ), b = fixef( lme.int.Bacteroidaceae ), Terms = 3 )$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.2390284 1.0000000 0.6249087

``` r
wald.test( Sigma = vcov( lme.int.Bacteroidaceae ), b = fixef( lme.int.Bacteroidaceae ), Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.0103796 1.0000000 0.3148119

``` r
wald.test( Sigma = vcov( lme.int.Bacteroidaceae ), b = fixef( lme.int.Bacteroidaceae ), Terms = 5 )$result # donor
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.08342805 1.00000000 0.77270448

``` r
wald.test( Sigma = vcov( lme.int.Bacteroidaceae ), b = fixef( lme.int.Bacteroidaceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##       chi2         df          P 
    ## 7.72418680 3.00000000 0.05206948

# Bacteroidalesfam.incertaesedis

``` r
sum(df$Bacteroidalesfam.incertaesedis == 0) / length(df$Bacteroidalesfam.incertaesedis)
```

    ## [1] 0.3277778

``` r
xyplot(Bacteroidalesfam.incertaesedis ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Bacteroidalesfam.incertaesedis", ylab = "Bifidobacteriaceae", ylim = c(0, 0.05),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
hist(df$Bacteroidalesfam.incertaesedis)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
df$Bacteroidalesfam.incertaesedis.new = asin(sqrt(df$Bacteroidalesfam.incertaesedis))
hist(df$Bacteroidalesfam.incertaesedis.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
xyplot(Bacteroidalesfam.incertaesedis.new ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Bacteroidalesfam.incertaesedis", ylab = "Bifidobacteriaceae", ylim = c(0, 0.15),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
lme.int.Bacteroidalesfam.incertaesedis = lmer(Bacteroidalesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

#lme.slope.Bacteroidalesfam.incertaesedis = lmer(Bacteroidalesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

glmm.int.Bacteroidalesfam.incertaesedis = glmmTMB(Bacteroidalesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi= ~ clinical_outcome_wk14)

#glmm.slope.Bacteroidalesfam.incertaesedis = glmmTMB(Bacteroidalesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ clinical_outcome_wk14) 

AICc(lme.int.Bacteroidalesfam.incertaesedis)
```

    ## [1] -401.876

``` r
AICc(glmm.int.Bacteroidalesfam.incertaesedis) # lower
```

    ## [1] -472.4924

``` r
# diagnostics zigmm
glmm.Bacteroidalesfam.incertaesedis.diag = DHARMa::simulateResiduals(glmm.int.Bacteroidalesfam.incertaesedis)

plot(glmm.Bacteroidalesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
plotResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
plotResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
plotResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->

``` r
plotResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->

``` r
plotResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-7.png)<!-- -->

``` r
testZeroInflation(glmm.Bacteroidalesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = 54.029, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(glmm.Bacteroidalesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.99185, p-value = 0.952
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(glmm.Bacteroidalesfam.incertaesedis.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-4-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.0085, p-value = 0.9893
    ## alternative hypothesis: true autocorrelation is not 0

``` r
# diagnostics lme
lme.Bacteroidalesfam.incertaesedis.diag = DHARMa::simulateResiduals(lme.int.Bacteroidalesfam.incertaesedis)

plot(lme.Bacteroidalesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plotResiduals(lme.Bacteroidalesfam.incertaesedis.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
plotResiduals(lme.Bacteroidalesfam.incertaesedis.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
plotResiduals(lme.Bacteroidalesfam.incertaesedis.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->

``` r
plotResiduals(lme.Bacteroidalesfam.incertaesedis.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-5.png)<!-- -->

``` r
plotResiduals(lme.Bacteroidalesfam.incertaesedis.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-6.png)<!-- -->

``` r
plotResiduals(lme.Bacteroidalesfam.incertaesedis.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-7.png)<!-- -->

``` r
testZeroInflation(lme.Bacteroidalesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(lme.Bacteroidalesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.81644, p-value = 0.424
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(lme.Bacteroidalesfam.incertaesedis.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-5-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 1.9281, p-value = 0.9101
    ## alternative hypothesis: true autocorrelation is not 0

``` r
wald.test(Sigma = vcov(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], b = fixef(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], Terms = 2)$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.5729724 1.0000000 0.2097759

``` r
wald.test(Sigma = vcov(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], b = fixef(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], Terms = 3)$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.0224265 1.0000000 0.3119441

``` r
wald.test(Sigma = vcov(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], b = fixef(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], Terms = 4)$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 2.0403603 1.0000000 0.1531734

``` r
wald.test(Sigma = vcov(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], b = fixef(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], Terms = 5)$result # donor
```

    ## $chi2
    ##        chi2          df           P 
    ## 0.004073397 1.000000000 0.949111029

``` r
wald.test(Sigma = vcov(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], b = fixef(glmm.int.Bacteroidalesfam.incertaesedis)[["cond"]], Terms = c(6, 9, 10))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.8840161 3.0000000 0.5968249

# Bifidobacteriaceae

``` r
sum(df$Bifidobacteriaceae == 0) / length(df$Bifidobacteriaceae)
```

    ## [1] 0.1833333

``` r
xyplot(Bifidobacteriaceae ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Bifidobacteriaceae", ylim = c(0, 0.03),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
hist(df$Bifidobacteriaceae)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
df$Bifidobacteriaceae.new = asin(sqrt(df$Bifidobacteriaceae))
hist(df$Bifidobacteriaceae.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

``` r
xyplot(Bifidobacteriaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Bifidobacteriaceae", ylim = c(0, 0.15),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
lme.int.Bifidobacteriaceae = lmer(Bifidobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

#lme.slope.Bifidobacteriaceae = lmer(Bifidobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

glmm.int.Bifidobacteriaceae = glmmTMB(Bifidobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi= ~ 1)

#glmm.slope.Bifidobacteriaceae = glmmTMB(Bifidobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ 1) 

AICc(lme.int.Bifidobacteriaceae)
```

    ## [1] -423.3096

``` r
AICc(glmm.int.Bifidobacteriaceae) # lower
```

    ## [1] -488.4

``` r
# diagnostics zigmm
glmm.Bifidobacteriaceae.diag = DHARMa::simulateResiduals(glmm.int.Bifidobacteriaceae)

plot(glmm.Bifidobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
plotResiduals(glmm.Bifidobacteriaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

``` r
plotResiduals(glmm.Bifidobacteriaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->

``` r
plotResiduals(glmm.Bifidobacteriaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

``` r
plotResiduals(glmm.Bifidobacteriaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->

``` r
plotResiduals(glmm.Bifidobacteriaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-6.png)<!-- -->

``` r
plotResiduals(glmm.Bifidobacteriaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-7.png)<!-- -->

``` r
testZeroInflation(glmm.Bifidobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(glmm.Bifidobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 1.0273, p-value = 0.76
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(glmm.Bifidobacteriaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-10-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 1.5908, p-value = 0.5153
    ## alternative hypothesis: true autocorrelation is not 0

``` r
# diagnostics lme
lme.Bifidobacteriaceae.diag = DHARMa::simulateResiduals(lme.int.Bifidobacteriaceae)

plot(lme.Bifidobacteriaceae.diag)
```

    ## Warning in newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS, L =
    ## G$L, : Fitting terminated with step failure - check results carefully

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
plotResiduals(lme.Bifidobacteriaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

``` r
plotResiduals(lme.Bifidobacteriaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->

``` r
plotResiduals(lme.Bifidobacteriaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

``` r
plotResiduals(lme.Bifidobacteriaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

``` r
plotResiduals(lme.Bifidobacteriaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

``` r
plotResiduals(lme.Bifidobacteriaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->

``` r
testZeroInflation(lme.Bifidobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(lme.Bifidobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.90665, p-value = 0.576
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(lme.Bifidobacteriaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 1.5675, p-value = 0.4912
    ## alternative hypothesis: true autocorrelation is not 0

``` r
wald.test(Sigma = vcov(glmm.int.Bifidobacteriaceae)[["cond"]], b = fixef(glmm.int.Bifidobacteriaceae)[["cond"]], Terms = 2)$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.6082058 1.0000000 0.2047442

``` r
wald.test(Sigma = vcov(glmm.int.Bifidobacteriaceae)[["cond"]], b = fixef(glmm.int.Bifidobacteriaceae)[["cond"]], Terms = 3)$result # age
```

    ## $chi2
    ##     chi2       df        P 
    ## 2.697083 1.000000 0.100532

``` r
wald.test(Sigma = vcov(glmm.int.Bifidobacteriaceae)[["cond"]], b = fixef(glmm.int.Bifidobacteriaceae)[["cond"]], Terms = 4)$result # pretreatment
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.06257952 1.00000000 0.80246440

``` r
wald.test(Sigma = vcov(glmm.int.Bifidobacteriaceae)[["cond"]], b = fixef(glmm.int.Bifidobacteriaceae)[["cond"]], Terms = 5)$result # donor
```

    ## $chi2
    ##       chi2         df          P 
    ## 5.17560979 1.00000000 0.02290613

``` r
wald.test(Sigma = vcov(glmm.int.Bifidobacteriaceae)[["cond"]], b = fixef(glmm.int.Bifidobacteriaceae)[["cond"]], Terms = c(6, 9, 10))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 5.9368420 3.0000000 0.1147238

# Clostridiaceae

## Transform the data

``` r
hist( df$Clostridiaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Clostridiaceae-1.png)<!-- -->

``` r
df$Clostridiaceae.new = asin( sqrt( df$Clostridiaceae ))
hist( df$Clostridiaceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Clostridiaceae-2.png)<!-- -->

``` r
sum( df$Clostridiaceae.new == 0 ) / length( df$Clostridiaceae.new ) # proportion of zeros
```

    ## [1] 0.01666667

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Clostridiaceae

``` r
a = xyplot( Clostridiaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Clostridiaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Clostridiaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Clostridiaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Clostridiaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Clostridiaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Clostridiaceae.plot
```

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Clostridiaceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Clostridiaceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Clostridiaceae = lme( Clostridiaceae.new~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
               data = df, 
               random = ~ 1 | subject_id )

# lme.slope.Clostridiaceae = lme( Clostridiaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7), 
              # data = df, 
              # random = ~ ns( timepoint.new.num, knots = 7 ) | subject_id )

glmm.int.Clostridiaceae = glmmTMB( Clostridiaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               family = gaussian, 
               data = df, 
               zi= ~ 1 )

# glmm.slope.Clostridiaceae = glmmTMB( Clostridiaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns(timepoint.new.num, knots = 7 ) | subject_id ), 
               # family = gaussian, 
               # data = df, 
               # zi= ~ 1 )

AICc( lme.int.Clostridiaceae )
```

    ## [1] -373.0101

``` r
AICc( glmm.int.Clostridiaceae ) # lower
```

    ## [1] -431.8529

# Diagnostics Clostridiaceae

``` r
glmm.Clostridiaceae.diag = DHARMa::simulateResiduals( glmm.int.Clostridiaceae )
plot( glmm.Clostridiaceae.diag )

plotResiduals( glmm.Clostridiaceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( glmm.Clostridiaceae.diag, form = df$timepoint.new )
plotResiduals( glmm.Clostridiaceae.diag, form = df$treated_with_donor )
plotResiduals( glmm.Clostridiaceae.diag, form = df$sex ) 
plotResiduals( glmm.Clostridiaceae.diag, form = df$age )
plotResiduals( glmm.Clostridiaceae.diag, form = df$pretreatment )

testZeroInflation( glmm.Clostridiaceae.diag )
testDispersion( glmm.Clostridiaceae.diag )

output1 = recalculateResiduals( glmm.Clostridiaceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

ZIGMM model had a lower AICc and similar diagnostics to the LMM.

``` r
conf.int = data.frame( confint( glmm.int.Clostridiaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

clostridfam = broom.mixed::tidy( glmm.int.Clostridiaceae )[ , -c( 1 )]
clostridfam = clostridfam[ -c( 11:23 ), ]
clostridfam = clostridfam[ , -c( 1 )]
effect = "fixed"

clostridfam = data.frame( model = "Clostridiales incertaesedis", effect = effect, clostridfam, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
clostridfam = clostridfam[ , -c( 8 )]
clostridfam = as_tibble( clostridfam )

dwplot( clostridfam ) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw() + theme( legend.title = element_blank())
```

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Clostridiaceae-1.png)<!-- -->

# Clostridialesfam.incertaesedis

## Transform the data

``` r
hist( df$Clostridialesfam.incertaesedis )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Clostridialesfam.incertaesedis-1.png)<!-- -->

``` r
df$Clostridialesfam.incertaesedis.new = asin( sqrt( df$Clostridialesfam.incertaesedis ))
hist( df$Clostridialesfam.incertaesedis.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Clostridialesfam.incertaesedis-2.png)<!-- -->

``` r
sum( df$Clostridialesfam.incertaesedis.new == 0 ) / length( df$Clostridialesfam.incertaesedis.new ) # proportion of zeros
```

    ## [1] 0.01111111

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Clostridialesfam.incertaesedis

``` r
a = xyplot( Clostridialesfam.incertaesedis.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Clostridialesfam.incertaesedis", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Clostridialesfam.incertaesedis.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Clostridialesfam.incertaesedis.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Clostridialesfam.incertaesedis.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Clostridialesfam.incertaesedis.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Clostridialesfam.incertaesedis.plot
```

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Clostridialesfam.incertaesedis-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Clostridialesfam.incertaesedis

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Clostridialesfam.incertaesedis = lme( Clostridialesfam.incertaesedis.new~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
               data = df, 
               random = ~ 1 | subject_id )

# lme.slope.Clostridialesfam.incertaesedis = lme( Clostridialesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
#                data = df, 
#                random = ~ ns( timepoint.new.num, knots = 7 ) | subject_id )

# glmm.int.Clostridialesfam.incertaesedis = glmmTMB( Clostridialesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
#                family = gaussian, 
#                data = df, 
#                zi = ~ clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ),
#                REML = TRUE )

glmm.slope.Clostridialesfam.incertaesedis = glmmTMB( Clostridialesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), 
               family = gaussian, 
               data = df, 
               zi = ~ clinical_outcome_wk14 + ns(timepoint.new.num, knots = 7 ) + ( 1 | subject_id ),
               REML = TRUE )

AICc( lme.int.Clostridialesfam.incertaesedis ) 
```

    ## [1] -276.8366

``` r
AICc( glmm.slope.Clostridialesfam.incertaesedis ) # lower
```

    ## [1] -274.9811

# Diagnostics Clostridialesfam.incertaesedis

``` r
glmm.Clostridialesfam.incertaesedis.diag = DHARMa::simulateResiduals( glmm.slope.Clostridialesfam.incertaesedis )
plot( glmm.Clostridialesfam.incertaesedis.diag )

plotResiduals( glmm.Clostridialesfam.incertaesedis.diag, form = df$clinical_outcome_wk14 )
plotResiduals( glmm.Clostridialesfam.incertaesedis.diag, form = df$timepoint.new )
plotResiduals( glmm.Clostridialesfam.incertaesedis.diag, form = df$treated_with_donor )
plotResiduals( glmm.Clostridialesfam.incertaesedis.diag, form = df$sex ) 
plotResiduals( glmm.Clostridialesfam.incertaesedis.diag, form = df$age )
plotResiduals( glmm.Clostridialesfam.incertaesedis.diag, form = df$pretreatment )

testZeroInflation( glmm.Clostridialesfam.incertaesedis.diag )
testDispersion( glmm.Clostridialesfam.incertaesedis.diag )

output1 = recalculateResiduals( glmm.Clostridialesfam.incertaesedis.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

As the ZIGMM showed better diagnostics, it was selected as the better
model. AICc values were very similar between LMM and ZIGMM models.

``` r
conf.int = data.frame( confint( glmm.slope.Clostridialesfam.incertaesedis ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

clostridfam = broom.mixed::tidy( glmm.slope.Clostridialesfam.incertaesedis )[ , -c( 1 )]
clostridfam = clostridfam[ -c( 11:23 ), ]
clostridfam = clostridfam[ , -c( 1 )]
effect = "fixed"

clostridfam = data.frame( model = "Clostridiales incertaesedis", effect = effect, clostridfam, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
clostridfam = clostridfam[ , -c( 8 )]
clostridfam = as_tibble( clostridfam )

dwplot( clostridfam ) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw() + theme( legend.title = element_blank())
```

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Clostridialesfam.incertaesedis-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], b = fixef( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], Terms = 2 )$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.6668434 1.0000000 0.1966819

``` r
wald.test( Sigma = vcov( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], b = fixef( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], Terms = 3 )$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 1.1672326 1.0000000 0.2799706

``` r
wald.test( Sigma = vcov( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], b = fixef( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.05131206 1.00000000 0.82079572

``` r
wald.test( Sigma = vcov( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], b = fixef( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], Terms = 5 )$result # donor
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.2329128 1.0000000 0.6293724

``` r
wald.test( Sigma = vcov( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], b = fixef( glmm.slope.Clostridialesfam.incertaesedis )[[ "cond" ]], Terms = c( 6, 9, 10 ))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.5442520 3.0000000 0.9090673

# Coriobacteriaceae

``` r
sum(df$Coriobacteriaceae == 0) / length(df$Coriobacteriaceae) 
```

    ## [1] 0.1444444

``` r
xyplot(Coriobacteriaceae ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Coriobacteriaceae", ylim = c(0, 0.04),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
hist(df$Coriobacteriaceae)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
df$Coriobacteriaceae.new = asin(sqrt(df$Coriobacteriaceae))
hist(df$Coriobacteriaceae.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->

``` r
xyplot(Coriobacteriaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Coriobacteriaceae", ylim = c(0, 0.2),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->

``` r
lme.int.Coriobacteriaceae = lmer(Coriobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

#lme.slope.Coriobacteriaceae = lmer(Coriobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

glmm.int.Coriobacteriaceae = glmmTMB(Coriobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi= ~ 1)

#glmm.slope.Coriobacteriaceae = glmmTMB(Coriobacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ 1) 

AICc(lme.int.Coriobacteriaceae)
```

    ## [1] -448.087

``` r
AICc(glmm.int.Coriobacteriaceae) # lower
```

    ## [1] -515.3126

``` r
# diagnostics zigmm
glmm.Coriobacteriaceae.diag = DHARMa::simulateResiduals(glmm.int.Coriobacteriaceae)

plot(glmm.Coriobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
plotResiduals(glmm.Coriobacteriaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

``` r
plotResiduals(glmm.Coriobacteriaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->

``` r
plotResiduals(glmm.Coriobacteriaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->

``` r
plotResiduals(glmm.Coriobacteriaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-5.png)<!-- -->

``` r
plotResiduals(glmm.Coriobacteriaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-6.png)<!-- -->

``` r
plotResiduals(glmm.Coriobacteriaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-7.png)<!-- -->

``` r
testZeroInflation(glmm.Coriobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(glmm.Coriobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 1.0031, p-value = 0.92
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(glmm.Coriobacteriaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-16-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.7067, p-value = 0.2503
    ## alternative hypothesis: true autocorrelation is not 0

``` r
# diagnostics lme
lme.Coriobacteriaceae.diag = DHARMa::simulateResiduals(lme.int.Coriobacteriaceae)

plot(lme.Coriobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plotResiduals(lme.Coriobacteriaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
plotResiduals(lme.Coriobacteriaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->

``` r
plotResiduals(lme.Coriobacteriaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-4.png)<!-- -->

``` r
plotResiduals(lme.Coriobacteriaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-5.png)<!-- -->

``` r
plotResiduals(lme.Coriobacteriaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-6.png)<!-- -->

``` r
plotResiduals(lme.Coriobacteriaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-7.png)<!-- -->

``` r
testZeroInflation(lme.Coriobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(lme.Coriobacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.88804, p-value = 0.464
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(lme.Coriobacteriaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-17-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.6932, p-value = 0.2602
    ## alternative hypothesis: true autocorrelation is not 0

``` r
wald.test(Sigma = vcov(glmm.int.Coriobacteriaceae)[["cond"]], b = fixef(glmm.int.Coriobacteriaceae)[["cond"]], Terms = 2)$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.2755894 1.0000000 0.5996068

``` r
wald.test(Sigma = vcov(glmm.int.Coriobacteriaceae)[["cond"]], b = fixef(glmm.int.Coriobacteriaceae)[["cond"]], Terms = 3)$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 3.4463738 1.0000000 0.0633907

``` r
wald.test(Sigma = vcov(glmm.int.Coriobacteriaceae)[["cond"]], b = fixef(glmm.int.Coriobacteriaceae)[["cond"]], Terms = 4)$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.4753432 1.0000000 0.4905392

``` r
wald.test(Sigma = vcov(glmm.int.Coriobacteriaceae)[["cond"]], b = fixef(glmm.int.Coriobacteriaceae)[["cond"]], Terms = 5)$result # donor
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.06039145 1.00000000 0.80587831

``` r
wald.test(Sigma = vcov(glmm.int.Coriobacteriaceae)[["cond"]], b = fixef(glmm.int.Coriobacteriaceae)[["cond"]], Terms = c(6, 9, 10))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 5.2953155 3.0000000 0.1514067

# Eubacteriaceae

``` r
sum(df$Eubacteriaceae == 0) / length(df$Eubacteriaceae)
```

    ## [1] 0.03333333

``` r
xyplot(Eubacteriaceae ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Eubacteriaceae",
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
hist(df$Eubacteriaceae)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

``` r
df$Eubacteriaceae.new = asin(sqrt(df$Eubacteriaceae))
hist(df$Eubacteriaceae.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->

``` r
xyplot(Eubacteriaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Eubacteriaceae",
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->

``` r
library(nlme)

lme.int.Eubacteriaceae = lme(Eubacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7), 
               data = df, 
               random = ~ 1 | subject_id)

lme.slope.Eubacteriaceae = lme(Eubacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7), 
               data = df, 
               random = ~ ns(timepoint.new.num, knots = 7) | subject_id)

anova_results = anova(lme.int.Eubacteriaceae, lme.slope.Eubacteriaceae)
mixtureLRT = anova_results[["L.Ratio"]][2]
0.5 * pchisq(mixtureLRT, 1, lower.tail = FALSE) + 0.5 * pchisq(mixtureLRT, 2, lower.tail = FALSE) # use random-slopes
```

    ## [1] 0.00172291

``` r
glmm.int.Eubacteriaceae = glmmTMB(Eubacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi = ~ 1,
               REML = TRUE)

lme.slope.Eubacteriaceae = lmer(Eubacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

#glmm.slope.Eubacteriaceae = glmmTMB(Eubacteriaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi = ~ 1,
               #REML = TRUE)

AICc(lme.slope.Eubacteriaceae) # lower
```

    ## [1] -298.627

``` r
AICc(glmm.int.Eubacteriaceae)
```

    ## [1] -296.4785

``` r
# diagnostics lme
lme.Eubacteriaceae.diag = DHARMa::simulateResiduals(lme.slope.Eubacteriaceae)

plot(lme.Eubacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
plotResiduals(lme.Eubacteriaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

``` r
plotResiduals(lme.Eubacteriaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-3.png)<!-- -->

``` r
plotResiduals(lme.Eubacteriaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-4.png)<!-- -->

``` r
plotResiduals(lme.Eubacteriaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-5.png)<!-- -->

``` r
plotResiduals(lme.Eubacteriaceae.diag, form = df$age)
```

    ## qu = 0.75, log(sigma) = -2.568769 : outer Newton did not converge fully.

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-6.png)<!-- -->

``` r
plotResiduals(lme.Eubacteriaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-7.png)<!-- -->

``` r
testZeroInflation(lme.Eubacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(lme.Eubacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.8617, p-value = 0.528
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(lme.Eubacteriaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-22-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.0318, p-value = 0.9601
    ## alternative hypothesis: true autocorrelation is not 0

``` r
# diagnostics zigmm
glmm.Eubacteriaceae.diag = DHARMa::simulateResiduals(glmm.int.Eubacteriaceae)

plot(glmm.Eubacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
plotResiduals(glmm.Eubacteriaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

``` r
plotResiduals(glmm.Eubacteriaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->

``` r
plotResiduals(glmm.Eubacteriaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->

``` r
plotResiduals(glmm.Eubacteriaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->

``` r
plotResiduals(glmm.Eubacteriaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-6.png)<!-- -->

``` r
plotResiduals(glmm.Eubacteriaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-7.png)<!-- -->

``` r
testZeroInflation(glmm.Eubacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(glmm.Eubacteriaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.85777, p-value = 0.304
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(glmm.Eubacteriaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-23-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 1.8573, p-value = 0.8225
    ## alternative hypothesis: true autocorrelation is not 0

``` r
wald.test(Sigma = vcov(lme.slope.Eubacteriaceae), b = fixef(lme.slope.Eubacteriaceae), Terms = 2)$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.4357754 1.0000000 0.5091686

``` r
wald.test(Sigma = vcov(lme.slope.Eubacteriaceae), b = fixef(lme.slope.Eubacteriaceae), Terms = 3)$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.1476811 1.0000000 0.7007613

``` r
wald.test(Sigma = vcov(lme.slope.Eubacteriaceae), b = fixef(lme.slope.Eubacteriaceae), Terms = 4)$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.4566400 1.0000000 0.4991985

``` r
wald.test(Sigma = vcov(lme.slope.Eubacteriaceae), b = fixef(lme.slope.Eubacteriaceae), Terms = 5)$result # donor
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.9215501 1.0000000 0.3370683

``` r
wald.test(Sigma = vcov(lme.slope.Eubacteriaceae), b = fixef(lme.slope.Eubacteriaceae), Terms = c(6, 9, 10))$result
```

    ## $chi2
    ##     chi2       df        P 
    ## 1.593260 3.000000 0.660919

# Firmicutesfam.incertaesedis

``` r
sum(df$Firmicutesfam.incertaesedis == 0) / length(df$Firmicutesfam.incertaesedis)
```

    ## [1] 0.1388889

``` r
xyplot(Firmicutesfam.incertaesedis ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Firmicutesfam.incertaesedis", ylim = c(0, 0.04),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

``` r
hist(df$Firmicutesfam.incertaesedis)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

``` r
df$Firmicutesfam.incertaesedis.new = asin(sqrt(df$Firmicutesfam.incertaesedis))
hist(df$Firmicutesfam.incertaesedis.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->

``` r
xyplot(Firmicutesfam.incertaesedis ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Firmicutesfam.incertaesedis", ylim = c(0, 0.04), 
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->

``` r
lme.int.Firmicutesfam.incertaesedis = lmer(Firmicutesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

#lme.slope.Firmicutesfam.incertaesedis = lmer(Firmicutesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

glmm.int.Firmicutesfam.incertaesedis = glmmTMB(Firmicutesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi = ~ 1,
               REML = TRUE)

#glmm.slope.Firmicutesfam.incertaesedis = glmmTMB(Firmicutesfam.incertaesedis.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi = ~ clinical_outcome_wk14 + ( 1 | subject_id),
               #REML = TRUE)

AICc(lme.int.Firmicutesfam.incertaesedis) # lower
```

    ## [1] -487.2808

``` r
AICc(glmm.int.Firmicutesfam.incertaesedis)
```

    ## [1] -484.9563

``` r
# diagnostics lme
lme.Firmicutesfam.incertaesedis.diag = DHARMa::simulateResiduals(lme.int.Firmicutesfam.incertaesedis)

plot(lme.Firmicutesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
plotResiduals(lme.Firmicutesfam.incertaesedis.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
plotResiduals(lme.Firmicutesfam.incertaesedis.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->

``` r
plotResiduals(lme.Firmicutesfam.incertaesedis.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-4.png)<!-- -->

``` r
plotResiduals(lme.Firmicutesfam.incertaesedis.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-5.png)<!-- -->

``` r
plotResiduals(lme.Firmicutesfam.incertaesedis.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-6.png)<!-- -->

``` r
plotResiduals(lme.Firmicutesfam.incertaesedis.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-7.png)<!-- -->

``` r
testZeroInflation(lme.Firmicutesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(lme.Firmicutesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.88674, p-value = 0.576
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(lme.Firmicutesfam.incertaesedis.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-28-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.5984, p-value = 0.3353
    ## alternative hypothesis: true autocorrelation is not 0

``` r
# diagnostics zigmm
glmm.Firmicutesfam.incertaesedis.diag = DHARMa::simulateResiduals(glmm.int.Firmicutesfam.incertaesedis)

plot(glmm.Firmicutesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
plotResiduals(glmm.Firmicutesfam.incertaesedis.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

``` r
plotResiduals(glmm.Firmicutesfam.incertaesedis.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-3.png)<!-- -->

``` r
plotResiduals(glmm.Firmicutesfam.incertaesedis.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-4.png)<!-- -->

``` r
plotResiduals(glmm.Firmicutesfam.incertaesedis.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-5.png)<!-- -->

``` r
plotResiduals(glmm.Firmicutesfam.incertaesedis.diag, form = df$age)
```

    ## qu = 0.25, log(sigma) = -2.845019 : outer Newton did not converge fully.

    ## qu = 0.25, log(sigma) = -2.749399 : outer Newton did not converge fully.

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-6.png)<!-- -->

``` r
plotResiduals(glmm.Firmicutesfam.incertaesedis.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-7.png)<!-- -->

``` r
testZeroInflation(glmm.Firmicutesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(glmm.Firmicutesfam.incertaesedis.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.87834, p-value = 0.48
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(glmm.Firmicutesfam.incertaesedis.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-29-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 2.7181, p-value = 0.2422
    ## alternative hypothesis: true autocorrelation is not 0

``` r
wald.test(Sigma = vcov(lme.int.Firmicutesfam.incertaesedis), b = fixef(lme.int.Firmicutesfam.incertaesedis), Terms = 2)$result # sex
```

    ## $chi2
    ##        chi2          df           P 
    ## 8.497455405 1.000000000 0.003556435

``` r
wald.test(Sigma = vcov(lme.int.Firmicutesfam.incertaesedis), b = fixef(lme.int.Firmicutesfam.incertaesedis), Terms = 3)$result # age
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.09330206 1.00000000 0.76002051

``` r
wald.test(Sigma = vcov(lme.int.Firmicutesfam.incertaesedis), b = fixef(lme.int.Firmicutesfam.incertaesedis), Terms = 4)$result # pretreatment
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.02402493 1.00000000 0.87682149

``` r
wald.test(Sigma = vcov(lme.int.Firmicutesfam.incertaesedis), b = fixef(lme.int.Firmicutesfam.incertaesedis), Terms = 5)$result # donor
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.3843411 1.0000000 0.5352891

``` r
wald.test(Sigma = vcov(lme.int.Firmicutesfam.incertaesedis), b = fixef(lme.int.Firmicutesfam.incertaesedis), Terms = c(6, 9, 10))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 2.7735028 3.0000000 0.4278804

# Lachnospiraceae

## Transform the data

``` r
hist( df$Lachnospiraceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Lachnospiraceae-1.png)<!-- -->

``` r
df$Lachnospiraceae.new = asin( sqrt( df$Lachnospiraceae ))
hist( df$Lachnospiraceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Lachnospiraceae-2.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Lachnospiraceae-1.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Lachnospiraceae-1.png)<!-- -->

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

``` r
wald.test( Sigma = vcov( lme.int.Lachnospiraceae ), b = fixef( lme.int.Lachnospiraceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##        chi2          df           P 
    ## 10.66533649  3.00000000  0.01368026

# Oscillospiraceae

## Transform the data

``` r
hist( df$Oscillospiraceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Oscillospiraceae-1.png)<!-- -->

``` r
df$Oscillospiraceae.new = asin( sqrt( df$Oscillospiraceae ))
hist( df$Oscillospiraceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Oscillospiraceae-2.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Oscillospiraceae-1.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Oscillospiraceae-1.png)<!-- -->

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

``` r
wald.test( Sigma = vcov( lme.int.Oscillospiraceae ), b = fixef( lme.int.Oscillospiraceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##       chi2         df          P 
    ## 9.89128300 3.00000000 0.01951324

# Prevotellaceae

## Transform the data

``` r
hist( df$Prevotellaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Prevotellaceae-1.png)<!-- -->

``` r
df$Prevotellaceae.new = asin( sqrt( df$Prevotellaceae ))
hist( df$Prevotellaceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Prevotellaceae-2.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Prevotellaceae-1.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Prevotellaceae-1.png)<!-- -->

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

``` r
wald.test( Sigma = vcov( lme.int.Prevotellaceae ), b = fixef( lme.int.Prevotellaceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##        chi2          df           P 
    ## 13.39249451  3.00000000  0.00386031

# Rikenellaceae

## Transform the data

``` r
hist( df$Rikenellaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Rikenellaceae-1.png)<!-- -->

``` r
df$Rikenellaceae.new = asin( sqrt( df$Rikenellaceae ))
hist( df$Rikenellaceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Rikenellaceae-2.png)<!-- -->

``` r
sum( df$Rikenellaceae.new == 0 ) / length( df$Rikenellaceae.new ) # proportion of zeros
```

    ## [1] 0.1166667

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Rikenellaceae

``` r
a = xyplot( Rikenellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Rikenellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Rikenellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Rikenellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Rikenellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Rikenellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Rikenellaceae.plot
```

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Rikenellaceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Rikenellaceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Rikenellaceae = lme( Rikenellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
               data = df, 
               random = ~ 1 | subject_id )

# lme.slope.Rikenellaceae = lme( Rikenellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ), 
               #data = df, 
               #random = ~ ns( timepoint.new.num, knots = 7 ) | subject_id )

glmm.int.Rikenellaceae = glmmTMB( Rikenellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               family = gaussian, 
               data = df, 
               zi= ~ clinical_outcome_wk14 + ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ))

#glmm.slope.Rikenellaceae = glmmTMB(Rikenellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 + ns(timepoint.new.num, knots = 7) + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ clinical_outcome_wk14 + ns(timepoint.new.num, knots = 7) + ( 1 | subject_id))

AICc( lme.int.Rikenellaceae ) 
```

    ## [1] -335.9933

``` r
AICc( glmm.int.Rikenellaceae ) # lower
```

    ## [1] -380.3248

# Diagnostics Rikenellaceae

``` r
glmm.Rikenellaceae.diag = DHARMa::simulateResiduals( glmm.int.Rikenellaceae )
plot( glmm.Rikenellaceae.diag )

plotResiduals( glmm.Rikenellaceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( glmm.Rikenellaceae.diag, form = df$timepoint.new )
plotResiduals( glmm.Rikenellaceae.diag, form = df$treated_with_donor )
plotResiduals( glmm.Rikenellaceae.diag, form = df$sex ) 
plotResiduals( glmm.Rikenellaceae.diag, form = df$age )
plotResiduals( glmm.Rikenellaceae.diag, form = df$pretreatment )

testZeroInflation( glmm.Rikenellaceae.diag )
testDispersion( glmm.Rikenellaceae.diag )

output1 = recalculateResiduals( glmm.Rikenellaceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

ZIGMM was chosen due to a lower AICc value comapred to LMM. However,
there are problems with diagnostics with both models.

``` r
conf.int = data.frame( confint( glmm.int.Rikenellaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

clostridfam = broom.mixed::tidy( glmm.int.Rikenellaceae )[ , -c( 1 )]
clostridfam = clostridfam[ -c( 11:23 ), ]
clostridfam = clostridfam[ , -c( 1 )]
effect = "fixed"

clostridfam = data.frame( model = "Clostridiales incertaesedis", effect = effect, clostridfam, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
clostridfam = clostridfam[ , -c( 8 )]
clostridfam = as_tibble( clostridfam )

dwplot( clostridfam ) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw() + theme( legend.title = element_blank())
```

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Rikenellaceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( glmm.int.Rikenellaceae )[[ "cond" ]], b = fixef( glmm.int.Rikenellaceae )[[ "cond" ]], Terms = 2 )$result # sex
```

    ## $chi2
    ##         chi2           df            P 
    ## 1.314554e+01 1.000000e+00 2.882046e-04

``` r
wald.test( Sigma = vcov( glmm.int.Rikenellaceae )[[ "cond" ]], b = fixef( glmm.int.Rikenellaceae )[[ "cond" ]], Terms = 3 )$result # age
```

    ## $chi2
    ##       chi2         df          P 
    ## 3.48853548 1.00000000 0.06179523

``` r
wald.test( Sigma = vcov( glmm.int.Rikenellaceae )[[ "cond" ]], b = fixef( glmm.int.Rikenellaceae )[[ "cond" ]], Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##       chi2         df          P 
    ## 0.08056699 1.00000000 0.77653050

``` r
wald.test( Sigma = vcov( glmm.int.Rikenellaceae )[[ "cond" ]], b = fixef( glmm.int.Rikenellaceae )[[ "cond" ]], Terms = 5 )$result # donor
```

    ## $chi2
    ##       chi2         df          P 
    ## 4.31030455 1.00000000 0.03788218

``` r
wald.test( Sigma = vcov( glmm.int.Rikenellaceae )[[ "cond" ]], b = fixef( glmm.int.Rikenellaceae )[[ "cond" ]], Terms = c( 6, 9, 10 ))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 4.9318040 3.0000000 0.1768594

# Ruminococcaceae

## Transform the data

``` r
hist( df$Ruminococcaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Ruminococcaceae-1.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Ruminococcaceae-1.png)<!-- -->

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

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Ruminococcaceae-1.png)<!-- -->

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

``` r
wald.test( Sigma = vcov( lme.int.Ruminococcaceae ), b = fixef( lme.int.Ruminococcaceae ), Terms = c( 6, 9, 10 ))$result # clinical outcome wk14
```

    ## $chi2
    ##        chi2          df           P 
    ## 11.12211779  3.00000000  0.01108354

# Sutterellaceae

``` r
sum(df$Sutterellaceae == 0) / length(df$Sutterellaceae) 
```

    ## [1] 0.07222222

``` r
xyplot(Sutterellaceae ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Sutterellaceae", ylim = c(0, 0.05),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
hist(df$Sutterellaceae)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->

``` r
df$Sutterellaceae.new = asin(sqrt(df$Sutterellaceae))
hist(df$Sutterellaceae.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-32-3.png)<!-- -->

``` r
xyplot(Sutterellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df, xlab = "Timepoint", ylab = "Sutterellaceae", ylim = c(0, 0.2),
           panel = function(x, y) {
            panel.average(x, y, horizontal = FALSE, col = "plum", lwd = 4)
       })
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-32-4.png)<!-- -->

``` r
lme.int.Sutterellaceae = lmer(Sutterellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), data = df)

#lme.slope.Sutterellaceae = lmer(Sutterellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), data = df)

glmm.int.Sutterellaceae = glmmTMB(Sutterellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (1 | subject_id), 
               family = gaussian, 
               data = df, 
               zi= ~ 1)

#glmm.slope.Sutterellaceae = glmmTMB(Sutterellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7) + (ns(timepoint.new.num, knots = 7) | subject_id), 
               #family = gaussian, 
               #data = df, 
               #zi= ~ 1)

AICc(lme.int.Sutterellaceae)
```

    ## [1] -448.1483

``` r
AICc(glmm.int.Sutterellaceae) # lower
```

    ## [1] -516.8995

``` r
glmm.Sutterellaceae.diag = DHARMa::simulateResiduals(glmm.int.Sutterellaceae)
plot(glmm.Sutterellaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
plotResiduals(glmm.Sutterellaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

``` r
plotResiduals(glmm.Sutterellaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-3.png)<!-- -->

``` r
plotResiduals(glmm.Sutterellaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-4.png)<!-- -->

``` r
plotResiduals(glmm.Sutterellaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-5.png)<!-- -->

``` r
plotResiduals(glmm.Sutterellaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-6.png)<!-- -->

``` r
plotResiduals(glmm.Sutterellaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-7.png)<!-- -->

``` r
testZeroInflation(glmm.Sutterellaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(glmm.Sutterellaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.99248, p-value = 0.992
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(glmm.Sutterellaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-34-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 1.4524, p-value = 0.3798
    ## alternative hypothesis: true autocorrelation is not 0

``` r
lme.Sutterellaceae.diag = DHARMa::simulateResiduals(lme.int.Sutterellaceae)
plot(lme.Sutterellaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
plotResiduals(lme.Sutterellaceae.diag, form = df$clinical_outcome_wk14)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-2.png)<!-- -->

``` r
plotResiduals(lme.Sutterellaceae.diag, form = df$timepoint.new)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-3.png)<!-- -->

``` r
plotResiduals(lme.Sutterellaceae.diag, form = df$treated_with_donor)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-4.png)<!-- -->

``` r
plotResiduals(lme.Sutterellaceae.diag, form = df$sex)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-5.png)<!-- -->

``` r
plotResiduals(lme.Sutterellaceae.diag, form = df$age)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-6.png)<!-- -->

``` r
plotResiduals(lme.Sutterellaceae.diag, form = df$pretreatment)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-7.png)<!-- -->

``` r
testZeroInflation(lme.Sutterellaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-8.png)<!-- -->

    ## 
    ##  DHARMa zero-inflation test via comparison to expected zeros with
    ##  simulation under H0 = fitted model
    ## 
    ## data:  simulationOutput
    ## ratioObsSim = Inf, p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
testDispersion(lme.Sutterellaceae.diag)
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-9.png)<!-- -->

    ## 
    ##  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
    ##  simulated
    ## 
    ## data:  simulationOutput
    ## dispersion = 0.89253, p-value = 0.416
    ## alternative hypothesis: two.sided

``` r
output1 = recalculateResiduals(lme.Sutterellaceae.diag, group=df$timepoint.new)
testTemporalAutocorrelation(output1, time = unique(df$timepoint.new))
```

![](Top_15_Mixed_Models_files/figure-gfm/unnamed-chunk-35-10.png)<!-- -->

    ## 
    ##  Durbin-Watson test
    ## 
    ## data:  simulationOutput$scaledResiduals ~ 1
    ## DW = 1.4656, p-value = 0.3918
    ## alternative hypothesis: true autocorrelation is not 0

``` r
wald.test(Sigma = vcov(glmm.int.Sutterellaceae)[["cond"]], b = fixef(glmm.int.Sutterellaceae)[["cond"]], Terms = 2)$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.4550010 1.0000000 0.4999696

``` r
wald.test(Sigma = vcov(glmm.int.Sutterellaceae)[["cond"]], b = fixef(glmm.int.Sutterellaceae)[["cond"]], Terms = 3)$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.7871152 1.0000000 0.3749738

``` r
wald.test(Sigma = vcov(glmm.int.Sutterellaceae)[["cond"]], b = fixef(glmm.int.Sutterellaceae)[["cond"]], Terms = 4)$result # pretreatment
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.6540075 1.0000000 0.4186835

``` r
wald.test(Sigma = vcov(glmm.int.Sutterellaceae)[["cond"]], b = fixef(glmm.int.Sutterellaceae)[["cond"]], Terms = 5)$result # donor
```

    ## $chi2
    ##       chi2         df          P 
    ## 6.90647377 1.00000000 0.00858842

``` r
wald.test(Sigma = vcov(glmm.int.Sutterellaceae)[["cond"]], b = fixef(glmm.int.Sutterellaceae)[["cond"]], Terms = c(6, 9, 10))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 8.7938926 3.0000000 0.0321605

# Veillonellaceae

## Transform the data

``` r
hist( df$Veillonellaceae )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Veillonellaceae-1.png)<!-- -->

``` r
df$Veillonellaceae.new = asin( sqrt( df$Veillonellaceae ))
hist( df$Veillonellaceae.new )
```

![](Top_15_Mixed_Models_files/figure-gfm/Histogramm%20Veillonellaceae-2.png)<!-- -->

``` r
sum( df$Veillonellaceae.new == 0 ) / length( df$Veillonellaceae.new ) # proportion of zeros
```

    ## [1] 0.3333333

``` r
# None responders
df.none = subset( df, clinical_outcome_wk14 == "None" ) # changed filter to subset

# Good responders
df.good = subset( df, clinical_outcome_wk14 == "Good" )
```

## Plot Veillonellaceae

``` r
a = xyplot( Veillonellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, xlab = "Timepoint", ylab = "Relative abundance (transformed)", main="Veillonellaceae", 
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#EC7063", lwd = 4 )
           },
           key = list( space = "right", 
                      lines = list( col = c( "#EC7063", "#82E0AA" ), lty = c( 1, 1 ), lwd = 2 ), 
                      text = list( c( "Non-Responders", "Responders" ))))

b = xyplot( Veillonellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, xlab = "Timepoint",
           panel = function( x, y ) {
             panel.average( x, y, horizontal = FALSE, col = "#82E0AA", lwd = 4, type = "l", lty = 1 )
           })

c = xyplot( Veillonellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.none, type = "p", col = "#EC7063" )

d = xyplot( Veillonellaceae.new ~ timepoint.new | clinical_outcome_wk14, data = df.good, type = "p", col = "#82E0AA" )

Veillonellaceae.plot = a + as.layer( b ) + as.layer( c ) + as.layer( d )
Veillonellaceae.plot
```

![](Top_15_Mixed_Models_files/figure-gfm/Plot%20Veillonellaceae-1.png)<!-- -->

## Choose between the lmer and glmmTMB models Veillonellaceae

We have used time as a numeric value and added a spline at knots = 7.

``` r
lme.int.Veillonellaceae = lmer( Veillonellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               data = df )

lme.slope.Veillonellaceae = lmer( Veillonellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), 
               data = df )

# test random slopes vs random intercepts
anova_results = anova( lme.int.Veillonellaceae, lme.slope.Veillonellaceae )
```

    ## refitting model(s) with ML (instead of REML)

``` r
mixtureLRT = anova_results[[ "Chisq" ]][ 2 ]
0.5 * pchisq( mixtureLRT, 1, lower.tail = FALSE ) + 0.5 * pchisq( mixtureLRT, 2, lower.tail = FALSE )
```

    ## [1] 0.09061474

``` r
# use random intercepts

glmm.int.Veillonellaceae = glmmTMB( Veillonellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns(timepoint.new.num, knots = 7 ) + ( 1 | subject_id ), 
               family = gaussian, 
               data = df, 
               zi= ~ 1 ) # can't make this more complicated

glmm.slope.Veillonellaceae = glmmTMB( Veillonellaceae.new ~ sex + age + pretreatment + treated_with_donor + clinical_outcome_wk14 * ns( timepoint.new.num, knots = 7 ) + ( ns( timepoint.new.num, knots = 7 ) | subject_id ), 
               family = gaussian, 
               data = df, 
               zi= ~ 1 ) 

# test random slopes vs random intercepts
anova_results = anova( glmm.int.Veillonellaceae, glmm.slope.Veillonellaceae )
mixtureLRT = anova_results[[ "Chisq" ]][ 2 ]
0.5 * pchisq( mixtureLRT, 1, lower.tail = FALSE ) + 0.5 * pchisq( mixtureLRT, 2, lower.tail = FALSE )
```

    ## [1] 7.668495e-05

``` r
# use random slopes model

AICc( lme.int.Veillonellaceae )
```

    ## [1] -535.7482

``` r
AIC( lme.slope.Veillonellaceae )
```

    ## [1] -532.4421

``` r
AICc( glmm.int.Veillonellaceae )
```

    ## [1] -604.6374

``` r
AICc( glmm.slope.Veillonellaceae ) # lower
```

    ## [1] -610.4771

# Diagnostics Veillonellaceae

``` r
glmm.Veillonellaceae.diag = DHARMa::simulateResiduals( glmm.slope.Veillonellaceae )
plot( glmm.Veillonellaceae.diag )

plotResiduals( glmm.Veillonellaceae.diag, form = df$clinical_outcome_wk14 )
plotResiduals( glmm.Veillonellaceae.diag, form = df$timepoint.new )
plotResiduals( glmm.Veillonellaceae.diag, form = df$treated_with_donor )
plotResiduals( glmm.Veillonellaceae.diag, form = df$sex ) 
plotResiduals( glmm.Veillonellaceae.diag, form = df$age )
plotResiduals( glmm.Veillonellaceae.diag, form = df$pretreatment )

testZeroInflation( glmm.Veillonellaceae.diag )
testDispersion( glmm.Veillonellaceae.diag )

output1 = recalculateResiduals( glmm.Veillonellaceae.diag, group = df$timepoint.new )
testTemporalAutocorrelation( output1, time = unique( df$timepoint.new ))
```

ZIGMM was chosen due to a lower AICc value. Diagnostics between LMM and
ZIGMM were similar.

``` r
conf.int = data.frame( confint( glmm.slope.Veillonellaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

clostridfam = broom.mixed::tidy( glmm.slope.Veillonellaceae )[ , -c( 1 )]
clostridfam = clostridfam[ -c( 11:23 ), ]
clostridfam = clostridfam[ , -c( 1 )]
effect = "fixed"

clostridfam = data.frame( model = "Clostridiales incertaesedis", effect = effect, clostridfam, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
clostridfam = clostridfam[ , -c( 8 )]
clostridfam = as_tibble( clostridfam )

dwplot( clostridfam ) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw() + theme( legend.title = element_blank())
```

![](Top_15_Mixed_Models_files/figure-gfm/Results%20Veillonellaceae-1.png)<!-- -->

``` r
wald.test( Sigma = vcov( glmm.slope.Veillonellaceae )[[ "cond" ]], b = fixef( glmm.slope.Veillonellaceae )[[ "cond" ]], Terms = 2 )$result # sex
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.2926860 1.0000000 0.5885044

``` r
wald.test( Sigma = vcov( glmm.slope.Veillonellaceae )[[ "cond" ]], b = fixef( glmm.slope.Veillonellaceae )[[ "cond" ]], Terms = 3 )$result # age
```

    ## $chi2
    ##      chi2        df         P 
    ## 0.4492686 1.0000000 0.5026825

``` r
wald.test( Sigma = vcov( glmm.slope.Veillonellaceae )[[ "cond" ]], b = fixef( glmm.slope.Veillonellaceae )[[ "cond" ]], Terms = 4 )$result # pretreatment
```

    ## $chi2
    ##         chi2           df            P 
    ## 1.165241e+01 1.000000e+00 6.411939e-04

``` r
wald.test( Sigma = vcov( glmm.slope.Veillonellaceae )[[ "cond" ]], b = fixef( glmm.slope.Veillonellaceae )[[ "cond" ]], Terms = 5 )$result # donor
```

    ## $chi2
    ##       chi2         df          P 
    ## 3.99328027 1.00000000 0.04568205

``` r
wald.test( Sigma = vcov( glmm.slope.Veillonellaceae )[[ "cond" ]], b = fixef( glmm.slope.Veillonellaceae )[[ "cond" ]], Terms = c( 6, 9, 10 ))$result
```

    ## $chi2
    ##      chi2        df         P 
    ## 2.7292301 3.0000000 0.4352828

# Coefficients Plot - Comparison

``` r
# Bacteroidaceae
bacteroid = broom.mixed::tidy( lme.int.Bacteroidaceae, conf.int = TRUE )
bacteroid = bacteroid[ -c( 11:12 ), ]
model = "Bacteroidaceae"
bacteroid = data.frame( model = model, bacteroid )
bacteroid = as_tibble( bacteroid )

# Bacteroidalesfam.incertaesedis
conf.int = data.frame( confint( glmm.int.Bacteroidalesfam.incertaesedis ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

bactfam = broom.mixed::tidy( glmm.int.Bacteroidalesfam.incertaesedis )[ , -c( 1 )]
bactfam = bactfam[ -c( 11:23 ), ]
bactfam = bactfam[ , -c( 1 )]
effect = "fixed"

bactfam = data.frame( model = "Bacteroidalesfam.incertaesedis", effect = effect, bactfam, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
bactfam = bactfam[ , -c( 8 )]
bactfam = as_tibble( bactfam )

# Bifidobacteriaceae
conf.int = data.frame( confint( glmm.int.Bifidobacteriaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

bifidobac = broom.mixed::tidy( glmm.int.Bifidobacteriaceae )[ , -c( 1 )]
bifidobac = bifidobac[ -c( 11:23 ), ]
bifidobac = bifidobac[ , -c( 1 )]
effect = "fixed"

bifidobac = data.frame( model = "bifidobac", effect = effect, bifidobac, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
bifidobac = bifidobac[ , -c( 8 )]
bifidobac = as_tibble( bifidobac )

# Clostridiaceae
conf.int = data.frame( confint( glmm.int.Clostridiaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

clostrid = broom.mixed::tidy( glmm.int.Clostridiaceae )[ , -c( 1 )]
clostrid = clostrid[ -c( 11:23 ), ]
clostrid = clostrid[ , -c( 1 )]
effect = "fixed"

clostrid = data.frame( model = "Clostridiaceae", effect = effect, clostrid, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
clostrid = clostrid[ , -c( 8 )]
clostrid = as_tibble( clostrid )

# Clostridiales incertaesedis
conf.int = data.frame( confint( glmm.slope.Clostridialesfam.incertaesedis ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

clostridfam = broom.mixed::tidy( glmm.slope.Clostridialesfam.incertaesedis )[ , -c( 1 )]
clostridfam = clostridfam[ -c( 11:23 ), ]
clostridfam = clostridfam[ , -c( 1 )]
effect = "fixed"

clostridfam = data.frame( model = "Clostridiales incertaesedis", effect = effect, clostridfam, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
clostridfam = clostridfam[ , -c( 8 )]
clostridfam = as_tibble( clostridfam )

# Coriobacteriaceae
conf.int = data.frame( confint( glmm.int.Coriobacteriaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:21 ), ]

coriobac = broom.mixed::tidy( glmm.int.Coriobacteriaceae )[ , -c( 1 )]
coriobac = coriobac[ -c( 11:23 ), ]
coriobac = coriobac[ , -c( 1 )]
effect = "fixed"

coriobac = data.frame( model = "Coriobacteriaceae", effect = effect, coriobac, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
coriobac = coriobac[ , -c( 8 )]
coriobac = as_tibble( coriobac )

# Eubacteriaceae 
eubact = broom.mixed::tidy( lme.slope.Eubacteriaceae, conf.int = TRUE )
eubact = eubact[ -c( 11:12 ), -c( 7:8 )]
model = "Eubacteriaceae"
eubact = data.frame( model = model, eubact )
eubact = as_tibble( eubact )

# Firmicutesfam.incertaesedis
firmifam = broom.mixed::tidy( lme.int.Firmicutesfam.incertaesedis, conf.int = TRUE )
firmifam = firmifam[ -c( 11:12 ), -c( 7:8 )]
model = "Firmicutesfam.incertaesedis"
firmifam = data.frame( model = model, firmifam )
firmifam = as_tibble( firmifam )

# Lachnospiraceae
lachno = broom.mixed::tidy( lme.int.Lachnospiraceae, conf.int = TRUE )
lachno = lachno[ -c( 11:12 ), -c( 7:8 )]
model = "Lachnospiraceae"
lachno = data.frame( model = model, lachno )
lachno = as_tibble( lachno )

# Oscillospiraceae
oscillo = broom.mixed::tidy( lme.int.Oscillospiraceae, conf.int = TRUE )
oscillo = oscillo[ -c( 11:12 ), -c( 6,8 )]
model = "Oscillospiraceae"
oscillo = data.frame( model = model, oscillo )
oscillo = as_tibble( oscillo )

# Prevotellaceae
prevo = broom.mixed::tidy( lme.int.Prevotellaceae, conf.int = TRUE )
prevo = prevo[ -c( 11:12 ), -c( 7:8 )]
model = "Prevotellaceae"
prevo = data.frame( model = model, prevo )
prevo = as_tibble( prevo )

# Rikenellaceae
conf.int = data.frame( confint( glmm.int.Rikenellaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:16 ), ]

rikenell = broom.mixed::tidy( glmm.int.Rikenellaceae )[ , -c( 1 )]
rikenell = rikenell[ -c( 11:18 ), ]
rikenell = rikenell[ , -c( 1 )]
effect = "fixed"

rikenell = data.frame( model = "Rikenellaceae", effect = effect, rikenell, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
rikenell = rikenell[ , -c( 8 )]
rikenell = as_tibble( rikenell )

# Ruminococcaceae
rumino = broom.mixed::tidy( lme.int.Ruminococcaceae, conf.int = TRUE )
rumino = rumino[ -c( 11:12 ), -c( 6,8 )]
model = "Ruminococcaceae"
rumino = data.frame( model = model, rumino )
rumino = as_tibble( rumino )

# Sutterellaceae
conf.int = data.frame( confint( glmm.int.Sutterellaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:16 ), ]

sutterell = broom.mixed::tidy( glmm.int.Sutterellaceae )[ , -c( 1 )]
sutterell = sutterell[ -c( 11:18 ), ]
sutterell = sutterell[ , -c( 1 )]
effect = "fixed"

sutterell = data.frame( model = "Sutterellaceae", effect = effect, sutterell, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
sutterell = sutterell[ , -c( 8 )]
sutterell = as_tibble( sutterell )

# Veillonellaceae
conf.int = data.frame( confint( glmm.slope.Veillonellaceae ))
colnames( conf.int )[ 1 ] = "conf.low"
colnames( conf.int )[ 2 ] = "conf.high"
conf.int = conf.int[ -c( 11:17 ), ]

veillon = broom.mixed::tidy( glmm.slope.Veillonellaceae )[ , -c( 1 )]
veillon = veillon[ -c( 11:18 ), ]
veillon = veillon[ , -c( 1 )]
effect = "fixed"

veillon = data.frame( model = "Veillonellaceae", effect = effect, veillon, conf.low = conf.int$conf.low, conf.high = conf.int$conf.high )
veillon = veillon[ , -c( 8 )]
veillon = as_tibble( veillon )
```

``` r
# combine
common <- intersect( colnames( prevo ), colnames( rumino ))
all_bacteria = rbind( bacteroid[ common ], bactfam[ common ], bifidobac[ common ], clostrid[ common ], clostridfam[ common ], coriobac[ common ], eubact[ common ], firmifam[ common ], lachno[ common ], oscillo[ common ], prevo[ common ], rikenell[ common ], rumino[ common ], sutterell[ common ], veillon[ common ] )

all_bacteria <- all_bacteria %>%
  filter(effect=='fixed' )

# pdf( "coefficientsplot.pdf", width = 15, height = 11 )
# dwplot( all_bacteria, dodge_size = 0.8, dot_args = list( size = 3 ),
#   whisker_args = list( size = 2 )) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw( base_size = 4 ) + theme( legend.title = element_blank(), text = element_text( size = 20 ))

dwplot( all_bacteria ) + xlab( "Coefficient Estimate" ) + geom_vline( xintercept = 0, colour = "grey60", linetype = 2 ) + theme_bw() + theme( legend.title = element_blank())
```

![](Top_15_Mixed_Models_files/figure-gfm/Coefficient%20Plot-1.png)<!-- -->

``` r
# dev.off()
```
