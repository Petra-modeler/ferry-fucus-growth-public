#---------

#---------

# IMPACTS OF FERRY TRAFFIC WAVES ON GROWTH OF FUCUS

#---------

#---------

getwd()
setwd("/Users/pesaari/Documents/PhD/R/Growth")
library(dplyr)

# Bring data in
data <- read.csv("private/data/growth.csv", stringsAsFactors = FALSE)

# Create new column 'tile', that combines SITE.ID and original.tile.number
data$tile <- paste(data$SITE.ID, data$Original.tile.number, sep = "_")

# Remove broken tips (SH2M & OK3O)
data <- data[!data$Fucus.ID.key %in% c("SH2M", "OK3O"), ]

data$Genotype.origin <- factor(data$Genotype.origin)
data$I.or.C <- factor(data$I.or.C)
data$Genotype <- factor(data$Genotype)

#-----------

# IMPACTS MEASURED IN BIOMASS

#-----------

######

# DIFFERENCES IN ORIGINAL BIOMASS
  # Does genotype origin affect biomass in June?  
  # Test if results tell about difference in starting biomass or differences in growth

######

library(glmmTMB)
m1 <- glmmTMB(
  Tip.dry.weight.June.ww_dw_ratio..g. ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=gaussian)

summary(m1)
# Tips collected from impact areas are 0.125 g smaller on average compared to tips from control areas in June
# The difference is highly significant

library(emmeans)
emmeans(m1, ~ Genotype.origin, type="response")
# 95 % CI (emmeans):
# Control: 0.307–0.385 g
# Impact: 0.182–0.260 g
# → no overlap

# Since the initial weight varies between the origin areas, growth change cannot be used in the models
# The initial weight must therefore be accounted for in the models

#######

# MODELS

#######

# How does the tip dry mass in July depend on the initial size in June, the area's impact status, and the origin of the genotype?

# 1. model
m2 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype) +
    (0+I.or.C|Genotype), 
  data=data, 
  family=gaussian)

# 2. model -> USE THIS
m3 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

# 3. model
m4 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile),
  data=data, 
  family=gaussian)

summary(m2)
summary(m3) # best AIC, no overfitting, same results as in full model
summary(m4) # probably the worst model

# 1. July biomass depends very strongly on initial size
# 2. No difference between I and C in July biomass when June biomass is the same; no main effect of treatment
# 3. In Impact areas, tips with larger initial dry weight in June tended to grow slightly less compared to those in Control areas
    # almost statistically significant (0.0505)
# 4. Whether a genotype was collected near or far from a shipping lane did not significantly affect tip growth
# 5. There is no evidence that genotype origin affects tip dry weight in July differently in impact and control sites

anova(m2, m3) # m3 has smaller AIC -> m3 better
anova(m3, m4) # m3 has smaller AIC -> m3 better 

emmeans(m3, ~ I.or.C, type="response") 


######

# RESIDUALS

######

pearson_resid <- residuals(m3, type = "pearson")
plot(pearson_resid)
hist(pearson_resid)

library(DHARMa)
resid1 <- simulateResiduals(m2, n=1000)
resid2<- simulateResiduals(m3, n=1000)
resid3 <- simulateResiduals(m4, n=1000)

plot(resid1)
# KS p: 0.00376
# Dispersion p: 0.888
# Outlier p: 0.00467
plot(resid2)
# KS p:0.00366
# Dispersion p: 0.906
# Outlier p: 6e-05
plot(resid3)
# KS p: 0.00283
# Dispersion p: 0.986
# Outlier p: 6e-04

######

# TRENDS

# How much the tip's dry weight in July increases per unit of the dry weight in June, 
# separately for Impact and Control areas

######

emtr <- emtrends(
  m3,
  ~ I.or.C,
  var = "Tip.dry.weight.June.ww_dw_ratio..g."
)

df_trends <- as.data.frame(emtr)

# In Control areas, July biomass increases by ~2.13 g for each 1 g of June biomass
# In Impact areas, the increase is ~1.86 g per 1 g of June biomass
# The slope indicates how much July dry weight increases per unit of June dry weight, separately for Control and Impact areas
# The slope is slightly lower in Impact areas (p = 0.05)

####

# TEST THE EFFECT OF RANDOM EFFECTS

# leave one out and compare models

#####

testm1 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)
  

testm2 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|SITE.ID) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

testm3 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile),
  data=data, 
  family=gaussian)

anova(m3, testm1) # taking site away from m3 doesn't make the model worse statistically
anova(m3, testm2) # having tile in m3 can have purpose in the model but not very strong evidence of that
anova(m3, testm3) # Genotype is important for the explanatory power of the model
# -> Some genotypes exhibit different average growth, independent of whether they are in Control or Impact areas.

#####

# Plot 1
# Treatment (Impact vs Control)

#####

library(ggeffects)
library(ggplot2)

# Change I.or.C to factor
data$I.or.C <- factor(data$I.or.C, levels = c("Control", "Impact"))
data$Genotype.origin <- factor(data$Genotype.origin, levels = c("Control", "Impact"))

# Predictions for treatment
pred_treatment <- ggpredict(
  m3,
  terms = c("I.or.C")
)

# Change to df
pred_treatment_df <- as.data.frame(pred_treatment)
# Copy x scolumn group and change to factor
pred_treatment_df$group <- factor(pred_treatment_df$x, levels = c("Control", "Impact"))

# Plot
ggplot(pred_treatment_df, aes(x = group, y = predicted, fill = group)) +
  geom_col(width = 0.5, color = "black", alpha=0.5) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  scale_fill_manual(values = c("Control" = "#00BFC4", "Impact" = "#F8766D"), guide = "none") +  # legendan poisto
  labs(
    x = "",
    y = "Predicted tip dry biomass in July (g)",
    subtitle = ""
  ) +
  theme_bw(base_size = 14)


#####

# Plot 2
# June - July effect in I and C

#####

library(dplyr)
library(ggplot2)

ggplot() +
  # Predictions
  geom_line(
    data = newdata,
    aes(
      x = Tip.dry.weight.June.ww_dw_ratio..g.,
      y = fit,
      color = I.or.C
    ),
    size = 1.4
  ) +
  
  # Confidence
  geom_ribbon(
    data = newdata,
    aes(
      x = Tip.dry.weight.June.ww_dw_ratio..g.,
      ymin = lower,
      ymax = upper,
      fill = I.or.C
    ),
    alpha = 0.25,
    color = NA,
    show.legend = FALSE
  ) +
  
  # July observations as points
  geom_point(
    data = data,
    aes(
      x = Tip.dry.weight.June.ww_dw_ratio..g.,
      y = Tip.dry.weight.July..g.,
      color = I.or.C,
      shape = I.or.C
    ),
    size = 2,
    alpha = 0.6
  ) +
  
  # Colors
  scale_color_manual(
    name = "Predicted biomass",
    values = c("Control" = "#00BFC4", "Impact" = "#F8766D")
  ) +
  scale_fill_manual(
    values = c("Control" = "#00BFC4", "Impact" = "#F8766D")
  ) +
  
  # Legend for points
  scale_shape_manual(
    name = "Observed biomass",
    values = c("Control" = 16, "Impact" = 16)
  ) +
  
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(linetype = 1, shape = NA)
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        color = c("#00BFC4", "#F8766D")
      )
    )
  ) +
  
  labs(
    x = "Biomass in June (g)",
    y = "Biomass in July (g)"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )



#-----------

# IMPACTS MEASURED IN LENGHT

#----------

library(glmmTMB)

# Does genotype origin affect length in June? 
m5 <- glmmTMB(
  June.tip.length.tot..cm. ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=gaussian)

summary(m5)

library(emmeans)
emmeans(m5, ~ Genotype.origin, type="response")
# There is no statistically significant difference in June tip lengths between Genotype.origin areas
# It is okay to use change variable here, we can do both


######

# MODELS

######

# June length and July length
m6 <- glmmTMB(
  July.tip.length.tot..cm. ~ 
    June.tip.length.tot..cm. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    June.tip.length.tot..cm.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

# Change
m7 <- glmmTMB(
  Change.tip.length.tot..cm. ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    June.tip.length.tot..cm.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

test <- glmmTMB(
  Change.tip.length.tot..cm. ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

summary(m6)
# July tip length depends primarily on June tip length, and not on area (I.or.C) or origin (Genotype.origin).

summary(m7)
# The change in growth (June → July) does not differ statistically between Impact and Control areas or between Genotype.origin groups.

summary(test)
# no difference to previous models


pearson_resid_m6 <- residuals(m6, type = "pearson")
plot(pearson_resid_m6)
hist(pearson_resid_m6)

library(DHARMa)
resid_m6 <- simulateResiduals(m6, n=1000)

plot(resid_m6) # -> good fit
# KS p: 0.97265
# Dispersion p: 0.786
# Outlier p: 1


# Check effect of random variables
m8 <- glmmTMB(
  Change.tip.length.tot..cm. ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    June.tip.length.tot..cm.:I.or.C + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

m9 <- glmmTMB(
  Change.tip.length.tot..cm. ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    June.tip.length.tot..cm.:I.or.C + 
    (1|SITE.ID) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

m10 <- glmmTMB(
  Change.tip.length.tot..cm. ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    June.tip.length.tot..cm.:I.or.C + 
    (1|SITE.ID) + 
    (1|tile),
  data=data, 
  family=gaussian)

anova(m8, m7) 
# SITE.ID explains a portion of the variation in growth change. 
# This indicates that different study sites have a small but significant effect on tip growth rates. 
# Individual genotypes (Genotype) still account for the largest part of the variation
anova(m9, m7) 
# Tile does not explain a significant portion of the variation in growth change.
anova(m10, m7) 
# Differences among individual genotypes explain a large portion of the variation in growth change.


#-----------

# IMPACTS MEASURED AS NUMBER OF TIPS

#-----------
# Bring data in
data <- read.csv("private/data/growth.csv", stringsAsFactors = FALSE)
data_june <- read.csv("private/data/data_june.csv", stringsAsFactors = FALSE)

colnames(data)[colnames(data) == "Apikal.tips..amount."] <- "Apikal_tips_amount_July"

# Create new column 'tile', that combines SITE.ID and original.tile.number
data$tile <- paste(data$SITE.ID, data$Original.tile.number, sep = "_")

# Remove broken tips (SH2M & OK3O)
data <- data[!data$Fucus.ID.key %in% c("SH2M", "OK3O"), ]

# Combine June tip count to dataframe
data <- data %>%
  left_join(
    data_june %>% select(site_tile_genotype, apical_tip_number_June),
    by = c("Fucus.ID.key" = "site_tile_genotype")
  )

data$Genotype.origin <- factor(data$Genotype.origin)
data$I.or.C <- factor(data$I.or.C)
data$Genotype <- factor(data$Genotype)

# Check the distribution
hist(data$apical_tip_number_June)

# Does genotype origin affect tip count in June? 

# Poisson
m11 <- glmmTMB(
  apical_tip_number_June ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=poisson)

# Negative binomial 1
m12 <- glmmTMB(
  apical_tip_number_June ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=nbinom1)

# Negative binomial 2
m13 <- glmmTMB(
  apical_tip_number_June ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=nbinom2)


summary(m11)
summary(m12)
summary(m13)
# Negbin1 best and genotype origin doesn't have impact on apical tip number in June
# I could model change in apical tip count or do as previously and use July number as response variable

########

# MODELS

########

m14 <- glmmTMB(
  Apikal_tips_amount_July ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=nbinom1)

m15 <- glmmTMB(
  Apikal_tips_amount_July ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=nbinom2)

m16 <- glmmTMB(
  Apikal_tips_amount_July ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=poisson)

summary(m14)
summary(m15)
summary(m16)

anova(m14, m15)
# nbinom1 (m14) and nbinom2 (m15) fit the data similarly, with slightly lower AIC for nbinom2
# extremely high dispersion estimate in nbinom2 (111) is unrealistic, while nbinom1 provides a reasonable estimate
# nbinom1 is the best
# tip count in June is the only one explaining variation in fixed terms

anova(m14, m16)
# both models fit the data similarly


pearson_resid_m14 <- residuals(m14, type = "pearson")
pearson_resid_m16 <- residuals(m16, type = "pearson")
plot(pearson_resid_m14)
plot(pearson_resid_m16)
hist(pearson_resid_m14)
hist(pearson_resid_m16)

library(DHARMa)
resid_m14 <- simulateResiduals(m14, n=1000)
resid_m16<- simulateResiduals(m16, n=1000)

plot(resid_m14)
# KS p: 0.2469
# Dispersion p: 0.968
# Outlier p: 0.28
plot(resid_m16)
# KS p:0.40255
# Dispersion p: 0.888
# Outlier p: 0.42

# model assumptions are reasonably met in both models
# both have good fit

#######

# RANDOM TERMS

######

summary(m16)$varcor


#######

# CHANGE IN TIP COUNT

######

# Test model with change in tip count as response variable
data$change_in_tip_number <- data$Apikal_tips_amount_July - data$apical_tip_number_June

hist(data$change_in_tip_number)

m17 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=nbinom2)

m18 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=nbinom1)

m19 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=poisson)

m20 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data, 
  family=gaussian)

summary(m20)
pearson_resid_m20 <- residuals(m20, type = "pearson")

plot(pearson_resid_m20)

hist(pearson_resid_m20)


library(DHARMa)
resid_m20 <- simulateResiduals(m20, n=1000)


plot(resid_m20)
# KS p: 0.0284
# Dispersion p: 0.848
# Outlier p: 7e-05

#######

# ONLY POSITIVE CHANGE

######

# There is many observations where change was negative, which is not biologically reasonable
data_positive <- data[data$change_in_tip_number >= 0, ]
hist(data_positive$change_in_tip_number)

m21 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data_positive, 
  family=nbinom2)

m22 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data_positive, 
  family=nbinom1)

m23 <- glmmTMB(
  change_in_tip_number ~ 
    apical_tip_number_June + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    apical_tip_number_June:I.or.C + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype),
  data=data_positive, 
  family=poisson)

summary(m21)
summary(m22)
summary(m23)

# ALL PREVIOUS IN DISTANCE FROM FAIRWAY AND WITH WAVES
