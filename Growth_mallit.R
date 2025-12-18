#---------

#---------

# GROWTH

#---------

#---------

getwd()
setwd("/Users/pesaari/Documents/PhD/R/Growth")
library(dplyr)

# Bring data in
data <- read.csv("data/growth.csv", stringsAsFactors = FALSE)

# Create new column 'tile', that combines SITE.ID and original.tile.number
data$tile <- paste(data$SITE.ID, data$Original.tile.number, sep = "_")

# Remove broken tips (SH2M & OK3O)
data <- data[!data$Fucus.ID.key %in% c("SH2M", "OK3O"), ]

data$Genotype.origin <- factor(data$Genotype.origin)
data$I.or.C <- factor(data$I.or.C)
data$Genotype <- factor(data$Genotype)

#-----------

# BIOMASS

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


# I or C Ei eroa July biomassissa, kun June biomass on sama, käsittelyn päävaikutusta ei ole
# Geneettinen alkuperä ei selitä kasvua, kun lähtökoko huomioidaan
# Ei näyttöä geneettisen alkuperän ja käsittelyn yhteisvaikutuksesta kasvussa
# viitteitä siitä, että Impact-alueilla kasvu skaalautuu hieman heikommin lähtökoon kanssa

# 1. Growth depends very strongly on initial size
# 2. In Impact areas, tips with larger initial dry weight in June tended to grow slightly less by July compared to those in Control areas
    # almost statistically significant (0.0505)
# 3. No difference between I and C in July biomass when June biomass is the same; no main effect of treatment
# 4. Whether a genotype was collected near or far from a shipping lane did not significantly affect tip growth once initial size and site treatment were accounted for.
# 5. There is no evidence that the effect of genotype origin differs between areas.

emmeans(m3, ~ I.or.C, type="response") # Piirrä tästä kuva kun mallin sopivuus on testattu
# interaktiodien takia tulokset voi olla misleading

######

# residuals

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

# Trendit mallista
emtr <- emtrends(
  m3,
  ~ I.or.C,
  var = "Tip.dry.weight.June.ww_dw_ratio..g."
)

df_trends <- as.data.frame(emtr)

# Control-alueilla July biomass kasvaa ~2.13 g jokaista 1 g June biomassaa kohden
# Impact-alueilla kasvu on ~1.86 g per 1 g June biomassaa
# Kasvunopeus (slope) on pienempi Impact-alueilla, p=0.05

####

# testaa satunnaismuuttujien vaikutusta

#####



#####

# Plot 1
# Aallokkoisuus päävaikutus (Impact vs Control)

#####

library(ggeffects)
library(ggplot2)

# Muutetaan I.or.C faktoriksi
data$I.or.C <- factor(data$I.or.C, levels = c("Control", "Impact"))
data$Genotype.origin <- factor(data$Genotype.origin, levels = c("Control", "Impact"))

# Marginaaliset ennusteet pelkälle treatmentille (I.or.C)
pred_treatment <- ggpredict(
  m3,
  terms = c("I.or.C")
)

# Muutetaan data.frameksi ggplotille
pred_treatment_df <- as.data.frame(pred_treatment)
# Kopioidaan x sarakkeesta group ja muutetaan faktoriksi
pred_treatment_df$group <- factor(pred_treatment_df$x, levels = c("Control", "Impact"))

# Piirretään
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
# Alkukoko-loppukoko vaikutus (aaltovaikutus, x=June, y=July)

#####

library(dplyr)
library(ggplot2)

ggplot() +
  # Ennustetut viivat
  geom_line(
    data = newdata,
    aes(
      x = Tip.dry.weight.June.ww_dw_ratio..g.,
      y = fit,
      color = I.or.C
    ),
    size = 1.4
  ) +
  
  # Luottamusvälit (ei legendaan)
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
  
  # Havainnot pisteinä
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
  
  # Värit viivoille ja pisteille
  scale_color_manual(
    name = "Predicted biomass",
    values = c("Control" = "#00BFC4", "Impact" = "#F8766D")
  ) +
  
  # Täytöt ribbonille
  scale_fill_manual(
    values = c("Control" = "#00BFC4", "Impact" = "#F8766D")
  ) +
  
  # Pisteiden legenda
  scale_shape_manual(
    name = "Observed biomass",
    values = c("Control" = 16, "Impact" = 16)
  ) +
  
  # 🔴 KORJAUS TÄSSÄ
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


#####

# Plot 3
# Genotyyppivaihtelu - pylväs

######

mean_June <- mean(data$Tip.dry.weight.June.ww_dw_ratio..g., na.rm = TRUE)

newdata3 <- expand.grid(
  Genotype = levels(data$Genotype),
  I.or.C = levels(data$I.or.C),
  Genotype.origin = levels(data$Genotype.origin)[1],  # referenssitaso
  Tip.dry.weight.June.ww_dw_ratio..g. = mean_June
)

newdata3$predicted <- predict(
  m3,
  newdata = newdata3,
  re.form = NA
)



library(ggplot2)

# Piirretään pylväät
ggplot(newdata3, aes(x = Genotype, y = predicted, fill = I.or.C)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.7) +
  labs(
    x = "Genotype",
    y = "Predicted tip dry mass in July (g)",
    fill = "Treatment"
  ) +
  theme_bw(base_size = 14)

# Halutessasi voit lisätä myös alkuperäiset datapisteet päälle
ggplot(newdata, aes(x = Genotype, y = predicted, fill = I.or.C)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.7) +
  geom_jitter(data = data,
              aes(x = Genotype, y = Tip.dry.weight.July..g., color = I.or.C),
              width = 0.2, size = 1.5, alpha = 0.7,
              inherit.aes = FALSE) +
  labs(
    x = "Genotype",
    y = "Predicted tip dry mass in July (g)",
    fill = "Treatment",
    color = "Treatment"
  ) +
  theme_bw(base_size = 14)





library(ggeffects)
library(ggplot2)

data$Genotype <- factor(data$Genotype)
data$I.or.C <- factor(data$I.or.C)
# Ennusteet Genotype ja käsittely (I.or.C) mukaan
pred_genotype_treatment <- ggpredict(
  m3,
  terms = c("Genotype", "I.or.C"),
  type = "fixed"  # ei huomioida satunnaistekijät
)

# Muutetaan data.frameksi ggplotille
pred_df <- as.data.frame(pred_genotype_treatment)

# Piirretään pylväskaavio käsittelyllä värikoodattuna
ggplot(pred_df, aes(x = x, y = predicted, fill = group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.7) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(
    x = "Genotype",
    y = "Predicted tip biomass in July (g)",
    fill = "Treatment"
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))




ggplot() +
  # pylväät ennusteista
  geom_col(data = pred_df,
           aes(x = x, y = predicted, fill = group),
           position = position_dodge(width = 0.8),
           width = 0.7, color = "black", alpha = 0.7) +
  geom_errorbar(data = pred_df,
                aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  # alkuperäiset datapisteet
  geom_jitter(data = data,
              aes(x = Genotype, y = Tip.dry.weight.July..g., color = I.or.C),
              width = 0.2, size = 1.5, alpha = 0.7,
              inherit.aes = FALSE) +
  labs(x = "Genotype", y = "Tip dry mass in July (g)", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14)




library(dplyr)

summary_df <- data %>%
  group_by(Genotype, I.or.C) %>%
  summarise(
    mean_July = mean(Tip.dry.weight.July..g.),
    sd_July = sd(Tip.dry.weight.July..g.),
    n = n(),
    se = sd_July / sqrt(n)
  )

ggplot(summary_df, aes(x = Genotype, y = mean_July, fill = I.or.C)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean_July - se, ymax = mean_July + se),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(x = "Genotype", y = "Mean tip dry mass in July (g)", fill = "Treatment") +
  theme_bw(base_size = 14)



library(ggeffects)
library(ggplot2)

data$Genotype <- factor(data$Genotype)

# sovita malli uudelleen
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

# Ennusteet Genotype-satunnaistekijän mukaan
pred_genotype <- ggpredict(
  m3,
  terms = c("Genotype"),
  type = "random"# huomioidaan satunnaistekijät
)

# Muutetaan data.frameksi ggplotille
pred_genotype_df <- as.data.frame(pred_genotype)

# Piirretään ennusteet
ggplot(pred_genotype_df, aes(x = x, y = predicted)) +
  geom_col(width = 0.7, fill = "lightblue", color = "black", alpha = 0.6) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(
    x = "Genotype",
    y = "Predicted tip biomass in July (g)",
  ) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

#####

# Genotyyppivaihtel - Viiva

######

library(ggplot2)
library(ggeffects)

# Ennustetaan genotyyppikohtaiset arvot huomioiden mallin satunnaistekijät
pred_genotype <- ggpredict(
  m3,
  terms = c("Tip.dry.weight.June.ww_dw_ratio..g.", "Genotype"),
  type = "random"  # huomioidaan satunnaistekijät
)

# Muutetaan data.frameksi ggplotille
pred_genotype_df <- as.data.frame(pred_genotype)

# Piirretään viivakuva genotyypeistä
ggplot(pred_genotype_df, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  labs(
    x = "Tip dry weight in June (g)",
    y = "Predicted tip dry weight in July (g)",
    color = "Genotype"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )


#-----------

# LENGHT

#----------

# testaa vaikuttaako genotype origin alkupituuteen

#-----------

# AMOUNT OF TIPS

#----------

# testaa vaikuttaako genotype origin tipsien määrään



