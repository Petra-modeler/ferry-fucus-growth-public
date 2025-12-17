getwd()
setwd("/Users/pesaari/Documents/PhD/R/Growth")
library(dplyr)

# Bring data in
data <- read.csv("dry_weight_plot_numbers_corrected.csv", stringsAsFactors = FALSE)

# Luo uusi sarake 'tile', joka yhdistää SITE.ID ja original.tile.number
data$tile <- paste(data$SITE.ID, data$Original.tile.number, sep = "_")

# Remove broken tips
  # SH2M: paino pienentynyt ja pala puuttuu -> poista taulukosta
  # OK3O: selkeästi puuttuu pala -> poista taulukosta
data <- data[!data$Fucus.ID.key %in% c("SH2M", "OK3O"), ]


#-----------

# BIOMASS

#----------

######

# vaikuttaako genotype origin alkubiomassaan (että kyse ei välttämättä ole kasvun eroissa vaan lähtötilanteen erossa)

######

library(glmmTMB)
m_start_mass2 <- glmmTMB(
  Tip.dry.weight.June.ww_dw_ratio..g. ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=gaussian)

summary(m_start_mass2)

library(emmeans)
emmeans(m_start_mass2, ~ Genotype.origin, type="response") # Piirrä tästä kuva kun mallin sopivuus on testattu

# Impact-genotyyppien alkukoko kesäkuussa on keskimäärin 0.125 g pienempi
# kuin Control-genotyyppien
# Ero on erittäin merkitsevä
# 95 % CI (emmeans):
# Control: 0.307–0.385 g
# Impact: 0.182–0.260 g
# → ei päällekkäisyyttä

# ei voi käyttää change muuttujaa eli kasvun määrää vaan pitää huomioida alkupaino
# mutta haittaako se koska impactilta kerätyt laitettiin sekä impactille että controllille? -> ei minusta

#######

# Mallit

#######

# Miten tipin kuiva massa heinäkuussa riippuu alkukoosta kesäkuussa, alueen impact-statuksesta ja genotyypin alkuperästä?

# 1. 
m_biomass_2 <- glmmTMB(
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

# 2. USE THIS
m_biomass_test_yhdysvaikutus_puuttuu <- glmmTMB(
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

#3. 
m_biomass_test_yhdysvaikutus_ja_geno_puuttuu <- glmmTMB(
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

summary(m_biomass_2)
summary(m_biomass_test_yhdysvaikutus_puuttuu) # paras AIC, ei ylisovittava, sama johtopäätös kuin täydessä mallissa
summary(m_biomass_test_yhdysvaikutus_ja_geno_puuttuu) #kaikista huonoin ehkä

# Kasvu riippuu erittäin voimakkaasti lähtökoosta
# Impact-yksilöiden kasvu lisääntyy hitaammin lähtökoon kasvaessa ja vaikutus on tilastollisesti lähes merkitsevä
# I or C Ei eroa July biomassissa, kun June biomass on sama, käsittelyn päävaikutusta ei ole
# Geneettinen alkuperä ei selitä kasvua, kun lähtökoko huomioidaan
# Ei näyttöä geneettisen alkuperän ja käsittelyn yhteisvaikutuksesta kasvussa
# viitteitä siitä, että Impact-alueilla kasvu skaalautuu hieman heikommin lähtökoon kanssa

emmeans(m_biomass_test_yhdysvaikutus_puuttuu, ~ I.or.C, type="response") # Piirrä tästä kuva kun mallin sopivuus on testattu
# interaktiodien takia tulokset voi olla misleading

######

# residuals

######

pearson_resid_1 <- residuals(m_biomass_test_yhdysvaikutus_puuttuu, type = "pearson")
plot(pearson_resid_2)
hist(pearson_resid_2)

library(DHARMa)
resid1 <- simulateResiduals(m_biomass_2, n=1000)
resid2<- simulateResiduals(m_biomass_test_yhdysvaikutus_puuttuu, n=1000)
resid3 <- simulateResiduals(m_biomass_test_yhdysvaikutus_ja_geno_puuttuu, n=1000)

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
  m_biomass_test_yhdysvaikutus_puuttuu,
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
  m_biomass_test_yhdysvaikutus_puuttuu,
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
  scale_fill_manual(values = c("Control" = "blue", "Impact" = "red"), guide = "none") +  # legendan poisto
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
    values = c("Control" = "blue", "Impact" = "red")
  ) +
  
  # Täytöt ribbonille
  scale_fill_manual(
    values = c("Control" = "blue", "Impact" = "red")
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
        color = c("blue", "red")
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
library(ggeffects)
library(ggplot2)

data$Genotype <- factor(data$Genotype)

# sovita malli uudelleen
m_biomass_test_yhdysvaikutus_puuttuu <- glmmTMB(
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
  m_biomass_test_yhdysvaikutus_puuttuu,
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
  m_biomass_test_yhdysvaikutus_puuttuu,
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



