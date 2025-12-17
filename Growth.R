# TO DO: 
# tuota viivamalli luotettavuusrajoilla ja alkuperäisella datalla (y=loppukoko, x=alkukoko)
# pudota geno-aalto interactio pois kun tekee kuvaajan
# tee kuvaaja jossa genotyypit erikseen viivoina ja näkee miten eri genottypit kasvaa (genotyyppiestimaatit mallista)
# 3 kuvaa: aallokkoisuus päävaikutus (emmeansin luvut), alkukoko-loppoukoko alltovaikutus, genotyyppivaihtelu kolmas
# sama harjoitus pituudelle (change voi käyttää jos ei eroa lähtötilanteessa, eli testaa onko i vs c alkutilanne eri) ja tipsien määrälle
# biomassassa alkutilanne vahvemmin vaikuttaa muutokseen kuin pituuden muutoksessa koska jos tipsejä on 
# enemmän niin biomassa muuttuu enemmän kuin jos tipsejä on vähemmän, vaikka pituuskasvu olisi sama 
# (eli alkumassa on massa-analyysissä merkittävämpi vaikutus tod. näk)
# tipseissa count data, poisson jakauma varmasti

getwd()
setwd("/Users/pesaari/Documents/PhD/R/Growth")
library(dplyr)

# Bring data in
data <- read.csv("dry_weight_plot_numbers_corrected.csv", stringsAsFactors = FALSE)
# Luo uusi sarake 'tile', joka yhdistää SITE.ID ja original.tile.number

# Remove one broken tips SH2N, SH2M tai OK3O (still some left that are below -0.2)
data <- data[!data$Fucus.ID.key %in% c("SH2N", "SH2M", "OK3O"), ]
data$tile <- paste(data$SITE.ID, data$Original.tile.number, sep = "_")
# data_July_noNA <- data[!is.na(data$Tip.dry.weight.July..g.), ]
# Change tule number to factor
#data$Tile.number <- as.factor(data$Tile.number)

# vaikuttaako genotype origin alkukokoon (että kyse ei välttämättä ole kasvun eroissa vaan lähtötilanteen erossa)
m_start_mass <- lmer(
  Tip.dry.weight.June.ww_dw_ratio..g. ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data)

summary(m_start_mass) # Genotype.originImpact -0.12489    0.02949  -4.235 -> jos lisätään 0.6 -0.12489:n niin ei mene nollaan niin on merkitsevä
# jos poikkeaa nollasta niin eroaa toisistaan

emmeans(m_start_mass, ~ Genotype.origin, type="response") # Piirrä tästä kuva kun mallin sopivuus on testattu
# kontrollilta kerätyt painavampia
# eli koska tämä on näin niin ei käyttetä change muuttujaa vaan alkupäinoa ja loppupainoa

m_start_mass2 <- glmmTMB(
  Tip.dry.weight.June.ww_dw_ratio..g. ~ 
    Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (1|Genotype), 
  data=data,
  family=gaussian)
# merkitsevä



summary(m_start_mass2)
# lmer4 voisi kokeilla

# Fit model (whole model and all variables) -> LOPULLINEN MALLI!!!
library(glmmTMB)
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


#####

#-------

#####
# -> LOPULLINEN MALLI!!!
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

####

#-----

#####

# testaa mallien ero
anova(m_biomass_2, m_biomass_test_yhdysvaikutus_puuttuu) # jälkimmäisellä pienempi AIC niin parempi mutta merkitsevä

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

anova(m_biomass_test_yhdysvaikutus_puuttuu, m_biomass_test_yhdysvaikutus_ja_geno_puuttuu) # genotyyppi pitää olla mutta yhdysvaikutusta ei tarvita, slopet ei poikkea, genotyypeillä eri kasvu mutta aaltovaikutus ei vaikuta


summary(m_biomass_2)
summary(m_biomass_test_yhdysvaikutus_puuttuu) # paras AIC, ei ylisovittava, sama johtopäätös kuin täydessä mallissa
summary(m_biomass_test_yhdysvaikutus_ja_geno_puuttuu) #kaikista huonoin ehkä

library(emmeans)
emmeans(m_biomass_2, ~ I.or.C, type="response") # Piirrä tästä kuva kun mallin sopivuus on testattu

library(DHARMa)
pearson_resid_m_biomass_2 <- residuals(m_biomass_2, type = "pearson")
plot(pearson_resid_m_biomass_2)
hist(pearson_resid_m_biomass_2) # ei ihan normaali mutta ehkä hyväksyttävä
resid_m_biomass_2 <- simulateResiduals(m_biomass_2, n=1000)
plot(resid_m_biomass_2) # katso mitä tarkoittaa KS test ym.


m_biomass_3 <- glmmTMB(
  Tip.dry.weight.change ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (I.or.C|Genotype), 
  data=data_mass_change_noNA, 
  family=gaussian)

summary(m_biomass_3)
diagnose(m_biomass_3)

data_mass_change_noNA <- data[!is.na(data$Tip.dry.weight.change), ]
#sum(is.na(data_mass_change_noNA$Tip.dry.weight.change))

library(lme4)
m_biomass_3 <- lmer(
  Tip.dry.weight.change ~ 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    (1|SITE.ID) + 
    (1|tile) + 
    (I.or.C|Genotype), 
  data=data_mass_change_noNA)

summary(m_biomass_3)
pearson_resid_m_biomass_3 <- residuals(m_biomass_3, type = "pearson")
plot(pearson_resid_m_biomass_3)
hist(pearson_resid_m_biomass_3) # ei ihan normaali mutta ehkä hyväksyttävä
resid_m_biomass_3 <- simulateResiduals(m_biomass_3, n=1000)
plot(resid_m_biomass_3) # katso mitä tarkoittaa KS test ym., näyttää hyvältä



emmeans(m_biomass_3, ~ I.or.C, type="response") 
emmeans(m_biomass_3, ~ Genotype.origin, type="response") # jos toisen keskiarvo menee toisen confidence limitien väliin niin ei ole merkitsevä
emmeans(m_biomass_3, ~  I.or.C:Genotype.origin, type="response") # 

# genottypin alkuperä vaikuttaa kasvuun change-muuttujalla


#####
######
#######




#diagnose(m_biomass_2)

summary(m_biomass_2)$varcor
# (I.or.C | Genotype):n standardipoikkeama on 0.0008, eli käytännössä nolla. 
# Lisäksi korrelaatio interceptin kanssa on lähes 1 (0.995), 
# mikä viittaa siihen, että slopea ei pystytty kunnolla erottamaan interceptistä.
# Tämä on tyypillistä merkkiä siitä, että satunnaisslope ei tunnistu datasta, 
# eli data ei tue sen estimointia. Se voi aiheuttaa epäluotettavia p-arvoja ja epävarmoja estimaatteja.

# Multicollinearity can't be a problem as I have only one numeric variable

# Laske havaintojen määrä jokaisessa yhdistelmässä
table(data$I.or.C, data$Genotype.origin) # not a problem as there are many samples in each combination

# Removed tile
m_biomass_3 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. + 
    I.or.C + 
    Genotype.origin + 
    I.or.C:Genotype.origin + 
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C + 
    (1|SITE.ID) + 
    (I.or.C|Genotype), 
  data=data, 
  family=gaussian)
# Not convergence problem
diagnose(m_biomass_3)


# Can we keep tile, if we take random slope away (only intercept for genotype)
m_simple1 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. +
    I.or.C +
    Genotype.origin +
    I.or.C:Genotype.origin +
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C +
    (1 | SITE.ID) +
    (1 | Tile.number) +
    (1 | Genotype),
  data = data,
  family = gaussian
) # Not convergence problem

diagnose(m_simple1)
# Satunnaisslope (I.or.C|Genotype) ei tunnistunut, joten poistamalla sen mallin rakenne ei muuttunut merkittävästi.
# diagnose() antaa edelleen saman varoituksen suurista Z-statistiikoista, koska ongelma ei ole vain satunnaisslopessa vaan myös 
# joissain kiinteissä efekti-parametreissa (esim. Tip.dry.weight.June.ww_dw_ratio..g.)
# ja pienissä hajonnoissa (theta_1|Genotype.1).

# Might be scaling problem in numeric variable
summary(data$Tip.dry.weight.June.ww_dw_ratio..g.) # not a scaling problem

# Model without Genotype.origin
m_test1 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. +
    I.or.C +
    I.or.C:Genotype.origin +
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C +
    (1 | SITE.ID) +
    (1 | Tile.number) +
    (1 | Genotype),
  data = data,
  family = gaussian
)

# Compare to your original model
anova(m_simple1, m_test1)
# The likelihood ratio test (anova(m_simple1, m_test1)) shows a Chi-squared = 0 and p-value = 1.
# Removing Genotype.origin from the model does not reduce the fit at all
# This suggests Genotype.origin is likely not needed in the model — it could be safely removed, 
# which might help stabilize the parameter estimates and make the Z-statistics more reliable.

# 1. Remove I.or.C
m_test_IorC <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. +
    Genotype.origin +
    I.or.C:Genotype.origin +
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C +
    (1 | SITE.ID) +
    (1 | Tile.number) +
    (1 | Genotype),
  data = data,
  family = gaussian
)
anova(m_simple1, m_test_IorC)

# This strongly suggests that the main effect I.or.C is not contributing independently 
# to explaining variation in July tip weight — the information is likely captured by interactions 
# like I.or.C:Genotype.origin

# 2. Remove the interaction I.or.C:Genotype.origin
m_test_interaction1 <- glmmTMB(
  Tip.dry.weight.July..g. ~ 
    Tip.dry.weight.June.ww_dw_ratio..g. +
    I.or.C +
    Genotype.origin +
    Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C +
    (1 | SITE.ID) +
    (1 | Tile.number) +
    (1 | Genotype),
  data = data,
  family = gaussian
)
anova(m_simple1, m_test_interaction1)
# That means removing the I.or.C:Genotype.origin interaction does not significantly worsen model fit.
# You can consider simplifying the model by dropping I.or.C:Genotype.origin, 
# while keeping the main effects and the other interaction (Tip.dry.weight.June.ww_dw_ratio..g.:I.or.C)






#####
####
##

# RESIDUALS

# Extract residuals
pearson_resid_2 <- residuals(m_biomass_2, type = "pearson")
plot(pearson_resid_2)

# Save the data model used in a data frame
model_data_2 <- model.frame(m_biomass_2)

# Combine residuals with the data model used
model_data_2$id <- 1:nrow(model_data_2)   # yksilöivä tunniste
model_data_2$pearson_resid_2 <- residuals(m_biomass_2, type = "pearson")

# Plot residuals with codes
plot(model_data_2$pearson_resid_2, pch = 19)
text(x = 1:nrow(model_data_2),
     y = model_data_2$pearson_resid_2,
     labels = model_data_2$id,
     pos = 3, cex = 0.7)

# Check residuals that have value under -0.2 from pictures, if they can be removed (broken)
model_data_2[model_data_2$pearson_resid_2 < -0.2, ]
# FL2A: onko irrallinen pala otettu kuivamiseen mukaan?
# KK2A: jotain on revennyt, mutta ei voi sanoa puuttuuko palaa
# KK2B: ehkä vähän rikki alhaalta
# AL1B: vaikea hahmottaa kuvaa
# FL3G: ei mitään rikkinäisen näköistä
# SH2M: paino pienentynyt ja pala puuttuu -> poista taulukosta
# OK3O: selkeästi puuttuu pala -> poista taulukosta




