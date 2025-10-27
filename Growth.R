getwd()
setwd("/Users/pesaari/Documents/PhD/R/Growth")
library(dplyr)

df <- read.csv("Growth_siistitty.csv", stringsAsFactors = FALSE)

#######

# Investigate growth/day and wet mass/day

#######

# Draw box plots
library(ggplot2)

# Tip length
ggplot(df, aes(x = I.or.C, y = Change.in.tip.length.day..cm.)) +
  geom_boxplot(fill = "gray80", color = "black", width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +  # add data points
  labs(
    title = "",
    x = "",
    y = "Growth rate (cm/day)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Tip weight
ggplot(df, aes(x = I.or.C, y = Change.in.tip.weigth.day..g.)) +
  geom_boxplot(fill = "gray80", color = "black", width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  labs(
    title = "",
    x = "",
    y = "Growth rate (g/day)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Draw histograms 
ggplot(df, aes(x = Change.in.tip.length.day..cm.)) +
  geom_histogram(fill = "gray60", bins = 20, color = "black") +
  facet_wrap(~ I.or.C) +
  labs(
    title = "",
    x = "Growth rate (cm/day)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggplot(df, aes(x = Change.in.tip.weigth.day..g.)) +
  geom_histogram(fill = "gray60", bins = 20, color = "black") +
  facet_wrap(~ I.or.C) +
  labs(
    title = "",
    x = "Growth rate (g/day)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Data contains negative growth so I remove those from the data
df_clean <- df[
  df$Change.in.tip.length.day..cm. >= 0 & 
    df$Change.in.tip.weigth.day..g. >= 0, 
]

# Tip length
ggplot(df_clean, aes(x = I.or.C, y = Change.in.tip.length.day..cm.)) +
  geom_boxplot(fill = "gray80", color = "black", width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +  # add data points
  labs(
    title = "",
    x = "",
    y = "Growth rate (cm/day)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Tip weight
ggplot(df_clean, aes(x = I.or.C, y = Change.in.tip.weigth.day..g.)) +
  geom_boxplot(fill = "gray80", color = "black", width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.6, size = 2) +
  labs(
    title = "",
    x = "",
    y = "Growth rate (g/day)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Draw histograms 
ggplot(df_clean, aes(x = Change.in.tip.length.day..cm.)) +
  geom_histogram(fill = "gray60", bins = 20, color = "black") +
  facet_wrap(~ I.or.C) +
  labs(
    title = "",
    x = "Growth rate (cm/day)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggplot(df_clean, aes(x = Change.in.tip.weigth.day..g.)) +
  geom_histogram(fill = "gray60", bins = 20, color = "black") +
  facet_wrap(~ I.or.C) +
  labs(
    title = "",
    x = "Growth rate (g/day)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Check data normality
shapiro.test(df_clean$Change.in.tip.length.day..cm.[df$I.or.C == "Control"]) # Normally distributed
shapiro.test(df_clean$Change.in.tip.length.day..cm.[df$I.or.C == "Impact"]) # Normally distributed

shapiro.test(df_clean$Change.in.tip.weigth.day..g.[df$I.or.C == "Control"]) # Non-normal distribution
shapiro.test(df_clean$Change.in.tip.weigth.day..g.[df$I.or.C == "Impact"]) # Non-normal distribution

######

# Test if impact sites differ from control sites significantly

#####

# Independent samples t-test
t.test(Change.in.tip.length.day..cm. ~ I.or.C, data = df_clean)
# Wilcoxon rank-sum test
wilcox.test(Change.in.tip.weigth.day..g. ~ I.or.C, data = df_clean)


#######

# Genotypes

#######

# Does the growth rate differ between different genotypes
df2 <- read.csv("Genotyyppitarkasteluun.csv", stringsAsFactors = FALSE)

# Remove negative values
df2_clean <- df2[df2$`Change.tip.length.tot..cm..day` >= 0, ]

# Make SITE.ID, I vs C and Genotype as factors instead of characters
df2_clean$Genotype <- as.factor(df2_clean$Genotype)
df2_clean$I.or.C <- as.factor(df2_clean$I.or.C)


# Two way ANOVA to test if genotype, I vs C or their interaction has an effect on growth.
model <- aov(Change.tip.length.tot..cm..day ~ Genotype * I.or.C, data = df2_clean)
summary(model)

# The analysis showed that genotype (p < 0.001) and environment (I or C and site) (p < 0.001) had a significant effect on growth. 
# However, there was no significant interaction between genotype and environment (p = 0.509 and p = 0.591), 
# suggesting that genotypes responded similarly across different environmental conditions.

# Check assumptions
plot(model, which = 2)  # QQ-plot
plot(model, which = 1)  # Residuals vs Fitted


# Plot different genotypes in impact and control sites
library(dplyr)
library(ggplot2)

# Create a new column to group genotypes A–H and I–P
df2_clean <- df2_clean %>%
  mutate(GenotypeGroup = ifelse(Genotype %in% c("A", "B", "C", "D", "E", "F", "G", "H"),
                                "A–H", "I–P")) %>%
  mutate(GenotypeGroup = as.factor(GenotypeGroup))

model_group <- aov(Change.tip.length.tot..cm..day ~ GenotypeGroup, data = df2_clean)
summary(model_group)

#P-value = 0.272
# This is not statistically significant (above the common threshold of 0.05).
# So, there's no evidence that genotype group (A–H vs I–P) has an effect on daily growth rate.
# Belonging to genotype group A–H vs I–P does not significantly affect growth in all sites

# Rajaa data vain Impact-alueen havaintoihin
impact_data <- df2_clean %>% filter(I.or.C == "Impact")

# Testaa vaikutus GenotypeGroup:iin pelkästään Impact-saitilla
model_impact <- aov(Change.tip.length.tot..cm..day ~ GenotypeGroup, data = impact_data)
summary(model_impact)

# p-value is 0.843
# This means genotype group membership does not influence daily growth in the Impact environment in your data.

# Varmistetaan, että GenotypeGroup ja I.or.C ovat faktoreita
df2_clean <- df2_clean %>%
  mutate(
    GenotypeGroup = factor(GenotypeGroup, levels = c("A–H", "I–P")),
    I.or.C = factor(I.or.C, levels = c("Impact", "Control"))
  )

# Piirretään boxplot
ggplot(df2_clean, aes(x = I.or.C, y = Change.tip.length.tot..cm..day, fill = GenotypeGroup)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("A–H" = "#1f77b4", "I–P" = "#ff7f0e")) +
  labs(
    title = "Daily Growth by Genotype Group and Site Type",
    x = "Site Type",
    y = "Daily Growth (cm/day)",
    fill = "Genotype Group"
  ) +
  theme_minimal()




# Calculate mean growth per genotype
mean_growths <- df2_clean %>%
  group_by(Genotype, Group) %>%
  summarise(mean_growth = mean(Change.tip.length.tot..cm..day, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(mean_growth))

# Plot
ggplot(mean_growths, aes(x = reorder(Genotype, -mean_growth), y = mean_growth, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = c("A–H" = "#1f77b4", "I–P" = "#ff7f0e")) +  # Custom colors
  labs(title = "Average Growth by Genotype",
       x = "Genotype",
       y = "Daily Growth (cm/day)",
       fill = "Genotype Group") +
  theme_minimal()




# Plot different genotypes in impact and control sites
ggplot(df2_clean, aes(x = Genotype, y = Change.tip.length.tot..cm..day, fill = I.or.C)) +
  geom_boxplot() +
  labs(title = "Kasvu per genotyyppi (Impact vs Control)", y = "Päivittäinen kasvu (cm/day)") +
  theme_minimal()



# Lasketaan keskiarvot genotyypin ja site-tyypin mukaan
mean_growths_by_site <- df2_clean %>%
  group_by(Genotype, Group, I.or.C) %>%
  summarise(mean_growth = mean(Change.tip.length.tot..cm..day, na.rm = TRUE)) %>%
  ungroup()

# Piirretään kuvaaja: vierekkäiset palkit Control ja Impact
ggplot(mean_growths_by_site, aes(x = Genotype, y = mean_growth, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), aes(group = I.or.C), width = 0.7) +
  facet_wrap(~ I.or.C) +  # Voit halutessasi käyttää facetointia tai vierekkäisiä palkkeja
  scale_fill_manual(values = c("A–H" = "#1f77b4", "I–P" = "#ff7f0e")) +
  labs(
    title = "Average Daily Growth by Genotype and Site Type",
    x = "Genotype",
    y = "Daily Growth (cm/day)",
    fill = "Genotype Group"
  ) +
  theme_minimal()




# Plot different genotypes in all sites

mean_growths <- df2_clean %>%
  group_by(Genotype) %>%
  summarise(mean_growth = mean(Change.tip.length.tot..cm..day, na.rm = TRUE)) %>%
  arrange(desc(mean_growth))  # Järjestetään suurimmasta pienimpään

# Piirrä kaavio
ggplot(mean_growths, aes(x = reorder(Genotype, -mean_growth), y = mean_growth)) +
  geom_col(fill = "steelblue") +
  labs(title = "Genotyyppien keskimääräinen kasvu",
       x = "Genotyyppi", y = "Päivittäinen kasvu (cm/day)") +
  theme_minimal()


########
########
########

# Modelling with GLMM

# It is particularly useful when you need:
# Both fixed and random effects (e.g., hierarchical or grouped data),
# To model non-normal response variables (e.g., binomial, Poisson, negative binomial),

# Let's use weight as it takes into account both growth in length and width
########
########
########

# Bring data in
df <- read.csv("Genotyyppitarkasteluun_dry_weight.csv", stringsAsFactors = FALSE)
df2 <- read.csv("Tips_dry_weigth.csv", stringsAsFactors = FALSE)

# Clean data and combine information of dry weights in the main dataframe
library(dplyr)
library(stringr)

df <- df %>%
  mutate(Fucus.ID.key = paste0(SITE.ID, Tile.number, Genotype))


# 1. Luo df2:een Fucus.ID.key, eli site + tile number + genotype
df2 <- df2 %>%
  mutate(Fucus.ID = chartr("ÅÖ", "AO", Fucus.ID))

df2 <- df2 %>%
  mutate(
    site = str_extract(Fucus.ID, "^[A-Z]+"),
    tile_number = str_extract(Fucus.ID, "\\d"),  # ensimmäinen numero
    genotype = str_extract(Fucus.ID, "\\(([A-Z])\\)") %>% str_remove_all("[()]"),
    Fucus.ID.key = paste0(site, tile_number, genotype)
  )



# 2. Yhdistä df2:n Tip.weight df:ään Fucus.ID.key:n mukaan
df <- df %>%
  left_join(df2 %>% select(Fucus.ID.key, Tip.weight), by = "Fucus.ID.key") %>%
  mutate(Tip.dry.weight.July..g. = Tip.weight) %>%
  select(-Tip.weight)  # Poistetaan ylimääräinen sarake

# Save the new dataframe into file
install.packages("writexl")
library(writexl)
# write_xlsx(df, "Dryweight.xlsx")

# fit generalized linear mixed model (GLMM), Random intercept and slope model
#Growth = Wave_environment + Genotype_origin + Wave_environment*Genotype_origin+ Random Site +Plot +Genotype +Site*Genotype +Wave_environment*Genotype

# create a column where it is said is the origin of genotypes impact or control
df <- df %>%
  mutate(Genotype.origin = ifelse(Genotype %in% LETTERS[1:8], "Control", "Impact"))

# Save to file again
# write_xlsx(df, "Dryweight2.xlsx")

# Box plots of dry weight change
df$Tip.dry.weight.change <- df$Tip.dry.weight.July..g. - df$Tip.dry.weight.June.ww_dw_ratio..g.

library(ggplot2)
ggplot(df, aes(x = I.or.C, y = Tip.dry.weight.change, fill = I.or.C)) +
  geom_boxplot() +
  labs(
    title = "Change in Tip Dry Weight (July - June)",
    x = "Site Type",
    y = "Tip Dry Weight Change (g)"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "steelblue", "Impact" = "firebrick"))

# Box plots alku- ja lopputilanteesta
library(tidyr)
df_long <- df %>%
  select(I.or.C, Tip.dry.weight.June = Tip.dry.weight.June.ww_dw_ratio..g., 
         Tip.dry.weight.July = Tip.dry.weight.July..g.) %>%
  pivot_longer(cols = starts_with("Tip.dry.weight"), names_to = "Month", values_to = "Weight")

ggplot(df_long, aes(x = Month, y = Weight, fill = I.or.C)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Tip Dry Weight in June and July",
    x = "Month",
    y = "Weight (g)",
    fill = "Site Type"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "steelblue", "Impact" = "firebrick"))

# Viiva plotti
df_long <- df %>%
  select(ID = Fucus.ID.key, I.or.C,
         June = Tip.dry.weight.June.ww_dw_ratio..g.,
         July = Tip.dry.weight.July..g.) %>%
  pivot_longer(cols = c(June, July), names_to = "Month", values_to = "Weight")

df_long$Month <- factor(df_long$Month, levels = c("June", "July"))

ggplot(df_long, aes(x = Month, y = Weight, group = ID, color = I.or.C)) +
  geom_line(alpha = 0.6) +
  labs(
    title = "Change in Tip Dry Weight (June to July)",
    x = "Month",
    y = "Tip Dry Weight (g)",
    color = "Site Type"
  ) +
  theme_minimal()

ggplot(df, aes(x = Tip.wet.weight.July..g., y = Tip.dry.weight.July..g.)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  labs(
    title = "Wet vs. Dry Tip Weight in July",
    x = "Wet Tip Weight (g, July)",
    y = "Dry Tip Weight (g, July)"
  ) +
  theme_minimal()
# Alkupaino SH2 N alkupaino väärä -> tarkista paperi

# Viiva plotti pituudelle
df$June.tip.length.tot..cm.[df$Fucus.ID.key == "SH2M"] <- NA

df_long <- df %>%
  select(ID = Fucus.ID.key, I.or.C,
         June = June.tip.length.tot..cm.,
         July = July.tip.length.tot..cm.) %>%
  pivot_longer(cols = c(June, July), names_to = "Month", values_to = "Length")

df_long$Month <- factor(df_long$Month, levels = c("June", "July"))

ggplot(df_long, aes(x = Month, y = Length, group = ID, color = I.or.C)) +
  geom_line(alpha = 0.6) +
  labs(
    title = "Change in Tip Dry Length (June to July)",
    x = "Month",
    y = "Tip Dry Length",
    color = "Site Type"
  ) +
  theme_minimal()

ggplot(df, aes(x = Tip.wet.weight.July..g., y = Tip.dry.weight.July..g.)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  labs(
    title = "Wet vs. Dry Tip Weight in July",
    x = "Wet Tip Weight (g, July)",
    y = "Dry Tip Weight (g, July)"
  ) +
  theme_minimal()

#SH2M miinus




library(glmmTMB)
# linear mixed model sen takia, että voidaan hävittää informaatiota, joka aiheutuu datan rakenteesta jos käytetään lm
# jos muuttuja saa saman arvon klusterin sisällä niin on fixed

# y = massa lopussa

###
# FIXED FACTORS: 
###
# (riippumaton muuttuja)
# aaltoympäristö (I vs C)
# donor I vs C
# initial biomass

###
# RANDOM FACTORS:
###
# (sellainen ryhmä, jonka sisällä vaihtelu ei randomia suhteutettuna kokonais kuvaan, ja lähtötason erot ryhmien välillä voi vaikuttaa tuloksiin vaikka
# ei vaikuta tutkittavaan kohteeseen, kuten tässä kasvuun, siis random effect on sellainen joka voi vaikuttaa lopputulokseen niin ettei tiedetä
# johtuuko tulokset oikeista eroista vai siitä että jotain muuttujaa ei ole huomioitu niinkuin pitäisi, eli laittamalla ne random effectiksi eikä fixed effectiksi)

# site
# Tile number 
# Genotype


head(df)
df$Tip.dry.weight.June.ww_dw_ratio..g.[df$Fucus.ID.key == "SH2N"] <- NA

# Veijon ohjeen mukaan
m <- glmmTMB(Tip.dry.weight.July..g. ~ Tip.dry.weight.June.ww_dw_ratio..g. + I.or.C + Genotype.origin + (1|SITE.ID) + (1|Tile.number) + (I.or.C|Genotype), data=df, family=gaussian)
# (I.or.C|Genotype) -> testaa tämän merkitsevyys
m_pituus <- glmmTMB(July.tip.length.tot..cm. ~ June.tip.length.tot..cm. + I.or.C + Genotype.origin + (1|SITE.ID) + (1|Tile.number) + (I.or.C|Genotype), data=df, family=gaussian)
m_pituus_ilman_geno <- glmmTMB(July.tip.length.tot..cm. ~ June.tip.length.tot..cm. + I.or.C + Genotype.origin + (1|SITE.ID) + (1|Tile.number), data=df, family=gaussian)
# voi pituuden kokeilla niin että y= pituuden muutos alkumitta muuttujana ja loppupituus pois selittävänä muuttujana, 



model_data <- m$model
head(model_data)

# Mitä tämä malli selvittää? 
# Selitettävä muuttuja: kuivapaino heinäkuussa.
# Selittäjät:
  # Alkupaino (kesäkuussa): huomioidaan, koska se vaikuttaa kasvuun.
  # I.or.C (impact vs. control): tutkii, vaikuttaako ihmistoiminta kasvuun.
  # Genotype.origin (geneettinen alkuperä): tutkitaan, vaikuttaako alkuperä kasvuun.
  # I.or.C * Genotype.origin: onko vaikutus erilainen eri alkuperän genotyypeillä.
# Ristikkäistermit (interaktiot):
  # SITE.ID * Genotype ja I.or.C * Genotype: monimutkaisia interaktioita, jotka voivat olla vaikeita tulkita ja lisätä mallin monimutkaisuutta merkittävästi.
# Satunnaiset efektit:
  # Satunnaisvaihtelu paikan (SITE.ID), tiilen (Tile.number) ja genotyypin mukaan.

# -> Ylispesifinen: sisältää paljon interaktioita ja päällekkäisiä termejä (esim. sekä SITE.ID * Genotype että (1 | SITE.ID)).
# -> I.or.C * Genotype on sama asia kuin I.or.C + Genotype + I.or.C:Genotype
# -> Sisältää molemmat päävaikutukset sekä niiden interaktion

# Cleaned version
m_clean <- glmmTMB(
  Tip.dry.weight.July..g. ~ Tip.dry.weight.June.ww_dw_ratio..g. +
    I.or.C * Genotype.origin +
    (1 | SITE.ID) + (1 | Tile.number) + (1 | Genotype),
  data = df, family = gaussian)
# Mitä tämä malli selvittää?
# Kasvun (kuivapaino heinäkuussa) erot selittyvät:
# Alkupainolla (kesäkuussa)
# Vaikutusalueella (I.or.C): Onko kasvu erilaista ihmistoiminnan vaikutusalueilla?
# Genotyypin alkuperällä (Genotype.origin)
# Näiden kahden interaktiolla: Onko ihmistoiminnan vaikutus erilainen genotyypin alkuperän mukaan?
# Satunnaiset efektit huomioivat:
# Eri paikkojen, tiilien ja genotyyppien satunnaisvaihtelun.

summary(m)
summary(m_pituus)
summary(m_pituus_ilman_geno)
anova(m_pituus, m_pituus_ilman_geno)
# mitä pienempi AIC sen parempi, 2-6 pienempi on jo parempi, tässä paljon parempi kun geno mukana -> genotyypin merkitsevyystesti oli tämä, oli merkitsevästi vaikutusta
# genotyypin kasvutaso eroaa, mutta ei testaa I ja C vaikutusta yhdessä eli slopea saada

# mallin keskiarvot otetaan
library(emmeans)
emmeans(m, ~ I.or.C, type="response")

# Katsotaan miten malli istuu
library(DHARMa)

pearson_resid <- residuals(m, type = "pearson")
plot(pearson_resid)
hist(pearson_resid)


r <- simulateResiduals(m, n=1000)
plot(r)
hist(residuals(m),
     main = "Histogram of Residuals",
     xlab = "Residuals",
     col = "lightblue",
     border = "white") # normal

length(pearson_resid)
nrow(m$model)  # tämä pitäisi täsmätä

pearson_df <- data.frame(
  ID = names(pearson_resid),
  pearson_resid = as.numeric(pearson_resid)
)

resid_df_fucus <- model_data pearson_resid


resid_df_fucus <- data.frame(
  Fucus.ID.key = model_data$Fucus.ID.key,
  pearson_resid = pearson_resid)




# Tarkista multicollinearity
