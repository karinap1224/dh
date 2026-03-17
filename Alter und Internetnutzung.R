
# Simulation für Missing-Data-Mechanismen MCAR, MAR, MNAR

# ESS-Datensatz 2023 (Runde 11)

# Variablen und Rahmen: Alter (yrbrn) und Internetnutzung (netustm) in Deutschland 

install.packages("mice")
library(conflicted)
library(haven)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(broom)
library(mice)
library(tidyverse)
library(dplyr)

# Vorgehensweise: 

# Schritt 1: Referenzdatensatz erstellen "Ground Truth"
# Pakete laden: mice, haven, tidyverse, broom, etc
# Datenimport & Breinigung: 
#   - ESS 11 Datensatz laden.
#   - Alter berechnen (dplyr::mutate()) -> alter
#   - Missin-Codes von netustm in NA umwandeln (dplyr::filter()) oder ifelse()
# Grundgesamtheit definieren -> Referenzdatensatz erstellen (df_true <- ess %>% select(alter, netustm) %>% drop_na())
# -> Modell für df_true berechnen -> lm_true <- lm(netustm ~ alter, data = df_true)

# Schritt 2: Simulationsstudie durch Amputation (mice::ampute)
# Seed setzen: set.seed(123) -> Ergebnisse müssen reproduzierbar bleiben
# 1. MCAR: - mcar_data <- ampute(df_true, prop = 0.2, mech = "MCAR")$amp
#          - 20 % der Werte werden gelöscht 
# 2. MAR:  - mar_data <- ampute(df_true, prop = 0.2, mech = "MAR", weights = c(1, 0))$amp
#          - Durch Gewichte: Wahrscheinlichkeit für ein fehlendes netustm hängt ausschlieslich vom Alter ab (Ältere haben eher ein NA)
# 3. MNAR: - mnar_data <- ampute(df_true, prop = 0.2, mech = "MNAR", weights = c(0, 1))$amp
#          - Wahrscheinlichkeit für ein fehlendes netustm hängt vom wahren Wert der Internetnutzung selbst ab. 

# Schritt 3: Analyse und Gegenmaßnahme 
# 1) Complete Case Analysis (CCA/Listwise Deletion)
#     - lm(netustm ~ alter) für alle drei Mechanismen (NAs werden automatisch gelöscht)
#     -> Speicherung der Koeffizienten. MAR und MNAR sollten nun vom wahren lm_true-Koeffizienten abweichen (Bias)
# 2) Inverse Probability Weighting (IPW) (packages: lmtest und sandwich für robuste Standardfehler)
#     - Für MAR und MNAR only: 
#     1. Missing-Indikator erstellen: Hilfsvariable, die angibt, ob der Wert fehlt oder nicht.
#       mar_data$R_obs <- ifelse(is.na(mar_data$netustm), 0, 1)
#     2. Propensity Modell schätzen: Rechnet eine logistische Reg., um die Beobachtungswahrscheinlichk. vorherzusagen
#       prop_model <- glm(R_obs ~ alter, family = binomial, data = mar_data)
#     3. Gewichte berechnen: Das Gewicht ist der Kehrwert der vorhergesagten Wahrsch. (Inverse Probability)
#       mar_data$weight <- 1 / fitted(prop_model)
#     4. Gewichtete Regression: Analyse nur mit vollständigen Fällen, aber mit Gewichten.
#       lm(netustm ~ alter), data = subset(mar_data, R_obs == 1), weights = weight)
#     -> Bei MAR: Gewichtung korregiert den Bias, da die Auswahlwahrscheinlichkeit
#         vollständig durch das Alter erklärt werden kann. Reg.koeff sollte wieder nah an df_true liegen.
#     -> Bei MNAR: IPW scheitert (genau wie 'mice'). Da bei MNAR der Ausfall von
#         Internetnutzung selsbt abhängt (die wir nicht kennen), bringt die Gewichtung
#         anhand des Alters nicht den gewünschten Erfolg. Der Bias bleibt bestehen.




# Schritt 1: Referenzdatensatz erstellen "Ground Truth"

library(readr)
ESS11e04_1 <- read_csv("C:/Naomi/Studium Master/1. Semester (WS25_26)/Datenerhebung und Fehlerquellen/ESS11e04_1.csv")
View(ESS11e04_1)

ess <- ESS11e04_1

# Beschränkung auf Deutschland

de <- ess %>%
  mice::filter(cntry == "DE")
View(de)

summary(de)
dim(de)
str(de)


# Deskriptive Statistik 

table(de$yrbrn)
de$gndr
table(de$netustm)

# Alter berechnen und bereinigen

gebj <- de$yrbrn

yrbrn_num <- as.numeric(gebj)

de_data1 <- de %>% mutate(
  gebj_num <- as.numeric(yrbrn_num),
  yrbrn_clean = ifelse(yrbrn_num %in% c(7777, 8888, 9999), NA, yrbrn_num),
  alter = 2026 - yrbrn_clean)


# Internetnutzung bereinigen

de_data2 <- de_data1 %>% mutate(
  netustm_num = as.numeric(netustm),
  netustm_clean = ifelse(netustm_num %in% c(6666, 7777, 8888, 9999), NA, netustm_num)
) %>%

# Nur unsere beiden relevanten Variablen behalten und alle NAs löschen
  select(alter, netustm = netustm_clean) %>% drop_na()

df_true <- de_data2 #Referenzdatensatz erstellt 
coef(df_true)
View(df_true)

# Modell berechnen für df_true
lm_true <- lm(netustm ~ alter, data = df_true)
View(lm_true)
summary(lm_true)


# Schritt 2: Simulationsstudie durch Amputation (mice::ampute)
# Bei jedem Mechanismus werden 20 % der Werte werden gelöscht.
set.seed(42) # Gleicher Zufall für Reproduzierbarkeit


# MCAR:20 % der Werte werden gelöscht. 
mcar_data <- mcar_data <- ampute(df_true, prop = 0.2, mech = "MCAR")$amp


# MAR: Durch Gewichte: Wahrscheinlichkeit für ein fehlendes netustm hängt 
# ausschlieslich vom Alter ab (Ältere haben eher ein NA)
set.seed(42)
mar_data <- ampute(df_true, prop = 0.2, mech = "MAR", weights = c(1, 0))$amp


# MNAR: Wahrscheinlichkeit für ein fehlendes netustm hängt vom wahren Wert der 
# Internetnutzung selbst ab. 
set.seed(42)
mnar_data <- ampute(df_true, prop = 0.2, mech = "MNAR", weights = c(0, 1))$amp


# Phase 3: Analyse und Korrekturmaßnahmen
# 1) Complete Case Analysis (CCA/Listwise Deletion)

lm_mcar <- lm(netustm ~ alter, data = mcar_data)
coef(lm_mcar)

lm_mar <- lm(netustm ~ alter, data = mar_data)
coef(lm_mar)

lm_mnar <- lm(netustm ~ alter, data = mnar_data)
coef(lm_mnar)


























