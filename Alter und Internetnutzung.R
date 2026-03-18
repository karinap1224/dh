
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
de$eduyrs
table(de$netustm)

# Kontrollvariablen bereinigen

de_data1 <- de %>% mutate(
  gndr_num = as.numeric(gndr),
  gndr_clean = ifelse(gndr_num == 9, NA, gndr_num),
  
  edu_num = as.numeric(eduyrs), 
  edu_clean = ifelse(edu_num %in% c(77, 88, 99), NA, edu_num)
)

# Alter berechnen und bereinigen

de_data2 <- de_data1 %>% mutate(
  yrbrn_num = as.numeric(yrbrn),
  yrbrn_clean = ifelse(yrbrn_num %in% c(7777, 8888, 9999), NA, yrbrn_num),
  alter = 2023 - yrbrn_clean)

# Internetnutzung bereinigen

de_data_clean <- de_data2 %>% mutate(
  netustm_num = as.numeric(netustm),
  netustm_clean = ifelse(netustm_num %in% c(6666, 7777, 8888, 9999), NA, netustm_num)
) %>%

# Nur unsere relevanten Variablen behalten und alle NAs löschen
select(
  alter, 
  netustm = netustm_clean, 
  gndr = gndr_clean, 
  eduyrs = edu_clean
) %>% drop_na()


df_true <- de_data_clean #Referenzdatensatz erstellt 
View(df_true)


# Modelle berechnen für df_true

lm_true1 <- lm(netustm ~ alter, data = df_true)
coef(lm_true1)
View(lm_true1)
summary(lm_true1)

lm_true2 <- lm(netustm ~ alter + gndr + eduyrs, data = df_true)
coef(lm_true2)
View(lm_true2)
summary(lm_true2)

# Schritt 2: Simulationsstudie durch Amputation (mice::ampute)
# Bei jedem Mechanismus werden 20 % der Werte werden gelöscht.
set.seed(42) # Gleicher Zufall für Reproduzierbarkeit


# MCAR:20 % der Werte werden gelöscht. 
mcar_data <- ampute(df_true, prop = 0.2, mech = "MCAR")$amp
summary(mcar_data)


# MAR: Durch Gewichte: Wahrscheinlichkeit für ein fehlendes netustm hängt 
# ausschlieslich vom Alter ab (Ältere haben eher ein NA)

ncol(df_true) == 4 # Ergebnis auf TRUE

# 2. PATTERNS: Sagen, WO Lücken entstehen sollen
# 1 = Behalten, 0 = Löschen. 
# Reihenfolge: alter(1), netustm(0), gndr(1), eduyrs(1)
mypattern <- c(1, 0, 1, 1)

# 3. WEIGHTS: Sagen, WARUM Lücken entstehen (Der MAR-Mechanismus)
# Wir geben nur dem "alter" (Spalte 1) ein Gewicht von 1, dem Rest eine 0.
# Dadurch steuert ausschließlich das Alter, ob netustm gelöscht wird.
myweights_mar <- c(1, 0, 0, 0)

# 4. Amputation durchführen
set.seed(42)
amp_mar <- ampute(df_true, 
                  prop = 0.2, 
                  mech = "MAR", 
                  patterns = mypattern, 
                  weights = myweights_mar)

mar_data <- amp_mar$amp
summary(mar_data$netustm)


# MNAR: Wahrscheinlichkeit für ein fehlendes netustm hängt vom wahren Wert der 
# Internetnutzung selbst ab. 
set.seed(42)
ncol(df_true) == 4 # Ergebnis auf TRUE

myweights_mnar <- c(0, 1, 0, 0)

set.seed(42)
amp_mnar <- ampute(df_true, 
                  prop = 0.2, 
                  mech = "MNAR", 
                  patterns = mypattern, 
                  weights = myweights_mnar)

mnar_data <- amp_mnar$amp
summary(mnar_data$netustm)


# Phase 3: Analyse und Korrekturmaßnahmen
# 1) Complete Case Analysis (CCA/Listwise Deletion)

lm_mcar <- lm(netustm ~ alter + gndr + eduyrs, data = mcar_data)
coef(lm_mcar)
summary(lm_mcar)

lm_mar <- lm(netustm ~ alter + gndr + eduyrs, data = mar_data)
coef(lm_mar)
summary(lm_mar)

lm_mnar <- lm(netustm ~ alter + gndr + eduyrs, data = mnar_data)
coef(lm_mnar)
summary(lm_mnar)

# Vergleich von Mechanismen und Referenzdatensatz
# MCAR
coef(lm_true2)
coef(lm_mcar)
summary(lm_true2)
summary(lm_mcar)

# MAR
coef(lm_true2)
coef(lm_mar)
summary(lm_true2)
summary(lm_mar)

#MNAR
coef(lm_true2)
coef(lm_mnar)
summary(lm_true2)
summary(lm_mnar)


# 2) Inverse Probability Weighting (IPW)

# MAR
# 1.	Missing-Indikator erstellen: Erstellt eine Hilfsvariable (Dummy), 
#     die angibt, ob der Wert fehlt oder nicht.
mar_data$R_obs <- ifelse(is.na(mar_data$netustm), 0, 1)

# 2.	Propensity Modell schätzen: Rechnet eine logistische Regression, 
#     um die Beobachtungswahrscheinlichkeit vorherzusagen.
prop_model <- glm(R_obs ~ alter + gndr + eduyrs, family = binomial, data = mar_data)
summary(prop_model)

# 3.	Gewichte berechnen: Das Gewicht ist der Kehrwert der vorhergesagten 
#     Wahrscheinlichkeit (Inverse Probability).
mar_data$weight <- 1 / fitted(prop_model)

# 4.	Gewichtete Regression: Führt die finale Analyse durch. Für die Analyse 
#     nutzt ihr nur die vollständigen Fälle, aber wendet die Gewichte an.
lm(netustm ~ alter + gndr + eduyrs, data = subset(mar_data, R_obs == 1), weights = weight)
























