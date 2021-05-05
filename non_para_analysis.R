# Data Reading and Tidying
library(tidyverse)
library(quantreg)
library(MASS)
library(aod)
library(survey)


plates <- read.csv("raw_data.csv")

plates_t <- reshape(plates, 
                    direction = "long",
                    varying = list(names(plates)[3:9]),
                    v.names = "Lysed_Wells",
                    idvar = c("Day", "Plate"),
                    timevar = "Column",
                    times = 1:7)


plates_t <- separate(plates_t, Plate, c("Ad", "Treatment", "Replicate"), sep = ";")

plates_t$AdT <- paste(plates_t$Ad, plates_t$Treatment)

plates_t$AdTCol <- paste(plates_t$AdT, plates_t$Column)

plates_t <- plates_t %>% mutate(Lysed_Wells_perc = (Lysed_Wells/ 8) * 100)

plates_t$Day <- as.character(plates_t$Day)
plates_t$Column <- as.character(plates_t$Column)
str(plates_t)


# Calculate Median

plates_med <- plates_t %>% 
              group_by(Day, Ad, Treatment, Column) %>% 
                summarise(Lysed_Perc_Med = median(Lysed_Wells_perc),
                          Lysed_Perc_IQR = IQR(Lysed_Wells_perc))

plates_med$AdT <- paste(plates_med$Ad, plates_med$Treatment)

# plates_med$Lysed_Perc_Med <- as.factor(plates_med$Lysed_Perc_Med)
str(plates_med)

plates_med <- plates_med %>% filter(Column != 7)
        

# Models

plates_med.rq <- rq(Lysed_Perc_Med ~ Treatment + Ad, data = plates_med)
summary(plates_med.rq)

plates_med.kw <- kruskal.test(Lysed_Perc_Med ~ Treatment, data = plates_med)

plates_med.f <- friedman.test(Lysed_Perc_Med ~ Treatment | Ad, data = plates_med)

plates_med.olr <- polr(Lysed_Perc_Med ~ Treatment + Ad + Day + Column, data = plates_med)
summary(plates_med.olr)

(ctable <- coef(summary(plates_med.olr)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
(ctable <- cbind(ctable, "p value" = p))

regTermTest(plates_med.olr, "Treatment", method = "Wald")



# OLR Pure Lysate

plates_med_pl <- plates_med %>% filter(Column == 2, Day == 14)

plates_med_pl.olr <- polr(Lysed_Perc_Med ~ Treatment + Ad + Day, data = plates_med_pl)
summary(plates_med_pl.olr)

(ctable_pl <- coef(summary(plates_med_pl.olr)))
p_pl <- pnorm(abs(ctable_pl[, "t value"]), lower.tail = FALSE) * 2
(ctable_pl <- cbind(ctable_pl, "p value" = p_pl))


plates_med_pl_2.olr <- polr(Lysed_Perc_Med ~ AdT, data = plates_med_pl)
summary(plates_med_pl_2.olr)

(ctable_pl_2 <- coef(summary(plates_med_pl_2.olr)))
p_pl_2 <- pnorm(abs(ctable_pl_2[, "t value"]), lower.tail = FALSE) * 2
(ctable_pl_2 <- cbind(ctable_pl_2, "p value" = p_pl_2))




(ggplot(plates_med, aes(fill = Column, x = AdT, y = Lysed_Perc_Med)) +
  geom_bar(stat = "identity", position = "dodge")) +
  facet_wrap(vars(Day))

plates_med_pl <- plates_med %>% filter(Column == 2)
(ggplot(plates_med_pl, aes(x = AdT, y = Lysed_Perc_Med)) +
    geom_bar(stat = "identity", position = "dodge")) +
  facet_wrap(vars(Day))



# First Attempt at Parametric Analysis ----

# Adsorbed Sediment Effect Day 10

ad.aov_data_10 <- filter(plates_t, Ad == "A", Day == "10")
ad.aov_10 <- aov(Lysed_Wells ~ Treatment, ad.aov_data_10)
summary(ad.aov_10)
turkey_ad_10 <- HSD.test(ad.aov_10, trt = "Treatment")

(Ad_10_bar <- ggplot(ad.aov_10, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))

# Adsorbed Sediment Effect Day 14

ad.aov_data_14 <- filter(plates_t, Ad == "A", Day == "14")
ad.aov_14 <- aov(Lysed_Wells ~ Treatment, ad.aov_data_14)
summary(ad.aov_14)
turkey_ad_14 <- HSD.test(ad.aov_14, trt = "Treatment")

(Ad_14_bar <- ggplot(ad.aov_14, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))


# Desorbed Sediment Effect Day 10

de.aov_data_10 <- filter(plates_t, Ad == "D", Day == "10")
de.aov_10 <- aov(Lysed_Wells ~ Treatment, de.aov_data_10)
summary(de.aov_10)
turkey_de_10 <- HSD.test(de.aov_10, trt = "Treatment")

(De_10_bar <- ggplot(de.aov_10, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))

# Desorbed Sediment Effect Day 14

de.aov_data_14 <- filter(plates_t, Ad == "D", Day == "14")
de.aov_14 <- aov(Lysed_Wells ~ Treatment, de.aov_data_14)
summary(de.aov_14)
turkey_de_14 <- HSD.test(de.aov_14, trt = "Treatment")

(De_14_bar <- ggplot(de.aov_14, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))


# Adsorped Sediment Effect Day 10 Pure Lysate

ad.aov_data_10_2 <- filter(plates_t, Ad == "A", Day == "10", Column == "2")
ad.aov_10_2 <- aov(Lysed_Wells ~ Treatment, ad.aov_data_10_2)
summary(ad.aov_10_2)
turkey_ad_10_2 <- HSD.test(ad.aov_10_2, trt = "Treatment")

(Ad_10_bar_2 <- ggplot(ad.aov_10_2, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))


# Adsorped Sediment Effect Day 14 Pure Lysate

ad.aov_data_14_2 <- filter(plates_t, Ad == "A", Day == "14", Column == "2")
ad.aov_14_2 <- aov(Lysed_Wells ~ Treatment, ad.aov_data_14_2)
summary(ad.aov_14_2)
turkey_ad_14_2 <- HSD.test(ad.aov_14_2, trt = "Treatment")

(Ad_14_bar_2 <- ggplot(ad.aov_14_2, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))


# Desorped Sediment Effect Day 10 Pure Lysate

de.aov_data_10_2 <- filter(plates_t, Ad == "D", Day == "10", Column == "2")
de.aov_10_2 <- aov(Lysed_Wells ~ Treatment, de.aov_data_10_2)
summary(de.aov_10_2)
turkey_de_10_2 <- HSD.test(de.aov_10_2, trt = "Treatment")

(de_10_bar_2 <- ggplot(de.aov_10_2, aes(x = Treatment, y = Lysed_Wells)) +
    geom_bar(stat = "summary"))


# Tidy Control - Sediment Size Classes ----

treatment_diff <- read.csv("treatment_diff.csv")

treatment_diff_t <- reshape(treatment_diff, 
                            direction = "long",
                            varying = list(names(treatment_diff)[2:8]),
                            v.names = "Lysed_Wells",
                            idvar = c("Plate"),
                            timevar = "Column",
                            times = 1:7)


treatment_diff_t <- separate(treatment_diff_t, Plate, c("Day", "Treatment_Diff", "Ad_De", "Replicate"), sep = ";")

treatment_diff_t$AdTCol <- paste(treatment_diff_t$Treatment_Diff, treatment_diff_t$Column)

treatment_diff_t$Column <- as.character(treatment_diff_t$Column)
str(treatment_diff_t)

treatment_diff_t <- filter(treatment_diff_t, Column != 7)


# Control - Sediment Size Class ANOVA

diff.aov_data_14 <- filter(treatment_diff_t, Day == "14", Column == "2")
diff.aov_14 <- aov(Lysed_Wells ~ Treatment_Diff * Ad_De, diff.aov_data_14)
summary(diff.aov_14)                                                          # Not significant, makes all values more closely related
turkey_diff_14 <- HSD.test(diff.aov_14, trt = c("Treatment_Diff", "Ad_De"))

diff.aov_data_17 <- filter(treatment_diff_t, Day == "17", Column == "2")
diff.aov_17 <- aov(Lysed_Wells ~ Treatment_Diff * Ad_De, diff.aov_data_17)
summary(diff.aov_17)                                                          # Not significant, makes all values more closely related
turkey_diff_17 <- HSD.test(diff.aov_17, trt = c("Treatment_Diff", "Ad_De"))

# Stepwise Linear Regression of Dilutions

#Get Data

plates_dils <- read.csv("raw_data_2.csv")

plates_dils_t <- reshape(plates_dils, 
                         direction = "long",
                         varying = list(names(plates_dils)[3:9]),
                         v.names = "Lysed_Wells",
                         idvar = c("Day", "Plate"),
                         timevar = "Dilution",
                         times = c("Control", "1", "0.1", "0.01", "0.001", "0.0001", "0.00001"))


plates_dils_t <- separate(plates_dils_t, Plate, c("Ad", "Treatment", "Replicate"), sep = ";")

plates_dils_t$AdT <- paste(plates_dils_t$Ad, plates_dils_t$Treatment)

plates_dils_t$AdTCol <- paste(plates_dils_t$AdT, plates_dils_t$Dilution)

plates_dils_t <- plates_dils_t %>% mutate(Lysed_Wells_perc = (Lysed_Wells/ 8) * 100)

plates_dils_t$Day <- as.character(plates_dils_t$Day)
plates_dils_t$Dilution <- as.character(plates_dils_t$Dilution)
plates_dils_t$Lysed_Wells_perc <- as.numeric(plates_dils_t$Lysed_Wells_perc)
str(plates_dils_t)

plates_dils_t <- plates_dils_t %>% filter(Dilution != "0.00001", Dilution != "Control")

plates_dils_t_14 <- plates_dils_t %>% filter(Day == 14)
plates_dils_t_17 <- plates_dils_t %>% filter(Day == 17)

dil.lm_14 <- lm(Lysed_Wells ~ Ad * Treatment * Dilution, plates_dils_t_14)  # Model with all variables
dil.step <- stepAIC(dil.lm_14, direction = "both")
summary(dil.lm)

dil.lm_17 <- lm(Lysed_Wells ~ Ad * Treatment * Dilution, plates_dils_t_17)  # Model with day 17 column 2
dil.step_17 <- stepAIC(dil.lm_17, direction = "both")
summary(dil.lm)


# Test if Ad has an effect on the model
dil_no_ad.lm <- lm(Lysed_Wells ~ Treatment * Column, plates_t)  # Model without Ad
summary(dil_no_ad.lm)

anova(dil.lm, dil_no_ad.lm)  # Ad does have an effect


