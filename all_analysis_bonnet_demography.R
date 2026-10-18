#encounter rate kruskal wallis

#analysis- encounter
er_data<-read.csv("D:/Demography paper/takanalysis.csv")
library(dplyr)

summary <- er_data %>%
  group_by(gridcode) %>%
  summarise(
    er_mean = mean(er, na.rm = TRUE),
    er_sd = sd(er, na.rm = TRUE),
    indkm_mean = mean(ind.km, na.rm = TRUE),
    indkm_sd = sd(ind.km, na.rm = TRUE)
  )

print(summary)
write.csv(summary, "takanalysis_results.csv", row.names = FALSE)


kruskal.test(er ~ gridcode, data = er_data)
kruskal.test(ind.km ~ gridcode, data = er_data)


datanew<-read.csv("D:/Demography paper/R/extract_values_bonnet_new.csv")
kruskal.test(total ~ factor(land_use), data = datanew)

#GROUP SIZE -error bars 
library(dplyr)
summary_gs <- datanew %>%
  group_by(land_use) %>%
  summarise(
    mean_gs = mean(total, na.rm = TRUE),
    sd = sd(total, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci = 1.96 * se   # 95% CI
  )

error_bar_gsize<-ggplot(summary_gs,aes(x=reorder(land_use, -mean_gs),y=mean_gs))+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_gs - ci,
                    ymax = mean_gs + ci),
                width = 0.2) +
  labs(
    x = "Land use",
    y = "Mean group size (± 95% CI)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, size = 10, colour ="black", hjust = 0.5),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid = element_blank()
  )

ggsave(
  filename = "group_sizeerrorbars_landuse.png",
  plot = error_bar_gsize,
  width = 6,
  height = 4,
  dpi = 300
)

#mean no of females and males per group
mean(datanew$female)
sd(datanew$female)

mean(datanew$male)
sd(datanew$male)

t.test(datanew$female, datanew$male, paired = TRUE)


#chi square test- proportion of males and females across all land use types
library(dplyr)
sex_table <- datanew %>%
  group_by(land_use) %>%
  summarise(
    male = sum(male, na.rm = TRUE),
    female = sum(female, na.rm = TRUE)
  )
sex_mat <- as.matrix(sex_table[, c("male", "female")])
rownames(sex_mat) <- sex_table$land_use
chisq.test(sex_mat)

#kruskal wallis demographics

datanew$immature<- datanew$infant + datanew$juvenile
vars <- c("male", "female", "mf", "fi", "ai", "immature")

results <- lapply(vars, function(v) {
  test <- kruskal.test(as.formula(paste(v, "~ land_use")), data = scaled_envlayers)
  data.frame(
    Variable = v,
    Chi_sq = round(test$statistic, 2),
    df = test$parameter,
    p_value = round(test$p.value, 3)
  )
})

results_table <- do.call(rbind, results)
results_table

library(FSA)
#INDIVIDUAL
kruskal.test(total ~ land_use, data = scaled_envlayers)

#post-hoc Dunn's Bonferroni
# Example for fi ratio
dunnTest(fi ~ land_use, data = datanew, method = "bonferroni")

# Example for ai ratio
dunnTest(ai ~ land_use, data = datanew, method = "bonferroni")

dunnTest(immature ~ land_use, data = datanew, method = "bonferroni")

#scaled
scaled_envlayers<-read.csv("D:/Demography paper/layers/scaled_output/scaled_extractionlayers.csv")
nb_model <- glm.nb(
  total ~ ndvi_scale + footprint_ + elevation_ +
    bio_9_scal + bio_3_scal + bio_2_scal +
    bio_19_sca + bio_17_sca + bio_14_sca + nightlight +  x + y,
  data = scaled_envlayers
)

summary(nb_model)

library(MASS)

step_model <- stepAIC(nb_model, direction = "both")
summary(step_model)

#male
mean(datanew$male)
var(datanew$male)
#mean=variance - using poisson
#scaled male
pois_male_scaled <- glm(male ~  ndvi_scale + footprint_ + elevation_ + 
                          bio_9_scal + bio_3_scal + bio_2_scal +
                          bio_19_sca + bio_17_sca + bio_14_sca + nightlight 
                        + x + y,
                        data = scaled_envlayers, 
                        family= poisson)
summary(pois_male_scaled)

step_model <- stepAIC(pois_male_scaled, direction = "both")
summary(step_model)

#female
mean(datanew$female)
var(datanew$female)

nb_female_scaled <- glm.nb(
  female ~ ndvi_scale + footprint_ + elevation_ + 
    bio_9_scal + bio_3_scal + bio_2_scal +
    bio_19_sca + bio_17_sca + bio_14_sca + nightlight 
  + x + y,
  data = scaled_envlayers
)

summary(nb_female_scaled)

step_model <- stepAIC(nb_female_scaled, direction = "both")
summary(step_model)

#immature
mean(datanew$immature)
var(datanew$immature)
datanew$immature

nb_immature_scaled <- glm.nb(
  immature ~ ndvi_scale + footprint_ + elevation_ + 
    bio_9_scal + bio_3_scal + bio_2_scal +
    bio_19_sca + bio_17_sca + bio_14_sca + nightlight 
  + x + y,
  data = scaled_envlayers
)

summary(nb_immature_scaled)

step_model <- stepAIC(nb_immature_scaled, direction = "both")
summary(step_model)

#offset method - females per male
model_nb <- glm.nb(
  female ~ ndvi_scale + footprint_ + elevation_ +
    bio_9_scal + bio_3_scal + bio_2_scal +
    bio_19_sca + bio_17_sca + bio_14_sca +
    nightlight + x + y +
    offset(log(male)),
  data = scaled_envlayers
)

summary(model_nb)

step_model <- stepAIC(
  model_nb,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#immatures per female
model_nb_if <- glm.nb(
  immature ~ ndvi_scale + footprint_ + elevation_ +
    bio_9_scal + bio_3_scal + bio_2_scal +
    bio_19_sca + bio_17_sca + bio_14_sca +
    nightlight + x + y +
    offset(log(female)),
  data = scaled_envlayers
)

step_model <- stepAIC(
  model_nb_if,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#immatures per adult
scaled_envlayers$adult <- 
  scaled_envlayers$male + scaled_envlayers$female


model_nb_ia <- glm.nb(
  immature ~ ndvi_scale + footprint_ + elevation_ +
    bio_9_scal + bio_3_scal + bio_2_scal +
    bio_19_sca + bio_17_sca + bio_14_sca +
    nightlight + x + y +
    offset(log(adult)),
  data = scaled_envlayers
)

summary(model_nb_ia)

step_model <- stepAIC(
  model_nb_ia,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)



#class level
install.packages("usdm")  
library(usdm)

predictors <- class[, c(
  
  # Rice
  "PLAND_Rice","PD_Rice","LPI_Rice","ED_Rice",
  "PAFRAC_Rice","TCA_Rice","CONNECT_Rice",
  
  # Forest
  "PLAND_Forest","PD_Forest","LPI_Forest","ED_Forest",
  "PAFRAC_Forest","TCA_Forest","CONNECT_Forest",
  
  # Urban
  "PLAND_Urban","PD_Urban","LPI_Urban","ED_Urban",
  "PAFRAC_Urban","TCA_Urban","CONNECT_Urban",
  
  # Teak
  "PLAND_Teak","PD_Teak","LPI_Teak","ED_Teak",
  "PAFRAC_Teak","TCA_Teak","CONNECT_Teak",
  
  # Tea
  "PLAND_Tea","PD_Tea","LPI_Tea","ED_Tea",
  "PAFRAC_Tea","TCA_Tea","CONNECT_Tea",
  
  # Coffee
  "PLAND_Coffee","PD_Coffee","LPI_Coffee","ED_Coffee",
  "PAFRAC_Coffee","TCA_Coffee","CONNECT_Coffee",
  
  # Rubber
  "PLAND_Rubber","PD_Rubber","LPI_Rubber","ED_Rubber",
  "PAFRAC_Rubber","TCA_Rubber","CONNECT_Rubber"
  
)]


tea_vars <- class[, c("PLAND_Tea",
                      "PD_Tea",
                      "LPI_Tea",
                      "ED_Tea",
                      "PAFRAC_Tea",
                      "TCA_Tea",
                      "CONNECT_Tea",
                      "x",
                      "y")]
tea_vars <- tea_vars[, sapply(tea_vars, is.numeric)]
vifcor_tea <- vifcor(tea_vars, th = 0.7)
vifcor_tea
tea_selected <- exclude(tea_vars, vifcor_tea)
colnames(tea_selected)
vifstep_tea <- vifstep(tea_selected, th = 10)
vifstep_tea
tea_final <- exclude(tea_selected, vifstep_tea)
colnames(tea_final)

rice_vars <- class[, c("PLAND_Rice",
                       "PD_Rice",
                       "LPI_Rice",
                       "ED_Rice",
                       "PAFRAC_Rice",
                       "TCA_Rice",
                       "CONNECT_Rice",
                       "x",
                       "y")]

rice_vars <- rice_vars[, sapply(rice_vars, is.numeric)]

vifcor_rice <- vifcor(rice_vars, th = 0.7)
vifcor_rice

rice_selected <- exclude(rice_vars, vifcor_rice)
colnames(rice_selected)

vifstep_rice <- vifstep(rice_selected, th = 10)
vifstep_rice

rice_final <- exclude(rice_selected, vifstep_rice)
colnames(rice_final)

forest_vars <- class[, c("PLAND_Forest",
                         "PD_Forest",
                         "LPI_Forest",
                         "ED_Forest",
                         "PAFRAC_Forest",
                         "TCA_Forest",
                         "CONNECT_Forest",
                         "x",
                         "y")]

forest_vars <- forest_vars[, sapply(forest_vars, is.numeric)]

vifcor_forest <- vifcor(forest_vars, th = 0.7)
vifcor_forest

forest_selected <- exclude(forest_vars, vifcor_forest)
colnames(forest_selected)

vifstep_forest <- vifstep(forest_selected, th = 10)
vifstep_forest

forest_final <- exclude(forest_selected, vifstep_forest)
colnames(forest_final)

urban_vars <- class[, c("PLAND_Urban",
                        "PD_Urban",
                        "LPI_Urban",
                        "ED_Urban",
                        "PAFRAC_Urban",
                        "TCA_Urban",
                        "CONNECT_Urban",
                        "x",
                        "y")]

urban_vars <- urban_vars[, sapply(urban_vars, is.numeric)]

vifcor_urban <- vifcor(urban_vars, th = 0.7)
vifcor_urban

urban_selected <- exclude(urban_vars, vifcor_urban)
colnames(urban_selected)

vifstep_urban <- vifstep(urban_selected, th = 10)
vifstep_urban

urban_final <- exclude(urban_selected, vifstep_urban)
colnames(urban_final)

teak_vars <- class[, c("PLAND_Teak",
                       "PD_Teak",
                       "LPI_Teak",
                       "ED_Teak",
                       "PAFRAC_Teak",
                       "TCA_Teak",
                       "CONNECT_Teak",
                       "x",
                       "y")]

teak_vars <- teak_vars[, sapply(teak_vars, is.numeric)]

vifcor_teak <- vifcor(teak_vars, th = 0.7)
vifcor_teak

teak_selected <- exclude(teak_vars, vifcor_teak)
colnames(teak_selected)

vifstep_teak <- vifstep(teak_selected, th = 10)
vifstep_teak

teak_final <- exclude(teak_selected, vifstep_teak)
colnames(teak_final)

coffee_vars <- class[, c("PLAND_Coffee",
                         "PD_Coffee",
                         "LPI_Coffee",
                         "ED_Coffee",
                         "PAFRAC_Coffee",
                         "TCA_Coffee",
                         "CONNECT_Coffee",
                         "x",
                         "y")]

coffee_vars <- coffee_vars[, sapply(coffee_vars, is.numeric)]

vifcor_coffee <- vifcor(coffee_vars, th = 0.7)
vifcor_coffee

coffee_selected <- exclude(coffee_vars, vifcor_coffee)
colnames(coffee_selected)

vifstep_coffee <- vifstep(coffee_selected, th = 10)
vifstep_coffee

coffee_final <- exclude(coffee_selected, vifstep_coffee)
colnames(coffee_final)

rubber_vars <- class[, c("PLAND_Rubber",
                         "PD_Rubber",
                         "LPI_Rubber",
                         "ED_Rubber",
                         "PAFRAC_Rubber",
                         "TCA_Rubber",
                         "CONNECT_Rubber",
                         "x",
                         "y")]

rubber_vars <- rubber_vars[, sapply(rubber_vars, is.numeric)]

vifcor_rubber <- vifcor(rubber_vars, th = 0.7)
vifcor_rubber

rubber_selected <- exclude(rubber_vars, vifcor_rubber)
colnames(rubber_selected)

vifstep_rubber <- vifstep(rubber_selected, th = 10)
vifstep_rubber

rubber_final <- exclude(rubber_selected, vifstep_rubber)
colnames(rubber_final)


#total
class_rice <-glm.nb(total~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                    +x +y,
                    data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm.nb(total~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                      +x +y,
                      data = class)
step_model_f<- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm.nb(total~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                       +x +y,
                     data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)

class_teak <-glm.nb(total~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                      +x +y,
                    data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm.nb(total~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                     +x +y,
                   data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm.nb(total~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                        +x +y,
                      data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm.nb(total~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                        +x +y,
                      data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)


#males
class_rice <-glm(male~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                 +x +y,
                 family = poisson,
                 data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm(male~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                   +x +y,
                   family = poisson,
                   data = class)
step_model_f <- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm(male~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                    +x +y,
                  family = poisson,
                  data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)

class_teak <-glm(male~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                   +x +y,
                 family = poisson,
                 data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm(male~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                  +x +y,
                family = poisson,
                data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm(male~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                     +x +y,
                   family = poisson,
                   data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm(male~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                     +x +y,
                   family = poisson,
                   data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)

# Extract AIC values
aic_table <- AIC(step_model_rub,
                 step_model_cof,
                 step_model_f,
                 step_model_ri,
                 step_model_tea,
                 step_model_teak,
                 step_model_u)

# Calculate Delta AIC
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)

# Sort by lowest AIC
aic_table <- aic_table[order(aic_table$AIC), ]

aic_table

#females
class_rice <-glm.nb(female~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                    +x +y,
                    data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm.nb(female~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                      +x +y,
                      data = class)
step_model_f<- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm.nb(female~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                       +x +y,
                     data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)


class_teak <-glm.nb(female~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                      +x +y,
                    data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm.nb(female~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                     +x +y,
                   data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm.nb(female~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                        +x +y,
                      data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm.nb(female~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                        +x +y,
                      data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)

# Extract AIC values
aic_table <- AIC(step_model_rub,
                 step_model_cof,
                 step_model_f,
                 step_model_ri,
                 step_model_tea,
                 step_model_teak,
                 step_model_u)

# Calculate Delta AIC
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)

# Sort by lowest AIC
aic_table <- aic_table[order(aic_table$AIC), ]

aic_table

#immature
class_rice <-glm.nb(immature~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                    +x +y,
                    data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm.nb(immature~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                      +x +y,
                      data = class)
step_model_f<- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm.nb(immature~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                       +x +y,
                     data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)


class_teak <-glm.nb(immature~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                      +x +y,
                    data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm.nb(immature~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                     +x +y,
                   data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm.nb(immature~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                        +x +y,
                      data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm.nb(immature~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                        +x +y,
                      data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)

# Extract AIC values
aic_table <- AIC(step_model_rub,
                 step_model_cof,
                 step_model_f,
                 step_model_ri,
                 step_model_tea,
                 step_model_teak,
                 step_model_u)

# Calculate Delta AIC
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)

# Sort by lowest AIC
aic_table <- aic_table[order(aic_table$AIC), ]

aic_table

#females per male
class_rice <-glm.nb(female~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                    +x +y
                    +offset(log(male)),
                    data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm.nb(female~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                      +x +y  +offset(log(male)),
                      data = class)
step_model_f<- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm.nb(female~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                       +x +y +offset(log(male)),
                     data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)


class_teak <-glm.nb(female~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                      +x +y +offset(log(male)),
                    data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm.nb(female~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                     +x +y +offset(log(male)),
                   data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm.nb(female~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                        +x +y +offset(log(male)),
                      data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm.nb(female~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                        +x +y +offset(log(male)),
                      data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)

# Extract AIC values
aic_table <- AIC(step_model_rub,
                 step_model_cof,
                 step_model_f,
                 step_model_ri,
                 step_model_tea,
                 step_model_teak,
                 step_model_u)

# Calculate Delta AIC
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)

# Sort by lowest AIC
aic_table <- aic_table[order(aic_table$AIC), ]

aic_table

#immatures per female
class_rice <-glm.nb(immature~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                    +x +y
                    +offset(log(female)),
                    data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm.nb(immature~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                      +x +y  +offset(log(female)),
                      data = class)
step_model_f<- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm.nb(immature~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                       +x +y +offset(log(female)),
                     data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)


class_teak <-glm.nb(immature~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                      +x +y +offset(log(female)),
                    data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm.nb(immature~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                     +x +y +offset(log(female)),
                   data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm.nb(immature~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                        +x +y +offset(log(female)),
                      data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm.nb(immature~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                        +x +y +offset(log(female)),
                      data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)

# Extract AIC values
aic_table <- AIC(step_model_rub,
                 step_model_cof,
                 step_model_f,
                 step_model_ri,
                 step_model_tea,
                 step_model_teak,
                 step_model_u)

# Calculate Delta AIC
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)

# Sort by lowest AIC
aic_table <- aic_table[order(aic_table$AIC), ]

aic_table

#immatures per adult
class$adult<-class$male+class$female
class_rice <-glm.nb(immature~LPI_Rice + PAFRAC_Rice + ED_Rice + CONNECT_Rice
                    +x +y
                    +offset(log(adult)),
                    data = class)
step_model_ri <- stepAIC(
  class_rice,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_ri)


class_forest <-glm.nb(immature~PD_Forest+ LPI_Forest+ PAFRAC_Forest+ CONNECT_Forest
                      +x +y  +offset(log(adult)),
                      data = class)
step_model_f<- stepAIC(
  class_forest,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_f)


class_urban <-glm.nb(immature~ED_Urban+ LPI_Urban+ PAFRAC_Urban+ CONNECT_Urban +
                       +x +y +offset(log(adult)),
                     data = class)
step_model_u <- stepAIC(
  class_urban,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_u)


class_teak <-glm.nb(immature~ED_Teak+ PD_Teak+PAFRAC_Teak+ CONNECT_Teak +
                      +x +y +offset(log(adult)),
                    data = class)
step_model_teak <- stepAIC(
  class_teak,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_teak)

class_tea <-glm.nb(immature~PD_Tea+ ED_Tea+ TCA_Tea+ PAFRAC_Tea+ CONNECT_Tea+
                     +x +y +offset(log(adult)),
                   data = class)
step_model_tea <- stepAIC(
  class_tea,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_tea)


class_coffee <-glm.nb(immature~LPI_Coffee+ PAFRAC_Coffee+ ED_Coffee + CONNECT_Coffee+
                        +x +y +offset(log(adult)),
                      data = class)
step_model_cof <- stepAIC(
  class_coffee,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_cof)

class_rubber <-glm.nb(immature~LPI_Rubber+ PAFRAC_Rubber + CONNECT_Rubber+ PD_Rubber + PLAND_Rubber+
                        +x +y +offset(log(adult)),
                      data = class)
step_model_rub <- stepAIC(
  class_rubber,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model_rub)

AIC(step_model_rub,step_model_cof, step_model_f, step_model_ri, 
    step_model_tea, step_model_teak, step_model_u)

# Extract AIC values
aic_table <- AIC(step_model_rub,
                 step_model_cof,
                 step_model_f,
                 step_model_ri,
                 step_model_tea,
                 step_model_teak,
                 step_model_u)

# Calculate Delta AIC
aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC)

# Sort by lowest AIC
aic_table <- aic_table[order(aic_table$AIC), ]

aic_table

#EUCLIDEAN DISTANCE
#euclidean distance
library(ggplot2)
euclidean<-read.csv("D:/Demography paper/R/euclidean10m/bonnetpointsfinal.csv")

#total
library(MASS)
nb_model <- glm.nb(
  total ~ EucDisturb +	EucDisttea	+ EucDistthe +	EucDistrub +
    EucDistpad +	EucDistfor +	EucDistcof +
    x + y,
  data = euclidean
)

summary(nb_model)

step_model <- stepAIC(
  nb_model,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#male
pois_male_scaled <- glm(male ~  EucDisturb +	EucDisttea	+ EucDistthe +
                          EucDistrub + EucDistpad +	EucDistfor +	
                          EucDistcof + x + y,
                        data = euclidean, 
                        family= poisson)
summary(pois_male_scaled)

step_model <- stepAIC(
  pois_male_scaled,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#female
nb_model_female <- glm.nb(
  female ~ EucDisturb +	EucDisttea	+ EucDistthe +	EucDistrub +
    EucDistpad +	EucDistfor +	EucDistcof +
    x + y,
  data = euclidean
)

summary(nb_model_female)

step_model <- stepAIC(
  nb_model_female,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#immature
nb_model_immature <- glm.nb(
  immature ~ EucDisturb +	EucDisttea	+ EucDistthe +	EucDistrub +
    EucDistpad +	EucDistfor +	EucDistcof +
    x + y,
  data = euclidean
)

summary(nb_model_immature)

step_model <- stepAIC(
  nb_model_immature,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#females per male
model_mf_offset<-glm.nb(female ~ EucDisturb +	EucDisttea	+ EucDistthe +
                          EucDistrub +  EucDistpad +	EucDistfor +	
                          EucDistcof + x + y +offset(log(male)),
                        data = euclidean)

summary(model_mf_offset)
step_model <- stepAIC(
  model_mf_offset,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

#immatures per female
model_fi_offset<-glm.nb(immature ~ EucDisturb +	EucDisttea	+ EucDistthe +
                          EucDistrub +  EucDistpad +	EucDistfor +	
                          EucDistcof + x + y +offset(log(female)),
                        data = euclidean)

summary(model_fi_offset)
step_model <- stepAIC(
  model_fi_offset,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)

euclidean$adult <- euclidean$male + euclidean$female
model_ai_offset<-glm.nb(immature ~ EucDisturb +	EucDisttea	+ EucDistthe +
                          EucDistrub +  EucDistpad +	EucDistfor +	
                          EucDistcof + x + y +offset(log(adult)),
                        data = euclidean)

summary(model_ai_offset)

step_model <- stepAIC(
  model_ai_offset,
  direction = "both",   # backward + forward
  trace = TRUE
)

summary(step_model)
