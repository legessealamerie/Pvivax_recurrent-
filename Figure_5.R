
##########################################################

############### Author: Dr. Jordache Ramjith, Radboud University Medical Center


#############################################################


library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(mgcv)
library(pROC)
library(boot)

cid <- read_excel("CID.xlsx")
ethi <- read_excel("data_adama.xlsx")
indie <- read_excel("data_recurrent.xlsx")
indie$date = as.Date(indie$date)

dim(indie)

cols = colnames(indie)

indiesub <- na.omit(indie %>% select(id,age, gender, sick_malaria_1_year_bl,sick_malaria_3_months_bl))

indie = indie %>% 
  select(-age, -gender, -sick_malaria_1_year_bl,-sick_malaria_3_months_bl) %>% left_join(indiesub, by="id")

indie = indie[, cols]


combdt <- bind_rows(cid %>% mutate(data="cid"),
                    ethi %>% mutate(data="ethiopia"),
                    indie %>% mutate(data="Indie"))

dim(combdt)

combdt$pvs25_copies_ul = ifelse(is.na(combdt$pvs25_copies_ul) & combdt$pvcopy_ul==0,0,combdt$pvs25_copies_ul)

# asym effect ------------------------------------------------------------
plot_asymdt = combdt %>% 
  #filter(mic_species %in% c("P. vivax","pv","Pv")) %>% 
  mutate(data=factor(data), fever=factor(`fever_>=37.5 °C`), asym=factor(ifelse(!(type_infection %in% c("symptomatic","symptomatic recruitment","symptomatic recurrence", "symptomatic persisting")) & !is.na(type_infection),"asymptomatic","symptomatic"))) %>%
  filter(!is.na(pvcopy_ul) & !is.na(percent_inf))

dim(plot_asymdt)


# Calculate sample size for each group
sample_sizes <- plot_asymdt %>% 
  group_by(data,asym) %>% 
  summarise(n = n()) %>%
  arrange(data)


sample_sizes$data = factor(sample_sizes$data, levels=levels(plot_asymdt$data))
sample_sizes


plot_asymdt$logpcr = log10(plot_asymdt$pvcopy_ul+0.01)
plot_asymdt$gender = factor(plot_asymdt$gender)



mod1 = gam(
  cbind(number_inf, number_dissect - number_inf) ~ logpcr + 
    s(data,bs="re"), family = binomial(), data = plot_asymdt)
summary(mod1)
plot(mod1)

newd = expand.grid(logpcr = seq(min(plot_asymdt$logpcr),max(plot_asymdt$logpcr), length.out=10000),data=factor(c("Indie","ethiopia","cid")))

newd$pvcopy_ul = 10^newd$logpcr
# Step 2: Generate predictions with standard errors
predicted <- predict(mod1, newdata = newd, type = "link", se.fit = TRUE, exclude=c("s(data)"))

fam <- family(mod1)
ilink <- fam$linkinv

newd$predicted_values <- ilink(predicted$fit)
newd$se_upper <- ilink(predicted$fit + predicted$se.fit)  # 95% CI upper bound
newd$se_lower <- ilink(predicted$fit - predicted$se.fit ) # 95% CI lower bound)

library(scales)

plot_asym=ggplot(data=plot_asymdt, aes(x=pvcopy_ul, y =percent_inf/100))+
  geom_point(size=2,alpha=0.5)+
  theme_prism()+
  scale_y_continuous(labels=scales::percent)+
  #scale_x_log10(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  scale_x_log10(limits=c(0.1, 10000000), breaks=c(-1,1, 10, 100, 1000, 10000, 100000, 1000000,10000000),labels = trans_format("log10", math_format(10^.x))) +
  
  xlab("Pv18S copies/uL")+
  ylab("Proportion of infected mosquitoes, %")+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "top")+
  geom_line(data=newd,aes(y = predicted_values), size = 1) +
  geom_ribbon(data=newd,aes(y = predicted_values,ymin = se_lower, ymax = se_upper), alpha = 0.1, color = NA) +
  coord_cartesian(ylim=c(0,1))+
  stat_cor(aes(x=pvcopy_ul+0.1, y =percent_inf/100),method = "spearman", p.digits = 4, p.accuracy = 0.0001, cor.coef.name = "rho") 

plot_asym


plot_asym <- ggplot(data=plot_asymdt, aes(x=pvcopy_ul, y=percent_inf/100)) +
  geom_point(size=2, alpha=0.5) +
  theme_prism() +
  scale_y_continuous(labels=scales::number) +
  #scale_x_log10(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
  scale_x_log10(limits=c(0.1, 10000000), breaks=c(-1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels=trans_format("log10", math_format(10^.x))) +
  xlab("Pv18S copies/uL") +
  ylab("Proportion of infected mosquitoes") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top") +
  geom_line(data=newd, aes(y=predicted_values), size=1) +
  geom_ribbon(data=newd, aes(y=predicted_values, ymin=se_lower, ymax=se_upper), alpha=0.1, color=NA) +
  coord_cartesian(ylim=c(0, 1)) +
  stat_cor(aes(x=pvcopy_ul + 0.1, y=percent_inf/100), method="spearman", p.digits=4, p.accuracy=0.0001, cor.coef.name="rho")

plot_asym

options(scipen =1000)
plot_asym <- ggplot(data=plot_asymdt, aes(x=pvcopy_ul, y=percent_inf)) +
  geom_point(size=2, alpha=0.5) +
  theme_prism() +
  scale_y_continuous(labels=scales::number, breaks=seq(0, 100, by=25)) +
  scale_x_log10(limits=c(0.1, 10000000), breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels=trans_format("log10", math_format(10^.x))) +
  xlab("Pv18S copies/uL") +
  ylab("Proportion of infected mosquitoes (%)") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top") +
  geom_line(data=newd, aes(y=predicted_values * 100), size=1) +
  geom_ribbon(data=newd, aes(y=predicted_values * 100, ymin=se_lower * 100, ymax=se_upper * 100), alpha=0.1, color=NA) +
  coord_cartesian(ylim=c(0, 100)) +
  stat_cor(aes(x=pvcopy_ul + 0.1, y=percent_inf), method="spearman", p.digits=4, p.accuracy=0.0001, cor.coef.name="rho")

plot_asym

plot_asym <- ggplot(data=plot_asymdt, aes(x=pvcopy_ul, y=percent_inf)) +
  geom_point(size=2, alpha=0.5) +
  theme_prism() +
  scale_y_continuous(labels=scales::number, breaks=seq(0, 100, by=25)) +
  scale_x_log10(limits=c(0.1, 10000000), breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), labels=trans_format("log10", math_format(10^.x))) +
  xlab("Pv18S copies/uL") +
  ylab("Proportion of infected mosquitoes (%)") +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme(legend.position="top") +
  geom_line(data=newd, aes(y=predicted_values * 100), size=1) +
  geom_ribbon(data=newd, aes(y=predicted_values * 100, ymin=se_lower * 100, ymax=se_upper * 100), alpha=0.1, color=NA) +
  coord_cartesian(ylim=c(0, 100)) +
  stat_cor(aes(x=pvcopy_ul + 0.1, y=percent_inf), method="spearman", p.digits=4, p.accuracy=0.0001, cor.coef.name="rho",size=4)

plot_asym



ggsave("Pcr_pardens_infe.pdf", plot = plot_asym, device = "pdf", dpi = 1200, width = 7, height = 5)

mod1 = gam(
  cbind(number_inf, number_dissect - number_inf) ~ logpcr*asym + 
    s(asym,data,bs="re"), family = binomial(), data = plot_asymdt)
summary(mod1)
plot(mod1)

newd = expand.grid(logpcr = seq(min(plot_asymdt$logpcr),max(plot_asymdt$logpcr), length.out=10000),data=factor(c("Indie","ethiopia","cid")),asym=levels(plot_asymdt$asym))

newd$pvcopy_ul = 10^newd$logpcr
# Step 2: Generate predictions with standard errors
predicted <- predict(mod1, newdata = newd, type = "link", se.fit = TRUE, exclude=c("s(asym,data)"))

fam <- family(mod1)
ilink <- fam$linkinv

newd$predicted_values <- ilink(predicted$fit)
newd$se_upper <- ilink(predicted$fit + predicted$se.fit)  # 95% CI upper bound
newd$se_lower <- ilink(predicted$fit - predicted$se.fit ) # 95% CI lower bound)



plot_asym=ggplot(data=plot_asymdt, aes(x=pvcopy_ul, y =percent_inf/100,col=factor(asym)))+
  geom_point(size=2,alpha=0.5)+
  theme_prism()+
  scale_y_continuous(labels=scales::percent)+
  scale_x_log10(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
  xlab("pvcopy_ul")+
  ylab("percent_inf")+
  scale_fill_brewer(palette = "Set1")+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "top")+
  geom_line(data=newd,aes(y = predicted_values, col=factor(asym)), size = 1) +
  geom_ribbon(data=newd,aes(y = predicted_values,ymin = se_lower, ymax = se_upper, fill = factor(asym)), alpha = 0.1, color = NA) +
  coord_cartesian(ylim=c(0,1))+
  stat_cor(aes(x=pvcopy_ul+0.1, y =percent_inf/100,col=factor(asym)),method = "spearman", p.digits = 4, p.accuracy = 0.0001, cor.coef.name = "rho") 

plot_asym


# Save plot
#ggsave("Pcr_pardens_infe_asym_comparison.pdf", plot = plot_asym, device = "pdf", dpi = 1200, width = 7, height = 5)






# Refine a prediction model -----------------------------------------------


mod1 = gam(
  cbind(number_inf, number_dissect - number_inf) ~  
    gender+logpcr+
    s(age,k=5)+
    s(data,bs="re"), family = binomial(), data = plot_asymdt,
  select = T)
summary(mod1)



# Predictive power  -----------------------------------------------

###### ALL data ######


# Expand data for binary format
plot_data <- plot_asymdt %>%
  # Select relevant columns
  select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
  na.omit() %>%
  # Create a list of rows for each original row, expanding based on number_dissect
  rowwise() %>%
  do(data.frame(
    logpcr = .$logpcr,
    asym = .$asym,
    gender = .$gender,
    age = .$age,
    data = .$data,
    outcome = c(rep(1, .$number_inf), rep(0, .$number_dissect - .$number_inf))
  )) %>%
  ungroup()



plot_data$predicted <- predict(mod1, newdata=plot_data, type = "response")

roc_obj <- roc(plot_data$outcome, plot_data$predicted)
plot(roc_obj, main="ROC Curve for Infectivity Model")

auc_value <- auc(roc_obj)  # Extract AUC

legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), box.lty = 1)



###### Asymptomatics ######

# Expand data for binary format
plot_data <- plot_asymdt %>%
  # Select relevant columns
  select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
  na.omit() %>%
  # Create a list of rows for each original row, expanding based on number_dissect
  rowwise() %>%
  do(data.frame(
    logpcr = .$logpcr,
    asym = .$asym,
    gender = .$gender,
    age = .$age,
    data = .$data,
    outcome = c(rep(1, .$number_inf), rep(0, .$number_dissect - .$number_inf))
  )) %>%
  ungroup() %>% 
  filter(asym=="asymptomatic")



plot_data$predicted <- predict(mod1, newdata=plot_data, type = "response")

roc_obj <- roc(plot_data$outcome, plot_data$predicted)
plot(roc_obj, main="ROC Curve for Infectivity Model")

auc_value <- auc(roc_obj)  # Extract AUC

legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), box.lty = 1)



# Calibration -------------------------------------------------------------

###### ALL data ######


# Calculate observed proportions in the original data (not expanded)
plot_data <- plot_asymdt %>%
  #select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
  #na.omit() %>%
  mutate(
    observed_prop = number_inf / number_dissect,
    predicted = predict(mod1, newdata = ., type = "response")
  )


summary(lm(observed_prop ~ predicted, data=plot_data))

# Calibration plot
ggplot(plot_data, aes(x = predicted, y = observed_prop)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Mean Predicted Probability", y = "Mean Observed Proportion") +
  ggtitle("Calibration plot for all observations") +
  theme_minimal()


set.seed(123)  # For reproducibility
n_boot <- 1000  # Number of bootstrap samples
n_rows <- nrow(plot_asymdt)  # Get the number of rows in plot_asymdt

bootstrap_results <- data.frame(intercept = numeric(n_boot), slope = numeric(n_boot))

for (i in 1:n_boot) {
  # Resample the data
  boot_sample <- plot_asymdt %>%
    slice_sample(n = n_rows, replace = TRUE) %>%
    mutate(
      observed_prop = number_inf / number_dissect,
      predicted = predict(mod1, newdata = ., type = "response")
    )
  

  # Fit a linear regression to get intercept and slope for calibration line
  calibration_model <- lm(observed_prop ~ predicted, data = boot_sample)
  bootstrap_results$intercept[i] <- coef(calibration_model)[1]
  bootstrap_results$slope[i] <- coef(calibration_model)[2]
  
}

# Summary statistics for intercept and slope
summary(bootstrap_results)

###### asymptomatics #####

# Calculate observed proportions in the original data (not expanded)
plot_data <- plot_asymdt %>%
  select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
  na.omit() %>%
  mutate(
    observed_prop = number_inf / number_dissect,
    predicted = predict(mod1, newdata = ., type = "response")
  )%>%
filter(asym=="asymptomatic")


m1 = lm(observed_prop ~ I(predicted), data=plot_data)
summary(m1)



# Calibration plot
ggplot(plot_data, aes(x = predicted, y = observed_prop)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Predicted proportion", y = "Observed proportion") +
  ggtitle("Calibration Plot for asymptomatics") +
  theme_minimal()


set.seed(123)  # For reproducibility
n_boot <- 1000  # Number of bootstrap samples
n_rows <- nrow(plot_asymdt)  # Get the number of rows in plot_asymdt

bootstrap_results <- data.frame(intercept = numeric(n_boot), slope = numeric(n_boot))

for (i in 1:n_boot) {
  # Resample the data
  boot_sample <- plot_asymdt %>%
    filter(asym=="asymptomatic")%>%
    slice_sample(n = n_rows, replace = TRUE) %>%
    mutate(
      observed_prop = number_inf / number_dissect,
      predicted = predict(mod1, newdata = ., type = "response")
    )
  
  
  # Fit a linear regression to get intercept and slope for calibration line
  calibration_model <- lm(observed_prop ~ predicted, data = boot_sample)
  bootstrap_results$intercept[i] <- coef(calibration_model)[1]
  bootstrap_results$slope[i] <- coef(calibration_model)[2]
  
}

# Summary statistics for intercept and slope
summary(bootstrap_results)


#  VALIDATION -------------------------------------------------------

###### ALL data ######

set.seed(050686)  # Ensure reproducibility

n_boot <- 1000  # Number of bootstrap iterations
bootstrap_results <- data.frame(
  binary_brier = numeric(n_boot),
  continuous_brier = numeric(n_boot),
  r2 = numeric(n_boot),
  auc = numeric(n_boot),
  intercept = numeric(n_boot),
  slope = numeric(n_boot)
)

bootstrap_results_asym <- data.frame(
  binary_brier = numeric(n_boot),
  continuous_brier = numeric(n_boot),
  r2 = numeric(n_boot),
  auc = numeric(n_boot),
  intercept = numeric(n_boot),
  slope = numeric(n_boot)
)

for (i in 1:n_boot) {
  # Step 1: Bootstrap sample from the original dataset
  boot_sample <- plot_asymdt %>% sample_frac(1, replace = TRUE)
  
  # Step 2: Re-fit the model on the bootstrap sample
  boot_model <- gam(
    cbind(number_inf, number_dissect - number_inf) ~  
      gender + logpcr + 
      s(age, k=5) + s(data, bs = "re"),
    family = binomial(), data = boot_sample
  )
  
  # Step 3: Predict on the original (complete) dataset
  original_data <- plot_asymdt %>%
    select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
    mutate(
      observed_prop = number_inf / number_dissect,
      predicted = predict(boot_model, newdata = ., type = "response")
    ) %>% na.omit()  # Ensure complete cases
  
  # Step 4a: Calculate Brier scores
  # Binary Brier Score using expanded binary data
  expanded_data <- plot_asymdt %>%
    select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
    na.omit() %>%
    rowwise() %>%
    do(data.frame(
      logpcr = .$logpcr,
      asym = .$asym,
      gender = .$gender,
      age = .$age,
      data = .$data,
      outcome = c(rep(1, .$number_inf), rep(0, .$number_dissect - .$number_inf))
    )) %>%
    ungroup() %>%
    mutate(predicted = predict(boot_model, newdata = ., type = "response"))
  
  binary_brier_boot <- mean((expanded_data$outcome - expanded_data$predicted)^2)
  
  # Continuous Brier Score using original data proportions
  continuous_brier_boot <- mean((original_data$observed_prop - original_data$predicted)^2)
  
  # Step 4b: Calculate R² from the original data predictions
  r2_boot <- summary(lm(observed_prop ~ predicted, data = original_data))$r.squared
  
  # Step 4c: Calculate AUC on expanded data
  auc_boot <- suppressMessages(auc(roc(expanded_data$outcome, expanded_data$predicted)))[[1]]
  
  # Step 4d: Calculate calibration intercept and slope
  calibration_model <- lm(observed_prop ~ predicted, data = original_data)
  intercept_boot <- coef(calibration_model)[1]
  slope_boot <- coef(calibration_model)[2]
  
  # Store results
  bootstrap_results[i, ] <- c(binary_brier_boot, continuous_brier_boot, r2_boot, auc_boot, intercept_boot, slope_boot)
  
  ###### asymptomatics #######
  # Step 3: Predict on the original (complete) dataset for asymptomatics
  original_data <- plot_asymdt %>%
    select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
    filter(asym=="asymptomatic") %>%
    mutate(
      observed_prop = number_inf / number_dissect,
      predicted = predict(boot_model, newdata = ., type = "response")
    ) %>% na.omit()  # Ensure complete cases
  
  # Step 4a: Calculate Brier scores
  # Binary Brier Score using expanded binary data
  expanded_data <- plot_asymdt %>%
    filter(asym=="asymptomatic") %>%
    select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
    na.omit() %>%
    rowwise() %>%
    do(data.frame(
      logpcr = .$logpcr,
      asym = .$asym,
      gender = .$gender,
      age = .$age,
      data = .$data,
      outcome = c(rep(1, .$number_inf), rep(0, .$number_dissect - .$number_inf))
    )) %>%
    ungroup() %>%
    mutate(predicted = predict(boot_model, newdata = ., type = "response"))
  
  binary_brier_boot <- mean((expanded_data$outcome - expanded_data$predicted)^2)
  
  # Continuous Brier Score using original data proportions
  continuous_brier_boot <- mean((original_data$observed_prop - original_data$predicted)^2)
  
  # Step 4b: Calculate R² from the original data predictions
  r2_boot <- summary(lm(observed_prop ~ predicted, data = original_data))$r.squared
  
  # Step 4c: Calculate AUC on expanded data
  auc_boot <- suppressMessages(auc(roc(expanded_data$outcome, expanded_data$predicted)))[[1]]
  
  # Step 4d: Calculate calibration intercept and slope for asymptomatics
  calibration_model <- lm(observed_prop ~ predicted, data = original_data)
  intercept_boot <- coef(calibration_model)[1]
  slope_boot <- coef(calibration_model)[2]
  
  # Store results
  bootstrap_results_asym[i, ] <- c(binary_brier_boot, continuous_brier_boot, r2_boot, auc_boot, intercept_boot, slope_boot)
  
}

# Output bootstrap summary statistics

# Summary for all data
bootstrap_summary_all <- bootstrap_results %>%
  summarise(
    binary_brier_median = median(binary_brier),
    binary_brier_lower = quantile(binary_brier, 0.025),
    binary_brier_upper = quantile(binary_brier, 0.975),
    continuous_brier_median = median(continuous_brier),
    continuous_brier_lower = quantile(continuous_brier, 0.025),
    continuous_brier_upper = quantile(continuous_brier, 0.975),
    r2_median = median(r2),
    r2_lower = quantile(r2, 0.025),
    r2_upper = quantile(r2, 0.975),
    auc_median = median(auc),
    auc_lower = quantile(auc, 0.025),
    auc_upper = quantile(auc, 0.975),
    intercept_median = median(intercept),
    intercept_lower = quantile(intercept, 0.025),
    intercept_upper = quantile(intercept, 0.975),
    slope_median = median(slope),
    slope_lower = quantile(slope, 0.025),
    slope_upper = quantile(slope, 0.975)
  )

# Summary for asymptomatics
bootstrap_summary_asym <- bootstrap_results_asym %>%
  summarise(
    binary_brier_median = median(binary_brier),
    binary_brier_lower = quantile(binary_brier, 0.025),
    binary_brier_upper = quantile(binary_brier, 0.975),
    continuous_brier_median = median(continuous_brier),
    continuous_brier_lower = quantile(continuous_brier, 0.025),
    continuous_brier_upper = quantile(continuous_brier, 0.975),
    r2_median = median(r2),
    r2_lower = quantile(r2, 0.025),
    r2_upper = quantile(r2, 0.975),
    auc_median = median(auc),
    auc_lower = quantile(auc, 0.025),
    auc_upper = quantile(auc, 0.975),
    intercept_median = median(intercept),
    intercept_lower = quantile(intercept, 0.025),
    intercept_upper = quantile(intercept, 0.975),
    slope_median = median(slope),
    slope_lower = quantile(slope, 0.025),
    slope_upper = quantile(slope, 0.975)
  )

# Display summaries
print(bootstrap_summary_all)
print(bootstrap_summary_asym)



plot_data <- plot_asymdt %>%
  select(number_inf, number_dissect, logpcr, asym, gender, age, data) %>%
  na.omit() %>%
  mutate(
    observed_prop = number_inf / number_dissect,
    predicted = predict(mod1, newdata = ., type = "response")
  )%>%
  filter(asym=="asymptomatic")

ggplot(plot_data, aes(x = predicted, y = observed_prop)) +
  geom_point(size = 2, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_abline(slope = 0.8945741, intercept = 0.007090831, linetype = "solid", color = "blue") +
  labs(x = "Predicted proportion", y = "Observed proportion") +
  ggtitle("Calibration Plot for asymptomatics") +
  theme_minimal()+
  coord_cartesian(ylim=c(0,1),xlim=c(0,1))


newdata = expand.grid(logpcr=c(1,3,5),gender="Male",age=c(10,20,60), data=c("Indie","ethiopia","cid"))
newdata$prop=predict(mod1, newdata=newdata, type="response")
print(newdata)

library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(mgcv)
library(pROC)
library(boot)
cid <- read_excel("CID.xlsx")
ethi <- read_excel("data_adama.xlsx")
indie <- read_excel("data_recurrent.xlsx")
indie$date = as.Date(indie$date)
cols = colnames(indie)

indiesub <- na.omit(indie %>%dplyr:: select(id,age, gender, sick_malaria_1_year_bl,sick_malaria_3_months_bl))

indie = indie %>% 
  dplyr::select(-age, -gender, -sick_malaria_1_year_bl,-sick_malaria_3_months_bl) %>% left_join(indiesub, by="id") %>%
  arrange(id,sno)

indiesum = indie %>% 
  group_by(id) %>%
  summarise(
    anysymrec = sum(I(type_infection=="symptomatic recurrence"),na.rm=T),
    anyasymrec = sum(I(type_infection%in%c("submicroscopic recurrence","asymptomatic recurrence")),na.rm=T),
    persasym = sum(ifelse(type_infection=="symptomatic recruitment" & (lead(type_infection)=="submicroscopic persisting" |lead(type_infection,2)=="submicroscopic persisting" |lead(type_infection,3)=="submicroscopic persisting"),1,0),na.rm=T),
    perssym = sum(ifelse(type_infection=="symptomatic recruitment" & (lead(type_infection)=="symptomatic persisting" ),1,0),na.rm=T)
    
  )

sum(indiesum$anysymrec>=1)
sum(indiesum$anyasymrec>=1)
sum(indiesum$anyasymrec)

sum(indiesum$persasym>=1,na.rm=T)
sum(indiesum$perssym>=1,na.rm=T)


indie = indie %>% 
  group_by(id) %>%
  mutate(type_infection = ifelse(type_infection=="submicroscopic recurrence","asymptomatic recurrence", type_infection),
         type_infection = ifelse(type_infection=="Pf",NA, type_infection),
         type_infection = ifelse(type_infection=="submicroscopic persisting","asymptomatic persisting", type_infection),
         new = ifelse((type_infection %in% c("asymptomatic recurrence","symptomatic recurrence","symptomatic recruitment")), 1,0),
         infid=cumsum(new),
         infected = ifelse(type_infection %in% c("asymptomatic recurrence","symptomatic recurrence","symptomatic recruitment", "asymptomatic persisting","symptomatic persisting"),1,0),
         infected = ifelse(((infid == lag(infid) & infid == lead(infid)) & ((infected==0 | is.na(infected)) & lag(infected)==1 & lead(infected)==1)) | ((infid == lag(infid) & infid == lead(infid,2)) & ((infected==0 | is.na(infected)) & lag(infected)==1 & lead(infected,2)==1))| ((infid == lag(infid,2) & infid == lead(infid)) & ((infected==0 | is.na(infected)) & lag(infected,2)==1 & lead(infected)==1)) ,1,infected),
         infected=ifelse(is.na(infected),0,infected))

indie = indie %>% filter(infected==1) %>% group_by(id, infid) %>% mutate(infected = sum(infected), infcat=first(type_infection)) %>% ungroup()

indie = indie[, cols] %>%
  mutate(infcat = ifelse(type_infection %in% c("symptomatic recruitment"), "Recruitment",ifelse(type_infection %in% c("symptomatic recurrence","symptomatic persisting"),"Symptomatic","Asymptomatic")), 
         time = 3*I(infcat %in% c("Recruitment","Symptomatic"))+14*I(infcat %in% "Asymptomatic"))


combdt <- bind_rows(cid %>% mutate(data="cid"),
                    ethi %>% mutate(data="ethiopia"),
                    indie[,colnames(cid)] %>% mutate(data="Indie"))


# asym effect ------------------------------------------------------------
plot_asymdt = combdt %>% 
  mutate(data=factor(data)) %>%
  filter(!is.na(pvcopy_ul) & !is.na(percent_inf))

table(plot_asymdt$data)

plot_asymdt$logpcr = log10(plot_asymdt$pvcopy_ul+0.001)
plot_asymdt$gender = factor(plot_asymdt$gender)

indie.c = indie
indie.c$logpcr = log10(indie.c$pvcopy_ul+0.001)
indie.c$gender = factor(indie.c$gender)
indie.c$data = factor("Indie", levels=c("cid", "ethiopia", "Indie"))

# indie.c$duration = ifelse(indie.c$infcat %in% c("Recruitment","Symptomatic"),3,
#                           ifelse(indie.c$type_infection=="asymptomatic persisting",28,14))

mod1 = gam(
  cbind(number_inf, number_dissect - number_inf) ~  
    gender+logpcr+
    s(age, k=5)+
    s(data,bs="re"), family = binomial(), data = plot_asymdt,
  select = T)
summary(mod1)

indie.c$prop=predict(mod1, newdata=indie.c, type="response")

for (i in 1:nrow(indie.c)){
  
  if(!is.na(indie.c$pvcopy_ul[i]==0) & indie.c$pvcopy_ul[i]==0){
    indie.c$prop[i]=0
  } else {
    if(!is.na(indie.c$percent_inf[i])) {
      indie.c$prop[i]=indie.c$percent_inf[i]/100
    }
  } 
}



mod2 = gam(prop ~ infcat, data=indie.c, family=binomial())
summary(mod2)
preds = predict(mod2, newdata=data.frame(infcat=c("Recruitment", "Symptomatic", "Asymptomatic")), type="response")

infres = indie.c %>%
  mutate(infcat = factor(infcat, levels=c("Recruitment", "Symptomatic", "Asymptomatic"), labels=c("Recruitment", "Symptomatic moments", "Asymptomatic moments"))) %>%
  group_by(infcat) %>% 
  arrange(infcat) %>%
  summarise(N=n(),duration = sum(time)) 

infres$prop = preds

infres2 = infres %>%
  mutate(duration.share = duration/sum(duration), PI = duration*prop/sum(duration*prop))

infres2$xmin = cumsum(lag(infres2$duration.share,default=0))
infres2$xmax = cumsum(infres2$duration.share)
infres2$xmid = (infres2$xmin+infres2$xmax)/2



infres3 = infres %>%
  mutate(N.share = N/sum(N), PI = N*prop/sum(N*prop))

infres3$xmin = cumsum(lag(infres3$N.share,default=0))
infres3$xmax = cumsum(infres3$N.share)
infres3$xmid = (infres3$xmin+infres3$xmax)/2





# Other calculations ------------------------------------------------------


table(indie.c$infcat)

quantile(indie.c$pvs25_copies_ul[indie.c$type_infection=="symptomatic recruitment"], probs=c(0.5, 0.25, 0.75), na.rm=T)

quantile(indie.c$pvs25_copies_ul[indie.c$type_infection=="symptomatic recurrence"], probs=c(0.5, 0.25, 0.75), na.rm=T)

quantile(indie.c$pvs25_copies_ul[indie.c$type_infection=="asymptomatic recurrence"|indie.c$type_infection=="submicroscopic recurrence"], probs=c(0.5, 0.25, 0.75), na.rm=T)



# Confidence intervals based on predictions using simulations -------------

N=1000

mod1ad = mod1
coefs = MASS::mvrnorm(n=N,mod1$coefficients,vcov(mod1))

PI.asym.rec2 = vector(length=N)
PI.sym.rec2 = vector(length=N)
PI.sym.recruit2 = vector(length=N)

PI.asym.rec2b = vector(length=N)
PI.sym.rec2b = vector(length=N)
PI.sym.recruit2b = vector(length=N)

indie.cis=indie.c

infres.cis = indie.cis %>%
  mutate(infcat = factor(infcat, levels=c("Recruitment", "Symptomatic", "Asymptomatic"), labels=c("Recruitment", "Symptomatic moments", "Asymptomatic moments"))) %>%
  group_by(infcat) %>% 
  arrange(infcat) %>%
  summarise(N=n(),duration = sum(time)) 

for(j in 1:N){
  mod1ad$coefficients = as.vector(coefs[j,])
  
  indie.cis$prop=predict(mod1ad, newdata=indie.cis, type="response")
  
  for (i in 1:nrow(indie.cis)){
    
    if(!is.na(indie.cis$pvcopy_ul[i]==0) & indie.cis$pvcopy_ul[i]==0){
      indie.cis$prop[i]=0
    } else {
      if(!is.na(indie.cis$percent_inf[i])) {
        indie.cis$prop[i]=indie.cis$percent_inf[i]/100
      }
    } 
  }
  
  
  mod2ad = gam(prop ~ infcat, data=indie.cis, family=binomial())
  summary(mod2ad)
  preds = predict(mod2ad, newdata=data.frame(infcat=c("Recruitment", "Symptomatic", "Asymptomatic")), type="response")
  
  infres.cis$prop = preds
  
  infres.cis = infres.cis %>%
    mutate(duration.share = duration/sum(duration), prop.share=prop/sum(prop), PI = duration.share*prop.share/sum(duration.share*prop.share), N.share = N/sum(N), PI2=N.share*prop.share/sum(N.share*prop.share))
  
  PI.sym.recruit2[j]=infres.cis$PI[1]
  PI.sym.rec2[j]=infres.cis$PI[2]
  PI.asym.rec2[j]=infres.cis$PI[3]
  
  PI.sym.recruit2b[j]=infres.cis$PI2[1]
  PI.sym.rec2b[j]=infres.cis$PI2[2]
  PI.asym.rec2b[j]=infres.cis$PI2[3]
  
  print(j)
}


mat2=rbind(quantile(PI.sym.recruit2, probs=c(0.025, 0.975),na.rm=T),
           quantile(PI.sym.rec2, probs=c(0.025, 0.975),na.rm=T),
           quantile(PI.asym.rec2, probs=c(0.025, 0.975),na.rm=T))

mat3=rbind(quantile(PI.sym.recruit2b, probs=c(0.025, 0.975),na.rm=T),
           quantile(PI.sym.rec2b, probs=c(0.025, 0.975),na.rm=T),
           quantile(PI.asym.rec2b, probs=c(0.025, 0.975),na.rm=T))

infres2$PIci = paste0(format(round(infres2$PI*100,1),nsmall=1)," (",format(round(mat2[,1]*100,1),nsmall=1),""," - ",format(round(mat2[,2]*100,1),nsmall=1),")")
infres3$PIci = paste0(format(round(infres3$PI*100,1),nsmall=1)," (",format(round(mat3[,1]*100,1),nsmall=1),""," - ",format(round(mat3[,2]*100,1),nsmall=1),")")


areaplotpi2=ggplot(infres2) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = prop,fill = infcat,color=infcat))+
  scale_y_continuous(labels=scales::percent)+
  scale_y_continuous(labels = scales::number_format(scale = 100), breaks = seq(0, .5, by = 0.1)) +
  
  #scale_x_continuous(labels=scales::percent)+
  scale_x_continuous(labels = scales::number_format(scale = 100), breaks = seq(0, 1, by = 0.25)) +
  
  theme_prism()+
  coord_cartesian(ylim = c(-0.001, 0.5), xlim = c(-0.01,1))+ 
  ggrepel::geom_label_repel(aes(xmid+0.1, prop+0.0235, label = PIci,fill=infcat),color="black",show.legend = FALSE)+
  scale_fill_manual(name="",values=c("#0072B2","#92C5DE","#D55E00"))+ 
  scale_color_manual(name="",values=c("#0072B2","#92C5DE","#D55E00"))+
  theme(legend.position = "top")+ guides(size=FALSE,color=FALSE,label=FALSE)+
  xlab("Share of infected days (%)")+
  ylab("Percent of infected mosquitoes (%)")

areaplotpi2



ggsave("PIplot_moments_w_duration.png",areaplotpi2,width = 11, height = 7, dpi = 300,device='png')


areaplotpi2b=ggplot(infres3) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = prop,fill = infcat,color=infcat))+
  #scale_y_continuous(labels=scales::percent)+
  # scale_x_continuous(labels=scales::percent)+
  
  scale_y_continuous(labels = scales::number_format(scale = 100), breaks = seq(0, .5, by = 0.1)) +
  
  #scale_x_continuous(labels=scales::percent)+
  scale_x_continuous(labels = scales::number_format(scale = 100), breaks = seq(0, 1, by = 0.25)) +
  
  theme_prism()+
  coord_cartesian(ylim = c(-0.001, 0.5), xlim = c(-0.01,1))+ 
  geom_label(aes(xmid, prop+0.02, label = PIci,fill=infcat),color="black",show.legend = FALSE)+
  scale_fill_manual(name="",values=c("#0072B2","#92C5DE","#E69F00"))+ 
  scale_color_manual(name="",values=c("#0072B2","#92C5DE","#E69F00"))+ 
  theme(legend.position = "top")+ guides(size=FALSE,color=FALSE,label=FALSE)+
  xlab("% of infections")+
  ylab("Percentage of mosquitoes infected (%)")

areaplotpi2b


