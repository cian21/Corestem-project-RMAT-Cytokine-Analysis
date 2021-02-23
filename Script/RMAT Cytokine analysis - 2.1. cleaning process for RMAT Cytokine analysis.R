# Corestem project - RMAT Cytokine Analysis.proj
# 20210201 Sehwan Chun at Corestem, Inc.
# 2.1. cleaning process for RMAT Cytokine analysis

#### 1. source Loading ####
load("./Data/RMAT Cytokine Analysis files.image")
source("./Script/Functions/RMAT Cytokine analysis - functions.R")

#### 2. remove unavailable samples ####
reg_file_cleaning(reg_file, impute = T, month = 6)

#### TMP ####
results = regsubsets(X0M.baseline. ~., data = GR4_reg_x_values_v5)
summary_results = summary(results)

summary_results
summary_results$bic
summary_results$adjr2

results = regsubsets(progression_rate_change_6m ~., data = reg_x_values_diff_percent)
summary_results = summary(results)

summary_results
summary_results$bic
summary_results$adjr2

results = lm(X0M.baseline. ~ TGFB1_V5 + IL6_V5, data = GR4_reg_x_values_v5)
summary_results = summary(results)
summary_results

results = lm(progression_rate_change_6m ~ TGFB2_diff + IL6_diff, data = reg_x_values_diff_percent)
summary_results = summary(results) 
summary_results

#### plot ####

plot_baseline = ggscatter(reg_x_values_v5,
          x = "TGFB2_V5",
          y = "X0M.baseline.",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "red", fill = "grey30"),
          size = 3,
          xlab = "Baseline TGF-B2 level",
          ylab = "Baseline ALSFRS-R score",
          ggtheme = theme_bw(base_size = 15)) + 
    stat_cor(
    r.digits = 3,
    p.digits = 3,
    label.y = 50
)

plot_prc = ggscatter(reg_x_values_diff_percent,
          x = "TGFB2_diff",
          y = "progression_rate_change_6m",
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "red", fill = "grey30"),
          size = 3,
          xlab = "TGF-B2 level changes(%)",
          ylab = "Progression rate change(6M)",
          ylim = c(0,6),
          ggtheme = theme_bw(base_size = 15)) +
    stat_cor(
    r.digits = 3,
    p.digits = 3,
    label.y = 6
)
plot_baseline
plot_prc

ggsave(filename = "./Output/TGFB2_ALSFRSR_nimputed_baseline.tiff", plot_baseline,
       width = 8, height = 8, units = "in")
ggsave(filename = "./Output/TGFB2_ALSFRSR_nimputed_prc.tiff", plot_prc,
       width = 8, height = 8, units = "in")


ggsave(filename = "./Output/TGFB2_ALSFRSR_nimputed.tiff", ggarrange(plot_baseline, plot_prc, nrow = 1, ncol = 2),
       width = 15, height = 8, units = "in")
ggsave(filename = "./Output/TGFB2_ALSFRSR_imputed.tiff", ggarrange(plot_baseline, plot_prc, nrow = 1, ncol = 2),
       width = 15, height = 8, units = "in")
