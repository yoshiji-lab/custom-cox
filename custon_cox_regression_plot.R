# toy_km_targetpattern_100y_Q5_adj.R
# Q1–Q5 KM over 100y + Cox adjusted for age, sex, centre (with adjusted curves)
pkgs <- c("survival","survminer","dplyr","ggplot2")
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)

library(survival); library(survminer); library(dplyr); library(ggplot2)

pal5 <- c("#3C5488FF", "#00A087FF", "#4DBBD5FF", "#F39B7FFF", "#DC0000FF")  # Q1..Q5
set.seed(42)

H <- 100
n_per_q <- 5000
lev <- paste0("Q", 1:5)

# output figure prefs
dir.create("output", FALSE, TRUE)
w <- 7; h <- 3.5; dpi <- 300

# Simulate covariates and quintiles
N <- n_per_q * 5
df <- tibble(
  quintile = factor(rep(lev, each = n_per_q), levels = lev),
  age      = runif(N, 40, 80),
  sex      = factor(rbinom(N, 1, 0.5), labels = c("female","male")),
  centre   = factor(sample(as.character(1:5), N, replace = TRUE), levels = as.character(1:5))
)

# Target cumulative incidence at 100y by quintile
target_inc_100y <- c(Q1=.02, Q2=.05, Q3=.08, Q4=.11, Q5=.15)
lambda_q <- -log(1 - unname(target_inc_100y)) / H

# Modest covariate effects (centered so per-quintile targets are preserved on average)
lp <- 0.03 * (df$age - mean(df$age)) +
  0.25 * (df$sex == "male") +
  0.05 * (as.integer(df$centre) - mean(as.integer(df$centre)))

rate <- lambda_q[as.integer(df$quintile)] * exp(lp)
rate <- pmax(rate, 1e-8)

T <- rexp(N, rate = rate)   # event time (years)
C <- rep(H, N)              # admin censor at 100y
df <- df %>% mutate(time_years = pmin(T, C), event = as.integer(T <= C))

# ----- Unadjusted KM (kept for RDS bundle) -----
fit_km <- survfit(Surv(time_years, event) ~ quintile, data = df)

# ----- Cox PH adjusted for age, sex, centre -----
df$quintile <- relevel(df$quintile, "Q1")
fit_cox <- coxph(Surv(time_years, event) ~ quintile + age + sex + centre, data = df)

# Adjusted curves at mean covariates (robust alternative to ggadjustedcurves)
nd <- expand.grid(
  quintile = factor(lev, levels = lev),
  age      = mean(df$age),
  sex      = factor("female", levels = levels(df$sex)),
  centre   = factor("3", levels = levels(df$centre))
)
fit_adj <- survfit(fit_cox, newdata = nd)

# Adjusted plot
adj <- ggsurvplot(
  fit_adj, data = nd, fun = "event",
  conf.int = TRUE, conf.int.style = "ribbon", conf.int.alpha = 0.15,
  censor = FALSE, risk.table = FALSE,
  palette = pal5, ggtheme = theme_classic(base_size = 11),
  xlab = "Year", ylab = "Adjusted cumulative incidence",
  legend.title = "", legend.labs = lev,
  xlim = c(0, H), break.time.by = 20
)
ggsave("output/betacell1_type2diabetes.pdf", adj$plot, width = w, height = h, dpi = dpi)

# HR table (optional)
su <- summary(fit_cox); ci <- su$conf.int
hr <- data.frame(term = rownames(ci),
                 HR   = ci[, "exp(coef)"],
                 LCI  = ci[, "lower .95"],
                 UCI  = ci[, "upper .95"],
                 P    = su$coef[, "Pr(>|z|)"])
write.csv(hr, "output/betacell1_type2diabetes.csv", row.names = FALSE)

# ----- Save everything needed for later re-plotting -----
bundle <- list(
  data = df,
  fit_km = fit_km,
  fit_cox = fit_cox,
  fit_adj = fit_adj,
  newdata = nd,
  palette = pal5,
  legend_labels = lev,
  horizon_years = H,
  km_df = survminer::surv_summary(fit_km, data = df),
  adj_df = survminer::surv_summary(fit_adj, data = nd),
  rec_width = w, rec_height = h, rec_dpi = dpi
)
saveRDS(bundle, "output/betacell1_type2diabetes.rds")

# #Example re-plot later (no refitting):
# b <- readRDS("output/betacell1_type2diabetes.rds")
# p <- ggsurvplot(
#   b$fit_adj, data = b$newdata, fun = "event",
#   conf.int = TRUE, conf.int.style = "ribbon", conf.int.alpha = 0.15,
#   censor = FALSE, palette = b$palette, ggtheme = theme_classic(base_size = 11),
#   xlab = "Year", ylab = "Adjusted cumulative incidence",
#   legend.title = "", legend.labs = b$legend_labels,
#   xlim = c(0, b$horizon_years), break.time.by = 20
# )
# ggsave("output/betacell1_type2diabetes.pdf", p$plot, width = b$rec_width, height = b$rec_height, dpi = b$rec_dpi)
