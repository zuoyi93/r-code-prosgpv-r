# -----------------
# REPLICATION CODE
# -----------------

library("ProSGPV")
library("BeSS")
library("SIS")
library("ggplot2")
library("ggpubr")

load("data/sim.data.RData")
source("utils.R")

# -----------------
# FIGURE 1
# -----------------

png("Figures/Figure 1.png", units = "in", width = 12, height = 8, res = 300)
get.fig.1()
dev.off()

# --------------------------
# SECTION: PACKAGE OVERVIEW
# --------------------------

library(ProSGPV)

set.seed(1)
sim.data <- gen.sim.data(
  n = 100, p = 10, s = 4, family = "gaussian",
  beta.min = 1, beta.max = 5, rho = 0.2, nu = 2
)

x <- sim.data[[1]]
y <- sim.data[[2]]
(true.index <- sim.data[[3]])
true.beta <- sim.data[[4]]

sgpv.out.2 <- pro.sgpv(x, y)
sgpv.out.2

png("Figures/Figure 2.png", units = "in", width = 9, height = 9, res = 300)
plot(sgpv.out.2)
dev.off()

summary(sgpv.out.2)
beta.hat <- coef(sgpv.out.2)
rbind(beta.hat, true.beta)

predict(sgpv.out.2)

sgpv.out.1 <- pro.sgpv(x, y, stage = 1)
sgpv.out.1

png("Figures/Figure 3.png", units = "in", width = 9, height = 9, res = 300)
plot(sgpv.out.1)
dev.off()

# --------------------------
# SIMULATION RESULTS
# --------------------------

# -----------------------------------------------------------------
# FIGURE 4 - suppport recovery + parameter estimation + prediction
# -----------------------------------------------------------------

sr.linear <- get.plot.sr(
  data = out.main.gaussian, n = 100, p = 100,
  title.p = "Linear regression", type = "hd"
)

sr.logistic <- get.plot.sr(
  data = out.main.logistic, n = 20, p = 100,
  title.p = "Logistic regression", type = "ld"
)

sr.poisson <- get.plot.sr(
  data = out.main.poisson, n = 40, p = 40,
  title.p = "Poisson regression", type = "hd"
)

sr.cox <- get.plot.sr(
  data = out.main.cox, n = 40, p = 40,
  title.p = "Cox regression", type = "ld"
)

sr.all <- ggarrange(sr.linear, sr.logistic, sr.poisson, sr.cox,
  ncol = 4, nrow = 1, legend = "none"
)


pe.linear <- get.plot.pe(
  data = out.main.gaussian, n = 100, p = 100,
  ylim = c(0, 0.08), ybreaks = seq(0, 0.08, 0.02), type = "hd",
  title.p = "Linear regression"
)

pe.logistic <- get.plot.pe(
  data = out.main.logistic, n = 20, p = 100, cap = 0.8,
  ylim = c(0, 0.8), ybreaks = seq(0, 0.8, 0.2), type = "ld",
  title.p = "Logistic regression"
)

pe.poisson <- get.plot.pe(
  data = out.main.poisson, n = 40, p = 40, type = "hd",
  ylim = c(0, 0.015), ybreaks = seq(0, 0.015, 0.003),
  title.p = "Poisson regression"
)

pe.cox <- get.plot.pe(
  data = out.main.cox, n = 40, p = 40,
  ylim = c(0, 0.8), ybreaks = seq(0, 0.8, 0.2), type = "ld",
  title.p = "Cox regression"
)

pe.all <- ggarrange(pe.linear, pe.logistic, pe.poisson, pe.cox,
  ncol = 4, nrow = 1, legend="none" 
)

pr.linear <- get.plot.pr(
  data = out.main.gaussian, n = 100, p = 100,
  ylab = "Prediction RMSE", ylim = c(3, 7),
  ybreaks = 3:7,
  title.p = "Linear regression", type = "hd"
)

pr.logistic <- get.plot.pr(
  data = out.main.logistic, n = 20, p = 100,
  ylab = "Prediction AUC", ylim = c(0.5, 1),
  ybreaks = seq(0.5, 1, 0.1),
  title.p = "Logistic regression", type = "ld"
)

pr.poisson <- get.plot.pr(
  data = out.main.poisson, n = 40, p = 40,
  ylab = "Prediction RMSE", ylim = c(0, 40),
  ybreaks = seq(0, 40, 10),
  title.p = "Poisson regression", type = "hd"
)

pr.all <- ggarrange(pr.linear, pr.logistic, pr.poisson,
  ggplot() + theme_minimal(),
  ncol = 4, legend = "none"
)

rt.linear <- get.plot.rt(
  data = out.run.gaussian, n = 100, p = 100,
  title.p = "Linear regression", np = "p"
)

rt.logistic <- get.plot.rt(
  data = out.run.logistic, n = 20, p = 100,
  title.p = "Logistic regression", np = "n"
)

rt.poisson <- get.plot.rt(
  data = out.run.poisson, n = 40, p = 40,
  title.p = "Poisson regression", np = "p"
)

rt.cox <- get.plot.rt(
  data = out.run.cox, n = 40, p = 40, 
  title.p = "Cox regression", np = "n"
)

rt.all <- ggarrange(rt.linear, rt.logistic, rt.poisson, rt.cox,
                    ncol = 4, legend = "bottom", common.legend = T
)


png("Figures/Figure 4.png", units = "in", width = 10, height = 10, res = 300)
ggarrange(sr.all, pe.all, pr.all, rt.all,
  ncol = 1, labels = LETTERS[1:4]
)
dev.off()




