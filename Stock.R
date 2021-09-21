# -------------
# STOCK MARKET
# -------------

# DESCRIPTION OF DATA IS AVAILABLE ON THE PAGE 27 OF THE PAPER (APPENDIX I)

library(ProSGPV)

dji.d <- read.csv("data/Processed_DJI.csv")
head(dji.d)

dim(dji.d) # 1984 by 84

dji.d <- dji.d[complete.cases(dji.d), ]
dji.d[, c(1, 59)] <- NULL # remove the date

dim(dji.d) # 1114 by 82

# ------------------
# GET TEST DATA
# ------------------

set.seed(1)
test.index <- sample(1:nrow(dji.d), 614)

test.x <- dji.d[test.index, -1]
test.y <- dji.d[test.index, 1]

train.all.index <- setdiff(1:nrow(dji.d), test.index)

# ------------------
# GET RESULTS
# ------------------

# output the size of the selection set and prediction rmse
f.one <- function(train.n = 50) {
  index.use <- sample(train.all.index, train.n)

  x <- dji.d[index.use, -1]
  y <- dji.d[index.use, 1]

  # pro.sgpv unscaled
  sgpv.m <- pro.sgpv(x, y)
  out.sgpv.1 <- length(sgpv.m$var.index)
  out.sgpv.2 <- sqrt(mean((predict(sgpv.m, newdata = test.x) - test.y)^2))

  # bess
  bess.m <- bess(x, y)
  out.bess.1 <- length(which(bess.m$beta != 0))
  out.bess.2 <- sqrt(mean((predict(bess.m, newx = test.x) - test.y)^2))

  # isis
  invisible(capture.output(sis.m <- SIS(as.matrix(x), y, tune = "ebic")))
  out.sis.1 <- length(sis.m$ix)
  out.sis.2 <- sqrt(mean((predict(sis.m, newx = as.matrix(test.x), type = "response") - test.y)^2))

  # lasso
  cv.lasso <- cv.glmnet(as.matrix(x), y)
  out.lasso.1 <- length(which(coef(cv.lasso, s = cv.lasso$lambda.min)[-1] != 0))
  out.lasso.2 <- sqrt(mean((predict(cv.lasso,
    s = cv.lasso$lambda.min,
    newx = as.matrix(test.x)
  ) - test.y)^2))


  return(c(
    out.sgpv.1, out.sgpv.2,
    out.bess.1, out.bess.2,
    out.sis.1, out.sis.2,
    out.lasso.1, out.lasso.2
  ))
}



f.rep <- function(train.n, num.sim = 1e3) {
  temp <- replicate(num.sim, suppressWarnings(suppressMessages(f.one(train.n))))
  return(c(
    median(temp[1, ]), # median model size: sgpv
    quantile(temp[1, ], 0.25), # 1st quartile model size: sgpv
    quantile(temp[1, ], 0.75), # 1st quartile model size: sgpv

    median(temp[2, ]), # median prediction rmse: sgpv
    quantile(temp[2, ], 0.25), # 1st quartile prediction rmse: sgpv
    quantile(temp[2, ], 0.75), # 1st quartile prediction rmse: sgpv

    median(temp[3, ]), # median model size: bess
    quantile(temp[3, ], 0.25), # 1st quartile model size: bess
    quantile(temp[3, ], 0.75), # 1st quartile model size: bess

    median(temp[4, ]), # median prediction rmse: bess
    quantile(temp[4, ], 0.25), # 1st quartile prediction rmse: bess
    quantile(temp[4, ], 0.75), # 1st quartile prediction rmse: bess

    median(temp[5, ]), # median model size: sis
    quantile(temp[5, ], 0.25), # 1st quartile model size: sis
    quantile(temp[5, ], 0.75), # 1st quartile model size: sis

    median(temp[6, ]), # median prediction rmse: sis
    quantile(temp[6, ], 0.25), # 1st quartile prediction rmse: sis
    quantile(temp[6, ], 0.75), # 1st quartile prediction rmse: sis

    median(temp[7, ]), # median model size: lasso
    quantile(temp[7, ], 0.25), # 1st quartile model size: lasso
    quantile(temp[7, ], 0.75), # 1st quartile model size: lasso

    median(temp[8, ]), # median prediction rmse: lasso
    quantile(temp[8, ], 0.25), # 1st quartile prediction rmse: lasso
    quantile(temp[8, ], 0.75) # 1st quartile prediction rmse: lasso
  ))
}


range.n <- seq(40, 300, 10)

# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
registerDoParallel(no_cores)

# get results
temp <- foreach(
  n = range.n,
  .combine = cbind,
  .multicombine = TRUE
) %dopar% {
  sapply(n, f.rep)
}

save(temp, file = "data/dji.data.RData")

# ---------------------------------
# GET GRAPH FOR MODEL SIZE AND AUC
# ---------------------------------

load("data/dji.data.RData")

plot.d <- data.frame(
  x = range.n,
  Method = rep(c(
    "ProSGPV", "BeSS", "ISIS", "Lasso"
  ), each = length(range.n)),
  size.m = c(temp[1, ], temp[7, ], temp[13, ], temp[19, ]),
  size.1 = c(temp[2, ], temp[8, ], temp[14, ], temp[20, ]),
  size.3 = c(temp[3, ], temp[9, ], temp[15, ], temp[21, ]),
  auc.m = c(temp[4, ], temp[10, ], temp[16, ], temp[22, ]),
  auc.1 = c(temp[5, ], temp[11, ], temp[17, ], temp[23, ]),
  auc.3 = c(temp[6, ], temp[12, ], temp[18, ], temp[24, ])
)

plot.d$Method <- factor(plot.d$Method,
  levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  )
)

dcols <- c("black", "springgreen3", "blue", "red")

p.auc <- ggplot(plot.d, aes(x = x, y = auc.m, color = Method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Prediction RMSE in a test set") +
  geom_errorbar(aes(ymin = auc.1, ymax = auc.3), width = 2) +
  scale_color_manual(values = dcols) +
  guides(fill = F) +
  scale_x_continuous(
    limits = c(40, 300),
    breaks = seq(40, 300, 40)
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

p.size <- ggplot(plot.d, aes(x = x, y = size.m, color = Method)) +
  geom_line() +
  geom_point() +
  labs(x = "n", y = "Model size in the training set") +
  geom_errorbar(aes(ymin = size.1, ymax = size.3), width = 2) +
  scale_color_manual(values = dcols) +
  guides(fill = F) +
  scale_x_continuous(
    limits = c(40, 300),
    breaks = seq(40, 300, 40)
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )

# ISIS needs more n to include variables (e.g., n=2000)

png("Figures/Figure 5.png", units = "in", width = 8, height = 4, res = 300)
ggarrange(p.size, p.auc, ncol = 2, common.legend = T, legend = "bottom",
          labels=LETTERS[1:2])
dev.off()


# pro.sgpv selected 6, 7, 8, 10 most often, sometimes it selected 6, 7, 8, 9, 10

# "ROC_5"  "ROC_10" "ROC_15" "EMA_10"
# "ROC_20"

# lasso selected
# 2  6  7  8 10


get.lasso <- function(train.n = 200) {
  while (T) {
    index.use <- sample(train.all.index, train.n)

    x <- dji.d[index.use, -1]
    y <- dji.d[index.use, 1]

    # lasso
    cv.lasso <- cv.glmnet(as.matrix(x), y)
    lasso.index <- which(coef(cv.lasso, s = cv.lasso$lambda.min)[-1] != 0)
    out.lasso.1 <- length(lasso.index)

    if (out.lasso.1 == 5) {
      return(paste(lasso.index, collapse = " "))
      break
    }
  }
}

temp.2 <- replicate(1e3, get.lasso())
sort(table(temp.2), decreasing = T)

# 6 7 8 10 63  2 6 7 8 10  6 7 8 9 10   2 6 7 10 63  6 7 9 10 63  2 6 7 9 10  2 6 8 10 63
# 453          432         102          4            2            1           1
# 6 7 8 10 40  6 7 8 10 50 6 7 8 10 52  6 7 8 10 65  6 7 8 10 69
# 1            1           1            1            1





get.sgpv <- function(train.n = 200) {
  while (T) {
    index.use <- sample(train.all.index, train.n)

    x <- dji.d[index.use, -1]
    y <- dji.d[index.use, 1]

    # sgpv
    sgpv.m <- pro.sgpv(x, y)
    out.sgpv.1 <- length(sgpv.m$var.index)


    if (out.sgpv.1 == 4) {
      return(paste(sgpv.m$var.index, collapse = " "))
      break
    }
  }
}


# sgpv
temp.3 <- replicate(1e3, get.sgpv())
sort(table(temp.3), decreasing = T)

# 6 7 8 10   6 7 9 10   2 6 7 10  6 7 10 63
# 935        50         8         7

# what are they?

names(x)[c(6, 7, 8, 10)]

# "ROC_5"  "ROC_10" "ROC_15" "EMA_10"
# ROC : RATE OF CHANGE
# EMA: EXPONENTIAL MOVING AVERAGE

summary(lm(y ~ ., data = data.frame(x[, c(6, 7, 8, 10)], y))) # adjusted R^2 = 0.9996

head(x[, c(6, 7, 8, 10)])

# --------------------------------------
# FIND OUT THE TIME PERIOD
# --------------------------------------

full.d <- read.csv("data/Processed_DJI.csv")
names(full.d)

head(full.d$Date)
tail(full.d$Date)

full.d$Date[1600:1620]
