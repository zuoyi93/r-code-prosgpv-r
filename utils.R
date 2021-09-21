# ------------------
# UTILITY FUNCTIONS
# ------------------

# ------------------
# FIGURE 1
# ------------------

get.fig.1 <- function() {

  # limits
  xlim <- c(-10, 10)
  ylim <- c(-0.5, 6.5)

  # cases
  t.x <- -8
  t.y <- c(6, 4, 2)
  t.text <- c(
    as.expression(bquote("0<" ~ p[delta] ~ "<1")),
    as.expression(bquote(p[delta] == 0)),
    as.expression(bquote(p[delta] == 1))
  )

  # confidence interval
  ci.l.x <- c(0, 6, -2)
  ci.l.y <- c(4, 2, 0)
  ci.u.x <- c(10, 9, 3)
  ci.u.y <- c(4, 2, 0)
  theta.hat.x <- c(5, 7.5, 0.5)
  theta.hat.y <- c(4, 2, 0)

  par(mar = c(5.1, 2.1, 4.1, 2.1))

  plot(0,
    type = "n", xlim = xlim, ylim = ylim, yaxt = "n", xaxt = "n",
    xlab = "Effect size", ylab = "", cex.lab = 1.5
  )
  axis(1, at = seq(xlim[1], xlim[2], 2), cex.axis = 1.5)
  text(t.x, t.y, t.text, cex = 2)
  # null interval
  rect(-4, 4.5, 4, 5, border = NA, col = rgb(34, 139, 34, 100, maxColorValue = 255))
  rect(-4, 2.5, 4, 3, border = NA, col = rgb(34, 139, 34, 100, maxColorValue = 255))
  rect(-4, 0.5, 4, 1, border = NA, col = rgb(34, 139, 34, 100, maxColorValue = 255))
  text(-4, c(5.25, 3.25, 1.25), bquote("(H"[0]^"-"),
    col = rgb(34, 139, 34, 255, maxColorValue = 255), cex = 1.5
  )
  text(4, c(5.25, 3.25, 1.25), bquote("H"[0]^"+" ~ ")"),
    col = rgb(34, 139, 34, 255, maxColorValue = 255), cex = 1.5
  )
  legend(4, 6.7,
    box.col = "white", box.lwd = 0,
    legend = c("Interval null", "Confidence interval"),
    fill = c(
      rgb(34, 139, 34, 100, maxColorValue = 255),
      rgb(0, 191, 255, 100, maxColorValue = 255)
    ), border = NA,
    cex = 1.5
  )
  col <- rgb(34, 139, 34, 255, maxColorValue = 255)
  # confidence interval
  rect(0, 4.5, 10, 5, border = NA, col = rgb(0, 191, 255, 100, maxColorValue = 255))
  rect(6, 2.5, 9, 3, border = NA, col = rgb(0, 191, 255, 100, maxColorValue = 255))
  rect(-2, 0.5, 3, 1, border = NA, col = rgb(0, 191, 255, 100, maxColorValue = 255))
  text(ci.l.x, ci.l.y, bquote("(" ~ theta[l]),
    cex = 1.5,
    col = rgb(0, 191, 255, 255, maxColorValue = 255)
  )
  text(ci.u.x, ci.u.y, bquote(theta[u] ~ ")"),
    cex = 1.5,
    col = rgb(0, 191, 255, 255, maxColorValue = 255)
  )
  text(theta.hat.x, theta.hat.y, bquote(hat(theta)),
    cex = 1.5,
    col = rgb(0, 191, 255, 255, maxColorValue = 255)
  )
  # theta hat
  segments(theta.hat.x, theta.hat.y + 0.5, theta.hat.x, theta.hat.y + 1, lty = 2)
  arrows(xlim[1], c(4.75, 2.75, 0.75), xlim[2], c(4.75, 2.75, 0.75))
  arrows(xlim[2], c(4.75, 2.75, 0.75), xlim[1], c(4.75, 2.75, 0.75))
  # null bound
  segments(4, c(4.5, 2.5, 0.5), 4, c(5, 3, 1),
    lwd = 4, lend = 1,
    col = rgb(34, 139, 34, 255, maxColorValue = 255)
  )
  segments(-4, c(4.5, 2.5, 0.5), -4, c(5, 3, 1),
    lwd = 4, lend = 1,
    col = rgb(34, 139, 34, 255, maxColorValue = 255)
  )
  # ci bound
  segments(c(0, 10), 4.5, c(0, 10), 5,
    lwd = 4, lend = 1,
    col = rgb(0, 191, 255, 255, maxColorValue = 255)
  )
  segments(c(6, 9), 2.5, c(6, 9), 3,
    lwd = 4, lend = 1,
    col = rgb(0, 191, 255, 255, maxColorValue = 255)
  )
  segments(c(-2, 3), 0.5, c(-2, 3), 1,
    lwd = 4, lend = 1,
    col = rgb(0, 191, 255, 255, maxColorValue = 255)
  )

  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

# ------------------
# SIMULATION
# ------------------


main.one.time <- function(n, p, s, rho = 0.3, sig = 2,
                          beta.min = 0.1, beta.max = 0.4,
                          intercept = 2, nu = 2, gvif = F,
                          family = c("gaussian", "binomial", "poisson", "cox"),
                          scale = 2, shape = 1, rateC = 0.2) {

  # simulate training data
  sim.data <- gen.sim.data(
    n = round(n * 5 / 3), p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max, nu = nu,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  true.index <- sim.data[[3]]
  true.beta <- sim.data[[4]]

  # generate training and testing indices
  train.index <- sample(1:round(n * 5 / 3), n, replace = F)
  test.index <- setdiff(1:round(n * 5 / 3), train.index)

  # -------------------
  # lasso
  # -------------------

  if (family != "cox") {
    cv.m <- cv.glmnet(x[train.index, ], y[train.index], family = family)
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)[-1]
  }
  else {
    cv.m <- cv.glmnet(x[train.index, ],
      Surv(y[train.index, 1], y[train.index, 2]),
      family = "cox"
    )
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)
  }
  lasso.min.index <- which(lasso.min.coef != 0)

  # out 1: support recovery
  out.lasso.1 <- setequal(true.index, lasso.min.index)

  # out 2: parameter estimation
  out.lasso.2 <- mean(abs(as.numeric(lasso.min.coef) - true.beta))

  # out 3: prediction performance
  if (family == "gaussian") {
    lasso.pred <- predict(cv.m,
      s = cv.m$lambda.min, newx = x[test.index, ]
    )

    out.lasso.3 <- sqrt(mean((lasso.pred - y[test.index])^2))
  } else if (family == "binomial") {
    lasso.pred <- predict(cv.m,
      s = cv.m$lambda.min, newx = x[test.index, ],
      type = "response"
    )

    out.lasso.3 <- roc(y[test.index] ~ lasso.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    lasso.pred <- predict(cv.m,
      s = cv.m$lambda.min, newx = x[test.index, ],
      type = "response"
    )

    out.lasso.3 <- sqrt(mean((lasso.pred - y[test.index])^2))
  } else {
    out.lasso.3 <- 0
  }

  # -------------------
  # pro.sgpv
  # -------------------

  if (family != "cox") {
    out.sgpv <- pro.sgpv(x[train.index, ], y[train.index], family = family, gvif = gvif)
  } else {
    out.sgpv <- pro.sgpv(x[train.index, ], y[train.index, ], family = family, gvif = gvif)
  }
  sgpv.index <- out.sgpv$var.index

  # out 1: support recovery
  out.sgpv.1 <- setequal(true.index, sgpv.index)

  # out 2: parameter estimation
  out.sgpv.2 <- mean(abs(coef(out.sgpv) - true.beta))

  # out 3: prediction performance
  if (family == "gaussian") {
    sgpv.pred <- predict(out.sgpv,
      newdata = x[test.index, ]
    )

    out.sgpv.3 <- sqrt(mean((sgpv.pred - y[test.index])^2))
  } else if (family == "binomial") {
    sgpv.pred <- predict(out.sgpv,
      newdata = x[test.index, ],
      type = "response"
    )

    out.sgpv.3 <- roc(y[test.index] ~ sgpv.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    sgpv.pred <- predict(out.sgpv,
      newdata = x[test.index, ]
    )

    out.sgpv.3 <- sqrt(mean((sgpv.pred - y[test.index])^2))
  } else {
    out.sgpv.3 <- 0
  }

  # -------------------
  # BeSS
  # -------------------

  if (family != "cox") {
    bess.fit <- bess(x[train.index, ], y[train.index], family = family)
    bess.index <- which(bess.fit$beta != 0)
  } else {
    bess.fit <- bess(x[train.index, ], y[train.index, ], family = family)
    bess.index <- which(bess.fit$beta != 0)
  }

  # out 1: support recovery
  out.bess.1 <- setequal(true.index, bess.index)

  # out 2: parameter estimation
  out.bess.2 <- mean(abs(bess.fit$beta - true.beta))

  # out 3: prediction performance
  if (family == "gaussian") {
    bess.pred <- predict(bess.fit, newx = x[test.index, ])
    out.bess.3 <- sqrt(mean((bess.pred - y[test.index])^2))
  } else if (family == "binomial") {
    bess.pred <- predict(bess.fit, newx = x[test.index, ], type = "response")

    out.bess.3 <- roc(y[test.index] ~ bess.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    if (length(bess.index) > 0) {
      d.bess <- data.frame(y = y, x[, bess.index])

      poisson.bess <- glm(y ~ ., data = d.bess[train.index, ], family = "poisson", maxit = 3e2)

      bess.pred <- predict(poisson.bess, newdata = d.bess[test.index, ], type = "response")
      rmse.bess <- sqrt(mean((bess.pred - y[test.index])^2))
    } else {
      d.train <- data.frame(y = y[train.index], x[train.index, ])
      poisson.bess <- glm(y ~ 1, data = d.train, family = "poisson", maxit = 3e2)
      d.test <- data.frame(y = y[test.index], x[test.index, ])
      bess.pred <- predict(poisson.bess, newdata = d.test, type = "response")
      rmse.bess <- sqrt(mean((bess.pred - y[test.index])^2))
    }

    out.bess.3 <- rmse.bess
  } else {
    out.bess.3 <- 0
  }

  # -------------------
  # SIS
  # -------------------

  if (family != "cox") {
    invisible(capture.output(sis.m <- SIS(x[train.index, ],
      y[train.index],
      family = family, tune = "ebic"
    )))
  } else {
    invisible(capture.output(sis.m <- SIS(x[train.index, ],
      Surv(
        y[train.index, 1],
        y[train.index, 2]
      ),
      family = "cox",
      tune = "ebic", penalty = "lasso"
    )))
  }
  sis.index <- sis.m$ix

  # out 1: support recovery
  out.sis.1 <- setequal(true.index, sis.index)

  # out 2: parameter estimation
  sis.coef <- integer(p)

  if (length(sis.index) > 0) {
    if (family != "cox") {
      out.sis.coef <- as.numeric(sis.m$coef.est[-1])
    } else {
      out.sis.coef <- as.numeric(sis.m$coef.est)
    }


    for (i in 1:length(sis.index)) {
      sis.coef[sis.index[i]] <- out.sis.coef[i]
    }
  }

  out.sis.2 <- mean(abs(sis.coef - true.beta))

  # out 3: prediction performance
  if (family == "gaussian") {
    sis.pred <- predict(sis.m, newx = x[test.index, ], type = "response")
    out.sis.3 <- sqrt(mean((sis.pred - y[test.index])^2))
  } else if (family == "binomial") {
    sis.pred <- predict(sis.m, newx = x[test.index, ], type = "response")
    out.sis.3 <- roc(y[test.index] ~ sis.pred, plot = F, print.auc = F)$auc
  } else if (family == "poisson") {
    sis.pred <- predict(sis.m, newx = x[test.index, ], type = "response")
    out.sis.3 <- sqrt(mean((sis.pred - y[test.index])^2))
  } else {
    out.sis.3 <- 0
  }


  return(c(
    out.sgpv.1, out.sgpv.2, out.sgpv.3,
    out.lasso.1, out.lasso.2, out.lasso.3,
    out.bess.1, out.bess.2, out.bess.3,
    out.sis.1, out.sis.2, out.sis.3
  ))
}



main.many <- function(num.sim = 1e3, n, p, s, rho = 0.3, sig = 2,
                      beta.min = 0.1, beta.max = 0.4,
                      intercept = 0, nu = 2, gvif = F,
                      family = c("binomial", "poisson", "cox", "gaussian"),
                      scale = 2, shape = 1, rateC = 0.2) {
  suppressMessages(
    out <- replicate(num.sim, main.one.time(
      n = n, p = p, s = s, rho = rho,
      sig = sig, beta.min = beta.min,
      beta.max = beta.max, nu = nu,
      intercept = intercept, gvif = gvif,
      family = family, scale = scale,
      shape = shape, rateC = rateC
    ))
  )

  return(c(
    mean(out[1, ]), # sgpv capture rate
    median(out[2, ]), # sgpv median parameter estimation
    quantile(out[2, ], 0.25), # sgpv parameter estimation first quartile
    quantile(out[2, ], 0.75), # sgpv parameter estimation third quartile
    median(out[3, ]), # sgpv median prediction
    quantile(out[3, ], 0.25), # sgpv prediction first quartile
    quantile(out[3, ], 0.75), # sgpv prediction third quartile

    mean(out[4, ]), # lasso capture rate
    median(out[5, ]), # lasso median parameter estimation
    quantile(out[5, ], 0.25), # lasso parameter estimation first quartile
    quantile(out[5, ], 0.75), # lasso parameter estimation third quartile
    median(out[6, ]), # lasso median prediction
    quantile(out[6, ], 0.25), # lasso prediction first quartile
    quantile(out[6, ], 0.75), # lasso prediction third quartile

    mean(out[7, ]), # bess capture rate
    median(out[8, ]), # bess median parameter estimation
    quantile(out[8, ], 0.25), # bess parameter estimation first quartile
    quantile(out[8, ], 0.75), # bess parameter estimation third quartile
    median(out[9, ]), # bess median prediction
    quantile(out[9, ], 0.25), # bess prediction first quartile
    quantile(out[9, ], 0.75), # bess prediction third quartile

    mean(out[10, ]), # sis capture rate
    median(out[11, ]), # sis median parameter estimation
    quantile(out[11, ], 0.25), # sis parameter estimation first quartile
    quantile(out[11, ], 0.75), # sis parameter estimation third quartile
    median(out[12, ]), # sis median prediction
    quantile(out[12, ], 0.25), # sis prediction first quartile
    quantile(out[12, ], 0.75) # sis prediction third quartile
  ))
}

# ------------------
# run time
# ------------------

run.sgpv <- function(x,y,family,gvif){
  
  out.sgpv <- try(pro.sgpv(x, y, family, gvif),silent=T)
  sgpv.index <- try(out.sgpv$var.index,silent=T)
  
}

run.lasso <- function(x,y,family){
  
  if (family != "cox") {
    cv.m <- cv.glmnet(x, y, family = family)
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)[-1]
  }
  else {
    cv.m <- cv.glmnet(x,
                      Surv(y[, 1], y[, 2]),
                      family = "cox"
    )
    lasso.min.coef <- coef(cv.m, s = cv.m$lambda.min)
  }
  lasso.min.index <- which(lasso.min.coef != 0)
  
}

run.bess <- function(x,y,family){
  
  bess.fit <- bess(x, y, family = family)
  bess.index <- which(bess.fit$beta != 0)

}

run.isis <- function(x,y,family){
  
  if (family != "cox") {
    invisible(capture.output(sis.m <- SIS(x,
                                          y,
                                          family = family, tune = "ebic"
    )))
  } else {
    invisible(capture.output(sis.m <- SIS(x,
                                          Surv(
                                            y[, 1],
                                            y[, 2]
                                          ),
                                          family = "cox",
                                          tune = "ebic", penalty = "lasso"
    )))
  }
  sis.index <- sis.m$ix
  
  
}

run.time.one <- function(n, p, s, rho = 0.3, sig = 2,
                         beta.min = 0.1, beta.max = 0.4,
                         intercept = 2, nu = 2, gvif = F, 
                         family = c("gaussian", "binomial", "poisson", "cox"),
                         scale = 2, shape = 1, rateC = 0.2) {
  
  # simulate training data
  sim.data <- gen.sim.data(
    n = n, p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max, nu = nu,
    intercept = intercept, family = family,
    scale = 2, shape = 1, rateC = 0.2
  )
  x <- sim.data[[1]]
  y <- sim.data[[2]]
  
  # method 1: sgpv
  time.sgpv <- as.numeric(system.time(run.sgpv(x, y, family, gvif))[3])
  
  # method 2: lasso
  time.lasso <- as.numeric(system.time(run.lasso(x, y, family))[3])
  
  # method 3: bess
  time.bess <- as.numeric(system.time(run.bess(x, y, family))[3])
  
  # method 4: isis
  time.isis <- as.numeric(system.time(run.isis(x, y, family))[3])
  
  return(c(
    time.sgpv,
    time.lasso,
    time.bess,
    time.isis
  ))
  
}



run.time.many <- function(num.sim = 300, # 1e3 
                          n, p, s, rho = 0.3, sig = 2,
                          beta.min = 0.1, beta.max = 0.4,
                          intercept = 2, nu = 2, gvif = F,
                          family = c("gaussian", "binomial", "poisson", "cox")) {

  out <- replicate(num.sim, run.time.one(
    n = n, p = p, s = s, rho = rho, sig = sig,
    beta.min = beta.min, beta.max = beta.max, nu = nu,
    intercept = intercept, family = family, gvif = gvif
  ))
  
  median.sgpv <- median(out[1, ])
  l.sgpv <- as.numeric(quantile(out[1, ], 0.25))
  u.sgpv <- as.numeric(quantile(out[1, ], 0.75))
  
  median.lasso <- median(out[2, ])
  l.lasso <- as.numeric(quantile(out[2, ], 0.25))
  u.lasso <- as.numeric(quantile(out[2, ], 0.75))
  
  median.bess <- median(out[3, ])
  l.bess <- as.numeric(quantile(out[3, ], 0.25))
  u.bess <- as.numeric(quantile(out[3, ], 0.75))
  
  median.isis <- median(out[4, ])
  l.isis <- as.numeric(quantile(out[4, ], 0.25))
  u.isis <- as.numeric(quantile(out[4, ], 0.75))
  
  return(c(
    median.sgpv, l.sgpv, u.sgpv,
    median.lasso, l.lasso, u.lasso,
    median.bess, l.bess, u.bess,
    median.isis, l.isis, u.isis
  ))
}


# ------------------
# support recovery
# ------------------


get.plot.sr <- function(data, n = 20, p = 100, title.p,
                           type = c("ld", "hd"), num.sim = 1e3) {

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type == "ld") {
    plot.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19),
      rate = c(data[1, ], data[8, ], data[15, ], data[22, ])
    )

    xlim <- c(1, 20) * n
    xbreaks <- seq(2, 20, 4) * n
    xlab <- "n"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19)
    )
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19),
      rate = c(data[1, ], data[8, ], data[15, ], data[22, ])
    )

    xlim <- c(p, p * 10)
    xbreaks <- seq(p, p * 10, p * 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19)
    )
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c(
      "ProSGPV", "Lasso", "BeSS", "ISIS"
    )
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(1, 8, 15, 22)) {
    pe <- as.numeric(data[i, ])
    lb.out <- c(lb.out, pe - 1.96 * sqrt(pe * (1 - pe) / num.sim))
    ub.out <- c(ub.out, pe + 1.96 * sqrt(pe * (1 - pe) / num.sim))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d$lb[ci.d$lb < 0] <- 0
  ci.d$ub[ci.d$ub > 1] <- 1

  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  ))

  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      x = xlab, y = "Support recovery rate", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# ---------------------
# parameter estimation
# ---------------------


get.plot.pe <- function(data, n, p, type = c("ld", "hd"), title.p, cap = NULL,
                           ylim, ybreaks) {

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hd") {
    plot.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19),
      rate = c(data[2, ], data[9, ], data[16, ], data[23, ])
    )

    xlim <- c(1, 20) * n
    xbreaks <- seq(2, 20, 4) * n
    xlab <- "n"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19)
    )
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19),
      rate = c(data[2, ], data[9, ], data[16, ], data[23, ])
    )

    xlim <- c(p, p * 10)
    xbreaks <- seq(p, p * 10, p * 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19)
    )
  }

  plot.d$Method <- factor(plot.d$Method,
    levels = c("ProSGPV", "Lasso", "BeSS", "ISIS")
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(3, 10, 17, 24)) {
    lb.out <- c(lb.out, as.numeric(data[i, ]))
    ub.out <- c(ub.out, as.numeric(data[i + 1, ]))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  ))

  if (!is.null(cap)) {
    plot.d$rate[plot.d$rate > cap] <- cap
    ci.d$lb[ci.d$lb > cap] <- cap
    ci.d$ub[ci.d$ub > cap] <- cap
  }

  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = rate, col = Method)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(
      limits = ylim,
      breaks = ybreaks
    ) +
    labs(
      x = xlab, y = "Estimation MAE", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(x = x, ymin = lb, ymax = ub, fill = method), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}


# ------------------
# prediction
# ------------------

get.plot.pr <- function(data, n, p = 200, title.p, ylab, ylim, ybreaks, cap = NULL,
                           type = c("ld", "hd")) {

  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")

  if (type != "hd") {
    plot.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19),
      rate = c(data[5, ], data[12, ], data[19, ], data[26, ])
    )

    xlim <- c(1, 20) * n
    xbreaks <- seq(2, 20, 4) * n
    xlab <- "n"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19)
    )
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), 4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19),
      rate = c(data[5, ], data[12, ], data[19, ], data[26, ])
    )

    xlim <- c(p, p * 10)
    xbreaks <- seq(p, p * 10, p * 2)
    xlab <- "p"

    # add confidence interval
    ci.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), 4),
      method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), each = 19)
    )
  }


  plot.d$Method <- factor(plot.d$Method,
    levels = c("ProSGPV", "Lasso", "BeSS", "ISIS")
  )

  lb.out <- NULL
  ub.out <- NULL

  for (i in c(6, 13, 20, 27)) {
    lb.out <- c(lb.out, as.numeric(data[i, ]))
    ub.out <- c(ub.out, as.numeric(data[i + 1, ]))
  }

  ci.d$lb <- lb.out
  ci.d$ub <- ub.out

  ci.d$method <- factor(ci.d$method, levels = c(
    "ProSGPV", "Lasso", "BeSS", "ISIS"
  ))

  if (!is.null(cap)) {
    plot.d$rate[plot.d$rate > cap] <- cap
    ci.d$lb[ci.d$lb > cap] <- cap
    ci.d$ub[ci.d$ub > cap] <- cap
  }

  ggplot() +
    geom_line(
      data = plot.d,
      aes(x = x, y = rate, col = Method)
    ) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = ylim, breaks = ybreaks) +
    labs(
      x = xlab, y = ylab, col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = ci.d, aes(
      x = x, ymin = lb, ymax = ub, group = method,
      fill = method
    ), alpha = .4) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

# ------------------
# run time
# ------------------

get.plot.rt <- function(data, n, p = 200, np = c("n", "p"), title.p, y.cap = 0.8) {
  
  data[data > y.cap] <- y.cap
  
  # color scheme
  dcols <- c("black", "springgreen3", "blue", "red")
  
  if (np == "n") {
    plot.d <- data.frame(
      x = rep(seq(2, 20, 1) * n, each=4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 19),
      pe = c(data[c(1, 4, 7, 10), ]),
      lb = c(data[c(2, 5, 8, 11), ]),
      ub = c(data[c(3, 6, 9, 12), ])
    )
    
    xlim <- c(1, 20) * n
    xbreaks <- seq(2, 20, 4) * n
    xlab <- "n"
    
  } else {
    plot.d <- data.frame(
      x = rep(seq(p, p * 10, p / 2), each=4),
      Method = rep(c(
        "ProSGPV", "Lasso", "BeSS", "ISIS"
      ), 19),
      pe = c(data[c(1, 4, 7, 10), ]),
      lb = c(data[c(2, 5, 8, 11), ]),
      ub = c(data[c(3, 6, 9, 12), ])
    )
    
    xlim <- c(p, p * 10)
    xbreaks <- seq(p, p * 10, p * 2)
    xlab <- "p"
  }
  
  
  plot.d$Method <- factor(plot.d$Method,
                          levels = c(
                            "ProSGPV", "Lasso", "BeSS", "ISIS"
                          )
  )
  
  ggplot() +
    geom_line(data = plot.d, aes(x = x, y = pe, col = Method)) +
    scale_color_manual(values = dcols) +
    scale_x_continuous(limits = xlim, breaks = xbreaks) +
    scale_y_continuous(limits = c(0, y.cap), breaks = seq(0, y.cap, 0.2)) +
    labs(
      x = xlab, y = "Running time", col = "Method",
      title = title.p
    ) +
    geom_ribbon(data = plot.d, aes(
      x = x, ymin = lb,
      ymax = ub, fill = Method
    ), alpha = .3) +
    scale_fill_manual(values = dcols) +
    guides(fill = F) +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black")
    )
}
