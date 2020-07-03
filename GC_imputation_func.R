impute_binary_tau <- function(Y1, delta, X1, T1) {
  size <- length(Y1)
  
  df.pair <- matrix(nrow = choose(n = size, k = 2), ncol = 8)
  k <- 1
  for(i in seq(1, size - 1)) {
    for(j in seq(i + 1, size)) {
      df.pair[k, ] <- c(X1[i], Y1[i], delta[i], T1[i], X1[j], Y1[j], delta[j], T1[j])
      k <- k + 1
    }
  }
  
  G.fit1 <- survfit(Surv(time = Y1[X1 == 1], event = 1 - delta[X1 == 1]) ~ 1)
  G.fit0 <- survfit(Surv(time = Y1[X1 == 0], event = 1 - delta[X1 == 0]) ~ 1)
  G_func1 <- stepfun(x = G.fit1$time, y = c(1, G.fit1$surv))
  G_func0 <- stepfun(x = G.fit0$time, y = c(1, G.fit0$surv))
  
  H10.hat <- function(x, y, X1.obs = X1, T1.obs = Y1) {
    if(x == 0) {
      if(G_func0(y) == 0) return(0)
      
      return(mean(T1.obs[X1.obs == 0] > y) / G_func0(y))
      
    } else {
      if(G_func1(y) == 0) return(0)
      
      return(mean(T1.obs[X1.obs == 1] > y) / G_func1(y))
    }
  }
  
  H11.hat <- function(x, y, X1.obs = X1, T1.obs = Y1) {
    left <- H10.hat(x, y - 1e-3)
    right <- H10.hat(x, y)
    
    return(left - right)
  }
  
  ## imputation
  c.hat <- numeric(nrow(df.pair))
  d.hat <- numeric(nrow(df.pair))
  
  for(k in seq(1, nrow(df.pair))) {
    if(df.pair[k, 1] > df.pair[k, 5]) {
      ## xi > xi'
      if(df.pair[k, 2] > df.pair[k, 6]) {
        ## yi > yi'
        H10.hat.i <- H10.hat(0, df.pair[k, 2])
        H10.hat.ip <- H10.hat(0, df.pair[k, 6])
        
        if(df.pair[k, 7] == 1) {
          c.hat[k] <- 1
          d.hat[k] <- 0
        } else if(df.pair[k, 3] == 1 & df.pair[k, 7] == 0) {
          d.hat[k] <- ifelse(H10.hat.ip == 0, 0, H10.hat.i / H10.hat.ip)
          c.hat[k] <- 1 - d.hat[k]
        } else {
          t <- Y1[Y1 > df.pair[k, 2] & X1 == 1]
          
          if(length(t) == 0) {
            c.hat[k] <- 0
            d.hat[k] <- 0
          } else {
            integral.c <- numeric(length(t))
            integral.d <- numeric(length(t))
            
            for(m in seq(1, length(t))) {
              integral.c[m] <- (H10.hat.ip - H10.hat(0, t[m])) * H11.hat(1, t[m])
              integral.d[m] <- H10.hat(0, t[m]) * H11.hat(1, t[m])
            }
            
            c.hat[k] <- sum(integral.c) / H10.hat(1, df.pair[k, 2]) / H10.hat.ip
            d.hat[k] <- sum(integral.d) / H10.hat(1, df.pair[k, 2]) / H10.hat.ip
            
            c.hat[k] <- ifelse(H10.hat.ip * H10.hat(1, df.pair[k, 2]) == 0, 0, c.hat[k])
            d.hat[k] <- ifelse(H10.hat.ip * H10.hat(1, df.pair[k, 2]) == 0, 0, d.hat[k])
          }
        }
        
      } else {
        ## yi < yi'
        H10.hat.i <- H10.hat(1, df.pair[k, 2])
        H10.hat.ip <- H10.hat(1, df.pair[k, 6])
        
        if(df.pair[k, 3] == 1) {
          c.hat[k] <- 0
          d.hat[k] <- 1
        } else if(df.pair[k, 3] == 0 & df.pair[k, 7] == 1) {
          c.hat[k] <- ifelse(H10.hat.i == 0, 0, H10.hat.ip / H10.hat.i)
          d.hat[k] <- 1 - c.hat[k]
        } else {
          s <- Y1[Y1 > df.pair[k, 6] & X1 == 0]
          
          if(length(s) == 0) {
            c.hat[k] <- 0
            d.hat[k] <- 0
          } else {
            integral.c <- numeric(length(s))
            integral.d <- numeric(length(s))
            
            for(m in seq(1, length(s))) {
              integral.c[m] <- H10.hat(1, s[m]) * H11.hat(0, s[m])
              integral.d[m] <- (H10.hat.i - H10.hat(1, s[m])) * H11.hat(0, s[m])
            }
            c.hat[k] <- sum(integral.c) / H10.hat.i / H10.hat(0, df.pair[k, 6])
            d.hat[k] <- sum(integral.d) / H10.hat.i / H10.hat(0, df.pair[k, 6])
            
            c.hat[k] <- ifelse(H10.hat.i * H10.hat(0, df.pair[k, 6]) == 0, 0, c.hat[k])
            d.hat[k] <- ifelse(H10.hat.i * H10.hat(0, df.pair[k, 6]) == 0, 0, d.hat[k])
          }
        }
      }
    }
    
    if(df.pair[k, 1] < df.pair[k, 5]) {
      ## xi < xi'
      
      if(df.pair[k, 2] > df.pair[k, 6]) {
        ## yi > yi'
        H10.hat.i <- H10.hat(1, df.pair[k, 2])
        H10.hat.ip <- H10.hat(1, df.pair[k, 6])
        
        if(df.pair[k, 7] == 1) {
          c.hat[k] <- 0
          d.hat[k] <- 1
        } else if(df.pair[k, 3] == 1 & df.pair[k, 7] == 0) {
          c.hat[k] <- ifelse(H10.hat.ip == 0, 0, H10.hat.i / H10.hat.ip)
          d.hat[k] <- 1 - c.hat[k]
        } else {
          t <- Y1[Y1 > df.pair[k, 2] & X1 == 0]
          
          if(length(t) == 0) {
            c.hat[k] <- 0
            d.hat[k] <- 0
          } else {
            integral.c <- numeric(length(t))
            integral.d <- numeric(length(t))
            
            for(m in seq(1, length(t))) {
              integral.c[m] <- H10.hat(1, t[m]) * H11.hat(0, t[m])
              integral.d[m] <- (H10.hat.ip - H10.hat(1, t[m])) * H11.hat(0, t[m])
            }
            c.hat[k] <- sum(integral.c) / H10.hat(0, df.pair[k, 2]) / H10.hat.ip
            d.hat[k] <- sum(integral.d) / H10.hat(0, df.pair[k, 2]) / H10.hat.ip
            
            c.hat[k] <- ifelse(H10.hat.ip * H10.hat(0, df.pair[k, 2]) == 0, 0, c.hat[k])
            d.hat[k] <- ifelse(H10.hat.ip * H10.hat(0, df.pair[k, 2]) == 0, 0, d.hat[k])
          }
        }
        
      } else {
        ## yi < yi'
        H10.hat.i <- H10.hat(0, df.pair[k, 2])
        H10.hat.ip <- H10.hat(0, df.pair[k, 6])
        
        if(df.pair[k, 3] == 1) {
          c.hat[k] <- 1
          d.hat[k] <- 0
        } else if(df.pair[k, 3] == 0 & df.pair[k, 7] == 1) {
          d.hat[k] <- ifelse(H10.hat.i == 0, 0, H10.hat.ip / H10.hat.i)
          c.hat[k] <- 1 - d.hat[k]
        } else {
          s <- Y1[Y1 > df.pair[k, 6] & X1 == 1]
          
          if(length(s) == 0) {
            c.hat[k] <- 0
            d.hat[k] <- 0
          } else {
            integral.c <- numeric(length(s))
            integral.d <- numeric(length(s))
            
            for(m in seq(1, length(s))) {
              integral.c[m] <- (H10.hat.i - H10.hat(0, s[m])) * H11.hat(1, s[m])
              integral.d[m] <- H10.hat(0, s[m]) * H11.hat(1, s[m])
            }
            c.hat[k] <- sum(integral.c) / H10.hat(1, df.pair[k, 6]) / H10.hat.i
            d.hat[k] <- sum(integral.d) / H10.hat(1, df.pair[k, 6]) / H10.hat.i
            
            c.hat[k] <- ifelse(H10.hat.i * H10.hat(1, df.pair[k, 6]) == 0, 0, c.hat[k])
            d.hat[k] <- ifelse(H10.hat.i * H10.hat(1, df.pair[k, 6]) == 0, 0, d.hat[k])
          }
        }
        
      }
    }
    
    if(is.na(c.hat[k]) | is.na(d.hat[k])) {
      print(df.pair[k, ])
    }
    
    if(abs(c.hat[k]) > 1 | abs(d.hat[k]) > 1) {
      # print(c(k, c.hat[k], d.hat[k]))
    }
  }
  
  tau.imputed <- sum(c.hat) - sum(d.hat)
  tau.true <- sum((df.pair[, 1] - df.pair[, 5]) * ifelse(df.pair[, 4] > df.pair[, 8], 1, -1))
  
  c.true <- sum((df.pair[, 1] - df.pair[, 5]) * ifelse(df.pair[, 4] > df.pair[, 8], 1, -1) > 0)
  d.true <- sum((df.pair[, 1] - df.pair[, 5]) * ifelse(df.pair[, 4] > df.pair[, 8], 1, -1) < 0)
  
  
  n0 <- sum(X1 == 0)
  n1 <- sum(X1 == 1)
  
  logrank <- survdiff(Surv(time = Y1, event = delta) ~ X1)
  gehan <- survdiff(Surv(time = Y1, event = delta) ~ X1, rho = 1)
  gehan2 <- Gehan_func(Y1, delta, X1)
  
  n1 <- sum(X1)
  n0 <- size - n1
  mu.1 <- n1 * (n1 + n0 + 1) * 0.5
  sigma.1 <- sqrt(n1 * n0 * (n1 + n0 + 1) / 12)
  wilcox <- (sum(rank(T1)[X1 == 1]) - mu.1) / sigma.1
  
  return(list(bias = (tau.imputed - tau.true) / n0 / n1,
              C.rate = mean(delta),
              tau.true = tau.true / n0 / n1,
              c.hat = sum(c.hat) / n0 / n1,
              d.hat = sum(d.hat) / n0 / n1,
              c.true = c.true / n0 / n1,
              d.true = d.true / n0 / n1,
              tau.hat = tau.imputed / n0 / n1,
              logrank = logrank$chisq,
              wilcox = wilcox ^ 2,
              gehan = gehan$chisq,
              gehan2 = gehan2))
}

Gehan_func <- function(Y1, delta, X1) {
  m <- sum(X1 == 1)
  n <- sum(X1 == 0)
  Gehan.score <- numeric(length(Y1))
  
  for(i in seq(1, length(Y1))) {
    Gehan.score[i] <- sum(Y1[i] > Y1 & (delta == 1)) - sum(Y1[i] < Y1 & (delta[i] == 1))
  }
  
  ts <- sum(Gehan.score[X1 == 1])
  sigma2 <- sum(Gehan.score ^ 2) / (m + n)
  v <- m * n * sigma2 / (m + n - 1)
  
  return(list(z = ts / sqrt(v), score = Gehan.score))
}

MannWhitney_func <- function(Y1, delta, X1, T1) {
  m <- sum(X1 == 1)
  n <- sum(X1 == 0)
  MW.score <- numeric(length(T1))
  
  for(i in seq(1, length(T1))) {
    MW.score[i] <- sum(T1[i] > T1) - sum(T1[i] < T1)
  }
  
  ts <- sum(MW.score[X1 == 1])
  sigma2 <- sum(MW.score ^ 2) / (m + n)
  v <- m * n * sigma2 / (m + n - 1)
  
  return(list(z = sum(ts) / sqrt(v), score = MW.score, tau = sum(ts) / m / n))
}