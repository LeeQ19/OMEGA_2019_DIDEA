# Load library
library(DJL)
library(abind)

# Dynamic Inverse DEA
didea <- function(xdata, ydata, target, rts = "crs", orientation) {
  
  # Initial checks
  if (length(unique(dim(xdata)[1], dim(ydata)[1], dim(target)[1])) != 1) 
    stop("Data must be balanced.")
  if (length(unique(rev(dim(xdata))[1], rev(dim(ydata))[1], rev(dim(target))[1])) != 1) 
    stop("Data must be balanced.")
  if (!(rts %in% c("crs", "vrs", "irs", "drs"))) 
    stop("rts must be \"crs\", \"vrs\", \"irs\", or \"drs\".")
  if (!(orientation %in% c("i", "o"))) 
    stop("orientation must be \"i\", \"o\".")
  if(!((orientation == "i" && dim(ydata)[2] == dim(target)[2]) || (orientation == "o" && dim(xdata)[2] == dim(target)[2])))
    stop("Data must be balanced.")
  
  names <- if (is.null(rownames(xdata))) 
    1:dim(xdata)[1]
  else rownames(xdata)
  xdata <- if (length(dim(xdata)) != 3) 
    array(xdata, c(dim(xdata)[1], 1, dim(xdata)[2]))
  else as.array(xdata)
  ydata <- if (length(dim(ydata)) != 3) 
    array(ydata, c(dim(ydata)[1], 1, dim(ydata)[2]))
  else as.array(ydata)
  target <- if (length(dim(target)) != 3) 
    array(target, c(dim(target)[1], 1, dim(target)[2]))
  else as.array(target)
  n <- dim(xdata)[1]
  m <- dim(xdata)[2]
  s <- dim(ydata)[2]
  t <- dim(xdata)[3]
  results.lambda <- array(NA, dim = c(n, n, t), dimnames = list(names, 
                                                                names))
  results.alpha  <- array(NA, dim = c(n, m, t), dimnames = list(names, 
                                                                paste0("alpha", 1:m)))
  results.beta   <- array(NA, dim = c(n, s, t), dimnames = list(names, 
                                                                paste0("beta", 1:s)))
  results.xslack <- array(NA, dim = c(n, m, t), dimnames = list(names, 
                                                                paste0("xslack", 1:m)))
  results.yslack <- array(NA, dim = c(n, s, t), dimnames = list(names, 
                                                                paste0("yslack", 1:s)))
  efficiency     <- apply(abind(xdata, ydata, along = 2), 3, function(data) {dm.dea(data[, 1:m], data[, (m + 1):(m + s)], rts = rts, orientation = orientation)$eff})
  
  p.alp <- n * t + 1
  p.bet <- n * t + m * t + 1
  p.xsl <- n * t + m * t + s * t + 1
  p.ysl <- n * t + m * t + s * t + m * t + 1
  p.end <- n * t + m * t + s * t + m * t + s * t + 1
  
  for (j in 1: n) {
    lp.dea <- make.lp(0, n * t + m * t + s * t + m * t + s * t)
    
    if (orientation == "i") 
      set.objfn(lp.dea, c(rep(1, m * t)), indices = c(p.alp:(p.bet - 1)))
    else
      set.objfn(lp.dea, c(rep(-1, s * t)), indices = c(p.bet:(p.xsl - 1)))
    
    if (rts == "crs") 
      for (k in 1:t) {
        set.constr.type(lp.dea, 0, k)
      }
    
    for (k in 1:t) {
      for (i in 1:m) {
        if (orientation == "i") {
          add.constraint(lp.dea, c(xdata[, i, k], -efficiency[j, k], 1), 
                         indices = c((n * (k - 1) + 1):(n * (k - 1) + n), p.alp - 1 + m * (k - 1) + i, p.xsl - 1 + m * (k - 1) + i), "=", 0)
          add.constraint(lp.dea, c(1), indices = c(p.alp - 1 + m * (k - 1) + i), ">=", xdata[j, i, k])
        }
        else {
          add.constraint(lp.dea, c(xdata[, i, k], -1, 1), 
                         indices = c((n * (k - 1) + 1):(n * (k - 1) + n), p.alp - 1 + m * (k - 1) + i, p.xsl - 1 + m * (k - 1) + i), "=", 0)
          add.constraint(lp.dea, c(1), indices = c(p.alp - 1 + m * (k - 1) + i), "<=", xdata[j, i, k])
        }
      }
      
      for (r in 1:s) {
        if (orientation == "i") {
          add.constraint(lp.dea, c(ydata[, r, k], -1, -1), 
                         indices = c((n * (k - 1) + 1):(n * (k - 1) + n), p.bet - 1 + s * (k - 1) + r, p.ysl - 1 + r * (k - 1) + r), "=", 0)
          add.constraint(lp.dea, c(1), indices = c(p.bet - 1 + s * (k - 1) + r), ">=", ydata[j, r, k])
        }
        else {
          add.constraint(lp.dea, c(ydata[, r, k], -efficiency[j, k], -1), 
                         indices = c((n * (k - 1) + 1):(n * (k - 1) + n), p.bet - 1 + s * (k - 1) + r, p.ysl - 1 + r * (k - 1) + r), "=", 0)
          add.constraint(lp.dea, c(1), indices = c(p.bet - 1 + s * (k - 1) + r), "<=", ydata[j, r, k])
        }
      }
    }
    
    if (orientation == "i") {
      for (r in 1:s) {
        add.constraint(lp.dea, c(rep(1, t)), 
                       indices = c(p.bet - 1 + s * 0:(t - 1) + r), "=", sum(target[j, r, ]))
      }
    }
    else {
      for (i in 1:m) {
        add.constraint(lp.dea, c(rep(1, t)), 
                       indices = c(p.alp - 1 + m * 0:(t - 1) + i), "=", sum(target[j, i, ]))
      }
    }
    
    set.bounds(lp.dea, lower = c(rep(0, p.end - 1)))
    
    solve.lpExtPtr(lp.dea)
    temp.p <- get.variables(lp.dea)
    results.lambda[j, , ] <- array(temp.p[1:(p.alp - 1)], c(n, t))
    results.alpha [j, , ] <- array(temp.p[p.alp:(p.bet - 1)], c(m, t))
    results.beta  [j, , ] <- array(temp.p[p.bet:(p.xsl - 1)], c(s, t))
    results.xslack[j, , ] <- array(temp.p[p.xsl:(p.ysl - 1)], c(m, t))
    results.yslack[j, , ] <- array(temp.p[p.ysl:(p.end - 1)], c(s, t))
    
    for (k in 1:t) {
      for (i in 1:m) {
        add.constraint(lp.dea, c(1), indices = c(p.alp - 1 + m * (k - 1) + i), "=", results.alpha[j, i, k])
      }
      
      for (r in 1:s) {
        add.constraint(lp.dea, c(1), indices = c(p.bet - 1 + s * (k - 1) + r), "=", results.beta[j, r, k])
      }
    }
    
    set.objfn(lp.dea, c(rep(-1, p.end - p.xsl)), indices = c(p.xsl:(p.end - 1)))
    
    solve.lpExtPtr(lp.dea)
    temp.s <- get.variables(lp.dea)
    results.lambda[j, , ] <- array(temp.s[1:(p.alp - 1)], c(n, t))
    results.xslack[j, , ] <- array(temp.s[p.xsl:(p.ysl - 1)], c(m, t))
    results.yslack[j, , ] <- array(temp.s[p.ysl:(p.end - 1)], c(s, t))
  }
  results <- list(lambda = results.lambda, eff.t = efficiency, 
                  alpha = results.alpha, beta = results.beta, 
                  xslack = results.xslack, yslack = results.yslack)
  return(results)
}
