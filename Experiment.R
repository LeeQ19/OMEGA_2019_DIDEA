source("didea.R")

# Test 1
data.test1 <- array(c(1, 1, 2, 2, 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 2), c(5, 3))
data.test1
id.x <- c(1:2)
id.y <- c(3)
rts  <- "crs"
ori  <- "i"
res.test1 <- dm.dea(data.test1[, id.x], data.test1[, id.y], rts = rts, orientation = ori)
res.test1$eff
res.test1$xslack

# Test 2
data.test2 <- array(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1), c(2, 3, 3))
data.test2
id.x <- c(1:2)
id.y <- c(3)
rts  <- "crs"
ori  <- "i"

res.test2 <- didea(data.test2[, id.x, ], data.test2[, id.y, ], data.test2[, id.y, ] * 2, rts = rts, orientation = ori)
abind(res.test2$alpha, res.test2$beta, along = 2)
