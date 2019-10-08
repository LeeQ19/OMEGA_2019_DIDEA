source("didea.R")

df.f.2d <- read.csv(url("http://bit.ly/RevisedFire4data"),   header = T)
df.f.3d <- simplify2array(by(df.f.2d[, -c(1, 11)], df.f.2d$Year, as.matrix))
id.t <- c(1)
id.x <- c(2:4)
id.y <- c(5:6)
rts   <- "crs"
ori <- "i"

results <- didea(df.f.3d[, id.x, ], df.f.3d[, id.y, ], df.f.3d[, id.y, ] * 1.1, rts = rts, orientation = ori)
