# maximizing multinomial likelihood
y <- c(1997, 906, 904, 32)
m <- sum(y)
# functions: log-likelihood, 1st derivative, 2nd derivative, and expected info
f.l <- function(theta, y) {
  temp <- y[1] * log(2 + theta) +
    (y[2] + y[3]) * log(1 - theta) +
    y[4] * log(theta)
  return(temp)
}
f.dl <- function(theta, y) {
  temp <- y[1] / (2+theta) +
    - (y[2] + y[3]) / (1 - theta) +
    y[4] / theta
  return(temp)
}
f.ddl <- function(theta, y) {
  temp <- - (y[1] / (2 + theta)^2 +
               (y[2] + y[3]) / (1-theta)^2 +
               y[4] / theta^2
  )
  return(temp)
}
f.info <- function(theta, y) {
  temp <- 0.25 * sum(y) * (1 / (2 + theta) +
                             2 / (1 - theta) +
                             1 / theta )
  return(temp)
}
# plot functions
library(ggplot2)
p1 <- ggplot(data.frame(theta = c(0.0001, 0.4)), aes(theta))
p1 <- p1 + stat_function(fun = f.l, args = list(y))
p1 <- p1 + labs(title = "log-likelihood")
#print(p1)
p2 <- ggplot(data.frame(theta = c(0.01, 0.4)), aes(theta))
p2 <- p2 + geom_hline(yintercept = 0, alpha = 0.5)
p2 <- p2 + stat_function(fun = f.dl, args = list(y))
p2 <- p2 + labs(title = "1st derivative")
#print(p2)
p3 <- ggplot(data.frame(theta = c(0.01, 0.4)), aes(theta))
p3 <- p3 + geom_hline(yintercept = 0, alpha = 0.5)
p3 <- p3 + stat_function(fun = f.ddl, args = list(y))
p3 <- p3 + labs(title = "2nd derivative")
#print(p3)
p4 <- ggplot(data.frame(theta = c(0.01, 0.4)), aes(theta))
p4 <- p4 + stat_function(fun = f.info, args = list(y))
p4 <- p4 + labs(title = "expected info")
#print(p4)
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol=2)







## 
###f.NR() function.

# NR routine for finding root of g(x) = 0.
# Requires predefined g(x) and gp(x) = deriv of g(x)
# The iteration is controlled by:
# eps = absolute convergence criterion
# maxit = maximum allowable number of iterations
# Input: xnew = user prompted starting value
# Output: number of root, steps, and note
f.NR <- function(g, gp, xnew = 1, eps = 0.001, maxit = 35, y = c(1,1,1,1)) {
  xold <- -Inf # needed so argument in while() loop is defined
  i <- 1; # initial iteration index
  NR.hist <- data.frame(i, xnew, diff = abs(xnew - xold)) # iteration history
  while ((i <= maxit) & (abs(xnew - xold) > eps)) {
    i <- i + 1 # increment iteration
    xold <- xnew # old guess is current guess
    xnew <- xold - g(xold, y) / gp(xold, y) # new guess
    NR.hist <- rbind(NR.hist, c(i, xnew, abs(xnew - xold))) # iteration history
  }
  out <- list()
  out$root <- xnew
  out$iter <- i
  out$hist <- NR.hist
  if (abs(xnew - xold) <= eps) {
    out$note <- paste("Absolute convergence of", eps, "satisfied")
  }
  if (i > maxit) {
    out$note <- paste("Exceeded max iterations of ", maxit)
  }
  return(out)
}



out0.01 <- f.NR(f.dl, f.ddl, xnew = 0.01, y = y)
out0.01
## $root
## [1] 0.03571
##
## $iter
## [1] 6
##
## $hist
## i xnew diff
## 1 1 0.01000 Inf
## 2 2 0.01734 0.0073377
## 3 3 0.02647 0.0091313
## 4 4 0.03344 0.0069732
## 5 5 0.03558 0.0021373
## 6 6 0.03571 0.0001323
##
## $note
## [1] "Absolute convergence of 0.001 satisfied"
out0.05 <- f.NR(f.dl, f.ddl, xnew = 0.05, y = y)
out0.05
## $root
## [1] 0.0357
##
## $iter
## [1] 4
##
## $hist
## i xnew diff
## 1 1 0.05000 Inf
## 2 2 0.03095 0.0190512
## 3 3 0.03512 0.0041720
## 4 4 0.03570 0.0005826
##
## $note
## [1] "Absolute convergence of 0.001 satisfied"
out0.20 <- f.NR(f.dl, f.ddl, xnew = 0.20, y = y)
out0.20
## $root
## [1] -0.4668
##
## $iter
## [1] 6
##
## $hist
## i xnew diff
## 1 1 0.20000 Inf
## 2 2 -0.09568 0.2956825
## 3 3 -0.26453 0.1688450
## 4 4 -0.44285 0.1783252
## 5 5 -0.46669 0.0238361
## 6 6 -0.46681 0.0001253
##
## $note
## [1] "Absolute convergence of 0.001 satisfied"
out0.40 <- f.NR(f.dl, f.ddl, xnew = 0.40, y = y)
out0.40
## $root
## [1] 0.0357
##
## $iter
## [1] 5
##
## $hist
## i xnew diff
## 1 1 0.40000 Inf
## 2 2 0.02246 0.3775390
## 3 3 0.03098 0.0085169
## 4 4 0.03513 0.0041502
## 5 5 0.03570 0.0005755
##
## $note
## [1] "Absolute convergence of 0.001 satisfied"
out0.50 <- f.NR(f.dl, f.ddl, xnew = 0.50, y = y)
out0.50
## $root
## [1] -0.4668
##
## $iter
## [1] 7
##
## $hist
## i xnew diff
## 1 1 0.5000 Inf
## 2 2 0.1413 0.3586592
## 3 3 -0.0699 0.2112391
## 4 4 -0.1985 0.1286382
## 5 5 -0.4080 0.2094407
## 6 6 -0.4659 0.0578853
## 7 7 -0.4668 0.0009514
##
## $note
## [1] "Absolute convergence of 0.001 satisfied"
# estimated standard deviation via Fisher's information
sqrt(1/f.info(out0.05$root, y))
## [1] 0.005838
# estimated standard deviation via second derivative
sqrt(-1/f.ddl(out0.05$root, y))
## [1] 0.006027
# plot functions
library(ggplot2)
p <- ggplot(data.frame(theta = c(0.0001, 0.55)), aes(theta))
p <- p + geom_hline(yintercept = 0, alpha = 0.5)
p <- p + stat_function(fun = function(theta, y)
{theta - f.dl(theta, y) / f.ddl(theta, y)}
, args = list(y))
p <- p + labs(title = "theta - f.dl(theta, y) / f.ddl(theta, y)")
print(p)

