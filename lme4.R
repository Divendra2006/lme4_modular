library(lme4)
library(nloptr)
library(glmmTMB)

data(cbpp)

lf <- lFormula(incidence / size ~ period + (period | herd), data = cbpp,
               control = lmerControl(check.nobs.vs.nRE = "ignore"))

lower_indices <- lf$reTrms$lower

devfun <- mkLmerDevfun(lf$fr, lf$X, lf$reTrms)

diagonal_wrapper <- function(theta_diag) {
  theta <- numeric(length(lf$reTrms$lower))
  theta[lower_indices == 0] <- theta_diag
  devfun(theta)
}

theta_diag_init <- rep(0, sum(lower_indices == 0))  
lower_bounds <- rep(0, length(theta_diag_init))    
upper_bounds <- rep(Inf, length(theta_diag_init))  

opt <- nloptwrap(
  par = theta_diag_init,
  fn = diagonal_wrapper,
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 1000) 
)

lmer_fit <- mkMerMod(
  rho = environment(devfun),  
  opt = opt,
  reTrms = lf$reTrms,
  fr = lf$fr,
  mc = match.call()
)

print(lmer_fit)

glmmTMB_fit <- glmmTMB(incidence / size ~ period + diag(period | herd), data = cbpp, REML = TRUE)
print(summary(glmmTMB_fit))

cat("lme4 fixed effects:\n")
print(fixef(lmer_fit))

cat("glmmTMB fixed effects:\n")
print(fixef(glmmTMB_fit)$cond)

cat("lme4 random effects standard deviations:\n")
print(VarCorr(lmer_fit)$herd)

cat("glmmTMB random effects standard deviations:\n")
print(VarCorr(glmmTMB_fit)$cond$herd)