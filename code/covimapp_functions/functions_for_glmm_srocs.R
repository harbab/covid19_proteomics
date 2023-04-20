##### Functions for SROC #####
##### FUNCTIONS #####
logit <- function(x) {
  y = log(x/(1-x))
  return(y)
}

#------------Inverse logit---------

inv.logit <- function(x){
  y = exp(x)/(1+exp(x))
  return(y)
}

# Reshaping data in long format for glmer input
reshape.data <- function(data){
  
  groupD <- data.frame(study=data$study, 
                       wellclassified=data$TP, 
                       misclassified=data$FN, 
                       group="disease",
                       model = 1:dim(data)[1])
  # We excluded the model_id variable, because we model sensitivity and specificity based only on the study
  
  groupH <- data.frame(study=data$study, 
                       wellclassified=data$TN, 
                       misclassified=data$FP, 
                       group="healthy",
                       model = 1:dim(data)[1])
  
  data_long <- bind_rows(groupD, groupH) %>%
    as.data.frame()
  
  return(data_long)
  
}


# Glmer output
tab.glmer <- function(m) {
  
  est <- m@beta
  conf <- confint(m, method = 'Wald')[c('groupdisease', 'grouphealthy'),]
  std_study <- attr(VarCorr(m)$study, 'stddev')
  # We excluded the model variable, because we model sensitivity and specificity based only on the study
  corr_study <- attr(VarCorr(m)$study, 'correlation')[1,2]
  n_study <- length(unique(factor(m@flist$study)))
  
  res <- cbind(est, conf) %>%
    as.data.frame() %>%
    dplyr::mutate_all(inv.logit)
  
  tab <- cbind(res, std_study, corr_study, n_study)
  rownames(tab) <-  c('Sensitivity', 'Specificity')
  tab
}

#------------------- SROC functions---------------------------

# Extract results from glmer output
fit.glmer <- function(model, data) {
  
  fit <- list(coefficients = fixef(model),
              Psi = VarCorr(model)$study, 
              # We excluded the Omega variable, because we model sensitivity and specificity based only on the study
              alphasens = 1,
              alphafpr = 1,
              logLik = logLik(model),
              freqdata = data[c('TP','FP','TN','FN')],
              vcov = matrix(vcov(model)@x,nrow=2),
              sens = data$TP/(data$TP+data$FN),
              fpr = data$FP/(data$FP+data$TN))
  
  # Change the sign, so that logit FPR is obtained
  fit$coefficients[2] <- - fit$coefficients[2] 
  # Change the sign 
  fit$Psi[1,2] <- fit$Psi[2,1] <- - fit$Psi[2,1]
  
  fit$logLik <- NA
  attr(fit$logLik, "df") <- 5
  colnames(fit$vcov) <- c('groupdisease', 'grouphealthy')
  
  return(fit)
  
}

# Functions from mada
# Logit transformation
trafo <- function(alpha, x){return(talpha(alpha)$linkfun(x))}
# Inverse logit transformation
inv.trafo <- function(alpha, x){return(talpha(alpha)$linkinv(x))}

# Modified mada functions for glmer output
calc_hsroc_coef <- function(fit){
  coef <- fit$coefficients
  coef <- as.numeric(coef)
  
  Psi <- fit$Psi  
  ran.sd <- sqrt(diag(Psi))
  
  if(attr(fit$logLik,"df")== 5){
    # HSROC parameters (formulae from Harbord et al. 2008)
    Theta <- 0.5*(sqrt(ran.sd[2]/ran.sd[1])*coef[1] + sqrt(ran.sd[1]/ran.sd[2])*coef[2]) 
    Lambda <- sqrt(ran.sd[2]/ran.sd[1])*coef[1] - sqrt(ran.sd[1]/ran.sd[2])*coef[2] 
    sigma2theta <- 0.5*(ran.sd[1]*ran.sd[2] + Psi[1,2]) 
    sigma2alpha <- 2*(ran.sd[1]*ran.sd[2] - Psi[1,2])  
    beta <- log(ran.sd[2]/ran.sd[1])             
    coef_hsroc <- list(Theta = Theta,
                       Lambda = Lambda,
                       beta = beta,
                       sigma2theta = sigma2theta,
                       sigma2alpha = sigma2alpha)
    coef_hsroc <- lapply(coef_hsroc, function(x){attr(x, "names") <- NULL; x})
  }else{
    coef_hsroc = NULL
  }
  if(is.null(coef_hsroc)){
    warning("Can only compute coefficients for SROC curves without covariates. Returning NULL.")
  }
  return(coef_hsroc)
}

# SROC function for glmer output
sroc.glmer <- function(fit, fpr = 1:99/100, type = "ruttergatsonis",
                       return_function = FALSE, ...){
  stopifnot(is.logical(return_function))
  stopifnot(type %in% c("ruttergatsonis", "naive"))
  if(type == "ruttergatsonis"){
    coef_hsroc <- calc_hsroc_coef(fit)
    Lambda <- coef_hsroc$Lambda    
    Beta <- coef_hsroc$beta
    f <- function(x){
      return(inv.trafo(fit$alphasens, (Lambda*exp(-Beta/2) + exp(-Beta)*trafo(fit$alphafpr, x))))
    }
    sens <- f(fpr)
    if(!return_function){
      return(cbind(fpr, sens))
    }else{
      return(f)
    }
  }
}


# ROC ellipse function for glmer output
ROC.ellipse2 <- function(fit, conf.level = 0.95, pch = 1, add = TRUE, 
                         predict = TRUE, predlty = 3, predlwd = 1, predcol = 1, ...)
{
  alpha.sens <- fit$alphasens
  alpha.fpr <- fit$alphafpr
  mu <- fit$coefficients
  Sigma <- fit$vcov
  
  vcov <- fit$vcov
  Psi <- fit$Psi
  Omega <- fit$Omega
  
  talphaellipse <- ellipse::ellipse(Sigma, centre = mu, level = conf.level)
  ROCellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse))
  ROCellipse[,1] <- inv.trafo(alpha.fpr, talphaellipse[,2])
  ROCellipse[,2] <- inv.trafo(alpha.sens, talphaellipse[,1])
  
  if(predict) {
    if(is.null(Omega)) {
      Sigma_pred <- Psi + vcov
    } else {
      Sigma_pred <- Psi + Omega + vcov
    }
    talphaellipse_pred <- ellipse::ellipse(Sigma_pred, centre = mu, level = conf.level)
    predellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse_pred))
    predellipse[,1] <- inv.trafo(alpha.fpr, talphaellipse_pred[,2])
    predellipse[,2] <- inv.trafo(alpha.sens, talphaellipse_pred[,1])
    
  }
  
  if(add){
    lines(ROCellipse, ...)
    points(inv.trafo(alpha.fpr, mu[2]), 
           inv.trafo(alpha.sens, mu[1]), pch = pch, ...)
    lines(predellipse, lty = predlty, lwd = predlwd, col=predcol) 
    return(invisible(NULL))
  }
  if(!add){
    return(list(ROCellipse = ROCellipse, 
                fprsens = matrix(c(inv.trafo(alpha.fpr, mu[2]), 
                                   inv.trafo(alpha.sens, mu[1])),nrow = 1)))
  }
  
}



ROCellipse.glmer <- function(x, level = 0.95, add = FALSE, pch = 1,
                             predict = FALSE, predlty = 3, predlwd = 1, predcol = 1, ...){
  ROC.ellipse2(x, conf.level = level, add = add, pch = pch,  
               predict = predict,  predlty = predlty, predlwd = predlwd, predcol = predcol, ...)
  
}

# Plot SROC only on observed data
sroc.glmer.range <- function(fit) {
  
  min_s <- min(fit$sens)
  max_f <- max(fit$fpr)
  
  s <- sroc.glmer(fit) %>%
    as.data.frame() %>%
    filter(fpr <= max_f & sens >= min_s)
}

##### Modified mada functions by us #####
### Calculate DOR
DOR_calc = function(sens, fpr){
  sens*(1-fpr)/((1-sens)*fpr)
}

### Zwindermann & Bossuyt (2008) MCMC procedure to generate summary points (positive and negative likelihood ratio, diagnostic odds ratio) for the Reitsma et al. (2005) bivariate model
SummaryPts.glmer <- function(object, n.iter = 10^6, FUN = NULL, ...){
  SummaryPts.default <- function(object, mu,Sigma,alphasens = 1, alphafpr = 1,
                                 n.iter = 10^6, FUN, ...){
    samples <- rmvnorm(n.iter, mu, Sigma)
    sens <- inv.trafo(alphasens,samples[,1])
    fpr <- inv.trafo(alphafpr,samples[,2])
    out <- lapply(FUN, function(x){x(sens, fpr)})
    class(out) <- "SummaryPts"
    out
  }
  
  fit <- object
  if(length(coef(fit)) > 2){
    stop("SummaryPts is not be used for meta-regression!")}
  if(is.null(FUN)){FUN <- list(posLR = function(sens,fpr){sens/fpr},
                               negLR = function(sens,fpr){(1-sens)/(1-fpr)},
                               invnegLR = function(sens, fpr){(1-fpr)/(1-sens)},
                               DOR = function(sens, fpr){sens*(1-fpr)/((1-sens)*fpr)})}
  SummaryPts.default(mu = coef(fit)[1:2], Sigma = fit$vcov, # vcov(fit) didn't work for our output
                     alphasens = fit$alphasens, 
                     alphafpr = fit$alphafpr,
                     n.iter = n.iter, FUN = FUN)
}

### AUC calculation
AUC_calc = function(tpr) {
  n = length(tpr)
  (tpr[1]/2 + sum(tpr[2:(n-1)]) + tpr[n]/2)/n
}