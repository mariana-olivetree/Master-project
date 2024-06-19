################### Screen gene from Cho et al.

screen.gene <- function(data, avg.detect = 1, gene.col = NULL) {
  
  # data: each column = sample, each row = genes
  # gene.col: the column where gene names are provided. If null, ignored.
  
  gene.name = data[,gene.col]
  if (!is.null(gene.col)) data = data[,-gene.col]
  n.genes = dim(data)[1]
  n.sample = dim(data)[2]
  
  useF = which(rowSums(data) > avg.detect * n.sample)
  print(paste(length(useF), "out of ", n.genes, " will be used."))
  list(sample.size = n.sample, 
       stat = c(use = length(useF), out.of = n.genes, percentage = round(length(useF)/n.genes,2)*100),
       feature = list(index = useF, names = if (is.null(gene.col)) NA else gene.name[useF]))
}
##################


################### LB Function from Cho et al.

LB <- function (data, fin.n = 1000, fin.z = FALSE) {
  # fin.n: minimum nonzero counts for nonzero test
  # fin.z: logical. if true, correct beta hat of logistic regression
  require(gamlss)
  n = dim(data)[1]
  data$y.prop <- data$y / data$sampleSum #sampleSum[j] = sum_g y_g,j (g:gene, j:sample)
  out <- matrix(NA, 3, 2, dimnames = list(c("LB.nonz", "LB.zero", "LB.glob"), c("Estimate", "pval")))
  drop.nonzero <- TRUE       # sum(data$y > 0, na.rm = TRUE) < fin.
    
    bereg = try(glm(ifelse(y.prop > 0, 1, 0) ~ phenotype, family = binomial, data = data))
    
    
    if (any(class(bereg) %in% "try-error")) return(out)
    tab.tmp <- summary(bereg)$coef
    
    if (fin.z) { #finite sample correction for zero models
      pred <- bereg$fitted.values
      X <- model.matrix(bereg)
      XDX.inv <- solve(t(X) %*% (X * pred))
      Q <- X %*% XDX.inv %*% t(X)
      xi <- diag(Q) * (pred - 0.5)
      bias <- XDX.inv %*% t(X) %*% (pred * (1-pred) * xi)
      
      tab.tmp[, "Estimate"] <- tab.tmp[, "Estimate"] - bias
      tab.tmp[, "Std. Error"] <- tab.tmp[, "Std. Error"] * n/(2 * n - bereg$df.residual)  # df = 3 = intercept + Pheno + batch (n + k = n + (n - df.res))
      tab.tmp[, "Pr(>|z|)"] <- 2 * pt(-abs(tab.tmp[, "Estimate"] / tab.tmp[, "Std. Error"]), df = bereg$df.residual)
    }
    
    out["LB.zero", ] <- tab.tmp["phenotype1", c("Estimate", "Pr(>|z|)")]
    out["LB.glob", ] <- out["LB.zero", ]
    
  }
###############################




############################## LB without batch

LB_nobatch <- function (data, fin.n = 1000, fin.z = FALSE) {
  # fin.n: minimum nonzero counts for nonzero test
  # fin.z: logical. if true, correct beta hat of logistic regression
  require(gamlss)
  n = dim(data)[1]
  data$y.prop <- data$y / data$sampleSum #sampleSum[j] = sum_g y_g,j (g:gene, j:sample)
  out <- matrix(NA, 3, 2, dimnames = list(c("LB.nonz", "LB.zero", "LB.glob"), c("Estimate", "pval")))
  drop.nonzero <- TRUE       # sum(data$y > 0, na.rm = TRUE) < fin.
  
  bereg = try(glm(ifelse(y.prop > 0, 1, 0) ~ 1, family = binomial, data = data))
  
  
  if (any(class(bereg) %in% "try-error")) return(out)
  tab.tmp <- summary(bereg)$coef
  
  if (fin.z) { #finite sample correction for zero models
    pred <- bereg$fitted.values
    X <- model.matrix(bereg)
    XDX.inv <- solve(t(X) %*% (X * pred))
    Q <- X %*% XDX.inv %*% t(X)
    xi <- diag(Q) * (pred - 0.5)
    bias <- XDX.inv %*% t(X) %*% (pred * (1-pred) * xi)
    
    tab.tmp[, "Estimate"] <- tab.tmp[, "Estimate"] - bias
    tab.tmp[, "Std. Error"] <- tab.tmp[, "Std. Error"] * n/(2 * n - bereg$df.residual)  # df = 3 = intercept + Pheno + batch (n + k = n + (n - df.res))
    tab.tmp[, "Pr(>|z|)"] <- 2 * pt(-abs(tab.tmp[, "Estimate"] / tab.tmp[, "Std. Error"]), df = bereg$df.residual)
  }
  
  out["LB.zero", ] <- tab.tmp[, c("Estimate", "Pr(>|z|)")]
  out["LB.glob", ] <- out["LB.zero", ]
  
}
