library(ggplot2)

# This function simulates an association z-score matrix,
# masks a specified percentage of z-scores, and
# imputes them using either matrix completion or KOMPUTE
#' @param n.genes The number of genes to simulate
#' @param pheno.cor The m by m phenotypic correlation matrix of the m phenotypes considered, functioning as Sigma in the MVN(0, Sigma) simulation
#' @param mask.prop The proportion of z-scores that should be masked for the imputation
#' @param method The imputation method to use; either "mc" for matrix completion or "kompute"
#' @param info.cutoff The minimum imputation information score necessary to be included; 0 indicates all values are kept regardless of info score (all are always kept for the matrix completion method)
#' @param seed The seed number for random number generation
#' @return A list containing two objects: the plot of imputed vs original z-scores, and the correlation coefficient between imputed and original z-scores


simulation <- function(n.genes, pheno.cor, mask.prop, method="mc", info.cutoff=0, seed=123){

  # Set the seed for reproducibility in simulation
  set.seed(seed)

  # Define the number of phenotypes based on the dimensions of the correlation matrix
  npheno <- nrow(pheno.cor)

  # Simulate association Z-score matrix
  sim.z <- mvrnorm(n=n.genes, mu=rep(0, npheno), Sigma=pheno.cor)
  sim.z <- t(sim.z)

  # Mask a specified percentage of measured z-scores
  nimp <- nrow(sim.z)*ncol(sim.z)*mask.prop
  all.i <- 1:(nrow(sim.z)*ncol(sim.z))

  mask.i <- sort(sample(all.i, nimp))
  org.z = as.matrix(sim.z)[mask.i]
  zvec <- as.vector(as.matrix(sim.z))
  zvec[mask.i] <- NA
  zmat.imp <- matrix(zvec, nrow=npheno)
  rownames(zmat.imp) <- rownames(sim.z)


  # Apply the specified imputation method
  r <- 6
  if(method=="mc"){ # Matrix complete method
    mc.res <- svd.impute(zmat.imp, r)
    org.z <- sim.z[mask.i]
    imp.z <- mc.res[mask.i]

  } else if(method=="kompute"){ # KOMPUTE method
    kompute.res <- kompute(t(zmat.imp), pheno.cor, 0.01)

    imp.z <- as.matrix(t(kompute.res$zmat))[mask.i]
    imp.info <- as.matrix(t(kompute.res$infomat))[mask.i]

    imp <- data.frame(org.z=org.z, imp.z=imp.z, info=imp.info)
    imp <- imp[complete.cases(imp),]
    imp <- subset(imp, info>=0 & info <= 1)

    # Apply imputation information cutoff
    imp.sub <- subset(imp, info>info.cutoff)
    org.z <- imp.sub$org.z
    imp.z <- imp.sub$imp.z
    info <- imp.sub$info
  }

  # Calculate correlation coefficient between imputed and original z-scores
  cor.val <- round(cor(imp.z, org.z), digits=3)
  type <- ifelse(method=="mc", "Matrix Completion", "KOMPUTE")

  # Create plot
  if(info.cutoff==0){
    plot <- ggplot() +
      geom_point(aes(x=imp.z, y=org.z), alpha=0.1) +
      labs(title=paste0(type, ", ", mask.prop*100, "% Removed, Cor=", cor.val),
           x="Imputed Z-scores", y = "Measured Z-scores") +
      theme_minimal()

  } else{

    plot <- ggplot() +
      geom_point(aes(x=imp.z, y=org.z, col=info), alpha=0.1) +
      labs(title=paste0(type, ", ", mask.prop*100, "% Removed, Cor=", cor.val),
           x="Imputed Z-scores", y = "Measured Z-scores") +
      theme_minimal()
  }

  return(list(plot=plot, cor=cor.val))
}
