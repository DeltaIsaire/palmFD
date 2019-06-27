#############################
# Edited version of FD::fdisp
#############################
#
# Removed the input checks, allowing the function to be run on presence/absence
# matrices with some colSums of zero. This does not affect results.

library(FD)

test <- FALSE

# Extract original source code
# ----------------------------
if (test) {
  capture.output(getAnywhere("fdisp"), file = "functions/source_fdisp.R")
}


# Function rewrite
# ----------------
CustomFdisp <- function(d, a, tol = 1e-07) {
  n <- attr(d, "Size")
  com <- nrow(a)

  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
      r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
      r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
      pres <- which(a[i, ] > 0)
      nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
      if (nb.sp >= 2) {
          w <- a[i, pres]
          centroid <- apply(vec, 2, weighted.mean, w = w)
          dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
          dist.pos <- rowSums(dist.pos^2)
          if (any(!pos)) {
              dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
              dist.neg <- rowSums(dist.neg^2)
          }
          else dist.neg <- 0
          zij <- sqrt(abs(dist.pos - dist.neg))
          avg.dist.cent[i] <- weighted.mean(zij, w)
      }
      else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}


