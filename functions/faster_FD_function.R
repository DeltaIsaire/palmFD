#########################################
# faster, streamlined version of FD::dbFD
#########################################
#
# The original function dbFD() from package 'FD' does a lot of things we do not
# need it to do. Runtime is a constraint for us, hence it is worthwile to rewrite
# this function to do only what it strictly needs to for our purposes.


library(FD)

source(file = "functions/edited_fdisp_function.R")


test <- FALSE

# Extract original source code
# ----------------------------
if (test) {
  capture.output(getAnywhere("dbFD"), file = "functions/source_dbFD.R")
}

# Input data to test with
# -----------------------
if (test) {
  trait.matrix <- data.frame(a.trait = 1:10,
                             b.trait = 15:6,
                             c.trait = c(7, 8, 9, 8, 6, 4, 9, 6, 7, 5),
                             row.names = paste0("species.", letters[1:10])
                             )
  trait.matrix <- as.matrix(trait.matrix)

  set.seed(123)
  pres.abs.matrix <- matrix(ncol = nrow(trait.matrix),
                            nrow = 10,
                            dimnames = list(row = paste0("area.", letters[1:10]),
                                            col = rownames(trait.matrix)
                                            ),
                            data = round(runif(10 * nrow(trait.matrix)))
                            )
  pres.abs.matrix[6, c(1, 6)] <- 1
  pres.abs.matrix[10, c(1, 6)] <- 1
}

# Function rewrite
# ----------------
# We only need two outputs: FRic and FDis

FastFD <- function(trait.matrix, pres.abs.matrix, pcoa.traits = NULL) {
# No checks. Assume input is valid.
# Trait.matrix: a matrix with trait values in columns, and rownames as species.
# pres.abs.matrix: a matrix with presence (1) and absense (0) of species
#                  (colnames) in tdwg3 units (rownames).
# nrow(trait.matrix) must equal ncol(pres.abs.matrix)
# rownames(trait.matrix) must equal colnames(pres.abs.matrix)
# Neither matrix should have NAs.
# No community should have zero species, and all species should occur somewhere.
# There shall be exactly 3 traits, all numeric.
# All traits will be standardized to mean 0 and unit variance.
# All communities must have > 3 species with unique trait combinations
#
# NOTE: The FRic convex hull calculation uses 'joggling' to handle precision errors,
#       by passing option 'QJ' to the qhull algorithm called in convhulln().
#       see http://www.qhull.org/html/qh-impre.htm for more information.
#
# pcoa.traits: precomputed PcoA-traits to use for FRic calculation.
  x.rn <- rownames(trait.matrix)

  c <- dim(pres.abs.matrix)[1]

  trait.matrix <- data.frame(trait.matrix)
  x.dist <- dist(trait.matrix)
  attr(x.dist, "Labels") <- x.rn

  if (is.null(pcoa.traits)) {
    x.pco <- dudi.pco(x.dist, scannf = FALSE, full = TRUE)
    traits <- round(x.pco$li, .Machine$double.exponent)
  } else {
    traits <- pcoa.traits
  }

  if (ncol(traits) == 1) {
    traits.FRic <- traits
    cat("FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume.\n")
  }
  if (ncol(traits) > 1) {
    traits.FRic <- traits
  }
  
  disp <- CustomFdisp(x.dist, pres.abs.matrix)
  FDis <- disp$FDis
  
  FRic <- rep(NA, c)
  names(FRic) <- row.names(pres.abs.matrix)
  
  for (i in seq_len(c)) {
    sppres <- which(pres.abs.matrix[i, ] > 0)
    S <- length(sppres)
    tr <- data.frame(traits[sppres, ])
    tr.FRic <- data.frame(traits.FRic[sppres, ])
    ab <- as.matrix(pres.abs.matrix[i, sppres])
    abundrel <- ab / sum(ab)
    if (dim(tr.FRic)[2] > 1) {
      thresh <- 3
      convhull <- convhulln(tr.FRic, options = c("FA", "QJ"))
      FRic[i] <- convhull$vol
    }
    if (dim(tr.FRic)[2] == 1) {
      tr.range <- range(tr.FRic[, 1])
      t.range <- tr.range[2] - tr.range[1]
      FRic[i] <- t.range
    }
  }

  res <- vector("list", length = 2)
  names(res) <- c("FRic", "FDis")
  res$FRic <- FRic
  res$FDis <- FDis

  invisible(res)
}


# Test new function: compare to original
# --------------------------------------
if (test) {
  system.time(new <- FastFD(trait.matrix, pres.abs.matrix))
  #   user  system elapsed 
  #  0.020   0.000   0.018 

  system.time( {
    old <- dbFD(x = trait.matrix,
                a = pres.abs.matrix,
                w.abun = FALSE,
                stand.x = FALSE,
                corr = "cailliez",
                calc.FRic = TRUE,
                m = "max",
                stand.FRic = FALSE,
                calc.CWM = FALSE,
                calc.FDiv = FALSE
                )
    old <- old[c("FRic", "FDis")]
  }
  )
  #   user  system elapsed 
  #  0.048   0.004   0.052 

  cat("New output matches old output EXACTLY:", identical(new, old), "\n")
  cat("New output matches old output NEARLY:", isTRUE(all.equal(new, old)), "\n")
  # FALSE
  # TRUE

  # Evaluation:
  # For this test dataset, the new function works as expected,
  # and is significantly faster. Mission Accomplished.
}

