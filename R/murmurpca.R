matrixCPU <- function(x, y, goal)
{
    p <- ncol(x)
    m <- matrix(ncol=p, nrow=p)

    for(i in 1:p) {
        for(j in 1:p) {
            m[i,j] <- goal(x[,i], x[,j], y)
        }
    }
}

goalSMIFE1 <- function(x1, x2, y) {
  t1 <- entropy(cbind(x1, x2, y))
  t2 <- entropy(x1)
  t3 <- entropy(x2)
  t4 <- entropy(y)
  t5 <- mutinformation(x1, y)
  t6 <- mutinformation(x2, y)
  t7 <- mutinformation(x1, x2)

  res <- t1 - t2 - t3 - t4 + t5 + t6 + t7
  res
}

goalSMIFE2 <- function(x1, x2, y) {
  t1 <- entropy(cbind(x1, x2, y))
  t2 <- entropy(x1)
  t3 <- entropy(x2)
  t4 <- entropy(y)
  t7 <- mutinformation(x1, x2)

  res <- t2 + t3 + t4 - t1 - t7
  res
}

goalMI <- function(x1, x2, y) {
  mutinformation(x1, x2)
}

goalmRRC <- function(x1, x2, y) {
  condentropy(y, x1) + condentropy(y, x2) - condentropy(y, cbind(x1, x2))
}

goalCI <- function(x1, x2, y) {
    condinformation(y, x1, x2) + condinformation(y, x2, x1)
}

genericPCA <- function(x, m) {
  e <- eigen(m)
  
  epos <- e$values + abs(min(e$values))
  # enorm <- epos / sum(epos)

  scores <- x %*% e$vectors

  res <- list(
    sdev = epos / sum(epos),
    loadings = e$vectors,
    n.obs = nrow(x),
    scale = rep(1, ncol(x)),
    scores = scores,
    call = match.call()
    
    )
  class(res) <- "princomp"
       
  res
}

SMIFE1PCA <- function(x, y) {
  m <- matrixCPU(x, y, goalSMIFE1)

  res <- genericPCA(x, m)

  res
}

SMIFE2PCA <- function(x, y) {
  m <- matrixCPU(x, y, goalSMIFE2)

  res <- genericPCA(x, m)

  res
}

mRRSupervisedPCA <- function(x, y) {
  m <- matrixCPU(x, y, goalCI)

  # populate the diagonal
  for(i in 1:ncol(x)) {
    m[i,i] <- 2 * mutinformation(x[,i], y)
  }
  
  # center the matrix around 0
  offset <- abs(min(m))
  m <- m + offset 
  m <- max(m) - m

  res <- genericPCA(x, m)

  res
}

mRRUnsupervisedPCA <- function(x) {
  m <- matrixCPU(x, c(), goalMI)
  
  # ensure that the matrix is symmetrical along
  # the diagonal
  m <- m + t(m) 

  res <- genericPCA(x, m)

  res
}

mrmrPCA <- function(x, y = c(), method="mRRSupervised") {
  res <- c()
  
  if(method == "SMIFE1") {
    res <- SMIFE1PCA(x, y)
  } else if(method == "SMIFE2") {
    res <- SMIFE2PCA(x, y)
  } else if(method == "MI") {
    res <- mRRUnsupervisedPCA(x)
  } else if(method == "CI") {
    res <- mRRSupervisedPCA(x, y, platform)
  }
  
  res
}
