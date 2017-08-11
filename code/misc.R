# Shorthand for machine precision.
eps <- .Machine$double.eps

# Return the n x n identity matrix.
eye <- function (n)
  diag(rep(1,n))

# Return a sparse n x n identity matrix.
speye <- function (n)
  .symDiagonal(n)

# Return a m x n sparse matrix of all zeros.
spzeros <- function (m, n)
  sparseMatrix(dims = c(m,n),i = NULL,j = NULL)

# Subtract b[i] from each column A[,i].
subtract.cols <- function (A, b)
  t(t(A) - b)

# Scale each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)
    
# Return the quadratic norm (L2-norm) of vector x.
norm2 <- function (x)
  sqrt(sum(x^2))

# Return a sparse matrix with diagonal entries x.
spdiag <- function (x)
  .symDiagonal(length(x),x)

