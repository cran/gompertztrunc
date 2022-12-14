#' Helper functions
#'
#' @keywords internal
#'

## parameterize gompertz with M (mode) and b (slope)

## we use flexsurv functions and modify just slightly

# convert gompertz b an M params to a param
bM2a <- function(b, M) {
  a <- b * exp(-b * M)
  a
}

# convert gompertz a and b params to M
ab2M <- function(a, b) {
  M <- -log(a / b) / b
  M
}


# cdf (-inf, q) of gompertz distribution with modal parameterization
pgompertz.M <- function(q, b, M, ...) {
  a <- bM2a(b, M)
  flexsurv::pgompertz(q, shape = b, rate = a, ...)
}

# pdf of gompertz distribution with modal parameterization
dgompertz.M <- function(x, b, M, ...) {
  a <- bM2a(b, M)
  flexsurv::dgompertz(x, shape = b, rate = a, ...)
}

# randomly generate n draws from gompertz distribution with modal parameterization
rgompertz.M <- function(n, b, M) {
  a <- bM2a(b, M)
  flexsurv::rgompertz(n, shape = b, rate = a)
}

hgompertz.M <- function(x, b, M) {
  a <- bM2a(b, M)
  q <- x
  p <- flexsurv::pgompertz(q, shape = b, rate = a)
  l <- 1 - p
  d <- flexsurv::dgompertz(x, shape = b, rate = a)
  h <- d / l
  return(h)
}

#########

## Old Gompertz Files


getMode <- function(alpha, beta) {
  M <- (1 / beta) * log(beta / alpha)
  return(M)
}

getAlpha <- function(M, beta) {
  alpha <- beta / exp(M * beta)
  return(alpha)
}


# computes lifetable hazards for gompertz distribution
get.trunc.mean.gomp <- function(alpha, b, l, u) {
  x <- 0:110
  hx <- alpha * exp(b * x)
  Hx <- cumsum(hx)
  lx <- c(1, exp(-Hx))
  sum(lx)
  dx <- -diff(lx)
  s <- x %in% l:u ## approximate truncation for these cohorts
  sum((dx * x)[s]) / sum(dx[s])
}

## hx calc
hx_calc <- function(x, b, M) {
  hx = b * exp(b*(x - M))
  return(hx)
}

## lx calc
lx_calc <- function(hx) {
  lx <- 1 * exp(-cumsum(hx))
  return(lx)
}

## dx calc
dx_calc <- function(n, lx) {
  lx <- n * lx
  dx <-  lx - dplyr::lead(lx)
  return(dx)
}

