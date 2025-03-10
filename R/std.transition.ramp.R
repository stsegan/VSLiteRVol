#' "Standard ramp" function for building growth response functions.
#'
#' \code{std.ramp.r} takes a lower and upper bounds, and creates a linear response
#' between 0 and 1 for values in between the bounds. Values below (above) the lower (upper)
#' bound are assigned a value of zero (one).
#'
#' @param x The value at which we want to evaluate the ramp function.
#' @param x1 The lower bound of the support of the nonzero part of the ramp.
#' @param x2 The lower bound of the range of the preimage of 1.
#'
#' @export

# Base VSLite Ramp
std.ramp <- function(x,x1,x2){return(
  apply(
    as.matrix(
      apply(
        (x-x1)/(x2-x1), 1:length(dim(x)), min, 1
      )
    ),
    1:length(dim(x)), max, 0
  )
)}

# Linear
std.ramp.lin <- function(x, x1, x2, k_lin){return(
  apply(
    as.matrix(
      apply(
        k * ((x-x1)/(x2-x1)), 1:length(dim(x)), min, 1
      )
    ),
    1:length(dim(x)), max, 0
  )
)}

# Quadratic
std.ramp.quad <- function(x, x1, x2, k, m) {
  linear_part <- (x - x1) / (x2 - x1)
  quadratic_part <- m*((linear_part)^k)
  return(
    apply(
      as.matrix(
        apply(
          quadratic_part, 1:length(dim(x)), min, 1
        )
      ),
      1:length(dim(x)), max, 0
    )
  )
}

# Sigmoid
std.ramp.sig <- function(x, x1, x2, k_sig, m_sig) {
  linear_part <- (x - x1) / (x2 - x1)
  sigmoid_part <- m_sig * (1 / (1 + exp(-k_sig * (x - x1) / (x2 - x1))))
  return(
    apply(
      as.matrix(
        apply(
          sigmoid_part, 1:length(dim(x)), min, 1
        )
      ),
      1:length(dim(x)), max, 0
    )
  )
}

# Transition Function: Linear to Sigmoid
std.ramp.transition <- function(x, x1, x2, year, year_transition, k_lin, k_sig, m_sig) {
  # Determine which function to use based on the year
  if (year < year_transition) {
    # Linear part
    response <- std.ramp.lin(x, x1, x2, k_lin)
  } else {
    # Sigmoid part
    response <- std.ramp.sig(x, x1, x2, k_sig, m_sig)
  }
  
  return(response)
}
