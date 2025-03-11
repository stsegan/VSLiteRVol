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

# Linear
std.ramp.lin <- function(x, x1, x2, k){return(
  apply(
    as.matrix(
      apply(
        k * ((x-x1)/(x2-x1)), 1:length(dim(x)), min, 1
      )
    ),
    1:length(dim(x)), max, 0
  )
)}


# Sigmoid
std.ramp.sig <- function(x, x1, x2, k, m) {
  linear_part <- (x - x1) / (x2 - x1)
  sigmoid_part <- m * (1 / (1 + exp(-k * (x - x1) / (x2 - x1))))
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
# year_transition: either as an index OR as "real" year, and additionally pass on syear:eyear vector
# version here: second version: years <- syear:eyear
std.ramp.transition <- function(x, x1, x2, years, year_transition, k_lin, k_sig, m_sig) {
  # compute both ramp functions
  sig <- std.ramp.sigstd.ramp.sig(x, x1, x2, k_sig, m_sig)
  lin <- std.ramp.lin(x, x1, x2, k_lin)
  # create output based on cutoff year
  cut_off_index <- which(years == year_transition)
  response <- numeric(length(years))
  response[1:cut_off_index] <- sig[1:cut_off_index]
  response[(cut_off_index + 1):length(years)] <- lin[(cut_off_index + 1):length(years)]
  return(response)
}
