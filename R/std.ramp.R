#' "Standard ramp" function for building growth response functions.
#' 
#' \code{std.ramp.r} takes a lower and upper bounds, and creates a linear response
#' between 0 and 1 for values in between the bounds.  Values below (above) the lower (upper) 
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


# Transition ramp
std.ramp.transition <- function(x, x1, x2, year, syear, eyear, k_lin, k_sig, m_sig) {
  
  # Convert x (temperature or moisture) to a corresponding year index
  nyrs <- length(syear:eyear)  # Total number of years
  year_index <- (year - syear + 1) / nyrs  # Normalize year to a 0-1 scale
  
  # Define the linear part
  linear_part <- k_lin * ((x - x1) / (x2 - x1))
  
  # Define the sigmoid part
  sigmoid_part <- m_sig * (1 / (1 + exp(-k_sig * ((x - x1) / (x2 - x1)))))
  
  # Define the weight for the transition
  # The weight smoothly transitions from 0 to 1 around the transition point
  transition_width <- 0.1  # Adjust this for smoother or sharper transitions
  weight <- 1 / (1 + exp(-10 * (year_index - year) / transition_width))
  
  # Combine the two parts using the weight
  combined_response <- (1 - weight) * linear_part + weight * sigmoid_part
  
  # Ensure the response is bounded between 0 and 1
  return(
    apply(
      as.matrix(
        apply(combined_response, 1:length(dim(x)), min, 1)
      ),
      1:length(dim(x)), max, 0
    )
  )
}

# Where k = gradient/steepness of slope, m = height of slope.
# Would have to change gT/gM to include k and m values. 
# Need to configure this to either build a 'master function' which includes
# all possible slopes (i.e. you select a 'slope type' when inputting your data),
# or just define a function for each shaped slope (current state).
# Could also add other slope types (maybe a step-wise graph), but not sure if 
# they are more relevant than the ones already included. 
