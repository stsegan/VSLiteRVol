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
  transition_width <- 0.001  # Adjust this for smoother or sharper transitions
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