

Shannon <- function (vector) {
  
  # Define corresponding abundance value for each score
  Abundance = c(0.5, 3, 7.5, 15, 25, 35, 45, 55, 65, 75, 85) / 100
  names(Abundance) = c('0.5', '1', '2', '3', '4', '5', '6', '7', '8','9','10')
  
  scores = as.character(vector)
  
  # Vector containing abundance of each species
  a = Abundance[scores]
  print(a)
  # Vector containing abundance of each species as a proportion of total
  p = a/sum(a)
  print(p)
  # Apply formula for Shannon's Diversity Index
  return(-sum(log(p)*p))
}

Shannon(c(1,2,3))

# Test function
dummy_data <- data.frame(
  'Species' = c('Species 1', 'Species 2', 'Species 3', 'Species 4'),
  'Scores' = c(0.5, 3, 2, 0.5)
)

input <- dummy_data$Scores

Shannon(input)
