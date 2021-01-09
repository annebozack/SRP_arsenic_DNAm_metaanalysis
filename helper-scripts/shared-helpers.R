### Shared Helper Functions -- DVP, DMP ###

## Lambda function 
lambda <- function(p){
  
  median(qchisq(p, df = 1, lower.tail = FALSE),
         na.rm = TRUE) / qchisq(0.5, df = 1)
  
}