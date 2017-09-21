

# This code was found at 
# https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4 - Thanks Brian Diggs! 
# And discussed and edited here: https://stackoverflow.com/a/24241954/8376488 - Thanks Jack Aidley! 




fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}