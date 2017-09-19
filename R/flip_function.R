
# function that generates sequence on opposite strand 

flip.function <- function(nucleotide){
  if(nucleotide=="A"){
    return("T")
  }
  if(nucleotide=="T"){
    return("A")
  }
  if(nucleotide=="C"){
    return("G")
  }
  if(nucleotide=="G"){
    return("C")
  }
  if(nucleotide=="N"){
    return("N")
  }
}





