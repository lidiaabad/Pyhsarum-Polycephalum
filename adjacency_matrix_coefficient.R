source("C:/Users/labad/Documents/5biotec/tfg/programming/la.functions.R")
metro= as.matrix(read.table("C:/Users/labad/Documents/5biotec/tfg/programming/adj_mat_metro.txt"))
renfe=as.matrix(read.table("C:/Users/labad/Documents/5biotec/tfg/programming/adj_mat_renfe.txt"))
metro_real = matrix(0, nrow=nrow(metro), ncol=ncol(metro))
metro_pp= matrix(0 , nrow=nrow(metro), ncol=ncol(metro))
renfe_real = matrix(0, nrow=nrow(renfe), ncol=ncol(renfe))
renfe_pp= matrix(0 , nrow=nrow(renfe), ncol=ncol(renfe))
#correlation matrix
la.dist(metro[1:nrow(metro), 1:ncol(metro)], method="pearson")
la.dist(renfe[1:nrow(renfe), 1:ncol(renfe)], method="pearson")

for (i in 1:nrow(metro)){ 
  for (j in 1:ncol(metro)){
    if (i<j){ 
      val=metro[i,j]
      metro_real[i,j]=val
      metro_real[j,i]=val
    }else{ 
      val=round(metro[i,j]) #dichotomyze
      metro_pp[i,j]=val
      metro_pp[j,i]=val
    }
  }
}

for (i in 1:nrow(renfe)){ 
  for (j in 1:ncol(renfe)){
    if (i<j){ 
      val=renfe[i,j]
      renfe_real[i,j]=val
      renfe_real[j,i]=val
    }else{ 
      val=round(renfe[i,j]) #dichotomyze
      renfe_pp[i,j]=val
      renfe_pp[j,i]=val
    }
  }
}

la.NetworkCompare(metro_real, metro_pp)
la.NetworkCompare(renfe_real, renfe_pp)

