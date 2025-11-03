library(dplyr)
library(tidyr)


# Função para resumir os resultados


resumindo <- function(arquivo, pars_true, ic = 0.05, n){
  
  tabela <- read.table(arquivo, header = TRUE)
  
  resu <- data.frame()
  
  media_pars <- colMeans(tabela, na.rm = T)
  var_pars <- apply(tabela, 2, var, na.rm = T)


  resu <- data.frame(
      Parametro = pars_true,
      Media = media_pars,
      RB = ((media_pars - pars_true) / pars_true) * 100, 
      Var = var_pars,
      EQM = ((media_pars - pars_true)^2 + var_pars),
      MAE = colMeans(abs(tabela - pars_true), na.rm = T)
      )
    
    resu <- resu |> mutate(across(where(is.numeric), round, 4))
    
  
  summary_filename <- paste0( "results_simu/", "summary_", arquivo)
  
  write.table(resu, file = summary_filename, row.names = T, quote = FALSE)
  
}


pars_true <- c(-1.2, 0.03, 0.35, 0.25, -0.1, -0.06, -0.008)

# nk 10

arquivo = "simu_2DMxARMA_nk10.txt"

tabela <- read.table(arquivo, header = T)
resumindo(arquivo, pars_true)


#nk 30

arquivo = "simu_2DMxARMA_nk30.txt"
resumindo(arquivo, pars_true)


#nk 50

arquivo = "simu_2DMxARMA_nk50.txt"
resumindo(arquivo, pars_true)


#nk 80

arquivo = "simu_2DMxARMA_nk80.txt"
resumindo(arquivo, pars_true)

