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
  
  return(resu)
  
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


#==========/==========/==========/==========/==========/==========/==========/==========/
# Tabelas para o latex


library(xtable)
library(tidyverse)
library(kableExtra)

n <- c(10, 30, 50, 80)

arquivos <- paste0("results_simu/summary_simu_2DMxARMA_nk", n ,".txt")


arquivo <- read.table(file = arquivos[1], header = T)

arquivo <- arquivo |> 
  select(Parametro, Media, RB, EQM)

arquivo <- arquivo |> mutate(
  Parametro = paste0("$", Parametro, "$"),
  Media    = paste0("$", format(Media, nsmall = 2), "$"),
  RB       = paste0("$", format(RB, nsmall = 2), "$"),
  EQM      = paste0("$", format(EQM, nsmall = 2), "$"),
)

arquivo$Measures <- recode(row.names(arquivo),
                            "alpha" = "$\\alpha$",
                            "phi1" = "$\\phi_{0,0}$",
                            "phi2" = "$\\phi_{1,0}$",
                            "phi3" = "$\\phi_{0,1}$",
                            "theta1" = "$\\theta_{0,0}$",
                            "theta2" = "$\\theta_{1,0}$",
                            "theta3" = "$\\theta_{0,1}$")

arquivo <- arquivo %>%
  relocate(Measures, .before = 1)

rownames(arquivo) <- NULL


tabela_latex <- arquivo |>
  kbl(format = "latex", booktabs = TRUE, escape = FALSE,
      align = "lrrrr",
      col.names = c("Measures", "Parameter", "Mean", "RB(\\%)", "MSE")) 


# Exibir o código LaTeX formatado
cat(tabela_latex)


class(arquivo)
