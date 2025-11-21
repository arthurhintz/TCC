library(tidyverse)

# Função para resumir os resultados

resumindo <- function(arquivo, pars_true, ic = 0.05, n){
  
  tabela <- read.table(arquivo, header = TRUE)
  
  resu <- data.frame()
  
  media_pars <- colMeans(tabela, na.rm = T)
  var_pars <- apply(tabela, 2, var, na.rm = T)


  resu <- data.frame(
      Parametro = pars_true,
      Media = media_pars,
      #RB = ((media_pars - pars_true) / pars_true) * 100,
      RMSE = sqrt(colMeans((tabela-pars_true)^2, na.rm = T)),
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


#nk 100


arquivo = "simu_2DMxARMA_nk100.txt"
resumindo(arquivo, pars_true)



#==========/==========/==========/==========/==========/==========/==========/==========/
# Tabelas para o latex


library(xtable)
library(tidyverse)
library(kableExtra)

n <- c(10, 30, 50, 100)

lista_tabelas <- list()

for(i in seq_along(n)){
  
  k <- n[i]
  
  arquivo_nome <- paste0("results_simu/summary_simu_2DMxARMA_nk", k ,".txt")
  
  arquivo <- read.table(file = arquivo_nome, header = TRUE)
  
  rn <- row.names(arquivo)
  
  # Vetor com os símbolos LaTeX sem o valor numérico
  latex_symbols <- c(
    alpha  = "\\alpha",
    phi1   = "\\phi_{0,0}",
    phi2   = "\\phi_{1,0}",
    phi3   = "\\phi_{0,1}",
    theta1 = "\\theta_{0,0}",
    theta2 = "\\theta_{1,0}",
    theta3 = "\\theta_{0,1}"
  )
  
  arquivo <- arquivo |>
    mutate(
      Mean = paste0("$", format(Media, nsmall = 2), "$"),
      RMSE  = paste0("$", format(RMSE, nsmall = 2), "$"),
      EQM   = paste0("$", format(EQM, nsmall = 2), "$"),
      Measures = paste0(
        "$", latex_symbols[rn], " = ", format(Parametro, nsmall = 2), "$"
      ),
      n = k
    ) |> 
    select(Measures, n, Mean, RMSE, EQM)
  
  lista_tabelas[[i]] <- arquivo
}


lista_tabelas

# Junta tudo em uma única tabela
tabela_final <- bind_rows(lista_tabelas)


tabela_final <- tabela_final %>%
  relocate(Measures, n,.before = 1) |> 
  arrange(Measures)


rownames(tabela_final) <- NULL


tabela_latex <- tabela_final |>
  kbl(format = "latex", booktabs = TRUE, escape = FALSE,
      align = "lrrrr",
      col.names = c("Measures", "N = M", "Mean", "RMSE", "MSE")) 


# Exibir o código LaTeX formatado
cat(tabela_latex)


#==========/==========/==========/==========/==========/==========/==========/==========/
# Gráfico do RMSE


n <- c(10, 30, 50, 80, 100)

lista_graph <- list()

for(i in seq_along(n)){
  
  k <- n[i]
  
  arquivo_nome <- paste0("results_simu/summary_simu_2DMxARMA_nk", k ,".txt")
  
  arquivo <- read.table(file = arquivo_nome, header = TRUE)
  
  rn <- row.names(arquivo)
  
  arquivo <- arquivo |> 
    mutate(
      Parameter = recode(rn,
                         "alpha"  = "alpha",
                         "phi1"   = "phi_{0,0}",
                         "phi2"   = "phi_{1,0}",
                         "phi3"   = "phi_{0,1}",
                         "theta1" = "theta_{0,0}",
                         "theta2" = "theta_{1,0}",
                         "theta3" = "theta_{0,1}")
    ) |> 
    select(Parameter, MAE) |> 
    mutate(n = k)
  
  lista_graph[[i]] <- arquivo
}

lista_graph

graph_final <- bind_rows(lista_graph)

rownames(graph_final) <- NULL

#graph_final$RMSE <- pmin(graph_final$RMSE, 100)



graph_final <- graph_final |>
  mutate(
    par_group = case_when(
      grepl("^alpha", Parameter) ~ "alpha",
      grepl("^phi",   Parameter) ~ "phi",
      grepl("^theta", Parameter) ~ "theta",
      TRUE ~ NA_character_
    )
  )




linetypes <- c("alpha" = "solid", "phi" = "dashed", "theta" = "dotdash")
colors <- c("alpha" = "darkblue", "phi" = "darkred", "theta" = "darkgreen")

graph_final |> 
  ggplot(aes(
    x = n,
    y = MAE,
    group = Parameter,        # legenda detalhada
    color = par_group,        # cor por família
    linetype = par_group      # linha por família
  )) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(x = "Sample size", y = "RB (%)", color = "Parameter", linetype = "Parameter") +
  theme_minimal()



