# -- bibliotecas
library(ggplot2)
library(dplyr)
library(purrr)
library(magrittr)
library(rootSolve)
set.seed(1242)

# -- Especificação do modelo AR1
# y(t) = rho * y(t-1) + eps, sendo y(1) = 0, eps ~ N(0, 1), 0 < rho < 1

# -- Simula série temporal
rho <- 0.7
N <- 100
y <- numeric(N)
for (i in 2:N){
    y[i] <- rho * y[i-1] + rnorm(1)
}

df_simulado <- data.frame(periodos=1:N, valores = y) 

# -- gráfico da série
ggplot(df_simulado, aes(x = periodos, y = valores)) +
    geom_line(size=1) +
    xlab("") +
    ylab("Serie 1")

# -- Função de verossimilhança aproximada (ignorando primeira observação)
# -- yt|yt-1 ~ N(rho*y[t-1], 1) para t > 1
gera_vero_aprox <- function(param, dados){
    N <- length(dados)
    sum(dnorm(dados[2:N], mean=param*dados[1:(N-1)], sd = 1, log = T))
}
gera_vero_aprox.V <- Vectorize(gera_vero_aprox, "param")

# -- Função deviance da verossimilhança aproximada
# -- D = 2*[l(est_param) - l(param)]
deviance_vero <- function(param, log_vero_fn, est_param, ...){
    2*(log_vero_fn(est_param, ...) - log_vero_fn(param, ...))
}

# -- obtendo rho
rho_emv <- optimize(gera_vero_aprox, interval = c(0, 1), dados = y, maximum = T)$max
rho_emv

# -- Gráfico da função de verossimilhança
valores_param <- seq(0, 1, 0.01)
valores_ll <- map_dbl(valores_param, gera_vero_aprox, y)
df_vero_appx <- data.frame(params = valores_param, ll = valores_ll, emv = rho_emv)

df_vero_appx %>% ggplot(aes(x = params, y = ll)) +
    geom_line() +
    geom_vline(xintercept = rho_emv) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    xlab(expression(rho)) +
    ylab(expression(l(rho)))

# -- Encontrando os intervalos de confiança usando a função deviance
# -- rho: r = 0.10, cD = -2 * log(r), deviance - cD = 0
# -- Explicação:
# -- O intervalo deviance vem da relação entre as verossimilhanças. Desse modo:
# -- Intervalo baseado em verossimilhança: r >= L(theta)/L(theta^hat)
# -- => -2log(r) = -2[l(theta) - l(theta^hat)]
# -- => -2log(r) = 2[l(theta^hat - l(theta))]
# -- => cD = 2[l(theta^hat - l(theta))] ~ X_d (distribuição assintótica)
ic_dev <- function(param, deviance_fn, cD, ...){
    deviance_fn(param, ...) - cD            
}

cD <- -2 * log(0.1)
vl_ic_dev_aprox <- uniroot.all(ic_dev, c(0, 1), 
                               deviance_fn = deviance_vero, 
                               cD = cD,  log_vero_fn = gera_vero_aprox.V,
                               est_param = rho_emv, 
                               dados = y)
valores_dev <- map_dbl(valores_param, 
                       deviance_vero, 
                       log_vero_fn = gera_vero_aprox,
                       est_param = rho_emv, 
                       dados = y)
df_vero_appx$dev <- valores_dev

df_vero_appx %>% ggplot(aes(x = params, y = dev)) +
    geom_line() +
    geom_vline(xintercept = vl_ic_dev_aprox[1]) +
    geom_vline(xintercept = vl_ic_dev_aprox[2]) +
    scale_x_continuous(limits = c(0.5, 0.9)) +
    scale_y_continuous(limits = c(0, 7)) +
    xlab(expression(rho)) +
    ylab(expression(D(rho)))

# -- Função de log-verossmilhança completa para o modelo AR1
# -- Y1 ~ N(0, 1 / (1 - rho ^ 2)) e Yt|Yt-1 ~ N(rho*y[t-1], 1) para t > 1
gera_vero_comp <- function(param, dados){
    N <- length(dados)
    dnorm(dados[1], mean = 0, sd = sqrt((1 / (1 - param^2))), log = T) + 
        sum(dnorm(dados[2:N], mean=param*dados[1:(N-1)], sd = 1, log = T))
}

gera_vero_comp.V <- Vectorize(gera_vero_comp, "param")
rho_emv_comp <- optimize(gera_vero_comp, 
                         c(0, 1),
                         dados = y, 
                         maximum = T)

# -- I.C com deviance para verossimilhança completa
vl_ic_dev_comp <- uniroot.all(ic_dev, c(0, 1), 
                               deviance_fn = deviance_vero, 
                               cD = cD,  log_vero_fn = gera_vero_comp.V,
                               est_param = rho_emv_comp$maximum, 
                               dados = y)
valores_dev <- map_dbl(valores_param, 
                       deviance_vero, 
                       log_vero_fn = gera_vero_aprox,
                       est_param = rho_emv, 
                       dados = y)
df_vero_appx$dev <- valores_dev

