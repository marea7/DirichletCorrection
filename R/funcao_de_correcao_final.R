#' Correction function in estimation in multivariate Dirichlet regressions.
#'
#' Function to calculate the improved estimate in multivariate Dirichlet regressions, based on Melo et al. (2020).
#'
#' @param X It is the n*k model matrix with the observations in k known covariates
#' @param Y It is a random p-vector that has a Dirichlet distribution(λ).
#' @return A list containing:
#' \item{Betas}{The values of the beta matrix without correction.}
#' \item{Erro Padrão}{The values of the standard error.}
#' \item{B(θ)}{The second-order bias vector of theta in matrix format.}
#' \item{Betas Corrigidos}{Matrix with estimates of the bias correction of θ.}
#' \item{Função Score Modificada}{The modified Score function, which is derived from the logarithm of the likelihood function.}
#'
#' @export

dirichletBiasCorrection <- function(X, Y) {
  set.seed(1992)
  library(pracma)

  ## Número de colunas da matriz X (matriz de regressores)
  k  = ncol(X)
  ## Número de colunas e linhas da matriz Y (núm. de componentes)
  p  = ncol(Y)
  n = nrow(Y)

  ## Dimensão de theta
  kp = k*p

  ## chutes iniciais

  Y_media <- colMeans(Y)  # Média de cada componente
  Y_var <- apply(Y, 2, var)  # Variância de cada componente

  # Escolher um valor inicial para phi (dispersão)
  phi_inicial <- mean(Y_media * (1 - Y_media) / Y_var - 1)

  # Calcular alpha inicial usando Y_mean
  alpha_inicial <- Y_media * phi_inicial  # Chute inicial para alpha
  alpha_inicial

  alpha_inicial <- matrix(rep(Y_media * phi_inicial, each = n), nrow = n, ncol = p, byrow = FALSE)

  ## chute inicial de beta
  beta_inicial <- matrix(0, nrow = k, ncol = p)  # Inicializar beta
  for (j in 1:p) {
    log_alpha_j <- log(alpha_inicial[,j]) # Log de alpha para o j-ésimo componente
    beta_inicial[, j] <- solve(t(X) %*% X) %*% t(X) %*% log_alpha_j
  }

  beta_vetor <- as.vector(beta_inicial)  # Transformar em vetor 1D para o optim

  ## função log-vero do modelo
  log_likelihood_dirichlet <- function(beta_vetor, Y, X) {
    beta_inicial <- matrix(beta_vetor, nrow = k, ncol = p)
    eta <- X %*% beta_inicial
    alfas_teste <- exp(eta)
    n <- nrow(Y)
    log_likelihood <- 0

    for (i in 1:n) {
      phi <- sum(alfas_teste[i, ])
      log_likelihood <- log_likelihood +
        (lgamma(phi) - sum(lgamma(alfas_teste[i, ]))) +
        sum((alfas_teste[i, ] - 1) * log(Y[i, ]))
    }
    return(-log_likelihood)
  }

  ## Calcular alfas e betas estimados:

  #Valor da função objetivo no ponto inicial
  f_x0 <- log_likelihood_dirichlet(beta_vetor, Y, X)

  #Precisão da máquina
  epsilon <- .Machine$double.eps  # Aproximadamente 1.1e-16

  #Calcular factr
  tolerancia_objetivo <- 1e-8
  factr <- tolerancia_objetivo / (epsilon * abs(f_x0))

  #Definir pgtol
  pgtol <- 1e-5

  #Definir maxit
  maxit <- 1000

  #Chamar optim com os parâmetros ajustados
  resultado <- optim(
    par = beta_vetor,
    fn = log_likelihood_dirichlet,
    Y = Y,
    X = X,
    method = "L-BFGS-B",
    control = list(factr = factr, pgtol = pgtol, maxit = maxit)
  )
  betas <- resultado$par
  betas <- matrix(betas, nrow = k, ncol = p)

  # Considerando que:
  galfa = X %*% betas
  alfas_estim = exp(galfa)
  # então:
  g_alfa <- galfa                 # ou, equivalentemente, log(alfas_estim)
  g_alfa_derivado <- exp(g_alfa)    # que é igual a alfas_estim
  g_alfa_segunda  <- exp(g_alfa)     # igual também a alfas_estim

  C <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
  for (i in 1:nrow(Y)) {
    for (j in 1:ncol(Y)) {
      C[i, j] <- g_alfa_derivado[i, j] * (psi(1, sum(alfas_estim[i,])) - psi(1, alfas_estim[i, j]) + log(Y[i, j]))
    }
  }
  C

  I_p <- diag(p)
  Kron_prod <- kronecker(I_p, X)

  k_teta <- matrix(0, nrow = k*p, ncol = k*p)
  w <- array(0, dim = c(p, p, n))
  W <- matrix(0, n*p, n*p)
  # vetor_ab <- w[ a, b, ]  # Retorna um vetor de tamanho ta

  for(u in 1:k) {
    for(v in 1:k) {
      for (a in 1:p) {
        for (b in 1:p) {
          r = (b-1)*k + v
          s = (a-1)*k + u
          for (i in 1:n) {
            if (a == b) {
              w[a, b,i] <- - (g_alfa_derivado[i, a]^2) *
                (psi(1, sum(alfas_estim[i,])) - psi(1, alfas_estim[i, a]))
            } else {
              w[a, b,i] <- -g_alfa_derivado[i, a] * g_alfa_derivado[i, b] *
                psi(1, sum(alfas_estim[i,]))
            }
          }
          k_teta[r,s] <- (X[,u] * w[a,b,]) %*% X[,v]
        }
      }
    }
  }

  for (a in 1:p) {
    for (b in 1:p) {
      # Criando a matriz diagonal para este bloco
      bloco <- diag(w[a, b, ])

      # Definindo os índices corretos dentro da matriz M
      linha_ini <- (a - 1) * n + 1
      linha_fim <- a * n
      col_ini <- (b - 1) * n + 1
      col_fim <- b * n

      # Inserindo a matriz no bloco correto de M
      W[linha_ini:linha_fim, col_ini:col_fim] <- bloco
    }
  }
  W

  k_theta <- t(Kron_prod) %*% W %*% Kron_prod ##igual a k_teta

  f <- array(0, dim = c(p, p, p, n))
  delta <- array(0, dim = c(p, p, p, n))
  mrts <- array(0, dim = c(k * p, k * p, k * p))

  for(u in 1:k) {
    for(v in 1:k) {
      for(z in 1:k) {
        for (a in 1:p) {
          for (b in 1:p) {
            for (c in 1:p) {
              t = (a-1)*k + u
              s = (b-1)*k + v
              r = (c-1)*k + z

              for (i in 1:n) {
                if (b == c && c == a) {
                  f[c, b, a, i] <- 3 * g_alfa_segunda[i, a] * g_alfa_derivado[i, a] *
                    (psi(1, sum(alfas_estim[i,])) - psi(1, alfas_estim[i, a])) +
                    g_alfa_derivado[i, a]^3 * (psi(2, sum(alfas_estim[i,])) - psi(2, alfas_estim[i, a]))
                }

                if (c == a && a != b) {
                  f[c, b, a, i] <- g_alfa_derivado[i, b] *
                    (g_alfa_segunda[i, a] * psi(1, sum(alfas_estim[i,])) + g_alfa_derivado[i, a]^2 * psi(2, sum(alfas_estim[i,])))
                }

                if (b == c && b != a) {
                  f[c, b, a, i] <- g_alfa_derivado[i, a] *
                    (g_alfa_segunda[i, b] * psi(1, sum(alfas_estim[i,])) + g_alfa_derivado[i, b]^2 * psi(2, sum(alfas_estim[i,])))
                }

                if (b == a && a != c) {
                  f[c, b, a, i] <- g_alfa_derivado[i, c] *
                    (g_alfa_segunda[i, a] * psi(1, sum(alfas_estim[i,])) + g_alfa_derivado[i, a]^2 * psi(2, sum(alfas_estim[i,])))
                }

                if (b != c && c != a && b != a) {
                  f[c, b, a, i] <- g_alfa_derivado[i, a] * g_alfa_derivado[i, b] *
                    g_alfa_derivado[i, c] * psi(2, sum(alfas_estim[i,]))
                }

                # deltas
                if (b == a && a == c) {
                  delta[b, a, c, i] <- 2 * g_alfa_segunda[i, c] * g_alfa_derivado[i, c] *
                    (psi(1, sum(alfas_estim[i,])) - psi(1, alfas_estim[i, c])) +
                    g_alfa_derivado[i, c]^3 * (psi(2, sum(alfas_estim[i,])) - psi(2, alfas_estim[i, c]))
                }

                if (b == c && c != a) {
                  delta[b, a, c, i] <- g_alfa_derivado[i, a] * g_alfa_derivado[i, c]^2 * psi(2, sum(alfas_estim[i,]))
                }

                if (b == a && a != c) {
                  delta[b, a, c, i] <- g_alfa_derivado[i, c] * g_alfa_derivado[i, a]^2 * psi(2, sum(alfas_estim[i,])) +
                    g_alfa_segunda[i, a] * g_alfa_derivado[i, c] * psi(1, sum(alfas_estim[i,]))
                }

                if (a == c && c != b) {
                  delta[b, a, c, i] <- g_alfa_derivado[i, b] * g_alfa_derivado[i, c]^2 * psi(2, sum(alfas_estim[i,])) +
                    g_alfa_segunda[i, c] * g_alfa_derivado[i, b] * psi(1, sum(alfas_estim[i,]))
                }

                if (a != b && b != c && a != c) {
                  delta[b, a, c, i] <- g_alfa_derivado[i, a] * g_alfa_derivado[i, b] *
                    g_alfa_derivado[i, c] * psi(2, sum(alfas_estim[i,]))
                }
              }
              mrts[r,s,t] = t(X[,v] * X[,u] * delta[b,a,c,]) %*% X[,z] - (0.5*t(X[,z] * X[,v] * f[c,b,a,]) %*% X[,u])
            }
          }
        }
      }
    }
  }
  f
  delta
  mrts

  dim(mrts)[1] <- k*p
  dim(mrts)[2] <- k*p
  dim(mrts)[3] <- k*p

  M_aux <- mrts[, , 1]
  if (kp > 1) {
    for (r in 2:kp) {
      M_aux <- cbind(M_aux, mrts[, , r])
    }
  }

  #Agora, W_aux é a concatenação de todas as fatias de W
  mrts <- M_aux

  inverta <- solve(k_teta) #K_teta e INVERSA ESTÃO CERTAS
  erro_padrao <- sqrt(diag(inverta))

  b_teta <- inverta %*% mrts %*% as.vector(inverta)

  teta_chapeu <- c(betas[,1], betas[,2], betas[,3], betas[,4])

  teta_bc = matrix(teta_chapeu - b_teta,nrow = k, ncol = p)

  U_theta <- as.vector(t(X) %*% C)

  #função score Modificada
  U_theta_f1 <- U_theta - (k_teta %*% b_teta)
  #ou
  U_theta_f <- U_theta - mrts %*% as.vector(inverta)
  #São iguais
  return(list(
    "Betas" = betas,
    "B(θ)" = matrix(b_teta,,nrow = k, ncol = p),
    "Betas Corrigidos" = teta_bc,
    "Erro Padrao" = matrix(erro_padrao,,nrow = k, ncol = p),
    "Funcao Score Modificada" = U_theta_f))
    }
