# Protocolo ecuac int: spectral_analyzer.R v1.0
# Analiza las propiedades espectrales de la matriz de corrección m = R - Z.

# Se asume que la librería 'numbers' está instalada para los caracteres de Dirichlet.
# install.packages("numbers")
library(numbers)

spectral_analyzer <- function(N) {
  # 1. Generar las matrices R y Z.
  rz_data <- rrz_generator(N)
  if (any(is.na(rz_data$R))) {
    return(list(error = "N no es primo, no se puede analizar el espectro."))
  }
  
  # 2. Construir la matriz de corrección m.
  m <- rz_data$R - rz_data$Z
  
  # 3. Calcular propiedades espectrales clave.
  # El rango se calcula contando los valores singulares mayores a una tolerancia.
  rank_m <- sum(svd(m)$d > 1e-9)
  
  # Traza de m y de m^2.
  trace_m <- sum(diag(m))
  trace_m_sq <- sum(diag(m %*% m))
  
  # 4. Verificación de la Dicótoma Espectral (Teorema #171, versión m).
  # Nota: El teorema original era para m', aquí verificamos las trazas de m.
  
  resultados <- list(
    N = N,
    N_mod_4 = N %% 4,
    rank_m_observado = rank_m,
    rank_m_teorico = N - 2,
    traza_m = trace_m,
    traza_m2 = trace_m_sq
  )
  
  return(resultados)
}

# Ejemplo de uso:
# print(spectral_analyzer(7))
# print(spectral_analyzer(13))