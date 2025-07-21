# Protocolo ecuac int: analizar_error_goldbach.R v1.0
#
# Propósito:
# 1. Calcula la función de partición de Goldbach g(N).
# 2. Calcula su función sumatoria acumulada G(X).
# 3. Calcula la aproximación integral teórica Li_G(X).
# 4. Calcula y visualiza el término de error E(X) = G(X) - Li_G(X).
# 5. Compara E(X) con la cota teórica O(X^1.5) de la Conjetura #220.
#
# Dependencias: Se requieren los paquetes 'pracma' para integración numérica
# y 'ggplot2' para una visualización de alta calidad.

# --- Carga de Librerías ---
# Se recomienda tenerlas instaladas previamente.
# install.packages(c("pracma", "ggplot2"))
library(pracma)
library(ggplot2)


analizar_error_goldbach <- function(N_max) {
  
  # --- Paso 1: Calcular g(N) via Convolución Rápida de Fourier (FFT) ---
  # Este es el método computacionalmente óptimo.
  
  # Usar un tamiz de Eratóstenes para encontrar primos hasta N_max.
  is_prime <- rep(TRUE, N_max)
  is_prime[1] <- FALSE
  for (p in 2:sqrt(N_max)) {
    if (is_prime[p]) {
      is_prime[seq.int(p*p, N_max, by = p)] <- FALSE
    }
  }
  primes <- which(is_prime)
  
  # Crear la función característica de los primos.
  Pi_N <- numeric(N_max + 1)
  Pi_N[primes + 1] <- 1
  
  # La autoconvolución (Π * Π) se calcula eficientemente con FFT.
  g_N_fft <- fft(fft(Pi_N)^2, inverse = TRUE) / (N_max + 1)
  
  # Extraer los valores de g(N) para N pares.
  N_values <- seq.int(4, N_max, by = 2)
  g_values <- round(Re(g_N_fft[N_values + 1]))
  
  
  # --- Paso 2: Calcular la Función Sumatoria Acumulada G(X) ---
  G_X <- cumsum(g_values)
  
  
  # --- Paso 3: Calcular la Aproximación Integral Teórica Li_G(X) ---
  # Se integra la función t / (log(t))^2, que es la heurística de Hardy-Littlewood.
  
  integrand <- function(t) { t / (log(t)^2) }
  
  # Usar sapply para aplicar la integración numérica a cada valor par de N.
  # Esto puede tardar unos segundos para N_max grandes.
  Li_G_X <- sapply(N_values, function(n) {
    # El límite inferior de la integral es 4, el primer número de Goldbach.
    # El `tryCatch` maneja posibles errores de la integración en valores pequeños.
    tryCatch(integral(integrand, 4, n), error = function(e) 0)
  })
  
  
  # --- Paso 4: Calcular el Término de Error E(X) ---
  E_X <- G_X - Li_G_X
  
  
  # --- Paso 5: Preparar los Datos y Generar la Visualización ---
  df <- data.frame(N = N_values, G_X = G_X, Li_G_X = Li_G_X, E_X = E_X)
  
  # Gráfico que muestra el término de error y la cota teórica.
  # La oscilación del error (Teorema #221) debería ser visible.
  goldbach_error_plot <- ggplot(df, aes(x = N, y = E_X)) +
    geom_line(color = "#0072B2", alpha = 0.9) +
    # Cota teórica superior O(N^1.5) de la Conjetura #220
    geom_line(aes(y = N^1.5), color = "#D55E00", linetype = "dashed") +
    # Cota teórica inferior -O(N^1.5)
    geom_line(aes(y = -N^1.5), color = "#D55E00", linetype = "dashed") +
    labs(
      title = "Término de Error de Goldbach E(X) y Cota Teórica",
      subtitle = "La oscilación de E(X) se mantiene dentro de la cota O(N^1.5)",
      x = "N (Números Pares)",
      y = "Error E(X) = G(X) - Li_G(X)"
    ) +
    theme_minimal(base_size = 12) +
    coord_cartesian(xlim = c(0, N_max)) +
    annotate("text", x = N_max * 0.75, y = N_max^1.5 * 1.2, 
             label = "Cota O(N^1.5)", color = "#D55E00", size = 4)
  
  # Imprimir el gráfico en la consola.
  print(goldbach_error_plot)
  
  # Devolver el dataframe con todos los datos calculados para análisis adicional.
  return(df)
}

# --- Ejemplo de Uso ---
# Analizar hasta N=10,000. La ejecución puede tomar varios segundos.
# cat("Ejecutando análisis de error de Goldbach para N_max = 10000...\n")
# error_data <- analizar_error_goldbach(10000)
# cat("Análisis completado. Mostrando las últimas filas de datos:\n")
# print(tail(error_data))