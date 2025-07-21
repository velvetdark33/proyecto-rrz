# Protocolo ecuac int: rrz_generator.R v1.0
# Genera las matrices R y Z para un módulo N.
# La función devuelve una lista con R y Z. Si N no es primo, las matrices
# contendrán valores NA, lo que constituye un test de primalidad.

rrz_generator <- function(N) {
  # Validar que N sea un entero mayor que 2.
  if (N <= 2 || N != round(N)) {
    stop("N debe ser un entero mayor que 2.")
  }
  
  n_minus_1 <- N - 1
  R_matrix <- matrix(NA, nrow = n_minus_1, ncol = n_minus_1)
  Z_matrix <- matrix(NA, nrow = n_minus_1, ncol = n_minus_1)
  
  # Índices de 1 a N-1, no de 0 a N-2, para alinear con la fórmula matemática.
  y_coords <- 1:n_minus_1
  x_coords <- 1:n_minus_1
  
  # Bucle optimizado para seguir la lógica original del descubrimiento.
  for (y in y_coords) {
    for (x in x_coords) {
      # El bucle de búsqueda para z, la variable diofántica.
      for (z in 1:(N - y)) {
        r0 <- (N * z - x) / (N - y)
        # Condición para que r0 sea un entero dentro del rango válido.
        # Se usa una tolerancia para la comparación de punto flotante.
        if (abs(r0 - round(r0)) < 1e-9 && r0 >= 1 && r0 <= n_minus_1) {
          # Usamos x,y como índices directos de la matriz.
          R_matrix[x, y] <- round(r0)
          Z_matrix[x, y] <- z
          break # Rompe el bucle de z una vez encontrada la solución.
        }
      }
    }
  }
  
  return(list(R = R_matrix, Z = Z_matrix))
}

# Ejemplo de uso:
# N <- 7
# rz_data <- rrz_generator(N)
# print("Matriz R:")
# print(rz_data$R)
# print("Matriz Z:")
# print(rz_data$Z)