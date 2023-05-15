#' Tabla ANOVA del DBCA.
#'
#' Obtiene la Tabla de Análisis de Varianza (ANOVA) para un Diseño en Bloques Completamente al Azar (DBCA).
#'
#' @param respuesta (string) nombre de la variable respuesta.
#' @param tratamiento (string) nombre de la variable que representa los tratamientos.
#' @param bloque (string) nombre de la variable que representa los bloques.
#' @param data (\code{data.frame}) Tabla de datos en formato largo con los datos
#'    de los tratamientos, bloques y de la variable respuesta.
#' @return Devuelve una tabla en formato \code{data.frame} con los cálculos correspondientes
#'    al análisis de varianza.
#' @export
#'
#' @examples
#' \dontrun{
#' ## 1 Limpiar la memoria de R
#' rm(list = ls())
#'
#' ## 2 Cargamos el paquete DBCA
#' library(DBCA)
#'
#' ## 3 Ruta del archivo y nombre del mismo
#' archivo <- "C:/Users/Jesus/Documents/datos.csv"
#'
#' ## 4 Lectura de archivo de datos
#' df <- read.csv(archivo)
#'
#' ## 5 Ejecutamos la función TablaAnova
#' TablaAnova(respuesta = "y", tratamiento = "trat", bloque = "bloque", data = df)
#' }
TablaAnova <- function(respuesta, tratamiento, bloque, data){

  # Defino la variable respuesta y los tratamientos y bloques como factores
  y <- data[,respuesta]
  trat <- factor(data[,tratamiento])
  bloque <- factor(data[,bloque])
  a <- nlevels(trat)
  b <- nlevels(bloque)

  # Corrección para la media
  suma_total <- sum(y)
  C <- suma_total^2 /(a*b)

  # SC Total
  sc_total <- sum(y^2) - C
  gl_total <- a*b-1
  cm_total <- sc_total / gl_total

  # SC tratamientos
  sumasxtrat <- tapply(y, INDEX = trat, FUN = sum)
  n_trat <- tapply(y, INDEX = trat, FUN = length)
  sc_trat <- sum(sumasxtrat^2 / n_trat) - C
  gl_trat <- a-1
  cm_trat <- sc_trat / gl_trat

  # SC bloques
  sumasxbloques <- tapply(y, INDEX = bloque, FUN = sum)
  n_bloques <- tapply(y, INDEX = bloque, FUN = length)
  sc_bloques <- sum(sumasxbloques^2 / n_bloques) - C
  gl_bloques <- b-1
  cm_bloques <- sc_bloques / gl_bloques

  #SC residuales
  sc_res <- sc_total - sc_trat - sc_bloques
  gl_res <- (a-1)*(b-1)
  cm_res <- sc_res / gl_res

  # Valores F
  F_trat <- cm_trat / cm_res
  F_bloques <- cm_bloques / cm_res

  # P-values
  p_value_trat <- pf(F_trat, gl_trat, gl_res, lower.tail = FALSE)
  p_value_bloques <- pf(F_bloques, gl_bloques, gl_res, lower.tail = FALSE)


  # Creamos el dataframe
  tabla <- data.frame(FV = c("Tratamientos", "Bloques", "Residuales", "Total"),
                      SC = c(sc_trat, sc_bloques, sc_res, sc_total),
                      GL = c(gl_trat, gl_bloques, gl_res, gl_total),
                      CM = c(cm_trat, cm_bloques, cm_res, NA),
                      F = c(F_trat, F_bloques, NA, NA),
                      `p-value` = c(p_value_trat, p_value_bloques, NA, NA),
                      check.names = FALSE)
  rownames(tabla) <- NULL
  anava <- format(tabla)
  anava[is.na(tabla)] <- ""

  return(anava)
}
