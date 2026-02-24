# Bioinform-tica-avanzada
#NIVEL 1 – Manipulación básica con dplyr

  #Para seleccionar el directorio de trabajo
  
setwd("C:/Users/Alumno 25-B/Desktop")
  
  #Para seleccionar la base de datos .csv
  
base=read.csv("genes2.csv")

  #Para leer las primeras 6 filas
  
head(base)

  #Para conocer la estructura de la base de datos
  
str(base)

  #Para cambiar una columna de la base de datos a factores (Columna Tejido)
  
base$Tejido=as.factor(base$Tejido) 

  #Para ordenar factores (bajo<medio<alto)
  
base$Estres=factor(base$Estres,levels=c("Bajo","Medio","Alto"),ordered=TRUE)


#Reto 1: ¿Cuál es el promedio de expresión por tratamiento?
  #Para agrupar datosstr
  
lista_tratamientos <- base %>% + group_by(Tratamiento) %>% + group_split()

  #Para  sacar promedio de cada dato agrupado
  
base %>% group_by(Tratamiento) %>% summarise(promedio = mean(Expresion_TPM, na.rm = TRUE))

<img width="184" height="119" alt="image" src="https://github.com/user-attachments/assets/071bc878-b53a-4923-8452-ec48d057b8d2" />


Reto 2: Filtrar solo: - Tejido = Higado - Estrés = Alto - Tiempo = 48h
Ordenar por mayor expresión.

  #Para filtrar los datos a Tejido = Higado - Estrés = Alto - Tiempo = 48h
  
HAT= filter(base, Tejido=="Higado", Estres=="Alto", Tiempo=="48h")







