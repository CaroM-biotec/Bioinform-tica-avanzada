# Bioinform-tica-avanzada

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


#Reto 2: Filtrar solo: - Tejido = Higado - Estrés = Alto - Tiempo = 48h Ordenar por mayor expresión.

  #Para filtrar los datos a Tejido = Higado - Estrés = Alto - Tiempo = 48h
  
HAT= filter(base, Tejido=="Higado", Estres=="Alto", Tiempo=="48h")

<img width="581" height="638" alt="image" src="https://github.com/user-attachments/assets/46c57aae-5a11-480a-ae92-959b7faf798b" />

#Reto 3: ¿Cuáles son los 10 genes más expresados en condiciones de CRISPR?

  #Para ordenar de mayor a menor
  
HAT %>% arrange(-Expresion_TPM)

<img width="595" height="181" alt="image" src="https://github.com/user-attachments/assets/13d8cac2-697e-44dd-8450-7d00a64bc71a" />

#Reto 4: Calcular media y desviación estándar por Gene y Tratamiento.

  #Para calcular desviación estándar
  
base %>% group_by(Tratamiento) %>% summarise(sd(Expresion_TPM, na.rm = TRUE))

<img width="387" height="117" alt="image" src="https://github.com/user-attachments/assets/c6a6fe50-d08c-4abc-b48e-286ccc9e1f58" />
<img width="184" height="119" alt="image" src="https://github.com/user-attachments/assets/071bc878-b53a-4923-8452-ec48d057b8d2" />

#Reto 5: ¿Cuál tratamiento aumenta más la expresión bajo estrés alto?



