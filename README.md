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

<img width="596" height="31" alt="image" src="https://github.com/user-attachments/assets/dfc79472-fed7-4d6c-86d0-3504a7897221" />

#Reto 6: Calcular fold change promedio respecto a Control.

media_con=8.17

media_CRIS=8.29

media_DrugA=8.31

media_DrugB=8.45

RNAi=8.35

fc_CRISPR=media_CRIS/media_con

fc_DrugA=media_DrugA/media_con

fc_DrugB=media_DrugB/media_con

fc_RNAi=RNAi/media_con

<img width="301" height="178" alt="image" src="https://github.com/user-attachments/assets/cd64ff50-642e-47a9-805f-772065e5a0e7" />

  #Reto 7: Crear tabla ancha donde: Filas = Gene Columnas = Tratamiento    Valores = promedio de expresión.

ancha= nueva %>% + pivot_wider(names_from = Gene, values_from = Tratamiento)

<img width="1280" height="192" alt="image" src="https://github.com/user-attachments/assets/f0f9deaa-b6e6-44de-a25c-73e1e4ba56ca" />

  #Reto 8: Regresar la tabla anterior a formato largo.

tabla_larga=tabla_ancha %>% pivot_longer(cols=-x,names_to="Gene", values_to="Tratamiento")
tabla_limpia <- na.omit(tabla_larga)
tabla_limpia1<- tabla_larga %>%
+     rename(
+         Gene = name,
+         Tratamiento = value)

<img width="203" height="790" alt="image" src="https://github.com/user-attachments/assets/eaf08935-68ec-4037-8eed-6e30e4614de1" />

  #Reto 9: Separar la variable Sample_ID en Tejido, Tiempo y Replica.

base$Gene=gsub("Gene_","",base$Gene)

<img width="516" height="236" alt="image" src="https://github.com/user-attachments/assets/43db909d-8f38-49de-b0dc-94aecd25381c" />

#Reto 10: Identificar genes con: - Expresión media > 15 TPM - Bajo estrés
alto - En al menos 2 tejidos

genes_filtrados <- base %>%
+     filter(Estres == "Bajo") %>%                
+     group_by(Gene, Tejido) %>%
+     summarise(media_TPM = mean(Expresion_TPM), .groups = "drop") %>%
+     filter(media_TPM > 5) %>%                   
+     group_by(Gene) %>%
+     summarise(n_tejidos = n_distinct(Tejido)) %>%
+     filter(n_tejidos >= 2)


#Reto 11: ¿Cuál combinación Tratamiento–Tejido produce mayor varianza?
  varianza_trat_teji <- base %>%
+     group_by(Tratamiento, Tejido) %>%
+     summarise(varianza = var(Expresion_TPM, na.rm = TRUE), .groups = "drop") %>%
+     arrange(desc(varianza))

  varianza_trat_teji %>%
+     slice_max(varianza, n = 1)


#Reto 12: Generar tabla con: - Media por Gene - Media global - Diferencia
respecto a la media global

  media_global <- mean(base$Expresion_TPM, na.rm = TRUE)
  
  tabla_genes <- base %>%
+     group_by(Gene) %>%
+     summarise(media_gene = mean(Expresion_TPM, na.rm = TRUE)) %>%
+     mutate(
+         media_global = media_global,
+         diferencia = media_gene - media_global)


#Reto 13: Graficar la distribución de expresión (histograma) general.
  grafica <- ggplot(base, aes(x = Expresion_TPM)) +
+     geom_histogram(bins = 30, na.rm = TRUE)


#Reto 14: Hacer boxplot de expresión por Tratamiento.
  ggplot(base, aes(x = Tratamiento, y = Expresion_TPM)) +
+     geom_boxplot(na.rm = TRUE) +
+     labs(
+         title = "Expresión génica por tratamiento",
+         x = "Tratamiento",
+         y = "Expresión (TPM)")


#Reto 15: Hacer boxplot de expresión por Estrés, coloreado por
Tratamiento.
  ggplot(base, aes(x = Estres, y = Expresion_TPM, fill = Tratamiento)) +
+     geom_boxplot(na.rm = TRUE) +
+     labs(
+         title = "Expresión génica por nivel de estrés",
+         x = "Estrés",
+         y = "Expresión (TPM)",
+         fill = "Tratamiento")


#Reto 16: Graficar media de expresión por Tratamiento (usar datos
resumidos).
  datos_resumen <- base %>%
+     group_by(Tratamiento) %>%
+     summarise(media_expresion = mean(Expresion_TPM, na.rm = TRUE))

  ggplot(datos_resumen, aes(x = Tratamiento, y = media_expresion)) +
+     geom_col() +
+     labs(
+         title = "Media de expresión génica por tratamiento",
+         x = "Tratamiento",
+         y = "Media de expresión (TPM)")



