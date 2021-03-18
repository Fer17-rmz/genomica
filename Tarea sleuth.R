#------------------------------------------------------------------------------------
# 17 / 03 / 21
##BiocManager :: install ("biomaRt")
library("sleuth") #libreria necesaria para actividad 
library(biomaRt)

#Hacer una funcion facilita el llevar a cabo la actividad de BiomaRt, de manera que una vez
#establecidos los argumentos solo se corre la funcion despues 
tx2gene <- function(){
  
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") 
  #biomRT= de doonde se mapea 
  #dataset= de donde tomanos los datos
  
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = "ensembl_transcript_id",
                       ens_gene = "ensembl_gene_id", ext_gene = "external_gene_name")
  return(t2g)
}

t2g <- tx2gene() #Conectarse a la base de datos de ensambled
t2g #cuando lo corres te sale una lista enorme 

#Directorio de la base de datos (cambia si lo haces de otra pcy asi, solo buscas )
base_dir <-"~/6 SEXTO SEMESTRE/Genómica/__MACOSX" 

#SAMPLES
#en un mundo ideal buscas aquellos que quieres analizar en este caso tomare los marcados, se pueden cambiar
#depende de lo que quieras hacer 
samples <- paste0("sample", c("2","3","4",
                                   "10","11","12")) 
samples #Nos aseguramos de escoger los que queriamos, para eso imprimimos el objeto 

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id)) 
kal_dirs
#En este objeto lo que buscamos es asignar la ID de donde esta el sample
#objeto<- sapply(archivos(sample), function(id)) de donde se toma

#CREACIÓN DE DATAFRAME
#Con las  direcciones
#samples escogidas previamente 
#y el nombre que se les van a colocar a las muestras (En un experimento o trabajo el nombre 
#dependerá del origen de la muestra )
#por ejemplo si tenemos Wild types y mutantes 
s2c <- data.frame(path=kal_dirs, sample=samples, muestras <- c("WT","WT","WT", "M",
                                                              "M", "M"), stringsAsFactors= FALSE)
s2c
#
so <- sleuth_prep (s2c, full_model=~muestras, target_mapping = t2g, extra_bootstrap_summary = TRUE) 
#Asignamos los datos a un objeto para manipularlos
#s2c es el dataframe que hicimos con las muestras 
#full_mode, target_mapping y extra_bootstrap forman parte de la formula 
#t2g funcion anterior 

#PRUEBA DE HIPOTESIS
so <- sleuth_fit(so) #ajustar el modelo
so <- sleuth_wt(so, which_beta = "wildtype-Mutante")   # realizar prueba de hipotesis            
sleuth_live(so) #visualizacion con shiny (hace más digerible la visualización)
#cuando abre la cargas 

#CARGAR TABLA DE RESULTADOS 
setwd("~/6 SEXTO SEMESTRE/Genómica/__MACOSX/") #dirección
resultados<-read.table("test_table.csv",sep=",",
                       header=TRUE) #se cargan a R

#Busqueda y selección de los resultados que son significativos 
significativos<-which(resultados$qval<0.1) #encontrar aquelos valores de resultados que cumplen la condicion
significativos<-resultados[significativos,] #seleccion 

#INCREMEMTO
upregulated<-which(significativos$b>0) #buscar
upregulated<-significativos[upregulated,] #seleccionar 

#DECRECIERON 
downregulated<-which(significativos$b<0) #busqueda
downregulated<-significativos[downregulated,] #seleccion

#GUARDAR RESULTADOS
#para ver bonito los colocamos en una tabla (asignamos los resultados finales) a cada tipo de 
#regulación, objeto a tomar, archivo (dirección donde se guarda)
write.table(upregulated,file="~/6 SEXTO SEMESTRE/Genómica/__MACOSX/Upreg.txt",sep="\t")
write.table(downregulated,file="~/6 SEXTO SEMESTRE/Genómica/__MACOSX/Downreg.txt",sep="\t")

#REGULACIÓN POSITIVA
positivos_como_io<-upregulated$ext_gene