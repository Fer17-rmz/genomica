#I.Análisis en  red de datos de expresión hepática de ratones hembra:
#   búsqueda de módulos  relacionados con el peso corporal

#==================== Entrada de datos y limpieza de datos =======================
# PRIMERO VER DIRECTORIO ACTUAL

#Es importante puesto que facilita el trabajo posterior, no andas perdido
getwd();
# al objeto le asignamos la dirección en donde se encuentran los datos que vamos
# a manejar ( es importante escoger la carpeta data y no el que se muestra 
# como archivo zip)
workingDir <- "~/R/GENOMICA/FemaleLiver-Data.zip";
#Esta funcion no la corremos porque ya estamos en la dirección que queremos
# de lo contrario nos permite movernos de la inicial ( que arrojo wget() a la
#indicada )
# setwd(workingDir); 

# Cargamos la libreria que nos permite manejar redes pesadas 
library(WGCNA);
# PASO IMPORTANTE 
options(stringsAsFactors = FALSE);
#Leemos uno de los archivos separados por comas de la carpeta que descargamos
femData = read.csv("~/R/GENOMICA/FemaleLiver-Data/LiverFemale3600.csv");
# Vemos un poco de la base para saber con que estamos trabajando, ya que no son 
#datos que sacamos nosotros es importante tener una aproximación

#Ver la dimensión de la base 
dim(femData);
#Ver los nombres 
names(femData);

#==== Segunda parte======
#Datos de expresión en una nueva data. frame
datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
names(datExpr0) = femData$substanceBXH; #Nombrar las columnas de la nueva data
# acorde a los nombres de la original 
rownames(datExpr0) = names(femData)[-c(1:8)]; #nombrar columnas

#======= Parte tres =======
# Verificar datos
# para un correcto analisis busca elementos que no cumplen con las caracteristicas
#apropiadas, por ejemplo estar debajo del umbral establecido 
# la salida da como resultado una lista de los genes que cumplen los criterios
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK # si es necesario se repite

#============= Parte cuatro  ============
if (!gsg$allOK)
{
        #Imprime la lista de genes anterior
        # y para um mayor control en el experimento dice cuales fueron removidos
        if (sum(!gsg$goodGenes)>0) 
                printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0) 
                printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
        
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# =========== Quinto paso ========
#Trazar el Árbol de muestra 
sampleTree = hclust(dist(datExpr0), method = "average");
#Tamaño de ventana, depende del experimento y la capacidad de computo
sizeGrWindow(12,9)
#parametros 
par(cex = 0.6);
par(mar = c(0,4,2,0))
#trazar como tal el arbol 
plot(sampleTree, main = "Agrupación de muestras para detectar valores atípicos", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#======== Parte seis =====
# El grafico anterior por si solo no explicaria mucho a personas que no llevan a cabo 
# el analisis, poner la linea de corte  
abline(h = 15, col = "darkorchid3");
# Ahora que se tiene el punto de corte  podemos determinar los grupos que estan 
# por debajo de esta 
# la funcion permite hacer cortes en el dendograma a una altura constante (h=15)
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#Ver el resultado del corte en forma de tabla 
table(clust)
# Ahora que tenemos los cluster, nos quedamos con el cluster 1 de 134 elementos
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ] # los datos de expresion quw usamos ahora son los 
# que estan en el cluster 1 
nGenes <- ncol(datExpr) # nombramos columnas
nSamples <- nrow(datExpr) # nombramos renglones 

#======= Parte siete =====
#cargamos el archivo desde la ubicación en el disco de la maquina
traitData<-  read.csv("~/R/GENOMICA/FemaleLiver-Data/ClinicalTraits.csv");
# hacemos un reconocimiento rapido de los datos disponibles 
dim(traitData) # diametro
names(traitData) #nombres 

# En la base hay datos que dado el filtrado que hicimos ya no son necesarios
# entonces los eliminamos
allTraits<-  traitData[, -c(31, 16)];
allTraits <-allTraits[, c(2, 11:36) ];
#vemos que quedaran los datos correctos
dim(allTraits)
names(allTraits)

#Forman un marco de datos análogo a los datos de expresión que 
#contendrán los rasgos clínicos.
 #Muestras de ratones hembra, se toma los nombres de renglones de la 
#base de datos de expresion 
femaleSamples <- rownames(datExpr);
traitRows <- match(femaleSamples, allTraits$Mice); #renglones de reasgos
datTraits <- allTraits[traitRows, -1];
rownames(datTraits) <-allTraits[traitRows, 1];

collectGarbage();
#========= Parte ocho ====== 
# Volver a formar los clusters 
sampleTree2 = hclust(dist(datExpr), method = "average")
# Para hacer más digerible el analisis cada rasgo estara dado por un color: 
#blanco significa lento
#rojo significa alto
#gris significa entrada faltante
# numbers2colors representa con color la entrada numerica
#signed determina la paleta que se usa 
traitColors = numbers2colors(datTraits, signed = FALSE);
# Realizar el dendograma y mostrar colores representativos
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Ejemplo de dendrograma y mapa de calor de rasgos")

#======= Parte diez ======
#Guardar el resultado de los clusters en un archivo, de manera automatica 
#se guarda en el directorio en el que estamos trabajando pero puede determinarse
#otra ruta
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")

#Construcción de redes y detección de módulos
#cargamos el archivo generado en el paso anterior  de lso datos de entrada

lnames = load(file = "~/R/GENOMICA/FemaleLiver-01-dataInput.RData");
#Ver nombres de las variables que hemos cargado 
lnames

# ========== parte dos ========

# Establecemos  elconjunto de poderes de umbral, para este experimento suaves 
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
## análisis de topología de red
#la función ayuda a escger  una potencia de umbral suave adecuada para la
#construcción de la red 
sft<- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Vemos la figura de resultados
sizeGrWindow(5, 2) # el valor de la ventana  puede cambiar acorde al poder computacional
# con que se cuente 
par(mfrow <- c(1,2));
cex1 <- 0.9;
# Índice de ajuste free.scale topology 
#en función de la potencia de umbral  escogido 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="darkseagreen");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="darksalmon")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


softPower = 6;
adjacency = adjacency(datExpr, power = softPower);


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")



#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the expression and trait data saved in the first part
lnames = load(file = "~/R/GENOMICA/FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "~/R/GENOMICA/FemaleLiver-02-networkConstruction-stepByStep.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


names(datExpr)[moduleColors=="brown"]


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


annot = read.csv(file = "~/R/GENOMICA/FemaleLiver-Data/GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]]);
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================


write.csv(geneInfo, file = "geneInfo.csv")

#RED
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================



# Load the expression and trait data saved in the first part
lnames = load(file = "~/R/GENOMICA/FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "~/R/GENOMICA/FemaleLiver-02-networkConstruction-stepByStep.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)



#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "~/R/GENOMICA/FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "~/R/GENOMICA/FemaleLiver-02-networkConstruction-stepByStep.RData");
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "~/R/GENOMICA/FemaleLiver-Data/GeneAnnotation.csv");
# Select module
module = "brown";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "~/R/GENOMICA/FemaleLiver-Data/GeneAnnotation.csv");
# Select modules
modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);


