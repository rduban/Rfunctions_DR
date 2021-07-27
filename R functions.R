if(!require(psych)){
  install.packages('psych')  #instalar paquete psych
  require(psych)}
if(!require(dplyr)){
  install.packages('dplyr')  #instalar paquete dplyr
  require(dplyr)}
if(!require(lavaan)){install.packages('lavaan')} 
require(lavaan)       #Activar lavaan package. Se usa para Analisis de ecuaciones estructurales
if(!require(lsr)){
  install.packages('lsr')  #instalar paquete psych
require(lsr)}  
  
################################
mycorr<-function(mydata,group){ #esta funcion recibe un conjunto de variables cuantitativas y arroja la matriz de correlaciones
  ## para obtener una matriz de correlaciones con cada diagonal segmentada por grupo la variable dicotomica debe ser la primera en el dataframe
  
  if(isFALSE(group)==TRUE){ 
    x1<-as.data.frame(round(cor(mydata),3)) #identificando las correlaciones significativas
  }
  if(isTRUE(group)==TRUE){
    
    borrar <- names(mydata[1])
    gr1<-subset(mydata,mydata[1]==unique(mydata[1])[1,]) #filtramos por grupo
    gr1<-gr1[ , !(names(gr1) %in% borrar)] #eliminamos variable grupo
    gr2<-subset(mydata,mydata[1]==unique(mydata[1])[2,])
    gr2<-gr2[ , !(names(gr2) %in% borrar)]
    x1<-as.data.frame(round(cor(gr1),3)) #matriz de correlaciones grupo 1
    x2<-as.data.frame(round(cor(gr2),3)) #matriz de correlaciones grupo 2
    y<-as.data.frame(matrix(NA, ncol=ncol(x1),nrow = nrow(x1)))
    for(i in 1:nrow(x1)){
      for(j in 1:nrow(x1)){
        if(i<j){
          y[j,i]<-x1[j,i]} 
        if(i>j){
          y[j,i]<-x2[j,i]}}}}
  z<-as.data.frame(matrix(NA, ncol=ncol(x1),nrow = nrow(x1)))
  if(isFALSE(group)==TRUE){
    for(i in 1:nrow(x1)){
      for(j in 1:nrow(x1)){
        if(cor.test(mydata[[i]],mydata[[j]])$p.value<0.05){
          z[i,j]<-paste(x1[i,j],"*")
          if(cor.test(mydata[[i]],mydata[[j]])$p.value<0.001){
            z[i,j]<-paste(x1[i,j],"**")}
        }else(z[i,j]<-x1[i,j])}}}
  if(isTRUE(group)==TRUE){
    for(i in 1:nrow(x1)){
      for(j in 1:nrow(x1)){
        if(i>j){
          dataf<-gr1}
        else(dataf<-gr2)
        if(cor.test(dataf[[i]],dataf[[j]])$p.value<0.05){
          z[i,j]<-paste(y[i,j],"*")
          if(cor.test(dataf[[i]],dataf[[j]])$p.value<0.001){
            z[i,j]<-paste(y[i,j],"**")}}else(z[i,j]<-y[i,j])}}}
  
  
  rownames(z)<-names(x1)
  if(isTRUE(group)==TRUE){
    message(paste("Diagonal izquierda:",unique(mydata[1])[1,])); 
    message(paste("Diagonal derecha:",unique(mydata[1])[2,]))}
  message("* = significancia <  0.05");
  message("** = significancia <  0.001")
  
  z
}

mycfa<-function(mymodels,mydata){  
  #Recibe una lista que incluye los modelos y el dataframe con las variables
  #arroja una lista con los modelos especificados, objetos lavaan, la tabla de indices de ajuste,
  #y la tabla donde se determina el numero de indices de ajuste dentro los umbrales aceptables para determinar el modelo ganador
  
  modelos<-list()
  for(i in 1:length(mymodels)){
    modelos[[i]] <-cfa(mymodels[[i]],data = mydata)  # computamos cfa y guardamos los objetos lavaan en una lista
  }
  names(modelos)<-names(mymodels)
  
  indextable0<-as.data.frame(matrix(NA,nrow=length(modelos),ncol = length(fitmeasures(modelos[[1]])))) 
  colnames(indextable0)<-names(fitmeasures(modelos[[1]]))
  row.names(indextable0)<-names(modelos)
  for(i in 1:nrow(indextable0)){
    indextable0[i,]<-round(as.vector(unlist(fitmeasures(modelos[[i]]))),3)  #obtemos los indices de ajuste de cada modelo y los guardamos en una tabla
  }
  indextable<-cbind(indextable0[3:5],indextable0[9:10],indextable0[13:14],indextable0[16],indextable0[23],indextable0[26],indextable0[29],indextable0[38:40],indextable0[42]) #seleccionamos indices de interes
  indextablevalues<-as.data.frame(matrix(NA,nrow = nrow(indextable), ncol = ncol(indextable)))
  names(indextablevalues)<-names(indextable)
  row.names(indextablevalues)<-row.names(indextable)
  indextablevalues <- rename(indextablevalues, chi2_gl = df)
  
  # a continuacion asignamos 1 cada vez que cada uno de los indices de ajuste de cada modelo se encuentran dentro de los umbrales aceptados. 
  #Se asigna 0 cuando no ocurre. Estos valores se guardan en una nueva tabla
  indextablevalues$chisq<- ifelse(indextable$chisq==min(indextable$chisq),1,0)
  indextablevalues$chi2_gl<- ifelse(nrow(mydata)>=500 & indextable$chisq/indextable$df<=2,1,ifelse(nrow(mydata)>=500 & indextable$chisq/indextable$df==min(indextable$chisq/indextable$df),1,0))
  indextablevalues$pvalue<-ifelse(indextable$pvalue<=0.05,1,0)
  indextablevalues$cfi<- ifelse(indextable$cfi>=0.9 & indextable$cfi<=1,1,0)
  indextablevalues$tli<- ifelse(indextable$tli>=0.9 & indextable$tli<=1,1,0)
  indextablevalues$nfi<- ifelse(indextable$nfi>=0.9,1,0)
  indextablevalues$pnfi<- ifelse(indextable$pnfi==max(indextable$pnfi),1,ifelse(max(indextable$pnfi)-indextable$pnfi<=0.09,1,0))
  indextablevalues$rni<-ifelse(nrow(mydata)<500 & indextable$rni==max(indextable$rni),1,0)
  indextablevalues$rmsea<- ifelse(nrow(mydata)>=100 & indextable$rmsea<=0.08,1,0)
  indextablevalues$rmsea.pvalue<-ifelse(indextable$rmsea.pvalue<=0.05,1,0)
  indextablevalues$srmr<- ifelse(nrow(mydata)<500 & indextable$srmr==min(indextable$srmr),1,0)
  indextablevalues$gfi<- ifelse(indextable$gfi==max(indextable$gfi),1,0)
  indextablevalues$agfi<- ifelse(indextable$agfi>=0.9,1,ifelse(indextable$agfi==max(indextable$agfi),1,0))
  
  indextablevalues$pgfi<- ifelse(indextable$pgfi==max(indextable$pgfi),1,ifelse(sqrt((max(indextable$pgfi)-indextable$pgfi)^2)<=0.09,1,0))
  indextablevalues$ecvi<- ifelse(nrow(mydata)<500 & indextable$ecvi==min(indextable$ecvi),1,0)
  
  indextablevalues$total <- rowSums(indextablevalues) #se suman el numero de indices de ajuste dentro de los umbrales aceptados
  
  output<-NULL
  output$models<-mymodels
  output$cfamodels<-modelos
  output$all_fit_index<-indextable0
  output$fit_index<-indextable
  output$values_fit_index<-indextablevalues
  
  cat('The comparison of the models and the scoring of the fit indices were based on the guidelines of Hair, et al. (2014). Cite: Hair, J., Black, W., Babin, B. & Anderson, R. (2014). Multivariate Data Analysis (7th ed.). USA. Pearson.' )
  output
}


sepWords<-function(x,y){  #x = vector, y = separador
  
  ##entra vector de caracteres
  # sale data.frame con valores binarios
  a<-strsplit(x, split = y) #separar elementos
  listednames<-as.data.frame(t(unique(unlist(a, use.names = FALSE)))) #tomar la unidad de cada elemento
  mitabla<-as.data.frame(matrix(NA, nrow = length(x), ncol = dim(listednames)[2])) #crear tabla para output
  names(mitabla)<-listednames[1,]  # poner los nombres a la tabla de output
  for(i in 1:length(x)){
    for(j in 1:dim(mitabla)[2]){
      x<-ifelse(listednames[j] %in% a[[i]],1,0)   # si la unidad se encuentra en el vector dentro de la lista dame 1 sino 0
      mitabla[i,j]<-x   #agregar valor a tabla output
    }
  }
  
  mitabla
}

my_anova<-function(x,n){
  #esta funcion recibe un objeto anova y arroja una tabla con "df", f value, p value y eta^2
  y<-summary(x)
  z<-etaSquared(x,type=2,anova = FALSE)
  tabla<-as.data.frame(matrix(NA, ncol=4,nrow = n)); names(tabla)<-c("Df","F value","p value","eta^2")
  tabla[,1]<- as.vector(y[[1]][["Df"]][1:n])
  tabla[,2]<-as.vector(round(y[[1]][["F value"]][1:n],2))
  tabla[,3]<-as.vector(round(y[[1]][["Pr(>F)"]][1:n],4))
  tabla[,4]<-as.vector(round(z[1:n,1],4))
  rownames(tabla)<- rownames(y[[1]])[1:n]
  tabla
}

print("Este set de funciones fue desarrallo por el investigador Duban Romero. Si detecta algún inconveniente al usar las funciones por favor escribir al correo: rduban@uninorte.edu.co")