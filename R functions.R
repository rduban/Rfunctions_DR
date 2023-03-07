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
 if(!require(car)){
  install.packages('car')  #instalar paquete psych
require(car)}
 if(!require(lmtest)){
  install.packages('lmtest')  #instalar paquete psych
require(lmtest)}
################################
mycorr<- function(mydata,group=FALSE){ #esta funcion recibe un conjunto de variables cuantitativas y arroja la matriz de correlaciones
  ## para obtener una matriz de correlaciones con cada diagonal segmentada por grupo la variable dicotomica debe ser la primera en el dataframe
  
  if(isFALSE(group)==TRUE){ 
    x1<-as.data.frame(round(cor(na.omit(mydata)),2)) #identificando las correlaciones significativas
  }
  if(isTRUE(group)==TRUE){
    
    borrar <- names(mydata[1])
    gr1<-subset(mydata,mydata[1]==unique(mydata[1])[1,]) #filtramos por grupo
    gr1<-gr1[ , !(names(gr1) %in% borrar)] #eliminamos variable grupo
    gr2<-subset(mydata,mydata[1]==unique(mydata[1])[2,])
    gr2<-gr2[ , !(names(gr2) %in% borrar)]
    x1<-as.data.frame(round(cor(na.omit(gr1)),2)) #matriz de correlaciones grupo 1
    x2<-as.data.frame(round(cor(na.omit(gr2)),2)) #matriz de correlaciones grupo 2
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
          if(cor.test(mydata[[i]],mydata[[j]])$p.value<0.01){
            z[i,j]<-paste(x1[i,j],"**")
            if(cor.test(mydata[[i]],mydata[[j]])$p.value<0.001){
              z[i,j]<-paste(x1[i,j],"***")}}
        }else(z[i,j]<-x1[i,j])}}}
  if(isTRUE(group)==TRUE){
    for(i in 1:nrow(x1)){
      for(j in 1:nrow(x1)){
        if(i>j){
          dataf<-gr1}
        else(dataf<-gr2)
        if(cor.test(dataf[[i]],dataf[[j]])$p.value<0.05){
          z[i,j]<-paste(y[i,j],"*")
          if(cor.test(dataf[[i]],dataf[[j]])$p.value<0.01){
            z[i,j]<-paste(y[i,j],"**")
            if(cor.test(dataf[[i]],dataf[[j]])$p.value<0.001){
              z[i,j]<-paste(y[i,j],"***")}}}else(z[i,j]<-y[i,j])}}}
  
  
  rownames(z)<-names(x1)
  if(isTRUE(group)==TRUE){
    message(paste("Diagonal izquierda:",unique(mydata[1])[1,])); 
    message(paste("Diagonal derecha:",unique(mydata[1])[2,]))}
  message("* = significancia <  0.05");
  message("** = significancia <  0.01")
  message("*** = significancia <  0.001")
  
  z
}

 mycfa<-function(mymodels,mydata,estimador = "ML"){  
  #Recibe una lista que incluye los modelos y el dataframe con las variables
  #arroja una lista con los modelos especificados, objetos lavaan, la tabla de indices de ajuste,
  #y la tabla donde se determina el numero de indices de ajuste dentro los umbrales aceptables para determinar el modelo ganador
modelos<-list()
  
  for(i in 1:length(mymodels)){
    modelos[[i]] <-cfa(mymodels[[i]],data = mydata, estimator = estimador)  # computamos cfa y guardamos los objetos lavaan en una lista
  }
  names(modelos)<-names(mymodels)
  
  indextable0<-as.data.frame(matrix(NA,nrow=length(modelos),ncol = length(fitmeasures(modelos[[1]])))) 
  colnames(indextable0)<-names(fitmeasures(modelos[[1]]))
  row.names(indextable0)<-names(modelos)
  for(i in 1:nrow(indextable0)){
    indextable0[i,]<-round(as.vector(unlist(fitmeasures(modelos[[i]]))),3)  #obtemos los indices de ajuste de cada modelo y los guardamos en una tabla
  }
  # seleccionamos indices de interes  
  indextable<-indextable0[c("chisq","df","pvalue","cfi","tli","nnfi","nfi","pnfi","rni","rmsea","rmsea.ci.lower","rmsea.ci.upper","srmr","gfi","agfi","pgfi","ecvi")]
  
  indextablevalues<-as.data.frame(matrix(NA,nrow = nrow(indextable), ncol = ncol(indextable)))
  names(indextablevalues)<-names(indextable)
  row.names(indextablevalues)<-row.names(indextable)
  names(indextablevalues)[2] <- "chi2_gl"
  
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

my_esem<-function(cargas,n.factors){

  esemmodel<-vector()

  for(i in 1:n.factors){
    esemmodel[i]<-paste0("f",i,"=~",paste0(c(cargas[,i]),"*",names(cargas[,1]),collapse = "+"))
  }

  esemmodel<-paste0(esemmodel,collapse = "\n") #Model Specification
  esemmodel  
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

mydcohen<-function(df,varcon,varchar){
  a<-as.data.frame(df %>% group_by({{varchar}}) %>% summarise(M = mean(na.omit({{varcon}})),SD = sd(na.omit({{varcon}}))))
  d.cohen<-round((a[1,2]-a[2,2])/((a[1,3]+a[2,3])/2),3)
  d.cohen
}

my_t.test<- function(varscon,varfact,paired,var.equal){
  datos<-na.omit(data.frame(varscon,varfact))
  varfact<-factor(datos$varfact)
  tabla<-as.data.frame(matrix(NA,ncol = 10, nrow = 1)); names(tabla)<-c("t value","df","p.value","d.cohen",paste("M",levels(varfact)[1]),paste("SD",levels(varfact)[1]),paste("n",levels(varfact)[1]),paste("M",levels(varfact)[2]),paste("SD",levels(varfact)[2]),paste("n",levels(varfact)[2]))
  a<-tapply(datos$varscon,datos$varfact,FUN = function(x){round(na.omit(mean(x)),2)})
  b<-tapply(datos$varscon,datos$varfact,FUN = function(x){round(na.omit(sd(x)),2)})
  n<-tapply(datos$varscon,datos$varfact,FUN = function(x){round(na.omit(length(x)),2)})
  
  if(paired==TRUE){
    x<-t.test(datos$varscon~datos$varfact, paired = TRUE, var.equal = ifelse(var.equal==TRUE,TRUE,FALSE))
    tabla[1,4]<-sqrt((round(x$estimate/min(tapply(datos$varscon,datos$varfact,sd)),3))^2)
  }
  if(paired==FALSE){
    x<-t.test(datos$varscon~datos$varfact, paired = FALSE,var.equal = ifelse(var.equal==TRUE,TRUE,FALSE))
    tabla[1,4]<-sqrt((round((max(x$estimate)-min(x$estimate))/mean(tapply(datos$varscon,datos$varfact,sd)),3))^2)
  }
  
  tabla[1,1]<- round(x$statistic,2)
  tabla[1,2]<- x$parameter
  tabla[1,3]<- round(x$p.value,3)
  tabla[1,5]<-a[1]
  tabla[1,6]<-b[1]
  tabla[1,7]<-n[1]
  tabla[1,8]<-a[2]
  tabla[1,9]<-b[2]
  tabla[1,10]<-n[2]
  
  tabla
}
  my_anova2x2<-function(modelos,level){
#a<-aov(lm(autoeficaciaTec ~fase*genero+fase*contexto+fase*niv.ensenanza+fase*area, data = b))
tabla<-as.data.frame(matrix(NA,ncol = 5, nrow = length(modelos)))
dimnames(tabla)<-list(names(modelos),c("df","df.denom","F value","p value","eta.sq"))
for(i in 1:length(modelos)){
  a<-summary(modelos[[i]])
  z<-etaSquared(modelos[[i]],type=2,anova = FALSE)
  
  if(level==3){
   b<-a[[1]][3,]  
   df_dn<-a[[1]][4,1]
   eta.sq<-z[3,1]  
  }
  if(level==2){
   b<-a[[1]][2,]  
   df_dn<-a[[1]][3,1]
   eta.sq<-z[2,1]
  }
  if(level==1){
   b<-a[[1]][1,]  
   df_dn<-a[[1]][2,1]
   eta.sq<-z[1,1]
  }
  tabla[i,1]<-b[1]; tabla[i,2]<-df_dn; tabla[i,3:4]<-round(b[4:5],3);tabla[i,5]<-round(eta.sq,3)
  }
tabla
}
 
supu_rlm<-function(modelo, studentize = TRUE){
  
  if(studentize){
    modelo$residuals<-rstudent(modelo)
    x<-bptest(modelo)
      }
  if(studentize==FALSE){
    modelo$residuals<-scale(modelo$residuals)
    x<-bptest(modelo, studentize = FALSE)}
  
  breush_pagan<-c(x[[1]],x[[4]])
  x<-durbinWatsonTest(modelo)
  durbin_wat<-c(x[[2]],x[[3]])
  x<-shapiro.test(modelo$residuals)
  shapiro<-c(x[[1]],x[[2]])
  
  
  y<-round(as.data.frame(rbind(durbin_wat,shapiro,breush_pagan)),3)
  cumple<-c("No","No","No")
  if(y[1,2]>0.05){cumple[1]<-"Si"};if(y[2,2]>0.05){cumple[2]<-"Si"};if(y[3,2]>0.05){cumple[3]<-"Si"}
  
  x<-cbind(y, cumple)
  names(x)<-c("Estadistico","p.valor","Cumple")
  if(length(modelo$terms)>3){
    cat("VIF:",vif(modelo))
  }

car::scatterplot(modelo$fitted.values,modelo$residuals, ylab = "Residuales",xlab = "Estimacion", main = "Homocedasticidad")
x
}
  
myFactorInvariance<-function(model,data,vec,estimator="ML"){

bf<-cfa(model,data = data ,estimator=estimator,group = vec) # configural metrical invariance
bf_metric<-cfa(model,data = data ,estimator=estimator,group = vec,group.equal=c("loadings")) #weak metrical invariance
bf_scalar<-cfa(model,data = data ,estimator=estimator,group = vec,group.equal=c("loadings","intercepts")) #strong scalar invariance
bf_factor_mean<-cfa(model,data = data ,estimator=estimator,group = vec,group.equal=c("loadings","intercepts","means")) # FACTOR MEAN INVARIANCE

a<-lavaan::anova(bf,bf_metric,bf_scalar,bf_factor_mean)


findex<-c("cfi","tli","rmsea")
fit_indices<-data.frame(rbind(lavaan::fitMeasures(bf,findex),
                              lavaan::fitMeasures(bf_metric,findex),
                              lavaan::fitMeasures(bf_scalar,findex),
                              lavaan::fitMeasures(bf_factor_mean,findex)))

delta_cfi<-round(abs(c(NA,(fit_indices[2:5,1]-fit_indices[1:4,1])[1:3])),3)
delta_tli<-round(abs(c(NA,(fit_indices[2:5,2]-fit_indices[1:4,2])[1:3])),3)
delta_rmsea<-round(abs(c(NA,(fit_indices[2:5,3]-fit_indices[1:4,3])[1:3])),3)

b<-cbind(a,fit_indices,delta_cfi,delta_tli,delta_rmsea)
row.names(b)<-c("CONFIGURAL","WEAK METRIC","STRONG SCALAR","FACTOR MEAN")
print(b)
b
}
  
my_alphaCI<-function(datos,reps=1000,n=(nrow(datos)/3)){
  id=1:nrow(datos)
  a=NULL
  for(i in 1:reps){
   x<-datos[sample(id,n, replace = F),]
    (a[i]<-as.numeric(psych::alpha(x)$total[1]))
  }
  a<-as.vector(quantile(a,probs = c(0.025,0.975)))
  a
}

print("Este set de funciones fue desarrallo por el investigador Duban Romero. Si detecta alg?n inconveniente al usar las funciones por favor escribir al correo: rduban@uninorte.edu.co")
