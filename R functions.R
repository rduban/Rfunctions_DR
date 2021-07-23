

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
