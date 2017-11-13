#' @title mmpc algorithm and additive noise model
#' @description Use the mmpc algorithm to learn a skeleton and use anm to recognize the direction
#' @param data The data
#' @export
mmpcAnm<-function(data){
  fitG=mmpc(data,test="mi-g-sh") #shrinkage estimator for the mutual information
  fitG=amat(fitG)
  n=ncol(fitG)
  pb <- txtProgressBar(0,sum(fitG==1)/2,style = 3)
  count=0
  for(i in 1:n){
    for(j in 1:n){
      if(fitG[i,j]==1&&fitG[j,i]==1){
        count=count+1
        setTxtProgressBar(pb,count)
        X=data[,c(i,j)]
        fit<-getParents(X=X,method = "bivariateANM")

        if(fit[1,2]==1){ #if i to j
          fitG[j,i]=0
        }else if(fit[2,1]==1){
          fitG[i,j]=0
        }else{
          fitG[j,i]=0
          fitG[i,j]=0
        }
      }
    }
  }
  return(fitG)
}
