#' @title Calculate the f1,precision,recall score of the graph
#' @description Calculate the f1,precision,recall score of the graph
#' @param pred Your predicted graph
#' @param real Your real graph
#' @return f1,precision,recall score.
#' @export
#' @examples
#' pred<-matrix(c(0,0,0,0,1,0,1,1,0),nrow=3,ncol=3)
#' real<-matrix(c(0,0,0,0,1,0,1,0,0),nrow=3,ncol=3)
#' indicators(pred,real)
indicators<-function(pred,real){
  TP=sum((real==1)&(pred==1))
  FN=sum((real==1)&(pred==0))
  FP=sum((real==0)&(pred==1))
  TN=sum((real==0)&(pred==0))
  recall=TP/(TP+FN)
  precision=TP/(TP+FP)
  F1=2*precision*recall/(precision+recall)
  data.frame(recall,precision,F1)
}
parents<-function(G,i){
  which(G[,i]==1)
}
compareG<-function(G1,G2){
  which(colSums(abs(G1-G2))>0)
}

#' @title Generate a random graph
#' @description Generate a random graph based on the given dimension size and indegree
#' @param dim The random graph dimension
#' @param indegree The average indegree of random graph for each nodes
#' @return Return a random graph
#' @export
#' @examples
#' randomGraph(dim=10,indegree=1)
randomGraph<-function(dim,indegree){
  G<-matrix(0,nrow=dim,ncol=dim)
  repeat{
    r=rbinom((dim*dim-dim)/2,1,indegree*2/(dim-1))
    if(sum(r)/dim==indegree){
      G[upper.tri(G)]<-r
      return(G)
    }
  }

  # for(i in 2:dim){
  #   if(i<=indegree){
  #     G[1:(i-1),i]=1
  #   }else{
  #     G[]
  #     G[sample.int(i-1,indegree),i]=1
  #   }
  # }
  #return(G)
}
#' @title synthetic nonlinear data base on the graph
#' @description synthetic nonlinear data base on the graph. The data generation mechanism is y=scale(a1b1x^2+a2b2x^3+a3b3x^4+a4b4sin(x)+a5b5sin(x^2)).
#' @param G An adjacency matrix.
#' @param sample_num The number of samples you want to synthetic
#' @param ratio The noise ratio. It will grow or shrink the value of the noise.
#' @param return_noise Whether return the noise of each nodes for further analysis.
#' @return Return a synthetic data
#' @export
#' @examples
#' G<-matrix(c(0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),nrow = 4,ncol = 4)
#' data=synthetic_data_nonlinear(G,100)
synthetic_data_nonlinear<-function(G,sample_num,ratio=1,return_noise=FALSE){
  f<-function(x){
    return(data.frame(x^2,x^3,x^4,sin(x),sin(x^2)))
  }
  if(return_noise){
    noise=data.frame()
  }
  if(is.null(colnames(G))){
    g<-empty.graph(paste0("V",1:ncol(G)))
  }else{
    g<-empty.graph(colnames(G))
  }
  amat(g)<-G
  order<-node.ordering(g)
  data<-matrix(0,nrow=sample_num,ncol=length(order))
  if(is.null(colnames(G))){
    colnames(data)<-paste0("V",1:ncol(G))
  }else{
    colnames(data)<-colnames(G)
  }

  #pb <- txtProgressBar(0,length(order),style = 3)
  for(i in 1:length(order)){
    #setTxtProgressBar(pb, i)
    if(length(g$nodes[[order[i]]]$parents)==0){
      data[,order[i]]<-runif(sample_num,-1,1)*ratio
      if(return_noise){
        if(ncol(noise)==0){
          noise=data.frame(data[,order[i]])
        }else{
          noise=cbind(noise,data[,order[i]])
        }
      }
    }else{
      parents=g$nodes[[order[i]]]$parents
      x=data[,parents,drop=F]
      for(p in parents){
        a=runif(5,-3,3)
        b=rbinom(5,1,0.5)
        b[round(runif(1,1,5))]<-1
        a=a[b==1]
        data[,order[i]]<-data[,order[i]]+scale(t(t(a)%*%t(as.matrix(f(x[,p])[,b==1]))))
      }
      ei=runif(sample_num,-1,1)*ratio
      data[,order[i]]<-data[,order[i]]+ei
      if(return_noise){
        if(ncol(noise)==0){
          noise=data.frame(ei)
        }else{
          noise=cbind(noise,data.frame(ei))
        }
      }
    }
  }
  if(return_noise){
    return(list(data=as.data.frame(data),noise=noise))
  }else{
    return(as.data.frame(data))
  }

}

#' @title synthetic linear data base on the graph
#' @description Synthetic linear data base on the graph. The noises are sampled from the super-gaussian distribution. The coefficients are sample from U(-1,-0.5),U(0.5,1)
#' @param G An adjacency matrix.
#' @param sample_num The number of samples you want to synthetic
#' @param ratio The noise ratio It will grow or shrink the value of the noise
#' @param return_noise Whether return the noise of each nodes for further analysis.
#' @return Return a synthetic data
#' @export
#' @examples
#' G<-matrix(c(0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0),nrow = 4,ncol = 4)
#' data=synthetic_data_linear(G,100)
synthetic_data_linear<-function(G,sample_num,ratio=1,return_noise=FALSE){
  if(is.null(colnames(G))){
    g<-empty.graph(paste0("V",1:ncol(G)))
  }else{
    g<-empty.graph(colnames(G))
  }
  if(return_noise){
    noise=data.frame()
  }
  amat(g)<-G
  order<-node.ordering(g)
  data<-matrix(0,nrow=sample_num,ncol=length(order))
  if(is.null(colnames(G))){
    colnames(data)<-paste0("V",1:ncol(G))
  }else{
    colnames(data)<-colnames(G)
  }


  #pb <- txtProgressBar(0,length(order),style = 3)
  for(i in 1:length(order)){
    #setTxtProgressBar(pb, i)
    if(length(g$nodes[[order[i]]]$parents)==0){
      #data[,order[i]]<-runif(sample_num,-1,1)
      data[,order[i]]<-rnorm(sample_num)*ratio
      #data[,order[i]]<-rlnorm(sample_num)
      if(return_noise){
        if(ncol(noise)==0){
          noise=data.frame(data[,order[i]])
        }else{
          noise=cbind(noise,data[,order[i]])
        }
      }
    }else{
      parents=sort(g$nodes[[order[i]]]$parents)
      x=data[,parents,drop=F]
      a=runif(length(parents)+1,0.5,1)
      b=rbinom(length(parents)+1,2,0.5)
      b[b==0]=-1
      a=a*b
      x=data.frame(x,c=1)
      #data[,order[i]]<-as.matrix(x)%*%a+runif(sample_num,-1,1)*rate
      e=rnorm(sample_num)
      s=sign(e)
      data[,order[i]]<-scale(as.matrix(x)%*%a+s*e^2*ratio)
      #data[,order[i]]<-as.matrix(x)%*%a+s*e^2*rate
      #data[,order[i]]<-as.matrix(x)%*%a+rnorm(sample_num)
      if(return_noise){
        if(ncol(noise)==0){
          noise=data.frame(s*e^2*ratio)
        }else{
          noise=cbind(noise,s*e^2*ratio)
        }
      }
    }
  }
  if(return_noise){
    return(list(data=as.data.frame(data),noise=noise))
  }else{
    return(as.data.frame(data))
  }
}
