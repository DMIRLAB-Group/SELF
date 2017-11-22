score_ml_cd<-function(G,D,bw,nrounds,
                      gamma,score_type,booster,...){
  NodeScore<-rep(0,ncol(G))
  NodeScore<-updateScore_ml_cd(G,D,NodeScore,1:ncol(G),bw =bw,nrounds=nrounds,
                               gamma=gamma,score_type=score_type,booster=booster,...)
  return(NodeScore)
}
updateScore_ml_cd<-function(G,D,NodeScore,nodes,score_type="bic",bw ="nrd0",booster="gbtree",gamma=10,nrounds=30,...){
  Vn<-ncol(G)
  Dn<-nrow(D)
  bicterm=log(Dn)/Dn/2
  for(i in nodes){
    pa<-parents(G,i)
    if(length(pa)==0){
      N<-D[,i]
      #N<-scale(N)
      if(is.factor(D[,i])){
        si=length(unique(D[,i]))
        d=si-1
        eta=unlist(table(N)/length(N))
        if(score_type=="log"){
          NodeScore[i]<-sum(eta*log(eta)) #log
        }else if(score_type=="bic"){
          NodeScore[i]<-sum(eta*log(eta))-d*bicterm #bic
        }else if(score_type=="aic"){
          NodeScore[i]<-sum(eta*log(eta))-d/Dn #aic
        }
        next
      }else{
        if(score_type=="log"){
          f <- stats::approxfun(stats::density(N,bw=bw,from=min(N),to=max(N)))
          NodeScore[i]<-sum(log(f(N)))/Dn
        }else if(score_type=="bic"){
          f <- stats::approxfun(stats::density(N,bw=bw,from=min(N),to=max(N)))
          NodeScore[i]<-sum(log(f(N)))/Dn
        }else if(score_type=="aic"){
          f <- stats::approxfun(stats::density(N,bw=bw,from=min(N),to=max(N)))
          NodeScore[i]<-sum(log(f(N)))/Dn-1
        }

        next
      }
    }
    x<-D[,pa,drop=F]
    y<-D[,i]
    #x=stats::predict(dummyVars(~.,data=x,fullRank = T),x)


    if(is.factor(y)){
      collevels<-stats::na.omit(levels(y))
      y<-factor(y,labels=0:(length(levels(y))-1))
    }
    x<-apply(x,2,as.numeric)
    x<-as.matrix(x)
    #x<-as.data.frame(x)
    if(is.factor(y)){
        num_class<-length(unique(y))
        model<-xgboost(data=x,label=as.integer(y)-1,
                       objective="multi:softmax",num_class=num_class,verbose=0,
                       nrounds=nrounds,gamma=gamma,booster=booster,save_period=NULL,...)
    }else{
      if(booster=="lm"){
        model<-.lm.fit(x,as.numeric(y))
      }else {
        model<-xgboost(data=x,label=as.numeric(y),
                       verbose=0,nrounds=nrounds,gamma=gamma,booster=booster,save_period=NULL,...)
      }
    }
    #xgb---end
    if(booster=="lm"){
      y_hat<-x%*%model$coefficients
    }else{
      y_hat<-stats::predict(model,x)
    }

    if(is.factor(y)){
      N<-as.numeric(y)-as.numeric(y_hat)-1
    }else{
      N<-as.numeric(y)-as.numeric(y_hat)
    }
    #N=scale(N)




    if(is.factor(y)){
      pa<-parents(G,i)
      pa_name<-names(D)[pa]
      dt<-as.data.table(D)
      ri<-nrow(dt[,.N,by=pa_name])
      #ri<-nrow(ddply(D,as.quoted(pa_name),length))
      si=length(unique(D[,i]))
      if(length(ri)==0){
        ri=1
      }
      d=ri*(si-1)
      eta=unlist(table(N)/length(N))
      if(score_type=="log"){
        NodeScore[i]<-sum(eta*log(eta)) #log
      }else if(score_type=="bic"){
        NodeScore[i]<-sum(eta*log(eta))-d*bicterm #bic
      }else if(score_type=="aic"){
        NodeScore[i]<-sum(eta*log(eta))-d/Dn #aic
      }

    }else{
      if(booster=="gbtree"){
        dump=xgb.dump(model)
        d=length(grep(dump,pattern = "leaf"))
        d=d+length(pa) #add new feature for penalty term
      }else{
        d=length(pa)+1  #gblinear
      }
      if(score_type=="log"){
        f <- stats::approxfun(stats::density(N,bw=bw,from=min(N),to=max(N)))
        NodeScore[i]<-sum(log(f(N)))/Dn
      }else if(score_type=="bic"){
        f <- stats::approxfun(stats::density(N,bw=bw,from=min(N),to=max(N)))
        NodeScore[i]<-sum(log(f(N)))/Dn-d*bicterm
      }else if(score_type=="aic"){
        f <- stats::approxfun(stats::density(N,bw=bw,from=min(N),to=max(N)))
        NodeScore[i]<-sum(log(f(N)))/Dn-bicterm
      }

    }
  }

  return(NodeScore)
}
#' @title Fast Hill-Climbing
#' @description The function for the causal structure learning.
#' @param D Input Data.
#' @param G An initial graph for hill climbing. Default: empty graph.
#' @param min_increase Minimum score increase for faster convergence.
#' @param score_type You can choose "bic","log","aic" score to learn the causal struture. Default: bic
#' @param file Specifies the output folder and its path to save the model at each iteration.
#' @param verbose Show the progress bar for each iteration.
#' @param save_model Save the meta data during the iteration so that you can easily restore progress and evaluate the model during iteration.
#' @param bw the smoothing bandwidth which is the parameter of the function stats::density(Kernel stats::density Estimation)
#' @param booster Choose the regression method, it could be "lm", "gbtree" and "gblinear". The "lm" and "gblinear" is the linear regression methods and "gbtree" is the nonlinear regression method. Default: gbtree
#' @param gamma The parameter in xgboost: minimum loss reduction required to make a further partition on a leaf node of the tree. the larger, the more conservative the algorithm will be.
#' @param nrounds the maximum number of trees for xgboost.Default:30.
#' @param ... other parameters for xgboost.see also: help(xgboost)
#' @return The adjacency matrix of the casual structure.
#' @export
#' @examples
#' \dontrun{
#' #x->y->z
#' set.seed(0)
#' x=rnorm(4000)
#' y=x^2+runif(4000,-1,1)*0.1
#' z=y^2+runif(4000,-1,1)*0.1
#' data=data.frame(x,y,z)
#' fhc(data,gamma=10,booster = "gbtree")
#'
#' #x->y->z linear data
#' set.seed(0)
#' x=rnorm(4000)
#' y=3*x+runif(4000,-1,1)*0.1
#' z=3*y+runif(4000,-1,1)*0.1
#' data=data.frame(x,y,z)
#' fhc(data,booster = "lm")
#'
#'#randomGraph with linear data
#'
#'set.seed(0)
#'G=randomGraph(dim=10,indegree=1.5)
#'data=synthetic_data_linear(G=G,sample_num=4000)
#'fitG=fhc(data,booster = "lm")
#'indicators(fitG,G)
#'}
#'

fhc<-function(D,G=NULL,min_increase=0.01,score_type="bic",
                      file="",verbose=TRUE,save_model=FALSE,
                      bw ="nrd0",booster="gbtree",gamma=10,nrounds=30,...){
  min_history_diff=0
  h=5
  #if(h<2)stop("h should be greater than 2")
  if(is.null(G)){
    G=matrix(0,nrow=ncol(D),ncol=ncol(D))
    initG=G
  }else{
    initG<-G
  }
  if(ncol(G)!=ncol(D)){
    stop("the number of nodes should be consistent to the data")
  }
  Vn<-ncol(G)
  Dn<-nrow(D)
  if(save_model){
    if(file!=""){
      if(!file.exists(file)){
        dir.create(file)
      }
      path=paste0(file,"/")
    }else{
      path=""
    }
    t=1

    file<-sprintf(paste0(path,"GD%s.RData"),t)
    while (file.exists(file)){
      t=t+1
      file<-sprintf(paste0(path,"GD%s.RData"),t)
    }
    if(t==1){
      if(file.exists(paste0(path,"initResult"))){
        load(paste0(path,"initResult"))
      }else{
        initResult<-list()
        history_score=rep(0,h)
        initResult$G<-G
        initResult$D<-D
        initResult$NodeScore<-score_ml_cd(G=initResult$G,D=initResult$D,bw =bw,nrounds=nrounds,
                                          gamma=gamma,score_type=score_type,booster=booster,...)
        initResult$score<-sum(initResult$NodeScore)
        if(save_model){
          save(initResult,file = paste0(path,"initResult.RData"))
        }
      }
      t=0
    }else{
      t=t-1
      file<-sprintf(paste0(path,"GD%s.RData"),t)
      load(file)
      if(h>length(history_score)){
        history_score<-c(history_score,rep(0,h-length(history_score)))
      }
    }
  }else{
      initResult<-list()
      history_score=rep(0,h)
      initResult$G<-G
      initResult$D<-D
      initResult$NodeScore<-score_ml_cd(initResult$G,initResult$D,bw =bw,nrounds=nrounds,
                                        gamma=gamma,score_type=score_type,booster=booster,...)
       initResult$score<-sum(initResult$NodeScore)
    t=0
  }

  if(!is.matrix(G)){
    G<-as.matrix(G)
  }
  #bestScore=initResult$score
  bestScore=-Inf
  result=initResult
  bestResult=initResult


  repeat{
    oldG<-G
    bestG<-G
    GList<-list()
    for (i in 1:ncol(G)) {
      for (j in 1:ncol(G)) {
        if(i!=j){
          # here i should choose the top ten j for each i in the graph,j the
          res<-AddDelReverseLine(G,i=i,j=j)
          GList<-append(GList,res)
        }
      }
    }
    GList<-unique(GList)
    # for(i in 1:length(GList)){
    #   if(all(GList[[i]]==initG)){
    #     GList[[i]]<-NULL
    #     break
    #   }
    # }
    if(verbose){
      pb <- utils::txtProgressBar(0,length(GList),style = 3)
    }
    bestIndex=0
    for(k in seq.int(length.out = length(GList))){
      if(verbose){
        #print(sprintf("k=%d/%d",k,length(GList)))
        #print(bestG)
        utils::setTxtProgressBar(pb, k)
      }
      nodes<-compareG(G,GList[[k]])

      result$G<-GList[[k]]
      result$NodeScore<-updateScore_ml_cd(GList[[k]],D,initResult$NodeScore,nodes,
                                          bw =bw,nrounds=nrounds,gamma=gamma,score_type=score_type,
                                          booster=booster,...)

      #print(result$NodeScore)
      result$score<-sum(result$NodeScore)
      score<-result$score
      if(score>bestScore){
        bestScore=score
        bestResult=result
        bestIndex=k
        oldG<-bestG
        bestG<-GList[[k]]
      }
    }

    G<-bestG
    len=length(history_score)
    history_score[2:len]<-history_score[1:(len-1)]
    history_score[1]=bestResult$score
    if(abs(initResult$score-bestResult$score)<min_increase)break
    if(diff(range(diff(history_score)))<min_history_diff)break
    initResult<-bestResult
    if(verbose){
      print(paste0(c("bestscore:",sum(bestScore))))
    }
    if(save_model){
      t=t+1
      output<-sprintf(paste0(path,"GD%s.RData"),t)
      save(G,D,initResult,t,history_score,file=output)
    }

    if(all(oldG==G))break
print(G)
  }
  return(bestG)
}
