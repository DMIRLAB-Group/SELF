#include <Rcpp.h>
#include <iostream>
using namespace std;
using namespace Rcpp;

bool* walked;
bool dfsCheckCircuit(IntegerMatrix G,int current)
{
  int i;
  int Node_num=G.ncol();
  if (walked[current])
    return true;
  walked[current] = true;
  for (i = 0;i <Node_num;i++)
  {
    if (G(current,i))
      if (dfsCheckCircuit(G,i))
        return true;
  }
  walked[current] = false;
  return false;
}


bool NoCheckCircuit(IntegerMatrix G ,int start)
{
  int Node_num=G.ncol();
  walked=new bool[Node_num];
  int i;
  for(i=0;i<Node_num;i++)
    walked[i]=false;
  if (dfsCheckCircuit(G,start)){
    delete[] walked;
    return false;
  }
  delete[] walked;
  return true;
}



bool isDegreeRight(IntegerMatrix G,int i=-1){

  int ncol=G.ncol();
  int indegree,outdegree;
  bool muti=false;
  if(i<0){
    i=0;
    muti=true;
  }
  if(!muti){
    indegree=0;
    outdegree=0;
    for(int j=0;j<ncol;++j){
      if(i!=j){
        outdegree=outdegree+G(i,j);
        indegree=indegree+G(j,i);
      }
    }
    //Rprintf("outdegree=%d,indegree=%d",outdegree,indegree);
    if(indegree>2||(indegree+outdegree)==0)
      return false;
  }else{
    for(;i<ncol;++i){
      indegree=0;
      outdegree=0;
      for(int j=0;j<ncol;++j){
        if(muti&&i!=j){
          outdegree=outdegree+G(i,j);
          indegree=indegree+G(j,i);
        }
      }
      if(indegree>2||(indegree+outdegree)==0)
        return false;
    }
  }

  return true;
}


// [[Rcpp::export]]
std::vector<IntegerMatrix> AddDelReverseLine(IntegerMatrix G,int i=0,int j=0){
  G=clone(G);
  int nrow=G.nrow();
  int ncol=G.ncol();


  vector<IntegerMatrix> GList;
  if(i>=1&&j>=1){
    i--;
    j--;
    if(G(i,j)==1){
      G(i,j)=0; //del
      //G(j,i)=0;
      GList.push_back(clone(G));
        G(i,j)=0; //reverse
        G(j,i)=1;
        if(NoCheckCircuit(G,i)){
          GList.push_back(clone(G));

        }
      }else{

        if(G(j,i)==0){//add without circle
          G(i,j)=1; //add
          if(NoCheckCircuit(G,j)){
            GList.push_back(clone(G));
          }
        }


    }

}
  return GList;
}

