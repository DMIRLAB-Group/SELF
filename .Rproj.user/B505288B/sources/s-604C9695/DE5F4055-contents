# SELF 

## Overview

Provides the SELF criteria to learn causal structure. 

Details of the algorithm can be found in "SELF: A Structural Equation Embedded Likelihood Framework for Causal Discovery" (AAAI2018).

## Installation

```{r, eval = FALSE}
install.packages("SELF")
```

### Quick Start

This package contain the data synthetic process and the casual structure learning algorithm. Here are some examples to make a quick start:

```{r example}
#x->y->z
set.seed(0)
x=rnorm(4000)
y=x^2+runif(4000,-1,1)*0.1
z=y^2+runif(4000,-1,1)*0.1
data=data.frame(x,y,z)
fhc(data,gamma=10,booster = "gbtree")

#x->y->z linear data
set.seed(0)
x=rnorm(4000)
y=3*x+runif(4000,-1,1)*0.1
z=3*y+runif(4000,-1,1)*0.1
data=data.frame(x,y,z)
fhc(data,booster = "lm")

#RandomGraph linear data
set.seed(0)
G=randomGraph(dim=10,indegree=1.5)
data=synthetic_data_linear(G=G,sample_num=4000)
fitG=fhc(data,booster = "lm")
indicators(fitG,G)
```

