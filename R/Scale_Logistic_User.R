############################################
############################################
######       Scale Algorithm         #######
######    Logistic Specification     #######
######         User Example          #######
######  Last Updated: 18/04/2018 MP  #######
############################################
############################################

user.logistic.example   <- function(){

############################################
############################################
### -1- Specification of Data Set
############################################
############################################

#############################################
#### 1.1 - Size, dimension, and parameters
#############################################

examp.design <- matrix(0,40,5); examp.design[,2:5] <- rnorm(160); examp.design[,1] <- 1
examp.data   <- rbinom(160,1,0.5)
beta.star    <- c(1,1,-1,2,-2) # What is the choosen centering value in the original space

if(exists("dsz")==FALSE){dsz <<- dim(examp.design)[1]}     # Size of data set
if(exists("dimen")==FALSE){dimen <<- dim(examp.design)[2]} # Dimensionality (>1)
if(exists("n.sigma")==FALSE){n.sigma <<- diag(rep(1,dimen))*dsz^(-1/2)} # Parameterise preconditioning matrix (inverse Fischers information)
if(exists("n.sigma.inv")==FALSE){n.sigma.inv <- solve(n.sigma)} # Specification of inverse preconditioning matrix

#############################################
#### 1.3 - Data Recaller
#############################################

datum     <- function(datum.idx){
    datum.x     <- examp.design[datum.idx,]
    datum.y     <- examp.data[datum.idx]
    list(idx=datum.idx,x=datum.x,y=datum.y)}


############################################
############################################
### -3- End Function Specification
############################################
############################################

}

