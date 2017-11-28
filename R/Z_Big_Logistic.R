############################################
############################################
######       Scale Algorithm         #######
######       Logistic Specification  #######
######         Exact Version         #######
######  Last Updated: 18/07/2017 MP  #######
############################################
############################################


    ### Set Seed
    #    run                 <- Sys.getenv(c("JOB_ID", "SGE_TASK_ID"))
    #    taskid              <- as.integer(Sys.getenv("SGE_TASK_ID"))
    #    set.seed(taskid)
    #    ncore               <- as.integer(Sys.getenv("NSLOTS"))
    #    options(cores=ncore)

    ### Load Data
    #    source("/home/stsmba/Big_Logistic_Data/Scale_R")
    #    source("/home/stsmba/Big_Logistic_Data/F_Big_Logistic_Data.R")
    #    load("/home/stsmba/Big_Logistic_Data/D_Big_Logistic.RData")

    ### Specification of ScaLE
        ss.on <- TRUE # Flag indicating whether to use subsampling (SCALE) or no subsampling (QSMC)
        scale_version <- scale_exact
        #spec.mat <- rbind(c(2^c(10:11)),c(10,10,20,20))

        p.num       <- 2^10 #spec.mat[1,taskid]      # Number of particles
        t.inc       <- 0.1
        T.extend    <- 0.1
        ss.size     <- 10 #spec.mat[2,taskid]

    ### File Name
        fnm 	<- paste("Big_Logistic_Data_P",p.num,"_S",ss.size,".RData",sep="")

    # Simulations
        # Initialise Algorithm
        scale_logistic_start <- function() {
            time.elapse <- system.time(simn <<- scale(p.num,t.inc,T.extend,ss.phi,ss.phiC,dimen,transform,un.transform,T.start=0,x.init=NULL,ss.size=ss.size,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale_zero.wei,prev.simn=NULL,progress.check=TRUE,phi.record=TRUE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew,scale_version=scale_version));
            list(time.elapse=time.elapse,simn=simn)
        }
        # Algorithm Loop
        scale_logistic_step<- function(time.elapse,simn) {
            time.inc <- system.time(simn <<- scale_extend(simn,t.inc,T.extend,scale_version))
            #cat("savingfile",fnm,"\n");
            #save.image(file=fnm)
            list(time.elapse=time.elapse,simn=simn)
        }
