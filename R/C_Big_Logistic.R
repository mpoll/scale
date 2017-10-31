############################################
############################################
######       Scale Algorithm         #######
######    Logistic Specification     #######
######  Exact Version Multiple DSZ   #######
######  Last Updated: 27/07/2017 MP  #######
############################################
############################################


    ### Set Seed
        run                 <- Sys.getenv(c("JOB_ID", "SGE_TASK_ID"))
        taskid              <- as.integer(Sys.getenv("SGE_TASK_ID"))
        set.seed(taskid)
        ncore               <- as.integer(Sys.getenv("NSLOTS"))
        options(cores=ncore)

    ### Load Data
        source("/home/stsmba/Big_Logistic_Data_Multiple/Scale.R")
        load("/home/stsmba/Big_Logistic_Data_Multiple/DCent_Big_Logistic.RData")

    ###
        reps        <- 100; reps.vec <- numeric(0); for(i in 1:reps){reps.vec <- c(reps.vec,rep(i,6))}
        spec.mat        <- rbind(rep(c(1:6),reps),reps.vec)

        alpha.cent      <- alpha.cent.mult[spec.mat[1,taskid],,drop=FALSE]
        alpha.p.cent    <- alpha.p.cent.mult[spec.mat[1,taskid]]
        alpha.cent.sq   <- alpha.cent.sq.mult[spec.mat[1,taskid]]
        phi.cent        <- phi.cent.mult[spec.mat[1,taskid]]

        dsz         <- dsz.mult[spec.mat[1,taskid]]
        replicate   <- spec.mat[2,taskid]

        source("/home/stsmba/Big_Logistic_Data_Multiple/F_Big_Logistic_Data.R")



    ### Specification of ScaLE
        ss.on <- TRUE # Flag indicating whether to use subsampling (SCALE) or no subsampling (QSMC)
        scale.version <- scale.exact

        p.num       <- 2^10      # Number of particles
        t.inc       <- 0.25
        T.extend    <- 1
        ss.size     <- 32

        base.vars <- c(3,3,3,3,3)
        x.init <- matrix(0,dimen,p.num); for(j in 1:dimen){x.init[j,] <- rnorm(p.num,beta.star[j]+1/(sqrt(dsz)),base.vars[j]*dsz^(-1/2))}

    ### File Name
        fnm 	<- paste("/scratch/stsmba/Big_Logistic_Data_Multiple/D_",spec.mat[1,taskid]+4,"_R",replicate,"_J",run[1],"_T",run[2],".RData",sep="")

        print(fnm); print(dsz); print(replicate)

    # Simulations
        # Initialise Algorithm
        time.elapse <- system.time(simn <- scale(p.num,t.inc,T.extend,ss.phi,ss.phiC,dimen,transform,un.transform,T.start=0,x.init=x.init,ss.size=ss.size,ess.thresh=0.5,resamp.method=resid.resamp,neg.wei.mech=scale.zero.wei,prev.simn=NULL,progress.check=TRUE,phi.record=TRUE,resamp.freq=p.num-1,theta=NULL,p.path.renew=p.path.renew,scale.version=scale.version))
        # Algorithm Loop
        repeat{
            time.inc <- system.time(simn <- scale.extend(simn,t.inc,T.extend,scale.version))
            time.elapse <- time.elapse+time.inc
            cat("savingfile",fnm,"\n");
            save.image(file=fnm)}