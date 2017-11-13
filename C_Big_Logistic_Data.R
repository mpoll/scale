############################################
############################################
######       Scale Algorithm         #######
######    Logistic Specification     #######
######         Exact Version         #######
######  Last Updated: 26/06/2017 MP  #######
############################################
############################################

	### Set Seed
		run     <- Sys.getenv(c("JOB_ID", "SGE_TASK_ID"))
		taskid	<- as.integer(Sys.getenv("SGE_TASK_ID"))
		set.seed(taskid)
		ncore 	<- as.integer(Sys.getenv("NSLOTS"))
		options(cores=ncore)
	### Load Data
        source("/home/stsmba/Big_Logistic_Data/F_Big_Logistic_Data.R")
	### Enviroment Settings
        initialise  <- (taskid-1)*10^5+1
        increment   <- 10^5
    ### Run
		centering.st <- system.time(centering <- data.init(initialise,increment))
        centering.alpha <- centering$grad.log.pi
        centering.alpha.p <- centering$lap.log.pi
    ### Save Simulation
		fnm 	<- paste("/scratch/stsmba/Big_Logistic_Data/Data_Files_Cent/D_",taskid,".RData",sep="")
		cat("savingfile",fnm,"\n")
		save(centering.st,centering,centering.alpha,centering.alpha.p,taskid,file=fnm)
