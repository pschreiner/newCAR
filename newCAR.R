#!/usr/bin/env Rscript
#
# Usage:
# Rscript newCAR.R config.txt

options(warn=-1)
args <- commandArgs(TRUE)
conf <- readLines(args[1])

source("newCAR_functions.R")

#####################
## Set Environment ##
#####################
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(MASS))				# Optional (needed for model building)

cores <- detectCores()
if(cores>6) { cores=6 }
register(MulticoreParam(workers = cores))

###############################
## Read in files from Config ##
###############################
cat(paste0("Reading files indicated in ", args[1], "...\n"))

# Disease (Required)
if(any(grepl("Disease:",conf))) {
	dis_file <- unlist(strsplit(conf[which(grepl("Disease:",conf))],"Disease:"))[2]
	dis_file <- gsub(" ","",dis_file)
	dis_file <- gsub("\t","",dis_file)
	dis <- read.table(dis_file,header=T,sep="\t") # Matrix of Tumor Sample Expression Data
	if(any(unlist(dis))<0) {
		stop(cat(paste0("ERROR: Negative values in", dis_file, "are not permitted. Quitting...")))
	}
} else {
	stop("ERROR: No \"Disease:\" file has been specified in your configuration file.")
}

# Control (Required)
if(any(grepl("Control:",conf))) {
	nodis_file <- unlist(strsplit(conf[which(grepl("Control:",conf))],"Control:"))[2]
	nodis_file <- gsub(" ","",nodis_file)
	nodis_file <- gsub("\t","",nodis_file)
	nodis <- read.table(nodis_file,header=T,sep="\t") # Matrix of Control Tissue Expression Data
	if(any(unlist(nodis))<0) {
		stop(cat(paste0("ERROR: Negative values in", nodis_file, "are not permitted. Quitting...")))
	}
} else {
	stop("No \"Control:\" file has been specified in your configuration file.")
}

# Microarray - Expression (Required)
if(any(grepl("Microarray:",conf))) {
	valueMatrix_file <- unlist(strsplit(conf[which(grepl("Microarray:",conf))],"Microarray:"))[2]
	valueMatrix_file <- gsub(" ","",valueMatrix_file)
	valueMatrix_file <- gsub("\t","",valueMatrix_file)
	valueMatrix <- read.table(valueMatrix_file,header=T,sep="\t") # Optional
	if(any(unlist(valueMatrix))<0) {
		stop(cat(paste0("ERROR: Negative values in", valueMatrix_file, "are not permitted. Quitting...")))
	}
} else {
	stop("No \"Microarray:\" file has been specified in your configuration file.")
}

# Microarray - Detection Calls
if(any(grepl("Detection:",conf))) {
	presenceMatrix_file <- unlist(strsplit(conf[which(grepl("Detection:",conf))],"Detection:"))[2]
	presenceMatrix_file <- gsub(" ","",presenceMatrix_file)
	presenceMatrix_file <- gsub("\t","",presenceMatrix_file)
	presenceMatrix 	<- read.table(presenceMatrix_file,header=T,sep="\t") # Optional
	predictModel=TRUE
} else { predictModel=FALSE }

# Membrane-Associated Definition (Required)
if(any(grepl("Membrane:",conf))) {
	mpbd_file <- unlist(strsplit(conf[which(grepl("Membrane:",conf))],"Membrane:"))[2]
	mpbd_file <- gsub(" ","",mpbd_file)
	mpbd_file <- gsub("\t","",mpbd_file)
	mpbd <- read.table(mpbd_file,header=TRUE,sep="\t")
} else {
	stop("No \"Membrane:\" file has been specified in your configuration file.")
}

# Protein Data
if(any(grepl("Protein:",conf))) {
	hpa_file <- unlist(strsplit(conf[which(grepl("Protein:",conf))],"Protein:"))[2]
	hpa_file <- gsub(" ","",hpa_file)
	hpa_file <- gsub("\t","",hpa_file)
	hpa <- read.table(hpa_file,header=TRUE,sep="\t")
	protein=TRUE
} else {
	stop("No \"Protein:\" file has been specified in your configuration file.")
}

# Cluster of Differentiation Definition (Optional)
if(any(grepl("CD:",conf))) {
	cd_ref_file <- unlist(strsplit(conf[which(grepl("CD:",conf))],"CD:"))[2]
	cd_ref_file <- gsub(" ","",cd_ref_file)
	cd_ref_file <- gsub("\t","",cd_ref_file)
	cd_ref <- read.table(cd_ref_file,header=TRUE,sep="\t")
	clustdiff=TRUE
} else { clustdiff=FALSE }

# Cancer-Testis Definition (Optional)
if(any(grepl("CancerTestis:",conf))) {
	ct_ref_file <- unlist(strsplit(conf[which(grepl("CancerTestis:",conf))],"CancerTestis:"))[2]
	ct_ref_file <- gsub(" ","",ct_ref_file)
	ct_ref_file <- gsub("\t","",ct_ref_file)
	ct_ref <- read.table(ct_ref_file,header=TRUE,sep="\t",skip=1)
	cancertesty=TRUE
} else { cancertesty=FALSE }

# Name of the Output (Optional)
if(any(grepl("Name:",conf))) {
	dname <- unlist(strsplit(conf[which(grepl("Name:",conf))],"Name:"))[2]
	dname <- gsub(" ","",dname)
	dname <- gsub("\t","",dname)
} else {
	dname <- format(Sys.Date(),"%b%d")
}
dir.create(as.character(dname))
setwd(as.character(dname))

############################
## Data Matrix Formatting ##
############################
cat("Formatting Data Matrices...\n")
dis <- dis[rownames(dis) %in% rownames(nodis),]
dis <- dis[rownames(dis) %in% rownames(valueMatrix),]
dis <- dis[order(rownames(dis)),]

nodis <- nodis[rownames(nodis) %in% rownames(dis),]
nodis <- nodis[order(rownames(nodis)),]

# Microarray Expression Value Matrix
valueMatrix <- valueMatrix[rownames(valueMatrix) %in% rownames(dis),]
valueMatrix <- valueMatrix[order(rownames(valueMatrix)),]

# Microarray Presence/Absence Matrix
if(isTRUE(predictModel)) {
	presenceMatrix <- presenceMatrix[rownames(presenceMatrix) %in% rownames(dis),]
	presenceMatrix <- presenceMatrix[order(rownames(presenceMatrix)),]
}

# (Qualitative) Protein Expression Data
if(isTRUE(protein)) {
	hpa <- hpa[rownames(hpa) %in% rownames(dis),]
	hpa <- hpa[order(rownames(hpa)),]
}

#####################################
## Expression Value Transformation ##
#####################################
cat("Executing Expression Value Transformation...\n")
# Step 1: weightedCumulPerc
dis 		<- weightedCumulPerc(dis)
nodis		<- weightedCumulPerc(nodis)
valueMatrix	<- weightedCumulPerc(valueMatrix)

# Step 2: Quantile Normalization (Transformed Microarray expression data used as target distribution)
valueMatrix	<- normalize.quantiles(as.matrix(valueMatrix),copy=F)
dis		<- normalize.quantiles.use.target(as.matrix(dis),target=valueMatrix[,1],copy=F)
nodis		<- normalize.quantiles.use.target(as.matrix(nodis),target=valueMatrix[,1],copy=F)

valueMatrix 	<- as.data.frame(valueMatrix)
dis 		<- as.data.frame(dis)
nodis 		<- as.data.frame(nodis)

write.table(dis,"Disease.TEVs.txt",quote=F,col.names=T,row.names=T,sep="\t")
write.table(nodis,"Control.TEVs.txt",quote=F,col.names=T,row.names=T,sep="\t")

res <- merge(nodis, dis, by=0)
rownames(res) <- res[,1]
res[,1] <- NULL
res <- res[order(rownames(res)),]

###############################################################
## Assess Gene Status using Microarray Presence/Absence Data ##
###############################################################
if(isTRUE(predictModel)) {
	cat("Assessing Gene Status in Tumor and Control...\n")
	# Build Model using Known Absent, Marginal, Present values from an array
	inputData <- unlist(res)
        inputData <- as.numeric(inputData[inputData!=0])

	if(abs(median(inputData)-median(unlist(valueMatrix)))>0.1) {
		cat("\n\nWARNING - Minor Severity: The median of the transformed values in your dataset vs median of the transformed values in the logistic regression model vary by more than 10% of the possible values (0-1). This warning is just a note to the user. \n\n")
	} else if(abs(median(inputData)-median(unlist(valueMatrix)))>0.2) {
		cat("\n\nWARNING - Moderate Severity: The median of the transformed values in your dataset vs median of the transformed values in the logistic regression model vary by more than 20% of the possible values (0-1).  Check the distributions of these values to ensure you are making meaningful predictions...\n\n")	
	} else if(abs(median(inputData)-median(unlist(valueMatrix)))>0.3) {
		cat("\n\nWARNING - Major Severity: The median of the transformed values in your dataset vs median of the transformed values in the logistic regression model vary by more than 30% of the possible values (0-1).  Model predictions would not be meaningful.  Check the distributions of these values for further details. Skipping Model Prediction...\n\n")
		break
	}

	ordLog <- cbind(as.character(unlist(presenceMatrix)),unlist(valueMatrix))
	colnames(ordLog) <- c("ProbeCall","Value")
	ordLog <- as.data.frame(ordLog)
	ordLog$ProbeCall <- factor(ordLog$ProbeCall,levels=c("A","P"),ordered=TRUE)
	ordLog$Value <- as.numeric(as.character(ordLog$Value))
	ordLog_model <- glm(ProbeCall ~ Value,data=ordLog,family="binomial")

	## Gene Status Assessment - Control
	# Prepare Transformed Values for "Absent" or "Present" Prediction
	pred_df_nodis <- cbind("?",unlist(res[,1:ncol(nodis)]))
	colnames(pred_df_nodis) <- c("ProbeCall","Value")
	pred_df_nodis <- as.data.frame(pred_df_nodis)
	pred_df_nodis$Value <- as.numeric(as.character(pred_df_nodis$Value))
	
	# Make Predictions Based on Transformed Values
	pred_nodis <- predict(ordLog_model,pred_df_nodis,type="response")
	pred_nodis <- ifelse(pred_nodis>=0.5,"P","A")
	pred_mat_nodis <- matrix(pred_nodis,ncol=ncol(nodis),nrow=nrow(res))	
	rownames(pred_mat_nodis) <- rownames(res)
	colnames(pred_mat_nodis) <- colnames(nodis)

	# Perform Calculations
	ap_nodis <- bplapply(1:nrow(pred_mat_nodis),function(x) {
		a_nodis <- length(which(pred_mat_nodis[x,]=="A"))
		p_nodis <- length(which(pred_mat_nodis[x,]=="P"))
		c(a_nodis,p_nodis,round( (p_nodis/ncol(nodis))*100 ,1))
	})
	ap_nodis <- t(as.data.frame(ap_nodis))
	rownames(ap_nodis) <- rownames(pred_mat_nodis)
	colnames(ap_nodis) <- c("Control_NumAbsent","Control_NumPresent","Control_Present%")

	## Gene Status Assessment - Disease
	pred_df_dis <- cbind("?",unlist(res[,(ncol(nodis)+1):ncol(res)]))
	colnames(pred_df_dis) <- c("ProbeCall","Value")
	pred_df_dis <- as.data.frame(pred_df_dis)
	pred_df_dis$Value <- as.numeric(as.character(pred_df_dis$Value))
	
	# Make Predictions Based on Transformed Values
	pred_dis <- predict(ordLog_model,pred_df_dis,type="response")
	pred_dis <- ifelse(pred_dis>=0.5,"P","A")
	pred_mat_dis <- matrix(pred_dis,ncol=ncol(dis),nrow=nrow(res))	
	rownames(pred_mat_dis) <- rownames(res)
	colnames(pred_mat_dis) <- colnames(dis)
	write.table(data.frame(Gene=rownames(pred_mat_dis),pred_mat_dis),
		"Disease.TEV_GeneStatusPredictions.xls",
		quote=F,col.names=TRUE,row.names=FALSE,
		sep="\t"
	)

	# Perform Calculations
	ap_dis <- bplapply(1:nrow(pred_mat_dis),function(x) {
		a_dis <- length(which(pred_mat_dis[x,]=="A"))
		p_dis <- length(which(pred_mat_dis[x,]=="P"))
		c(a_dis,p_dis,round( (p_dis/ncol(dis))*100 ,1))
	})
	ap_dis <- t(as.data.frame(ap_dis))
	rownames(ap_dis) <- rownames(pred_mat_dis)
	colnames(ap_dis) <- c("Disease_NumAbsent","Disease_NumPresent","Disease_Present%")

	ap <- cbind(ap_dis,ap_nodis)
	write.table(data.frame(Gene=rownames(ap),ap),
		paste0(dname,".Detection_Calling_Summary.txt"),
		quote=FALSE,row.names=FALSE,col.names=TRUE,
		sep="\t"
	)
}

# Set Outlier Values - by default, outlier is a value below the median value in control or above the median value in control
nd_med <- median(unlist(nodis))
d_med <- median(unlist(dis))

################################################################
## Perform Calculations to Assess Tumor vs Control Expression ##
################################################################
cat("Identifying CAR Targets...\n")
summary_res <- bplapply(1:nrow(res),function(i) {	
	nodis_min <- min(as.numeric(res[i,1:ncol(nodis)]))
	nodis_mean <- mean(as.numeric(res[i,1:ncol(nodis)]))
	nodis_med <- median(as.numeric(res[i,1:ncol(nodis)]))
	nodis_sd <- sd(as.numeric(res[i,1:ncol(nodis)]))
	nodis_max <- max(as.numeric(res[i,1:ncol(nodis)]))	

	if(	any(grepl("Aorta",colnames(nodis),fixed=TRUE)) | 
		any(grepl("Aorta",colnames(nodis),fixed=TRUE)) | 
		any(grepl("Aorta",colnames(nodis),fixed=TRUE)) |
		any(grepl("Aorta",colnames(nodis),fixed=TRUE))
	) {
		inds <- which(nodis[i,]>nd_med)
		if(length(inds) == 0) {
			nd_num_out <- 0
			perc_nd_out_most <- "NA"
			name_nd_out_most <- " "
		} else {
			tmp <- colnames(nodis)[inds]	# Above (in Control) or Below (in Disease) the Median of Values within the Cohort is designated threshold for outliers
			spl <- strsplit(tmp,"_")
			tmp <- sapply(1:length(spl),function(x) spl[[x]][2])
			tmp <- table(tmp)[order(-table(tmp))]
			anyVitalTissue <- ifelse(grepl("Aorta",tmp) |
						grepl("Coronary",tmp) |
						grepl("Brain",tmp) |
						grepl("Heart",tmp),
						TRUE,
						FALSE
			)
			nd_num_out <- sum(tmp)
			perc_nd_out_most <- paste(round((tmp[1]/nd_num_out)*100,0),"%",sep="")
			name_nd_out_most <- names(tmp)[1]
		}
	} else {
		inds <- which(nodis[i,]>nd_med)
		if(length(inds) == 0) {
			nd_num_out <- 0
			perc_nd_out_most <- "NA"
			name_nd_out_most <- " "
			anyVitalTissue <- "NA"
		} else {
			nd_num_out <- length(inds)
			perc_nd_out_most <- "NA"
			name_nd_out_most <- paste(colnames(nodis)[inds],collapse=",")
			anyVitalTissue <- "NA" 
		}
	}
	if(is.null(name_nd_out_most)) { name_nd_out_most = "NONE" }

	dis_min <- min(as.numeric(res[i,(ncol(nodis)+1):ncol(res)]))
	dis_mean <- mean(as.numeric(res[i,(ncol(nodis)+1):ncol(res)]))
	dis_med <- median(as.numeric(res[i,(ncol(nodis)+1):ncol(res)]))
	dis_sd <- sd(as.numeric(res[i,(ncol(nodis)+1):ncol(res)]))
	dis_max <- max(as.numeric(res[i,(ncol(nodis)+1):ncol(res)]))	

	tmp <- dis[i,dis[i,]<d_med]			# Above (in Control) or Below (in Disease) the Median of Values within the Cohort is designated threshold for outliers
	d_num_out <- length(tmp)
	d_num_out_perc <- round((d_num_out/ncol(dis))*100,1)
	d_num_out_perc <- paste0(d_num_out_perc,"%")
	if(d_num_out==0) {
		name_d_out_most = "NONE"
	} else { 
		name_d_out_most <- paste(colnames(tmp),collapse=",")
	}
	pval <- wilcox.test(as.numeric(res[i,1:ncol(nodis)]),as.numeric(res[i,(ncol(nodis)+1):ncol(res)]),alternative="less",exact=TRUE)$p.value

	# Add Human Protein Atlas Data
	if(isTRUE(protein)) {
		x <- hpa[i,]
		hpa_median <- median(as.numeric(x[!is.na(x)]))
		hpa_mean <- mean(as.numeric(x[!is.na(x)]))
		perc_na <- round((length(x[is.na(x)])/ncol(hpa)) * 100,1)
		perc_na <- paste(perc_na,"%",sep="")
		if(is.na(hpa_median)) {
			al_low <- "NA"
			al_med <- "NA"
			al_high <- "NA"
		} else {
			al_low <- length(x[!is.na(x) & x>=1]) # at least "low" (i.e. low or higher than low)
			al_med <- length(x[!is.na(x) & x>=2])
			al_high <- length(x[!is.na(x) & x>=3])
		}
	} else {
		hpa_median <- "NA"
		hpa_mean <- "NA"
		perc_na <- "NA"
		al_low <- "NA"
		al_med <- "NA"
		al_high <- "NA"
	}

	med_diff <- dis_med - nodis_med
	c(dis_min, dis_med, dis_mean, dis_sd, dis_max, d_num_out, d_num_out_perc, name_d_out_most, nodis_min, nodis_med, nodis_mean, nodis_sd, nodis_max, nd_num_out, perc_nd_out_most, name_nd_out_most, med_diff, pval, hpa_median, hpa_mean, perc_na, al_low, al_med, al_high, anyVitalTissue)
})

summary_res <- t(as.data.frame(summary_res))
rownames(summary_res) <- rownames(res)
if(isTRUE(predictModel)) {
	summary_res <- cbind(summary_res,ap[,c("Disease_Present%","Control_Present%")])
	colnames(summary_res) <- c("Dis_Min","Dis_Median","Dis_Mean","Dis_SD","Dis_Max","Dis_NumOutliers","Dis_PercOutliers","Dis_OutliersTopName","Out_Min","Out_Median","Out_Mean","Out_SD","Out_Max","Out_NumOutliers","Out_PercOutliersTop","Out_OutliersTopName","Difference in Median","Wilcoxon_Pvalue","HPA Median","HPA Mean","HPA % \"NA\"","Protein Expression >= \"Low\"","Protein Expression >= \"Medium\"","Protein Expression >= \"High\"","Potential Expression in Vital Tissue","Disease_Present%","Control_Present%")
} else {
	colnames(summary_res) <- c("Dis_Min","Dis_Median","Dis_Mean","Dis_SD","Dis_Max","Dis_NumOutliers","Dis_PercOutliers","Dis_OutliersTopName","Out_Min","Out_Median","Out_Mean","Out_SD","Out_Max","Out_NumOutliers","Out_PercOutliersTop","Out_OutliersTopName","Difference in Median","Wilcoxon_Pvalue","HPA Median","HPA Mean","HPA % \"NA\"","Protein Expression >= \"Low\"","Protein Expression >= \"Medium\"","Protein Expression >= \"High\"","Potential Expression in Vital Tissue")
}
summary_res <- as.data.frame(summary_res)
write.table(data.frame(Gene=rownames(summary_res),summary_res),
	paste(dname,".AllMetrics.txt",sep=""),
	quote=FALSE,row.names=FALSE,col.names=TRUE,
	sep="\t"
)

summary_res[,1] <- as.character(summary_res[,1]); summary_res[,2] <- as.character(summary_res[,2]);
summary_res[,3] <- as.character(summary_res[,3]); summary_res[,4] <- as.character(summary_res[,4]);

mpbd <- mpbd[as.character(mpbd[,1]) %in% rownames(summary_res),]
mpbd <- mpbd[order(as.character(mpbd[,1])),]

# Table only including genes that have potential membrane association
if(any(rownames(summary_res) %in% as.character(mpbd[,1]))) {
	membs <- summary_res[rownames(summary_res) %in% as.character(mpbd[,1]),]
} else {
	stop("Membrane gene names do not match data matrix. Quitting...")
}

membs <- data.frame(row.names=rownames(membs),MembraneAssociated=mpbd[,2],membs)
membs <- membs[membs$MembraneAssociated=="Y",]
write.table(data.frame(Gene=rownames(membs),membs),
	paste(dname,".MembraneGenes.txt",sep=""),
	quote=FALSE,row.names=FALSE,col.names=TRUE,
	sep="\t"
)

membs$Wilcoxon_Pvalue <- as.numeric(as.character(membs$Wilcoxon_Pvalue))
fin <- membs[membs$Wilcoxon_Pvalue<=0.05,]
fin$Wilcoxon_Pvalue <- as.numeric(as.character(fin$Wilcoxon_Pvalue))
fin$Dis_Median <- as.numeric(as.character(fin$Dis_Median)); fin$Out_Median <- as.numeric(as.character(fin$Out_Median))
fin$Dis_Mean <- as.numeric(as.character(fin$Dis_Mean)); fin$Out_Mean <- as.numeric(as.character(fin$Out_Mean))

# Membrane genes with significanctly different expression in tumor vs control (5% level) 
targets <- fin[order(fin$Wilcoxon_Pvalue,-abs(fin$Dis_Median-fin$Out_Median),-abs(fin$Dis_Mean-fin$Out_Mean)),]
write.table(data.frame(Gene=rownames(targets),targets),
	paste(dname,".MembraneGenes_p05.txt",sep=""),
	quote=FALSE,row.names=FALSE,col.names=TRUE,
	sep="\t"
)

##############################
## Format and Print Results ##
##############################
cat("Creating Graphics...\n")
# List of computationally predicted Tumor-Associated Antigens
toptargs <- membs[as.numeric(as.character(membs$Wilcoxon_Pvalue))<=0.05,]
toptargs <- toptargs[as.numeric(as.character(toptargs$Dis_Median))>=0.5,]
toptargs <- toptargs[as.numeric(as.character(toptargs$Out_Median))<=0.3,]
if(isTRUE(protein)) {
	toptargs <- toptargs[which(as.integer(as.character(toptargs$Protein.Expression.....Low.))<=7 | toptargs$Protein.Expression.....Low.=="NA"),]
}

# Establish ranking metric
rnk1 <- -log10(toptargs$Wilcoxon_Pvalue)
rnk2 <- rnk1 * (as.numeric(as.character(toptargs$Dis_Median))-as.numeric(as.character(toptargs$Out_Median)))
toptargs$Rank_Metric <- rnk2
toptargs <- toptargs[order(-toptargs$Rank_Metric),]

if(nrow(toptargs)>=1) {
	write.table(data.frame(Gene=rownames(toptargs),toptargs),
		paste(dname,".TopTargets_ProteinwithNA.txt",sep=""),
		quote=FALSE,row.names=FALSE,col.names=TRUE,
		sep="\t"
	)
	toptargs_nona <- toptargs[as.character(toptargs$Protein.Expression.....Low.)!="NA",]
	write.table(data.frame(Gene=rownames(toptargs_nona),toptargs_nona),
		paste(dname,".TopTargets.txt",sep=""),
		quote=FALSE,row.names=FALSE,col.names=TRUE,
		sep="\t"
	)
} else {
	print("No Tumor-Associated Antigens Identified.")
}

## Create Boxplots ##
invisible(targetsBoxplot(res,rownames(toptargs),10,paste(dname,".Targets.MembraneGenes.Boxplot.pdf",sep="")))

# Cluster of Differentiation
if(isTRUE(clustdiff)) {
	if(any(rownames(toptargs) %in% as.character(cd_ref[,2]))) {
		toptargs_cd <- toptargs[rownames(toptargs) %in% as.character(cd_ref[,2]),]
	} else {
		cat("No CD genes found. Check the names in your data matrices and reference files. Skipping...")
		break
	}
	if(nrow(toptargs_cd)>=1) {
		write.table(data.frame(Gene=rownames(toptargs_cd),toptargs_cd),
			paste(dname,".Targets.CD_Family.txt",sep=""),
			quote=FALSE,row.names=FALSE,col.names=TRUE,
			sep="\t"
		)
	}
	invisible(targetsBoxplot(res,rownames(toptargs_cd),10,paste(as.character(dname),".Targets.CD_Family.Boxplot.pdf",sep="")))
}

# Cancer-testis Antigens
if(isTRUE(cancertesty)) {
	if(any(rownames(toptargs) %in% as.character(ct_ref[,2]))) {
		toptargs_ct <- toptargs[rownames(toptargs) %in% as.character(ct_ref[,2]),]
	} else {
		cat("No Cancer-Testis genes found. Check the names in your data matrices and reference files. Skipping...")
		break
	}
	if(nrow(toptargs_ct)>=1) {
		write.table(data.frame(Gene=rownames(toptargs_ct),toptargs_ct),
			paste(dname,".Targets.CTgenes.txt",sep=""),
			quote=FALSE,row.names=FALSE,col.names=TRUE,
			sep="\t"
		)
	}
	invisible(targetsBoxplot(res,rownames(toptargs_ct),10,paste(as.character(dname),".Targets.CTgenes.Boxplot.pdf",sep="")))
}

cat("== Tumor-Associated Antigens have been successfully predicted! ==\n")
