#!/usr/bin/env Rscript
weightedCumulPerc <- function(data_mat) {
	require(BiocParallel)

	if(any(unlist(data_mat)<0)) { return("\nNegative Values are not permitted in the input matrix\n\n") }
	data_mat <- data_mat[order(rownames(data_mat)),]
	cp_matrix <- bplapply(1:ncol(data_mat), function(i) {
		curr <- data_mat[,i]
		names(curr) <- rownames(data_mat)
		curr <- curr[order(-curr)]
		tot <- sum(data_mat[!duplicated(data_mat[,i]),i])

		for(j in 1:length(curr)) {
			if(j==1) {
				recur <- curr[1]
				vect <- 1 - (recur/tot)
			} else {
				if(curr[j]==0) {
					vect <- append(vect,0)
				} else if(curr[j] == curr[j-1]) {
					vect <- append(vect,vect[j-1])
				} else {
					recur <- recur + curr[j]
					vect <- append(vect,1 - (recur/tot))
				}
			}
		}
		names(vect) <- names(curr)
		vect <- vect[order(names(vect))]
		return(vect)
	})
	cp_matrix <- as.data.frame(cp_matrix)
	colnames(cp_matrix) <- colnames(data_mat)
	rownames(cp_matrix) <- rownames(data_mat)
	cp_matrix <- round(cp_matrix,4)

	return(cp_matrix)
}
# weightedCumulPerc(data_matrix)

targetsBoxplot <- function(res,gene_list,num_res,out_name) {
	if(length(gene_list)==0) {
		return("No Gene targets to plot.")
	} else if (length(gene_list)>num_res) {
		gene_list <- gene_list[1:num_res]
	}

	pdf(out_name)
	par(mar=c(4,8,4,2))

	# Only displaying Top 10
	bx_nodis <- res[rownames(res) %in% gene_list,1:ncol(nodis)]
	bx_nodis <- bx_nodis[match(gene_list,rownames(bx_nodis)),]
	bx_nodis <- as.data.frame(bx_nodis)
	bx_dis <- res[rownames(res) %in% gene_list,(ncol(nodis)+1):ncol(res)]
	bx_dis <- bx_dis[match(gene_list,rownames(bx_nodis)),]
	bx_dis <- as.data.frame(bx_dis)

	for(i in 1:nrow(bx_dis)) {
		if(i==1) {
			tmp <- cbind(rownames(bx_dis)[i],"In",t(bx_dis[i,]))
			bx <- rbind(tmp,cbind(rownames(bx_nodis)[i],"Out",t(bx_nodis[i,])))
		} else {
			bx <- rbind(bx,cbind(rownames(bx_dis)[i],"In",t(bx_dis[i,])))
			bx <- rbind(bx,cbind(rownames(bx_nodis)[i],"Out",t(bx_nodis[i,])))
		}
	}
	bx <- as.data.frame(bx)
	colnames(bx) <- c("Gene","Condition","Value"); rownames(bx) <- NULL
	bx$Gene <- as.character(bx$Gene); bx$Condition <- as.character(bx$Condition); bx$Value <- as.numeric(as.character(bx$Value));
	attach(bx)

	if(grepl("CT",out_name)) {
		bords <- NULL
		for(x in 1:length(gene_list)) {
			if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C1") {
				bords <- append(bords,rep("red",2))
			} else if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C2") {
				bords <- append(bords,rep("blue",2))
			} else if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C3") {
				bords <- append(bords,rep("darkgreen",2))
			} else if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C4") {
				bords <- append(bords,rep("purple",2))
			} else if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C5") {
				bords <- append(bords,rep("orange",2))
			} else if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C6a") {
				bords <- append(bords,rep("brown",2))
			} else if(as.character(ct_ref[which(grepl(gene_list[x],ct_ref[,2])),10]) == "C6b") {
				bords <- append(bords,rep("lightgray",2))
			} else {
				bords <- append(bords,rep("yellow",2))
			}
		}

		y=invisible(	boxplot(Value ~ Condition*Gene,
				col=c("blue","red"),
				main=dname, las=1,
				ylim=c(0,1),
				xlab="Transformed Expression Value", ylab="",
				border=bords, horizontal=TRUE
		))
		mtext(paste("Num_Samples_In:", y$n[1], "              Num_Samples_Out:", y$n[2], sep=" "))	

		plot.new()
		legend("center",legend=c("C1: High-Confidence Testis-Specific Coding Genes","C2: High-Confidence Testis-Specific Non-Coding Genes","C3: Moderate-Confidence Testis-Specific Coding Genes","C4: Moderate-Confidence Testis-Specific Non-Coding Genes","C5: Low-Confidence Testis-Specific Genes","C6a: Genes with Testis-Specific Transcripts","C6b: Genes without Testis-Specific Transcripts"),fill=c("red","blue","darkgreen","purple","orange","brown","lightgray","yellow"))
		mtext("Boxplot (Border) Color Key:")
	} else {
		y=invisible(	boxplot(Value ~ Condition*Gene,
				col=c("blue","red"),
				main=dname, las=1,
				ylim=c(0,1),
				xlab="Transformed Expression Value", ylab="",
				horizontal=TRUE
		))
		mtext(paste("Num_Samples_In:", y$n[1], "              Num_Samples_Out:", y$n[2], sep=" "))
	}
	detach(bx)
	dev.off()
}
