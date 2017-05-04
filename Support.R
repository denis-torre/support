#################################################################
#################################################################
############### Pipeline Support
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####

##### 2. Other libraries #####


#######################################################
#######################################################
########## S1. PCA Analysis
#######################################################
#######################################################

#############################################
########## 1. runPCA
#############################################

runPCA <- function(expressionDataframe)
{
	# Calculate PCA results
	pcaResults <- prcomp(t(expressionDataframe), scale=FALSE, center=TRUE)

	# Get variance explained
	varExplained <- (pcaResults$sdev)^2/sum(pcaResults$sdev^2)

	# Add names
	names(varExplained) <- colnames(pcaResults$x)

	# Get labels
	varExplainedLabels <- sapply(names(varExplained), function(x) paste0(x, ' (', round(varExplained[x]*100, digits=1),'% var. explained)'))

	# Return results
	return(list(x=pcaResults$x, var=varExplained, varLabels=varExplainedLabels))
}

#############################################
########## 2. getPcaDataframe
#############################################

getPcaDataframe <- function(expressionDataframe, PCs=c('PC1', 'PC2', 'PC3'), sampleColumn=FALSE)
{
	# Calculate PCA results
	pca <- runPCA(expressionDataframe)

	# Get matrix
	pcaMatrix <- pca$x

	# Add variance explained
	pcaMatrix <- rbind(varExplained=pca$varLabels[colnames(pcaMatrix)], pcaMatrix)

	# Get subset
	pcaMatrix <- as.data.frame(pcaMatrix[,PCs])

	# Add sample ID
	if (sampleColumn != FALSE) {
		pcaMatrix <- as.data.frame(cbind(sample_id=rownames(pcaMatrix), pcaMatrix))
	}

	# Return matrix
	return(pcaMatrix)
}

#######################################################
#######################################################
########## S2. Dataframe Labeling
#######################################################
#######################################################

#############################################
########## 1. getLabels
#############################################

getLegend <- function(ids, labels, colors=c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')) {

	# Get label-color matching
	label2color <- setNames(colors[1:length(unique(labels))], as.character(unique(labels)))

	# Get id-label matching
	id2label <- setNames(as.character(labels), ids)

	# Get color vector
	colorVector <- sapply(ids, function(x) label2color[[id2label[[x]]]])

	# Return
	return(list(colors=colorVector, legend=label2color))
}

#############################################
########## 2. getDesign
#############################################

getDesign <- function(annotationDataframe, column) {

	# Create new dataframe
	designDataframeMelt <- data.frame(rows=rownames(annotationDataframe),
	                                  columns=annotationDataframe[,column],
	                                  values=1)

	# Cast
	designDataframe <- dcast(rows ~ columns, data=designDataframeMelt, value.var='values', fill=0)

	# Fix rownames
	rownames(designDataframe) <- designDataframe$rows
	designDataframe$rows <- NULL

	# Convert to matrix
	designMatrix <- as.matrix(designDataframe)

	# Return matrix
	return(designMatrix)
}

#######################################################
#######################################################
########## S3. Expression Normalization
#######################################################
#######################################################

#############################################
########## 1. runVST
#############################################

runVST <- function(rawcountDataframe) {

	# Load library
	require(DESeq2)

	# Convert to matrix
	rawcountMatrix <- as.matrix(rawcountDataframe)

	# Perform VST normalization
	vstMatrix <- varianceStabilizingTransformation(rawcountMatrix)

	# Return matrix
	return(vstMatrix)
}

#############################################
########## 2. runVoom
#############################################

runVoom <- function(rawcountDataframe, annotationDataframe=FALSE) {

	# Load library
	require(limma)

	# Convert to matrix
	rawcountMatrix <- as.matrix(rawcountDataframe)

	# Get design, if specified
	if (annotationDataframe != FALSE) {
		designMatrix <- getDesign(annotationDataframe, colnames(annotationDataframe))
	} else {
		designMatrix <- NULL
	}

	# Perform VST normalization
	voomMatrix <- voom(rawcountMatrix, design=designMatrix)$E

	# Return matrix
	return(voomMatrix)
}

#############################################
########## 3. runSizeFactorNormalization
#############################################

runSizeFactorNormalization <- function(rawcountDataframe, log=TRUE) {

	# Load library
	require(DESeq2)

	# Convert to matrix
	rawcountMatrix <- as.matrix(rawcountDataframe)

	# Estimate size factors
	sizeFactors <- estimateSizeFactorsForMatrix(rawcountMatrix)

	# Perform VST normalization
	normalizedCountMatrix <- t(t(rawcountMatrix)/sizeFactors)

	# Log transform, if specified
	if (log) {
		normalizedCountMatrix <- log10(normalizedCountMatrix + 1)
	}

	# Return matrix
	return(normalizedCountMatrix)
}

#############################################
########## 4. Quantile Normalization
#############################################

runQuantileNormalization <- function(expressionDataframe) {

	# Load library
	require(preprocessCore)

	# Convert to matrix
	expressionMatrix <- as.matrix(expressionDataframe)

	# Normalize
	normalizedExpressionMatrix <- normalize.quantiles(expressionMatrix)

	# Fix dimension names
	rownames(normalizedExpressionMatrix) <- rownames(expressionMatrix)
	colnames(normalizedExpressionMatrix) <- colnames(expressionMatrix)

	# Return matrix
	return(normalizedExpressionMatrix)
}

#############################################
########## 5. Batch Effect Removal
#############################################

runComBat <- function(expressionDataframe, annotationDataframe, covariateFormula=NULL, batchColumn='batch') {
	
    # Load library
    require(sva)

    # Get batch
    batch <- annotationDataframe[,batchColumn]

    # Get covariate, if specified
    if (!is.null(covariateFormula)) {
        covariateDesign <- model.matrix(as.formula(covariateFormula), data=annotationDataframe)
    } else {
        covariateDesign <- NULL
    }
    
    # Run Combat
    combatDataframe <- ComBat(dat=expressionDataframe, batch=batch, mod=covariateDesign, par.prior=TRUE, prior.plots=FALSE)

    # Return dataframe
    return(combatDataframe)
}

#######################################################
#######################################################
########## S4. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. runCharacteristicDirection
#############################################

runCharacteristicDirectionOld <- function(expressionDataframe, experimentColumns, controlColumns, minVariance=0, r=1) {

	# Load libraries
	source('/Users/denis/Documents/Projects/scripts/libraries/geoDE/chdir.R')
	source('/Users/denis/Documents/Projects/scripts/libraries/geoDE/nipals.R')

	# Convert to R character
	experimentColumns <- as.character(experimentColumns)
	controlColumns <- as.character(controlColumns)

	# Get dataframes
	experimentDataframe <- expressionDataframe[,experimentColumns]
	controlDataframe <- expressionDataframe[,controlColumns]

	# Get variance
	experimentGeneVariance <- apply(experimentDataframe, 1, var)
	controlGeneVariance <- apply(controlDataframe, 1, var)

	# Get common genes
	commonGenes <- intersect(names(experimentGeneVariance)[experimentGeneVariance > minVariance], names(controlGeneVariance)[controlGeneVariance > minVariance])

	# Filter
	experimentDataframeFiltered <- experimentDataframe[commonGenes,]
	controlDataframeFiltered <- controlDataframe[commonGenes,]

	# Run Characteristic Direction
	cdResults <- chdir(controlDataframeFiltered, experimentDataframeFiltered, commonGenes)

	# Convert to dataframe
	cdDataframe <- data.frame(CD=cdResults)

	# Return results
	return(cdDataframe)
}

##########
##### New Function - Updated filter.  Should return more complete, slightly different results
##########


runCharacteristicDirection <- function(expressionDataframe, experimentColumns, controlColumns, constantThreshold=1e-5) {

	# Load libraries
	source('/Users/denis/Documents/Projects/scripts/libraries/geoDE/chdir.R')
	source('/Users/denis/Documents/Projects/scripts/libraries/geoDE/nipals.R')

	# Convert to R character
	experimentColumns <- as.character(experimentColumns)
	controlColumns <- as.character(controlColumns)

	# Get dataframes
	experimentDataframe <- expressionDataframe[,experimentColumns]
	controlDataframe <- expressionDataframe[,controlColumns]

	# Get variable genes
	experimentVariableGenes <- names(which(diag(var(t(experimentDataframe))) > constantThreshold))
	controlVariableGenes <- names(which(diag(var(t(controlDataframe))) > constantThreshold))

	# Get common genes
	commonGenes <- intersect(experimentVariableGenes, controlVariableGenes)

	# Filter dataframes
	experimentDataframeFiltered <- experimentDataframe[commonGenes,]
	controlDataframeFiltered <- controlDataframe[commonGenes,]

	# Run Characteristic Direction
	cdResults <- chdir(controlDataframeFiltered, experimentDataframeFiltered, commonGenes)

	# Convert to dataframe
	cdDataframe <- data.frame(CD=cdResults)

	# Return results
	return(cdDataframe)
}

#############################################
########## 2. runDESeq2
#############################################

runDESeq2 <- function(countDataframe, annotationDataframe, design='~ sample_type', fixColNames=TRUE) {

	# Load library
	require(DESeq2)

	# Fix colnames
	if (fixColNames == TRUE){
		colnames(countDataframe) <- gsub('.', '-', colnames(countDataframe), fixed=TRUE)
	}


	# Prepare dds object
	dds <- DESeqDataSetFromMatrix(countData = countDataframe, colData = annotationDataframe, design = as.formula(design))

	# Filter
	dds <- dds[rowSums(counts(dds)) > 1,]

	# Run analysis
	dds <- DESeq(dds)

	# Get results
	res <- as.data.frame(results(dds))

	# Return results
	return(res)
}



#######################################################
#######################################################
########## S5. Gene ID Conversion
#######################################################
#######################################################

#############################################
########## 1. convertIDs
#############################################

convertIDs <- function(genes, orgDb='org.Hs.eg.db', originalIdType='ENTREZID', newIdType='SYMBOL') {

	# Load libraries
	require(AnnotationDbi)
	require(orgDb, character.only=TRUE)

	# Convert to character
	genes <- as.character(genes)

	# Remove NAs and get unique
	genes <- unique(genes[!is.na(genes)])

	# Get conversion dataframe
	conversionDataframe <- select(get(orgDb), keys=genes, keytype=originalIdType, columns=newIdType)

	# Return dataframe
	return(conversionDataframe)
}


#######################################################
#######################################################
########## S6. Statistical Analysis
#######################################################
#######################################################

#############################################
########## 1. Matrix correlation
#############################################

correlateMatrices <- function(matrix1, matrix2, method='pearson', use='everything') {
	return(cor(t(matrix1), matrix2, method=method, use=use))
}

#############################################
########## 2. Survival Association
#############################################

getSurvivalAssociation <- function(survivalDataframe, groupDataframe) {

	# Load library
	require(survival)

	# Create survival object
	survivalObj <- Surv(time=survivalDataframe$last_checked, event=survivalDataframe$event)

	# Get difference
	sdiff <- survdiff(survivalObj ~ groupDataframe$group)

	# Get p-value
	p <- 1 - pchisq(sdiff$chisq, df=length(sdiff$n) - 1)

	# Return
	return(p)
}

#######################################################
#######################################################
########## S7. L1000CDS2
#######################################################
#######################################################

#############################################
########## 1. Count small molecules
#############################################

getL1000CDS2Counts <- function(signatureDataframe, aggravate, group=FALSE, plot=TRUE, nDrugs=10, col=FALSE, main='', ylab='', order=FALSE, legendTitle=NULL) {
    
    # Filter
    signatureDataframeFiltered <- signatureDataframe[signatureDataframe$aggravate == aggravate,]

    # Set colors
    if (col == FALSE) {
    	if (aggravate == TRUE) {
    		col <- c('#ffffb2','#fecc5c','#fd8d3c','#f03b20','#bd0026')
    	} else if (aggravate == FALSE) {
    		col <- c('#f1eef6','#bdc9e1','#74a9cf','#2b8cbe','#045a8d')
    	}
    }
    
    # If group
    if (group != FALSE) {
        
        # Count
        drugCounts <- table(signatureDataframeFiltered$pert_desc, signatureDataframeFiltered[,group])
        
        # Sort
        drugCounts <- drugCounts[order(apply(drugCounts, 1, sum), decreasing=TRUE),]
        
        # Get total counts
        totalDrugCounts <- sort(apply(drugCounts, 1, sum), decreasing=TRUE)
        
        # Filter
        drugCounts <- drugCounts[totalDrugCounts > 0,]
        
        # Plot
        if (plot) {
            if (order == FALSE) {
                order <- 1:ncol(drugCounts)
            }
            topDrugs <- names(rev(totalDrugCounts[1:nDrugs]))
            barplot(t(drugCounts[topDrugs, order]), horiz=TRUE, las=2, col=col,
                    xlab='Times in top 50 signatures', ylab=ylab, main=main)
            legend('bottomright', fill=col, legend=colnames(drugCounts)[order], title=legendTitle)
        }
        
    } else {
        
        # Count
        drugCounts <- table(signatureDataframeFiltered$pert_desc)
                
        # Sort
        drugCounts <- sort(drugCounts, decreasing=TRUE)
        
        # Filter
        drugCounts <- drugCounts[drugCounts > 0]
        
        # Plot
        if (plot) {
            barplot(rev(drugCounts[1:nDrugs]), horiz=TRUE, las=2, col=col,
                    xlab='Times in top 50 signatures', ylab=ylab, main=main)
        }
    }
    # Return counts
    return(drugCounts)
}

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

