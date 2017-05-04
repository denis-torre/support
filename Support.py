#################################################################
#################################################################
############### Pipeline Support
#################################################################
#################################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. General support #####
import os, sklearn.decomposition, requests, json
import rpy2.robjects as robjects
import numpy as np
import pandas as pd
from plotly.graph_objs import Scatter3d, Layout

##### 2. Other libraries #####


#######################################################
#######################################################
########## S1. R Connection
#######################################################
#######################################################

#############################################
########## 1. Get R source
#############################################

def setupR(rSource):

	# Setup R
	r = robjects.r

	# Add source
	r.source(rSource)

	# Return object
	return r

#######################################################
#######################################################
########## S2. Plotting
#######################################################
#######################################################

#############################################
########## 1. 3D Plotly Scatter
#############################################

def plot3dScatter(dataframe, categoricalColumn, varLabelDict, annotationColumn=False, title='', PCs=['PC1', 'PC2', 'PC3'], colors=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'], width=900, height=600, size=5):

	# Define empty plot
	p = []

	# Get unique categories
	uniqueCategories = list(dataframe[categoricalColumn].unique())

	# Get dict with unique categories
	categoricalExpressionDict = {x:dataframe[dataframe[categoricalColumn] == x] for x in uniqueCategories}

	# Loop through categories
	for i in range(len(uniqueCategories)):

	    # Get category
	    category = uniqueCategories[i]

	    # Get plot dataframe
	    plotDataframe = categoricalExpressionDict[category]
	    
	    # Get annotation, if specified
	    annotation = plotDataframe[annotationColumn] if annotationColumn else ''

	    # Append trace
	    p.append(
	        Scatter3d(
	            x = plotDataframe[PCs[0]],
	            y = plotDataframe[PCs[1]],
	            z = plotDataframe[PCs[2]],
	            mode='markers',
	            text=annotation,
	            name=category,
	            marker=dict(
	                size=size,
	                color=colors[i],
	                opacity=0.9
	            ),
	        )
	    )
	    
	# Add layout
	layout = Layout(
	    title=title,
	    hovermode='closest',
	    width=width,
	    height=height,
	    scene=dict(
	        xaxis=dict(title=varLabelDict[PCs[0]]),
	        yaxis=dict(title=varLabelDict[PCs[1]]),
	        zaxis=dict(title=varLabelDict[PCs[2]]),
	    ),
	    margin=dict(
	        l=50,
	        r=50,
	        b=50,
	        t=50
	    )
	)

	# Prepare figure
	fig = dict(data=p, layout=layout)

	# Return figure
	return fig

#######################################################
#######################################################
########## S3. Enrichr
#######################################################
#######################################################

#############################################
########## 1. Add gene lists 
#############################################

def addGeneLists(geneset):
	ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
	genes_str = '\n'.join(geneset)
	payload = {
	    'list': (None, genes_str),
	}
	response = requests.post(ENRICHR_URL, files=payload)
	if not response.ok:
	    raise Exception('Error analyzing gene list')
	data = json.loads(response.text)
	return data

#############################################
########## 2. Get enrichment results
#############################################

def getEnrichmentResults(user_list_id, gene_set_library='GO_Biological_Process_2015', overlappingGenes=False):
	ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
	query_string = '?userListId=%s&backgroundType=%s'
	response = requests.get(
	    ENRICHR_URL + query_string % (user_list_id, gene_set_library)
	 )
	if not response.ok:
	    raise Exception('Error fetching enrichment results')

	data = json.loads(response.text)
	resultDataframe = pd.DataFrame(data[gene_set_library], columns=['rank', 'term_name', 'pvalue', 'zscore', 'combined_score', 'overlapping_genes', 'FDR', 'old_pvalue', 'old_FDR'])
	selectedColumns = ['term_name','zscore','combined_score','FDR'] if not overlappingGenes else ['term_name','zscore','combined_score','FDR', 'overlapping_genes']
	resultDataframe = resultDataframe.loc[:,selectedColumns]
	return resultDataframe

#############################################
########## 3. Get and add gene lists 
#############################################

def uploadToEnrichr(cdDataframe, column, nGenes=500):

	# Get sorted genes
	sortedGenesAscending = cdDataframe.sort_values(column, ascending=True).index.tolist()

	# Get genesets
	genesets = {'downregulated': sortedGenesAscending[:nGenes], 'upregulated': sortedGenesAscending[-nGenes:]}

	# Upload and get links
	genesetLinks = {x: addGeneLists(genesets[x]) for x in genesets.keys()}

	# Add URLs
	for key in genesetLinks.keys():
	    genesetLinks[key]['URL'] = 'http://amp.pharm.mssm.edu/Enrichr/enrich?dataset='+genesetLinks[key]['shortId']

    # Convert to dataframe
	enrichrLinkDataframe = pd.DataFrame(genesetLinks).T.reset_index().rename(columns={'index': 'geneset'})

	# Return links
	return enrichrLinkDataframe

#######################################################
#######################################################
########## S4. L1000CDS2
#######################################################
#######################################################

#############################################
########## 1. Run L1000CDS2
#############################################

def runL1000CDS2(cdDataframe, column, aggravate=False):
    # Set data
    data = {"genes": cdDataframe.index.tolist(), "vals":cdDataframe[column].tolist()}
    data['genes'] = [x.upper() for x in data['genes']]
    
    # Set configuration
    config = {"aggravate":aggravate, "searchMethod":"CD", "share":True, "combination":False, "db-version":"latest"}
    payload = {"data":data,"config":config}
    headers = {'content-type':'application/json'}
    
    # Perform request
    r = requests.post('http://amp.pharm.mssm.edu/L1000CDS2/query',data=json.dumps(payload),headers=headers)
    resCD= r.json()
    
    # Add URL
    resCD['URL'] = 'http://amp.pharm.mssm.edu/L1000CDS2/#/result/' + resCD['shareId']
    
    # Return result
    return resCD

#############################################
########## 2. Get L1000CDS2 results
#############################################

def getL1000CDS2Results(cdDataframe, column):

	# Define result dataframe and list
	resultSignatureDataframe = pd.DataFrame()
	linkList = []

	# Loop through aggravate
	for aggravate in [True, False]:
	    
	    # Run analysis
	    resCD = runL1000CDS2(cdDataframe, column, aggravate)
	    
	    # Get signature dataframe
	    signatureDataframe = pd.DataFrame(resCD['topMeta']).drop('overlap', axis=1).replace('-666', np.nan)
	    
	    # Add aggravate column
	    signatureDataframe['aggravate'] = aggravate
	    
	    # Concatenate
	    resultSignatureDataframe = pd.concat([resultSignatureDataframe, signatureDataframe])
	    
	    # Add link
	    linkList.append({'URL': resCD['URL'], 'aggravate': aggravate})
	    
	# Convert link list to dataframe
	linkDataframe = pd.DataFrame(linkList)

	# Create result dict
	resultDict = {'signatures': resultSignatureDataframe, 'links': linkDataframe}

	# Return dictionary
	return resultDict

#######################################################
#######################################################
########## S5. Statistical Analysis
#######################################################
#######################################################

#############################################
########## 1. PCA
#############################################

def runPCA(expressionDataframe, n_components=3):
    
    # Initialize PCA object
    pcaResults = sklearn.decomposition.PCA(n_components=n_components)
    
    # Run PCA
    pcaResults.fit(expressionDataframe)
    
    # Get PC dataframe
    pcaDataframe = pd.DataFrame(pcaResults.components_).T

    # Fix dimension names
    pcaDataframe.columns = ['PC'+str(x+1) % locals() for x in pcaDataframe.columns]
    pcaDataframe.index = expressionDataframe.columns
    
    # Get variance explained
    varExplainedDict = {'PC'+str(i+1): 'PC'+str(i+1)+' ('+str(round(pcaResults.explained_variance_[i], ndigits=1))+'% var. explained)'for i in range(len(pcaResults.explained_variance_))}

    # Return result
    return {'PCs': pcaDataframe, 'varExplained': varExplainedDict}

#############################################
########## 2. Plot 3D PCA
#############################################

def plot3DPCA(expressionDataframe, annotationDataframe, categoricalColumn, annotationColumn=False, title='PCA Analysis'):

	# Run PCA
	pcaResults = runPCA(expressionDataframe)

	# Create plot dataframe
	plotDataframe = pcaResults['PCs'].merge(annotationDataframe, left_index=True, right_index=True)

	# Create figure
	fig = plot3dScatter(plotDataframe, categoricalColumn, pcaResults['varExplained'], annotationColumn=annotationColumn, title=title)

	# Return fig
	return fig

#######################################################
#######################################################
########## S6. GEO2Enrichr
#######################################################
#######################################################

#############################################
########## 1. API
#############################################

def submitG2E(dataset, platform, A_cols, B_cols, is_geo=True, diffexp_method='chdir', numgenes=500, threshold=0.05, normalize=True):

	# Create dictionary
	data = dict(dataset = dataset,
				platform = platform,
				A_cols = A_cols,
				B_cols = B_cols,
				is_geo = is_geo,
				diffexp_method = diffexp_method,
				numgenes = numgenes,
				threshold = threshold,
				normalize = normalize)

	# Submit request
	r = requests.post('http://amp.pharm.mssm.edu/g2e/api/extract/geo', data=data)

	# Get extraction id
	try:
		extractionId = dict(r.json())['extraction_id']
		return extractionId
	except:
		return None

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

