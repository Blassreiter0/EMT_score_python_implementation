import numpy as np
import pandas as pd
from scipy.stats import kstest
import os
import matplotlib.pyplot as plt
from scipy.stats import mstats



def EMT76GS(log2Exp,input_gseid):
    print("Calculating EMT Score ....")

    EMTSignature = pd.read_excel("./Gene_signatures/76GS/EMT_signature_76GS.xlsx", header=0)
    EMTSignature = EMTSignature.dropna(subset=['Affymetrix Probe'])

    log2Exp = log2Exp.merge(EMTSignature[['Gene Symbol']], left_on='Gene Symbol', right_on='Gene Symbol')

    print(str(log2Exp.shape[0]) + "gene's expression values found")

    genes = log2Exp[['Gene Symbol']]

    log2Exp.to_csv("./geodata/76GS_genes.csv", index=False)
    
    log2Exp = log2Exp.drop(columns=['Gene Symbol'])

    ## Add sd
    print("adding noise to the data...")
    sdVal = np.random.normal(loc=0, scale=0.01, size=(log2Exp.shape[1]))
    log2Exp = log2Exp.astype(float)
    for i in range(log2Exp.shape[1]):
        log2Exp.iloc[:, i] += sdVal[i]


    ## get the weights for each gene
    if "CDH1" in genes.iloc[:,0].tolist():
        cdh1index = genes.iloc[:,0].tolist().index("CDH1")
        ecadhExpVal = log2Exp.iloc[cdh1index]
        weightVal = np.apply_along_axis(lambda x: np.corrcoef(x, ecadhExpVal)[0, 1], 1, log2Exp)
        EMTMatWt = weightVal.reshape(-1, 1) * log2Exp
        EMTScore = EMTMatWt.sum(axis=0)
        EMTScoreMean = EMTScore.mean()
        EMTScoreStd = EMTScore - EMTScoreMean

        
    else:
        print("CDH1 gene not found- 76 GS EMT score cannot be calculated")
        EMTScoreStd = np.zeros(log2Exp.shape)


    emtWrite = np.column_stack((log2Exp.columns.values, EMTScoreStd))
    emtWrite = pd.DataFrame(emtWrite, columns=["Sample_name", "GS76"])

    outfile = './output/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    emtWrite.to_csv(outfile + '/76GS_result.csv', index=False)
    print("Done")
    
    return emtWrite



def KSScore(expMat,input_gseid):
    print("Calculating EMT Score by KS score method ....")

    EMTSignature = pd.read_excel("./Gene_signatures/KS/EM_gene_signature_cellLine_KS.xlsx", header=None)
    EMTSignature.columns = ['Gene Symbol', 'state']
    expMat = expMat.merge(EMTSignature[['Gene Symbol','state']], left_on='Gene Symbol', right_on='Gene Symbol')

    mes_rows = expMat[expMat['state'] == 'Mes']
    epi_rows = expMat[expMat['state'] == 'Epi']
    sampleScore2 = []

    for i in range(1,expMat.shape[1]-1):
        # Two sided test
        ksTwoSided = kstest(mes_rows.iloc[:, i], epi_rows.iloc[:, i], alternative='two-sided')
        # One sided test: ecdf(Mes) > ecdf(Epi)
        ksResGrt = kstest(mes_rows.iloc[:, i], epi_rows.iloc[:, i], alternative="greater")
        ## One sided test: ecdf(Epi) > ecdf(Mes)
        ksResLess = kstest(epi_rows.iloc[:, i], mes_rows.iloc[:, i], alternative="greater")
        sampleScore2.append([ksTwoSided.statistic, ksTwoSided.pvalue, ksResGrt.statistic, ksResGrt.pvalue, ksResLess.statistic, ksResLess.pvalue])

    sampleScore2 = pd.DataFrame(sampleScore2)

    ## Assign signs to EMT score of sample based on test statistic and pvalue
    finalScore = np.zeros((sampleScore2.shape[0], 1))

    for i in range(sampleScore2.shape[0]):
        if sampleScore2.iloc[i, 3] < 0.05:
            finalScore[i] = -1 * sampleScore2.iloc[i, 2]
        elif sampleScore2.iloc[i, 5] < 0.05:
            finalScore[i] = sampleScore2.iloc[i, 4]
        else:
            if sampleScore2.iloc[i, 4] == max(sampleScore2.iloc[i, 2], sampleScore2.iloc[i, 4]):
                finalScore[i] = max(sampleScore2.iloc[i, 2], sampleScore2.iloc[i, 4])
            else:
                finalScore[i] = -1 * max(sampleScore2.iloc[i, 2], sampleScore2.iloc[i, 4])
    
    finalScore = pd.DataFrame(finalScore)
    finalScore.columns = ['KS']
    col_names = expMat.columns.tolist()
    col_names.pop(0)
    col_names.pop(-1)
    col_names = pd.DataFrame({'Sample_name': col_names})
    finalScore = pd.concat([col_names,finalScore],axis=1)
    outfile = './output/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    finalScore.to_csv(outfile + '/KS_result.csv', index=False) 
    print("Done")

    return finalScore



def plot_76GSvsKS(result_76GS, result_KS,input_gseid):
    plt.scatter(result_76GS['GS76'],result_KS['KS'], marker='*', color='black')
    plt.xlabel('EMT Score (76GS)')
    plt.ylabel('EMT Score (KS)')
    outfile = './graphs_generated/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    plt.savefig(outfile +'/76GSvsKS.png')
    plt.close()

    
    
    

def MLR3(A1, A2, A3, A4):
    a, b = A4.shape
    if A4.shape[0] == A2.shape[0]:
        A5 = mstats.mnrval(A1[:,0], A4.T, method='ordinal')
        A6 = np.zeros((b,1))
        for j in range(b):
            if np.array_equal(A5[j,:], np.nan*np.ones((1,A1.shape[0]-1))):
                A6[j,0] = np.nan
            else:
                A6[j,0] = np.argmax(A5[j,:]) + 1
    else:
        index = []
        for j in range(A2.shape[0]):
            temp = np.where(A3 == np.array(A2.iloc[j,0]))
            index.append(temp[0][0])
        A7 = np.zeros((A2.shape[0],b))
        for j in range(A2.shape[0]):
            A7[j,:] = A4.iloc[index[j],:]
        A5 = mnrval(A1.iloc[:,0], A7.T, method='ordinal')
        A6 = np.zeros((b,1))
        for j in range(b):
            if np.array_equal(A5[j,:], np.nan*np.ones((1,A1.shape[0]-1))):
                A6[j,0] = np.nan
            else:
                A6[j,0] = np.argmax(A5[j,:]) + 1
    return A5, A6



def MLRScore(geneExpMat,input_gseid,annot_table):

    #Predictor + Normalizer gene expression to run MLR model for EMT prediction
    geneList = pd.read_csv("./Gene_signatures/MLR/genes_for_EMT_score.txt", header=None, delim_whitespace=True)
    geneList = pd.DataFrame(geneList)
    geneList.columns = ['MLR_genes']

    geneExpMat = geneList.merge(geneExpMat, left_on='MLR_genes', right_on='Gene Symbol')

    temp1 = set(geneExpMat['Gene Symbol'].tolist())
    temp2 = set(geneList['MLR_genes'].tolist())

    # Find the genes that are not common
    notFound = temp1.symmetric_difference(temp2)
    for i in notFound:
        print("Gene not found in dataset: " + str(i))

    geneExpMat = geneExpMat.drop(columns=['MLR_genes'])
    outfile = './geodata/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    geneExpMat.to_csv(outfile + '_MLR_exp_levels.csv', index=False) 


    DataNCI60 = pd.read_csv('./default_data/DataNCI60.txt', header=None, delim_whitespace=True)
    
    
    # col = np.arange(0,DataNCI60.shape[1],1)
    # DataNCI60 = pd.DataFrame(columns=col)
    DataNCI60[0] = DataNCI60[0].str.strip("'")
    DataNCI60 = annot_table[['Gene Symbol','ID']].merge(DataNCI60, left_on='ID', right_on=0)
    DataNCI60 = DataNCI60.drop(columns = ['ID'])


    # replace nan and none values with 0.
    DataNCI60 = pd.DataFrame(DataNCI60)
    DataNCI60.fillna(0, inplace=True)

    # taking mean expression of all probes per gene.
    DataNCI60 = DataNCI60.groupby('Gene Symbol').mean(numeric_only=True)

    # merging MLR_genes with DataNCI60
    DataNCI60 = geneList.merge(DataNCI60, left_on='MLR_genes', right_on='Gene Symbol')

    NCI_data = DataNCI60.iloc[:,1:]
    data = geneExpMat.iloc[0:25,1:]


    NCI_norm = np.mean(np.mean(NCI_data.iloc[5:25, :], axis=0), axis=0)
    new_norm = np.mean(np.mean(data.iloc[5:25, :], axis=0), axis=0)

    d = new_norm - NCI_norm
    Norm_data = data - d

    a1,a2,a3 = 0,1,2    
    c1,c2,c3 = 0,1,2    


    ## Predictions
    # 2. (CLDN7, VIM/CDH1)

    NUM = data.shape[1]
    NormData = Norm_data

    B1 = pd.read_excel('./default_data/RelevantData.xlsx', sheet_name='B1')
    GeneList1 = pd.read_excel('./default_data/RelevantData.xlsx', sheet_name='GeneList1')
    


    #YhatDataGSEEMTGenes8pair1, PredictionsDataGSEEMTGenes8pair1 = MLR3(B1, GeneList1, [{'CLDN7'}, {'VIM/CDH1'}], pd.concat([data.iloc[a1,:], data.iloc[a2,:]/data.iloc[a3,:]],axis = 1))
    #YhatDataGSEEMTGenes8pair1Norm, PredictionsDataGSEEMTGenes8pair1Norm = MLR3(B1, GeneList1, [{'CLDN7'}, {'VIM/CDH1'}], [NormData.iloc[a1,:], NormData.iloc[a2,:]/NormData.iloc[a3,:]])


    return "MLR Skipped"
