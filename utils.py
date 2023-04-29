import numpy as np
import pandas as pd
from scipy.stats import kstest
import os
import matplotlib.pyplot as plt
import MLR_function
import matlab



def EMT76GS(log2Exp,input_gseid):
    print("Calculating 76GS Score ....")

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

    


def MLR3(A1, A2, A3, A4):
    my_MLR_function = MLR_function.initialize()

    b,a = A4.shape
    if a == A2.shape[0]:
        A1_linear_array = np.concatenate(A1.values.T)
        A4_linear_array = np.concatenate(A4.values.T)
        A1In = matlab.double(A1_linear_array, size=(4, 1))
        A4In = matlab.double(A4_linear_array, size=(b, a))
        A5 = my_MLR_function.MLR_function(A1In, A4In)
        # A5 = mstats.mnrval(A1[:,0], A4.T, method='ordinal')
        A5 = pd.DataFrame(A5)
        A6 = np.zeros((b,1))
        for j in range(b):
            if np.array_equal(A5.iloc[j,:], np.nan*np.ones((1,A1.shape[0]-1))):
                A6[j,0] = np.nan
            else:
                A6[j,0] = np.argmax(A5.iloc[j,:]) +1
    else:
        index = []
        for j in range(A2.shape[0]):
            temp = np.where(A3 == np.array(A2.iloc[j,0]))
            index.append(temp[0][0])
        A7 = np.zeros((b,A2.shape[0]))
        for j in range(A2.shape[0]):
            A7[:,j] = A4.iloc[:,index[j]]
        
        temp1,temp2 = A7.shape
        A1_linear_array = np.concatenate(A1.values.T)
        A7_linear_array = np.concatenate(A7.values.T)
        A1In = matlab.double(A1_linear_array, size=(4, 1))
        A7In = matlab.double(A7_linear_array, size=(temp1,temp2))
        A5 = my_MLR_function.MLR_function(A1In, A7In)
        # A5 = mnrval(A1.iloc[:,0], A7.T, method='ordinal')
        A5 = pd.DataFrame(A5)
        A6 = np.zeros((b,1))
        for j in range(b):
            if np.array_equal(A5.iloc[j,:], np.nan*np.ones((1,A1.shape[0]-1))):
                A6[j,0] = np.nan
            else:
                A6[j,0] = np.argmax(A5.iloc[j,:]) + 1

    my_MLR_function.terminate()
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


    ## Now the following code takes the NCI60 dataset and takes the mean of probes per gene and selects only the genes MLR genes.

    DataNCI60 = pd.read_csv('./default_data/DataNCI60.txt', header=None, delim_whitespace=True)
    DataNCI60[0] = DataNCI60[0].str.strip("'")
    DataNCI60 = pd.DataFrame(DataNCI60)



    # DataNCI60 = annot_table[['Gene Symbol','ID']].merge(DataNCI60, left_on='ID', right_on=0)
    # DataNCI60 = DataNCI60.drop(columns = ['ID'])

    # # replace nan and none values with 0.
    # DataNCI60.fillna(0, inplace=True)

    # # taking mean expression of all probes per gene.
    # DataNCI60 = DataNCI60.groupby('Gene Symbol').mean(numeric_only=True)

    # # merging MLR_genes with DataNCI60
    # DataNCI60 = geneList.merge(DataNCI60, left_on='MLR_genes', right_on='Gene Symbol')

    # NCI_data = DataNCI60.iloc[:,1:]
    # NCI_norm = np.mean(np.mean(NCI_data.iloc[5:25, :], axis=0), axis=0)




    NCI_data = pd.read_csv('./default_data/NCI_data.txt', header=None, delim_whitespace=True)
    NCI_norm = np.mean(np.mean(NCI_data.iloc[5:25, :], axis=0), axis=0)




    data = geneExpMat.iloc[0:25,1:]
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
    
    

    YHat1, PredictionsDataGSEEMTGenes8pair1 = MLR3(B1, GeneList1, ['CLDN7', 'VIM/CDH1'], pd.concat([data.iloc[a1,:], data.iloc[a2,:]/data.iloc[a3,:]],axis = 1))
    
    YHat1_norm, PredictionsDataGSEEMTGenes8pair1Norm = MLR3(B1, GeneList1, [{'CLDN7'}, {'VIM/CDH1'}], pd.concat([NormData.iloc[a1,:], NormData.iloc[a2,:]/NormData.iloc[a3,:]],axis = 1))


    ScoreEMT1 = np.full((NUM, 1), np.nan)
    for j in range(NUM):
        if YHat1.iloc[j, 0] > YHat1.iloc[j, 2]:
            ScoreEMT1[j] = YHat1.iloc[j, 1]
        else:
            ScoreEMT1[j] = 2.0 - YHat1.iloc[j, 1]

    ScoreEMT1_norm = np.full((NUM, 1), np.nan)
    for j in range(NUM):
        if YHat1_norm.iloc[j, 0] > YHat1_norm.iloc[j, 2]:
            ScoreEMT1_norm[j] = YHat1_norm.iloc[j, 1]
        else:
            ScoreEMT1_norm[j] = 2.0 - YHat1_norm.iloc[j, 1]



    plt.figure(1)
    plt.plot(0, 0, '.k')
    plt.plot(0, 0, '.b')
    plt.plot(0, 0, '.r')
    plt.scatter(DataNCI60.iloc[c1, 1:],np.log2((DataNCI60.iloc[c2, 1:].astype(np.float32) / DataNCI60.iloc[c3, 1:].astype(np.float32) + 1)), marker='*', color='black', label='NCI60')
    plt.scatter(data.iloc[a1, :],np.log2((data.iloc[a2, :].astype(np.float32) / data.iloc[a3, :].astype(np.float32) + 1)), marker='*', color='blue', label='Data')
    plt.scatter(Norm_data.iloc[a1, :],np.log2((Norm_data.iloc[a2, :].astype(np.float32) / Norm_data.iloc[a3, :].astype(np.float32) + 1)), marker='*', color='red', label='NormData')
    plt.xlabel('CLDN7 Expression')
    plt.ylabel('VIM/CDH1 Expression')
    plt.legend()
    outfile = './graphs_generated/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    plt.savefig(outfile +'/MLRresult.jpg')
    plt.close()

    print("Done")


    ScoreEMT1_norm = pd.DataFrame(ScoreEMT1_norm)
    ScoreEMT1_norm.columns = ['MLR']
    temp = pd.DataFrame(geneExpMat.iloc[:,1:].columns)
    # temp.pop(0)
    ScoreEMT1_norm = pd.concat([temp,ScoreEMT1_norm],axis = 1)
    outfile = './output/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    ScoreEMT1_norm.to_csv(outfile + '/MLR_result.csv', index=False) 

    return ScoreEMT1_norm






def plot_76GSvsKS(result_76GS, result_KS,result_MLR,input_gseid):

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))


    ax1.scatter(result_76GS['GS76'],result_KS['KS'], marker='*', color='black')
    ax1.set_xlabel('EMT Score (76GS)')
    ax1.set_ylabel('EMT Score (KS)')

    ax2.scatter(result_76GS['GS76'],result_MLR['MLR'], marker='*', color='black')
    ax2.set_xlabel('EMT Score (76GS)')
    ax2.set_ylabel('EMT Score (MLR)')

    ax3.scatter(result_MLR['MLR'],result_KS['KS'], marker='*', color='black')
    ax3.set_xlabel('EMT Score (MLR)')
    ax3.set_ylabel('EMT Score (KS)')

    fig.suptitle('EMT scores coorelation')

    outfile = './graphs_generated/' + input_gseid 
    if not os.path.exists(outfile):
        os.makedirs(outfile)
    plt.savefig(outfile +'/EMT_results_correlation.png')
    plt.close()

