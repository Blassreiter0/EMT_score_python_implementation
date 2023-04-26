import numpy as np
import pandas as pd
from scipy.stats import kstest
import os
import matplotlib.pyplot as plt



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
