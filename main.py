import GEOparse
import pandas as pd
import numpy as np
from utils import *
import sys
import os
import subprocess



# installing necessary python packages.
os.chdir('./default_data')
subprocess.call(['python', 'setup.py', 'install'])
os.chdir('..')




input_gseid = str(input("Enter the GSE id: "))
# input_gseid = "GSE21422"

gse = GEOparse.get_GEO(geo=input_gseid, destdir="./geodata")

# make a table of all gsm samples with their microarray data
gse_table = gse.pivot_samples('VALUE')


try:
    # download the platform from GEO
    gpl_id = gse.metadata['platform_id'][0]
    
    gpl = GEOparse.get_GEO(geo=gpl_id, destdir='./geodata')

    # make the annotation data
    annot_table = gpl.table

    # merge the two tables based on the 'ID_REF' column in the gsm_table and the 'ID' column in the annot_table
    merged_table = gse_table.merge(annot_table[['ID', 'Gene Symbol']], left_on='ID_REF', right_on='ID')


except ZeroDivisionError:
    print("Faced some problem in getting the GPL annotation table. Using the default annotation table of GPL570")
    gpl_id = "GPL570"
    gpl = GEOparse.get_GEO(geo=gpl_id, destdir='./geodata')

    # make the annotation data
    annot_table = gpl.table

    # merge the two tables based on the 'ID_REF' column in the gsm_table and the 'ID' column in the annot_table
    merged_table = gse_table.merge(annot_table[['ID', 'Gene Symbol']], left_on='ID_REF', right_on='ID')
    

# move last 2 columns to the front.
cols = merged_table.columns.tolist()
cols = cols[-2:] + cols[:-2]
merged_table = merged_table[cols]

# replace nan and none values with 0.
merged_table = pd.DataFrame(merged_table)
merged_table.fillna(0, inplace=True)


# save the merged table as a new CSV file
temp = "./geodata/" + input_gseid + "_gene_expression.csv"
merged_table.to_csv(temp, index=False)


genes = merged_table[['ID','Gene Symbol']]
merged_table = merged_table.drop(columns = ['ID'])

# add 1 to each value of the table ,take log2 of each value and finally get mean expression of all probes per gene.
merged_table[merged_table.columns[1:]] = merged_table[merged_table.columns[1:]].add(1)
merged_table.iloc[:, 1:] = np.log2(merged_table.iloc[:, 1:])
merged_table = merged_table.groupby('Gene Symbol').mean()

# reset the index to convert gene name from index to a column
merged_table = merged_table.reset_index()

# save the updated table to another file.
temp = "./geodata/" + input_gseid + "_gene_expression_log2_mean.csv"
merged_table.to_csv(temp, index=False)


# calculate 76GS , KS score, and MLR score.
score1 = EMT76GS(merged_table,input_gseid)
score2 = KSScore(merged_table,input_gseid)
score3 = MLRScore(merged_table,input_gseid,annot_table)

# plotting
plot_76GSvsKS(score1,score2,score3,input_gseid)
