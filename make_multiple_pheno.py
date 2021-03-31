import csv
import pandas as pd

df = pd.read_csv("Soma_inf_case_EUR_pheno.csv")
protein = list(df.columns.values)

outpath='/scratch/richards/sirui.zhou/GWAS_proteome/Somalogic/'
for protein in df:
	x=pd.DataFrame(df, columns=["GWAS_ID","GWAS_ID",protein])
	outfile=outpath + protein + '_inf_EUR' + '.csv'
	x.to_csv(outfile, index=False, header=False, sep=' ')