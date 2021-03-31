
#!/bin/bash
#PBS -N Somalogic_AFR
#PBS -o logs/
#PBS -e logs/
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=10G
#PBS -l vmem=10G
#PBS -l epilogue=/scratch/richards/sirui.zhou/GWAS_proteome/Afr_epilogue.sh

cd /scratch/richards/sirui.zhou/GWAS_proteome/
protein=($(awk '{ print $1}' Somalogic_list1.txt))
for ((i=0;i<${#protein[@]};++i)); do gcta64 \
  --mbfile BQC_filelist.txt \
  --fastGWA-lr \
  --grm AFR/Soma_inf_covidpos_AFR \
  --pheno Somalogic_AFR/${protein[i]}_inf_AFR.csv \
  --covar Soma_inf_case_AFR_cov.txt \
  --qcovar Soma_inf_case_AFR_qcov.txt \
  --maf 0.05 \
  --thread-num 2 \
  --out Somalogic_res_AFR/${protein[i]}_inf_AFR_fastGWA_lr; done
