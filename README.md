# PheWIS
Phenowide-wide interaction study (PheWIS) is an R computational tool to evaluate the association between phenotypes and gene-environmental interactions.

## Data preparation
Data resource: UK Biobank.
### Genotype data

1. bgen files are in /project/kylab/lab\_shared/UKB/imputationbgen\_v1.2_UKBsource.
2. Run scripts 1.convert\_bgen1.2\_to\_pgen.sh and 2.filter\_imputed_snps.sh
3. Use R to do the QC of participants by using script UKB\_QC-07082020.R and data file ukb48818\_rel\_s488282_output.dat
4. Run script 3.convert\_pfile_bgenUKB.sh to get the bgen files that will be used to generate the genotype data.
5. Run scripts 4.generate\_geno.sh and generate_geno.R to modify the genotype data into a specific format of PheWIS.


### Phenotype data

  1. Use software "createUKBphenome" to transfer the phenotype data into phecode. "createUKBphenome" is an R software that used to transfer data from UK Biobank into phecode.
  2. The path of file ukb34137.tab is /project/kylab/lab\_shared/UKB/pheno. If the others tables will be used, please make sure the table contains the following data: f.41270, f.40002, f.40006, f.40001, f.40013, f.41271, f.41201, f.41202, f.41203, f.41204, f.41205. These data are ICD9 and ICD10 used to generate phecode.

```bash
cp ukb34137.tab phew.tab
git clone https://github.com/umich-cphds/createUKBphenome
git clone https://github.com/PheWAS/PheWAS
Rscript ./scripts/function.createUKBphenome.r
```
The phenotpye\_data is in the results directory and the name should follow the pattern "UKB\_PHECODE_date.txt"

### Exposure
1. Find the field code of your exposure on the website of UKBthen and extract exposure information from the UKB table. For example, the field ID of sex is f.31. The following commands are using to extract the column number of sex.
```bash
head -n 1 ukb34137.tab | awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}' | grep -n "f.31"
```

2. Use the following commands to write the columns to a new file. For example, sex is in the second column of the ukb34137.tab.
```bash
cut -f 1,2 ukb34137.tab > sex.txt
```
### Covariates
The process of covariates data generation is same with exposure data generation. Please see exposure for details.

## PheWIS
### PheWIS results

use script dophewis.R.

### Manhattan plots
use script Manhattan.R

contact: Mengyuan Zhang, mz59443@uga.edu
