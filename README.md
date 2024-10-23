# NanoString nCounter data analysis
The aim of this programme is to process the gene expression data from NanoString nCounter platform. The analysis is mainly initiated by R Markdown file "nano_string_ncounter_process", with the assistant of self-defined functions in "nanostring_RUVSeq_functions".

The analysis pipeline requires:
1. RCC files (Gene expression reads from NanoString nCounter platform)
2. Phenotype data(Document all the information of the samples)

The analysis pipeline is consist of 7 steps:
1. packages and function import
2. data import 
3. data preprocessing
4. data quality check
5. normalisation
6. downstream data analysis
7. visualisation

#To utilise this analysis pipeline, you have to make sure:
1. The pipelines is initiated under R ver. 4.0.0
2. Change the directory information in the pipline