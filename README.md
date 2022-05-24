# Y4_synoptic_project_ambientDNA
## optimized ambient DNA profile identification tools (scRNA-seq)
scRNA-Seq is widely used to study the heterogeneity of cells in tissues. Droplet-based approach is one of most popular scRNA-Seq methods. However, the contamination from ambient RNAs has been discovered, which leads to inaccurate estimate of transcriptomes. Although several tools/methods have been developed to correct the data, the method of identifying the genes encoding ambient RNAs is not present. In this project, you will try to develop a new in-house method/algorithm to identify the top ambient RNA-coding genes from scRNA-Seq dataset without spike-in control. A scRNA- Seq dataset as well as a list of contamination-causing genes will be provided. You can use any statistic method and algorithms you have learned in your undergraduate study. The principle, the software/algorithm and the results from provided dataset should be presented.

## synoptic_project_data
data of raw and filtered matrices was given
Data was provided by chaochen wang. The two matrices were derived from GSM4228187: mIslets_I_DMSO; Mus musculus; RNA-Seq (https://www.ncbi.nlm.nih.gov/sra/SRX7426316[accn]) after analyzing with cell ranger.

## synoptic_project_Rfiles
### report_soupx_sole.R
run the original soupx (https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)
### report_soupx_dropletfilter.R
run optimized code
