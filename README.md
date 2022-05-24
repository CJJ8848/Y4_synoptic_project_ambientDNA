# Y4_synoptic_project_ambientDNA
## optimized ambient DNA profile identification tools (scRNA-seq)
scRNA-Seq is widely used to study the heterogeneity of cells in tissues. Droplet-based approach is one of most popular scRNA-Seq methods. However, the contamination from ambient RNAs has been discovered, which leads to inaccurate estimate of transcriptomes. Although several tools/methods have been developed to correct the data, the method of identifying the genes encoding ambient RNAs is not present. In this project, you will try to develop a new in-house method/algorithm to identify the top ambient RNA-coding genes from scRNA-Seq dataset without spike-in control. A scRNA- Seq dataset as well as a list of contamination-causing genes will be provided. You can use any statistic method and algorithms you have learned in your undergraduate study. The principle, the software/algorithm and the results from provided dataset should be presented.

## OS and R version
linux, R 4.1.3
## pre-requirement
R package: 
library(DropletUtils) 1.12.3
library(Seurat) 4.1.0
library(dplyr) 1.0.8
library(patchwork) 1.1.1
library(SoupX) 1.5.2
library(ggplot2) 3.3.5
library(scales) 1.1.1
## synoptic_project_data
data of raw and filtered matrices was given
Data was provided by chaochen wang. The two matrices were derived from GSM4228187: mIslets_I_DMSO; Mus musculus; RNA-Seq (https://www.ncbi.nlm.nih.gov/sra/SRX7426316[accn]) after analyzing with cell ranger.

Due to the file size limitation of GitHub, the full data can not be uploaded.
You can download the two matrix folders from the institute server by:
scp 3180110107bit@10.105.100.153:/public/workspace/3180110107bit/synoptic_project_ambientDNAdecontamination.zip ./to/your/path

with password:111111

## synoptic_project_Rfiles
### report_soupx_sole.R
run the original soupx (https://academic.oup.com/gigascience/article/9/12/giaa151/6049831)
### report_soupx_dropletfilter.R
run optimized code
## example command
Rscript report_soupx_dropletfilter.R path/to/raw path/to/filtered

Rscript report_soupx_dropletfilter.R ../synoptic_project_data/raw_feature_bc_matrix/ ../synoptic_project_data/filtered_feature_bc_matrix/

## output
figures of my report under working dir
txt file containing ambientDNAprofile (ambientDNAprofile.txt)
