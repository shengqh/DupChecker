dupchecker
==========
* [Introduction](#Introduction)
* [Download and install](#download)
* [Example](#example)
* [Usage](#usage)
* [Report](#report)
* [Reproduce the report](#reproduce)

<a name="Introduction"/>
# Introduction #
Meta-analysis has become a popular approach for high-throughput genomic data analysis because it often can signifi-cantly increase power to detect biological signals or patterns in datasets. However, when using public-available databases for meta-analysis, duplication of samples is an often encoun-tered problem, especially for gene expression data. Not remov-ing duplicates would make study results questionable. We de-veloped a Bioconductor package Dupchecker that efficiently identifies duplicated samples by generating MD5 fingerprints for raw data.

<a name="download"/>
# Download and install #
You can install DupChecker package in R from [github](https://github.com/shengqh/dupchecker/) by following codes:

	library(devtools)
	install_github("DupChecker", user="shengqh")

