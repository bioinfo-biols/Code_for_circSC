# Exploring the cellular landscape of circular RNAs using full-length single-cell RNA sequencing

This repository contains the analysis pipeline for our circSC manuscript

- Step1.gene_count.py: gene quantification using HISAT2 and StringTie and pipeline

- Step2.circRNA_detection.py: circRNA identification and quantification using bwa, CIRI2 and CIRI_AS pipeline

- Step3.cell_cluster.R: scRNA-seq data integration and cell-type annotation

- Step4.integrated_data.py: generate circRNA BSJ matrix for all samples
