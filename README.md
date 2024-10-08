# Data-Analysis-on-Senescence-with-Focus-on-Role-of-Mitochondria

this contains the code written for the MSc Bioinformatics course - Data Exploration and Interpretation for Bioinformatics. grade: A3:20

Used R to answer the questions: What happens to the fibroblast transcriptome when proliferating cells become senescent? What then happens to the fibroblast transcriptome when senescent cells have their mitochondria artificially depleted?

Input Required - a Transcriptomics dataset, in particular:
1) Expression Table with sample names as column names and gene IDs as row names
2) Sample Metadata
3) Annotations table
4) Differential Tables between all three conditions with log2foldchange, p-value, p adjusted (the three conditions being proliferating, senescent, and senescent with mitochondria depletion)

Note on the biology - a bulk RNA-seq dataset comparing proliferating cells (prolif) to replicative senescent cells (senes), and replicative senescent cells which have had their mitochondria depleted (senes_MtD). The cells are human IMR90 fibroblasts. There are three replicates of each group, and so nine samples in total. 

What the Code does:
  Uses these files to do a bunch of calculations, and produces plots. A figures file is attached to demonstrate the kinds of plots produced.

note - unfortunately the input data for this analysis is unpublished, and hence cannot be shared.

Note - this analysis was tailored to this particular dataset with the intent of investigating the above biological questions and hence can't be applied to any generic transcriptomics dataset. that said, the functions found in the script may be reused as is for other analyses.
