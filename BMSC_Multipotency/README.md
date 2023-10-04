## BMSC_Multipotency Section
In this section, we present the entire pipeline that illustrates the process of obtaining the results from the single-nuclei multiomics obtained from Day0, Day1, and Day3 (Figure 4, 5).

### Individual Data Analysis (Results Not Shown in the Article)
In this step, each individual dataset was analyzed separately. The analysis served the purpose of excluding low-quality cells and was complemented with clustering analysis to identify key populations.

### Merge All Datasets (Figure 4B. 4C, 4D)
The scripts, All_Days_ARC.py and All_Days_Signac_Seurat.R, encompass the entire codebase essential for merging individual datasets and generating conclusive clustering results. In addition, the All_Days_Signac_Seurat.R script provides a detailed pipeline explaining how to conduct Motif analysis and assess RNA levels based on ATAC.

### Advanced Analysis (Figure 4E, 4F, 4G, 5D)
This section also includes the code for performing advanced analyses, which encompass scVelo, StemFinder, Tricycle, and GRN. These analyses contribute to a comprehensive exploration of the data.





