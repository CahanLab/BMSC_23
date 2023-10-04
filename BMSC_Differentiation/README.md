## BMSC_Differentiation Section
In this section, we present the entire pipeline that illustrates the process of obtaining the results from the single-nuclei multiomics obtained from trilineage differentiation (Figure 2, 3).

### Individual Data Analysis (Results Not Shown in the Article)
In this step, each individual dataset was analyzed separately. The analysis served the purpose of excluding low-quality cells and was complemented with clustering analysis to identify key populations.

### Merge All Datasets (Figure 2B. 2C)
The scripts, AD+OS+CH_ARC.py and All_Days_Signac_Seurat.R, encompass the entire codebase essential for merging individual datasets and generating conclusive clustering results. In addition, the AD+OS+CH_Signac_Seurat.R script provides a detailed pipeline explaining how to conduct Motif analysis and assess RNA levels based on ATAC.

### Advanced Analysis (Figure 2D, 2E, 3)
This section also includes the code for performing advanced analyses, which encompass scVelo and GRN. These analyses contribute to a comprehensive exploration of the data.





