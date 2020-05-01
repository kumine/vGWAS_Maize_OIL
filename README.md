# vGWAS_Maize_OIL

the source code for [Genetic variants and underlying mechanisms influencing variance heterogeneity in maize](https://onlinelibrary.wiley.com/doi/10.1111/tpj.14786)


vGWAS analysis

The vGWAS was performed using a two-step approach. In the first step, to correct for population stratification, the trait was fit in a linear mixed model with kinship matrix, which calculated by the polygenic  function  in  the  R-package  GenABEL to  get  Grammar  + residuals.  In  the  second  step,  the  variance-heterogeneity  between  the  Grammar  +  residuals  and SNPs  were  tested  using  Levene’s  test.  And  the  Levene’s  test  is  based  on  an  ANOVA  of  the absolute deviation from the median and detailed information is described in previous studies 
