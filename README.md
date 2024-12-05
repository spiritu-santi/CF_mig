# Retrieving cloud forest elevation shifts with occurrence data
Files needed to run the results from:

Ramírez-Barahona S *et al.*, 2024. Upslope plant species shifts in Mesoamerican cloud forests driven by climate and land use change. *Under review* 

- The GBIF raw occurrence dataset used to produce the final results can be directly accessed here (https://doi.org/10.15468/dl.544ewq) and should be properly cited.  The GBIF raw occurrence dataset is also available directly in the Zenodo repository.  The use of this dataset is licensed under CC BY-NC 4.0 (Attribution-NonCommercial).  
    - Prior to processing the data, we reduced file size by keeping only essential information using bash:  
      - spiritusanti$ unzip 0188599-210914110416597.zip  
      - spiritusanti$ awk -F"\t" '{print $6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$13"\t"$22"\t"$23"\t"$33}' 0188599-210914110416597.csv > Traqueos_NeoTropics_COR.csv  
      - spiritusanti$ wc -l TRAQUEOS_raw.csv  
      - spiritusanti$ wc -l Traqueos_NeoTropics_COR.csv  
    - Also available in the Zenodo repository is the processed occurrence dataset as an R-binary file: 'Historical_data_FULL_LUL.v.1.R'.  

- NOTE: we urge the user to dowload and cite the climate or land-use data from the original sources.  
- NOTE: If using the processed data used here, we strongly advice to cite the original sources.
    - Climate series: 
      - D. N. Karger, O. Conrad, J. Böhner, T. Kawohl, H. Kreft, R. W. Soria-Auza, N. E. Zimmermann, H. P. Linder, M. Kessler, Climatologies at high resolution for the earth’s land surface areas. Sci. Data. 4, 170122 (2017).
    - Land-use series: 
       - H. Liu, P. Gong, J. Wang, N. Clinton, Y. Bai, S. Liang, Annual dynamics of global land cover and its long-term changes from 1982 to 2015. Earth Syst. Sci. Data. 12, 1217–1243 (2020).
       - D. García-Álvarez, J. Lara-Hinojosa, F. J. Jurado-Pérez, J. Quintero-Villaraso, in Land Use Cover Datasets and Validation Tools, D. García-álvarez, M. Teresa, C. Olmedo, Eds. (Springer, 2022).

[![DOI](https://zenodo.org/badge/402174680.svg)](https://zenodo.org/badge/latestdoi/402174680)
