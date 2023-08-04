# MPIM_Intern


##Information about the latest used datasets with variables T,RH,SH over full available time range: 
s. ./Spec_Hum.pptx p.2-3

##Project part 1: Analyse RH, T:
Evaluation in ./eval_rh, different visualizations in *.ipynb notebooks:
*OneSet_Eval* contains basic analysis (horizontal maps, profiles, ...)
*SurfaceVal* contains a comparison of calculated lowest-level pressure values and 2D-Values provided by 2D model output
*Trendanalysis* plots trends and fits in different heights and longitudes
*VerticalAnalysis* plots vertical distributions of RH and T in different variations
Basic helper functions in *helper.py, especially dataload* contains dataload-,renaming-, and elevation-extracting-routines
Uses Datasets
*MERRA2*3D*1980-2023*T,RH,SP (full T, RH isobaric data)
*JRA*isobaric*1959-2021*T,RH,SP (full T, RH isobaric data)
* (MERRA2*2D (full T, RH surface data))
*JRA*2D (full T, RH surface data, for surface value extraction)
*Topo* (Topographie)

##Project part 2: Analyse RH, T, SH:
Evaluation in ./eval_sh, different visualizations in *.ipynb notebooks:
*Compre_Spfh_Rh.ipynb contains basic analysis (profiles, ...)
*VerticalAnalysis* plots vertical distributions of RH and T in different variations
Basic helper functions in *helper.py, especially dataload* contains dataload-,renaming-, and elevation-extracting-routines
Uses Datasets
*MERRA2*3D*1980-2023*all (full T, RH isobaric data)
*JRA*isobaric*1959-2021*sh (full T, RH isobaric data)
* (MERRA2*2D (full T, RH surface data))
*JRA*2D (full T, RH surface data, for surface value extraction)
*Topo* (Topographie)

##Project part 3: Radiosonde ananlysis:
Evaluation in ./Radiosond, different visualizations in *.ipynb notebooks:
*Compre_Spfh_Rh.ipynb contains basic analysis (profiles, ...)
*VerticalAnalysis* plots vertical distributions of RH and T in different variations
Basic helper functions in *helper.py, especially dataload* contains dataload-,renaming-, and elevation-extracting-routines, phiscis* contains calculation attempts from SH 2 RH
Uses Datasets
*MERRA2*3D*1980-2023*all (full T, RH isobaric data)
*JRA*isobaric*1959-2021*sh (full T, RH isobaric data)
* (MERRA2*2D (full T, RH surface data))
*JRA*2D (full T, RH surface data, for surface value extraction)
