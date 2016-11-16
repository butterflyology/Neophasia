**Neophasia** *project code and data for Chris Hamm's part
Basic quesitons:

Are they the same or different? 
Compare sympatric (allochronic) flights
Compare the overall variation of the species with the variation between sympatric flights. 


Using data provided by KB, 12 Type I landmarks for each wing (appears to be right wing), males only. 

Based on discussions with KB I thought that a linear discriminant analysis and discriminant function would be good additions. 

Files in the repo:
** *Neophasia_code.R
** *Neophasia_raw_data.txt *Tab-delimited data
** *README.md 
** *Neophasia_data.R

Population Key and sample sizes:
dp = Donner Pass (univoltine), 23 
ge = Goat Mountain Early (early flight), 40
gl = Goat Mountain Late (late flight), 42
la = Lang (univoltine), 14
me = Mendocino Pass Early (early flight at two flight location), 40
ml = Mendocino Pass Late (late flight at two flight location), 20
or = Oregan (univoltine), 14
wo = Woodfords (univoltine), 29

History of commits:

1. 11 April 2015 - Initial commit of code and data
1. 25 April 2015 - Added data file for Procrustes ANOVA: nwo_43_re_measure.csv
Also, the covariates file is:
Neophasia_wings.csv
1. 27 April 2015 - Checked remeasured data for errors using MANOVA, nothing systematic found. 
1. 3 May 2015 - Updated code to fit new version of Geomorph.

