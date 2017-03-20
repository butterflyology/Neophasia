**Neophasia** project code and data for Chris Hamm's part
Basic questions:

Are they the same or different?
Compare sympatric (allochronic) flights
Compare the overall variation of the species with the variation between sympatric flights.

Using data provided by KB, 12 Type I landmarks for each wing (appears to be right wing), males only.

Based on discussions with KB I thought that a linear discriminant analysis and discriminant function would be good additions.

Files in the repo:

1. `.gitignore`
1. `LDA-ellipse.R` - Code to plot ellipses by factor around a cloud of points.
1. `Neophasia_code.R` - Code needed to reproduce our results.  
1. `Neophasia_data.RData` - The output and results saved as a `.RData` object.
1. `README.md`

  `/Data`:
    * `Neophasia_raw_data.txt` - Raw geometric morphometric data.
    * `Neophasia_wings.csv` - Raw covariate data.
    * `nwo_43_re_measure.csv` - Geometric morphometric data for MANOVA to examine systematic error.


Population Key and sample sizes:

* dp = Donner Pass (univoltine), 23
* ge = Goat Mountain Early (early flight), 40
* gl = Goat Mountain Late (late flight), 42
* la = Lang (univoltine), 14
* me = Mendocino Pass Early (early flight at two flight location), 40
* ml = Mendocino Pass Late (late flight at two flight location), 20
* or = Oregan (univoltine), 14
* wo = Woodfords (univoltine), 29


History of commits:

1. 2015-04-11 - Initial commit of code and data.
1. 2015-04-25 - Added data file for Procrustes ANOVA.
1. 2015-04-27 - Checked remeasured data for errors using MANOVA, nothing systematic found.
1. 2015-05-3 - Updated code to fit new version of Geomorph.
1. 2017-03-20 - Updated code and `README.md` and prepared files for ZENODO DOI. 
