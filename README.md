# ExpMetFilter
An R-package designed to remove the exposure compound, its' metabolites and other impurities from non-targeted data analysis of toxicological exposure studies. Useful for metabolomics multivariate analysis but could also be used purely for removing blank features not caught by other data processing methods.

### Installation
```
install.packages("devtools")
library("devtools")
install_github("parasitetwin/ExpMetFilter")
```

## Functionality of package

"ExpMetFilter"
- A function containing 5 different filter options for non-targeted data as well as the possiblity of internal standard correction of evaporation. The user controls which filters they want to apply to their dataset through the meta-data in the input file which can be found in "Example files".

"GapFilter"
- A function which utilizes "Gap Status", gathered by Compound Discoverer (Thermo, USA), to filter out features which are not found in a specific sample type (as specified by the user).

## Example files 

"InputFileDescription"
- Files which in-depth explains how datasets should be structured in order to be utilized by the two functions of this package.

"ExpMetFilter-ExampleFile"
- A file containing two datasets (Orbitrap data from both HILIC and C18 chromatography analysis) of non-exposed embryos which have both undergone "GapFilter" and batchCorr processing.
