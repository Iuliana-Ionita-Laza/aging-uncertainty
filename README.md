# Quantile regression for proteomic age prediction

This repository provides code and documentation for proteomic age prediction and predictive uncertainty quantification using quantile regression (QR), as described in our paper, "Beyond point estimates: quantifying predictive uncertainty reveals hidden dimensions of biological age acceleration and improves risk interpretation".



## Install dependent packages in R

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)

tidyr [https://CRAN.R-project.org/package=tidyr](https://CRAN.R-project.org/package=tidyr)

stringr [https://CRAN.R-project.org/package=stringr](https://CRAN.R-project.org/package=stringr)

splines2 [https://CRAN.R-project.org/package=splines2](https://CRAN.R-project.org/package=splines2)

fsQRPPA (for model training) [https://anonymous.4open.science/r/fsQRPPA-6764/](https://anonymous.4open.science/r/fsQRPPA-6764/)



## Organ aging models

The organ aging models developed using UK Biobank Olink plasma proteomics data are available in the [model directory](/model). 

The proteins included in the multi-organ and organ-specific models were defined according to the paper, ["Plasma proteomics links brain and immune system aging with healthspan and longevity"](https://www.nature.com/articles/s41591-025-03798-1).



## Example of organ age prediction

A simulated toy dataset is provided in the [example directory](/example) to demonstrate the prediction workflow. 

Users can apply the provided [models](model/qr.beta.proteomic.age.tsv) to predict organ age from the simulated dataset.

The example input data includes 

1) the [chronological age at blood sample collection](example/example.age.tsv) for 1,000 individuals (column "IID").
2) the [Olink plasma protein expression matrix](example/example.proteomics.tsv) containing simulated proteomic profiles for 1,458 proteins (remaining columns) across 1,000 individuals (column "IID").

Run the R script [apply_age_model.R](apply_age_model.R) to perform organ age prediction and predictive uncertainty quantification on the simulated dataset.



## Output of organ age prediction

The [script](apply_age_model.R) generates a tab-delimited text file containing organ age predictions and predictive uncertainty measures. 

Each row corresponds to one organ-individual pair.

The output includes the following columns:

- "IID": individual identifier from the simulated dataset.
- "Organ": organ model name.
- "Age.QR0.5": predicted organ age using the corresponding organ model indicated in the "Organ" column.
- "Gap.QR0.5": standardized age gap (measure of accelerated or decelerated aging) for the corresponding organ model.
- "Prob": tail probability (uncertainty quantification) for the corresponding organ model.
- "Interval.Length": length of prediction interval (uncertainty quantification) for the corresponding organ model.



## Model training

Use fsQRPPA [https://anonymous.4open.science/r/fsQRPPA-6764/](https://anonymous.4open.science/r/fsQRPPA-6764/) and simulated training dataset.

TBA.



## Citation

TBA.


