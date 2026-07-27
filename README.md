# Quantile regression for proteomic age prediction

This repository provides code and documentation for proteomic age prediction and predictive uncertainty quantification using quantile regression (QR), as described in our paper, "Beyond point estimates: quantifying predictive uncertainty reveals hidden dimensions of biological age acceleration and improves risk interpretation".



## Application of organ aging models trained on UK Biobank Olink proteomics data



### Install dependent packages in R

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)

tidyr [https://CRAN.R-project.org/package=tidyr](https://CRAN.R-project.org/package=tidyr)

stringr [https://CRAN.R-project.org/package=stringr](https://CRAN.R-project.org/package=stringr)

splines2 [https://CRAN.R-project.org/package=splines2](https://CRAN.R-project.org/package=splines2)

fsQRPPA (for model training) [https://anonymous.4open.science/r/fsQRPPA-6764/](https://anonymous.4open.science/r/fsQRPPA-6764/)



### Organ aging models

The organ aging models developed using UK Biobank Olink plasma proteomics data are available in the [model directory](/model). 

The proteins included in the multi-organ and organ-specific models were defined according to the paper, ["Plasma proteomics links brain and immune system aging with healthspan and longevity"](https://www.nature.com/articles/s41591-025-03798-1).



### Example of organ age prediction

A simulated toy dataset is provided in the [example directory](/example) to demonstrate the prediction workflow. 

Users can apply the provided [models](model/qr.beta.proteomic.age.tsv) to predict organ age from the simulated dataset.

The example input data includes 

1) the [chronological age at blood sample collection](example/example.age.tsv) for 1,000 individuals (column "IID").
2) the [Olink plasma protein expression matrix](example/example.proteomics.tsv) containing simulated proteomic profiles for 1,458 proteins (remaining columns) across 1,000 individuals (column "IID").

Run the R script [apply_age_model.R](apply_age_model.R) to perform organ age prediction, predictive uncertainty quantification, and downstream analyses on the simulated dataset.



### Output of organ age prediction

The [script](apply_age_model.R) generates a tab-delimited text file containing organ age predictions and predictive uncertainty measures. 

Each row corresponds to one organ-individual pair.

The output includes the following columns:

- "IID": individual identifier from the simulated dataset.
- "Organ": organ model name.
- "Age.QR0.5": predicted organ age using the corresponding organ model indicated in the "Organ" column.
- "Gap.QR0.5": standardized age gap (measure of accelerated or decelerated aging) for the corresponding organ model.
- "Tail.Prob": tail probability (uncertainty quantification) for the corresponding organ model.
- "Interval.Length": length of prediction interval (uncertainty quantification) for the corresponding organ model.



## Model training with user's own Olink proteomics data

Users can train their own models using the [fsQRPPA](https://arxiv.org/abs/2601.02826) method on a Olink proteomics dataset. The fsQRPPA R package, available at [https://anonymous.4open.science/r/fsQRPPA-6764/](https://anonymous.4open.science/r/fsQRPPA-6764/), is required for model training.

An example, [train_age_model.R](train_age_model.R), is provided to demonstrate the model develpment with a [simulated toy dataset](/example).

Although the provided code allows users to train organ aging models using their own proteomics data, the models rely on the predefined sets of organ-enriched proteins. These protein sets were derived from approximately 3,000 proteins measured by the Olink Explore 3072 platform, and were originally established in the study ["Plasma proteomics links brain and immune system aging with healthspan and longevity"](https://www.nature.com/articles/s41591-025-03798-1). Therefore, to achieve optimal prediction performance, users' datasets should cover as many of these organ-enriched proteins as possible. Limited overlap between the available proteins and the predefined organ-enriched protein sets may reduce model accuracy.



## Citation

TBA.


