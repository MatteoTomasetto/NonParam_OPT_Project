# Obstetrics-and-Periodontal-Therapy-Analysis

Several observational studies on pregnant women have suggested an association
between pre-term births and low birthweight with respect to periodontal
disease. 

The dataset is composed by 823 pregnant women that are affected by periodontal
disease; they are subdivided into two groups based on the severity
of their oral condition: 413 women are assigned to the treatment group and
410 to the control one. The patients in the treatment group received periodontal
treatment, oral hygiene instruction, and tooth polishing at their
follow-ups, while those assigned to control underwent only brief oral exams.
Regardless of assignment, all participants attend 4 follow-up visits where
periodontal data are collected.

The main goal is to assess the association between pre-term births and periodontal disease through nonparametric statistics tools in order to understand if treatment of
maternal periodontal disease can reduce risks of pre-term birth and low
birthweight.

## Data

Data are available at https://higgi13425.github.io/medicaldata/reference/opt.html

## Code

- `OPT.R` contains the code to import the dataset

- `DATA PREPROCESSING.R` performs data preprocessing

- `OUTLIERS.R` contains the code to detect outliers though non parametric tools

- `TEST.R` contains the code to find risk factors and to assess the differences between the treatment and control groups via non parametric permutational tests

- `BOOSTRAP.R` performs a bootstrap analysis for the mean and median of gestional age and birthweight in control and treatment groups 

- `CONFORMAL PREDICTION.R` contains the code to get a confidence interval for the mean and median of gestional age and birthweight in control and treatment groups 

- `REGRESSION.R` contains the code to model the probability of pre-term birth and low birthweight through some covaraites via non parametric B-Splines logistic regression

- `SURVIVAL ANALYSIS.R` performs survival analysis for gestional age via non parametric tools such as Kaplan-Meier estimator, log-rank test and Cox model

- `PLOTS.R` contains the code to get some plots of interest for data exploration and analysis

## Authors
* [Matteo Tomasetto](https://github.com/MatteoTomasetto)
* [Laura Gamba](https://github.com/lauragamba)
* [Asia Salpa](https://github.com/asiasalpa)
