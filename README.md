## Processing Steps

Create the results directories
```
mkdir -p results/raw
mkdir -p results/batch_corrected
mkdir -p results/harmonized
```

Download the data by using:
```
bash workflow/scripts/download_from_webpage_links.sh
```

## Multiomics & Modeling Steps
To run my current code it's best to work locally with VSCode (Python code) + Rstudio (R code). The scripts should be run in the following way: 

### 1. Preparing-Computable-Matrices.ipynb
```
Open Visual Studio + Jupyter
Install the necessary dependencies (pandas, numpy)
Run notebook
```

### 2. Run-Jive-On-Computable-Matrices.R
```
Open R Studio
Set the working directory to the project dir
Install the necessary dependencies (r.jive, dplyr)
Run code
```

### 3. Model-With-Jive-Factors.ipynb
```
Open Visual Studio + Jupyter
Install the necessary dependencies
Run notebook
```