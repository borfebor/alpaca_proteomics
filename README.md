![ALPACA](https://github.com/borfebor/alpaca_app/blob/62b6915c377ccc9af4bc85aed6e715ea27c583b3/ALPACA_LOGO2.png)

> An app by [B. Ferrero-Bordera](https://www.linkedin.com/in/borjaferrero/)

![Alpaca pipeline modules](https://github.com/borfebor/alpaca_proteomics/blob/main/Tutorial/alpaca_pipeline.png)

# Table of Contents
- [Documentation](#Documentation)
   - [Getting started](##Getting-started)
   - [GUI version](##GUI-version)
- [Usage Example](#Usage-Example)
   - [Step-by-step protocols](##Step-by-step-rotocols)
   - [Example datasets](##Example-datasets)
- [Library features](#Library-Features)
   - [Input data requirements](##Input-Data-Requirements)
   - [Data importation & pre-processing](##Data-Importation-&-Pre-processing)
   - [Protein quantification](##Protein-Quantification)
   - [Integration of sample preparation details ](##Data-Integration)


# Documentation

📘 [Documentation and Tutorial](https://borfebor.github.io/alpaca_proteomics/)

## Getting started

```Python
pip install alpaca-proteomics
```

```Python
from alpaca_proteomics import alpaca
```

## GUI version

Check our GUI version on the [🦙 **alpaca app**](https://alpaca.nube.uni-greifswald.de/) or run it locally using [Docker](www.docker.com/get-started):

1. Ensure Docker is Installed and Running
Make sure you have Docker installed on your system and that the Docker daemon is running.

2. Clone and Build the Application

```bash
# Clone the repository
git clone https://github.com/borfebor/alpaca_app.git

# Navigate to the folder containing the cloned repository
cd alpaca_app

# Build the Alpaca app (the -t flag specifies the name of the Docker image)
docker build -t alpaca-app . 

# Start the Alpaca app from the terminal
docker run -p 8501:8501 alpaca-app
```

3. Access the Web Interface:
Once the container is running, open the following link in your browser:
[http://localhost:8501/](http://localhost:8501/)


## Cite us

Stay tuned, the paper is submitted.

# Usage Example

A data analysis using Jupyter Notebooks is described [here](https://github.com/borfebor/alpaca_proteomics/blob/main/Tutorial/Tutorial.ipynb) using the dataset published in [Ferrero-Bordera et al. 2024. Microbiology Spectrum](https://doi.org/10.1128/spectrum.02616-23)

## Step-by-step protocols

1. [Pre-processing Proteomics Data with the Alpaca Pipeline](https://www.protocols.io/view/pre-processing-proteomics-data-with-the-alpaca-pip-dm6gp95j5vzp/v1)
2. [Absolute Quantification of Proteome Abundances with the Alpaca Pipeline](https://www.protocols.io/view/absolute-quantification-of-proteome-abundances-wit-5jyl8dx57g2w/v1)
3. [Integration of Sample Preparation and Subcellular Fraction Enrichment in Alpaca](https://www.protocols.io/view/integration-of-sample-preparation-and-subcellular-j8nlk9z46v5r/v1)
   
## Example datasets

Example datasets are available in the following [folder](https://github.com/borfebor/alpaca_app/tree/main/Datasets).
- **Enriched_example.txt**: Exoproteome dataset from [Ferrero-Bordera et al. 2024. Microbiology Spectrum](https://doi.org/10.1128/spectrum.02616-23)
- **Membrane_example.txt**: Membrane proteome dataset from [Antelo-Varela et al. 2019. Anal. Chem.](https://doi.org/10.1021/acs.analchem.9b02869)

# Library-Features

### Requirements

Required packages and last tested working versions:

- pandas 2.1.4
- numpy 1.24.3
- matplotlib 3.8.1
- scipy 1.11.3
- scikit-learn 1.4.1.post1
- sklearn 1.0.1
- seaborn 0.11.2
- thefuzz 0.20.0

## Input Data Requirements

- **Protein Groups** file containing the protein Groups and the quantified intensities as `.csv`, `.txt` or `.xlsx`.
- **Quantification Standards** file as `.csv`, `.txt` or `.xlsx`. It requires 3 columns for a proper execution (Accession, MW (kDa) & Amount (fmol). See Quantification for more details.
- **Enrichment Standards** (Optional) file as `.csv`, `.txt` or `.xlsx`. It requires 3 columns for a proper execution (Accession, MW (kDa) & StdConcentration (µg/µl). See Enrichment for more details.

#### Labwork input

Experimental details (in our example `params.txt`) can be added as txt, csv or xlsx formats. This file can include the columns described in the following table:

**Table 1.** Experimental parameters table. This example covers all possible columns. Nonetheless, not all columns are necessary. For example, Enrichment columns (EnrichmentDirection, StdDilution, StdVolume) are only used if any enrichment step was performed. More information about this is described in the Enrichment section.

| Condition   | SampleVolume | ProteinConcentration | AmountMS | CellsPerML | TotalCultureVolume | ProteinSRM | fmolSRM | Enrichment | EnrichmentMode | StdDilution | StdVolume |
|-------------|--------------|----------------------|----------|------------|--------------------|------------|---------|------------|---------------------|-------------|-----------|
| Cond1_t0    | 2.31         | 2.99                 | 9.67     | 4.54       | 7.54               | TNAMLN     | 4.44    | False      |                 | 3.96        | 1.22      |
| Cond2_t1    | 2.50         | 0.20                 | 4.10     | 5.13       | 2.62               | AJFVYC     | 4.85    | True       | Enrichment                | 2.43        | 1.51      |
| Cond3_t2    | 7.38         | 6.56                 | 2.77     | 3.66       | 3.80               | BYEKSC     | 9.71    | True       | Concentration                | 5.71        | 8.53      |

- **Condition**: Identifier for the condition or timepoint in which the parameters were applied.
- **SampleVolume**: Protein extract volume (µL) used for protein digestion.
- **ProteinConcentration**: Determined protein concentration (µg/µl) in the sample.
- **AmountMS**: Sample amount (in µg) injected into the mass spectrometer.
- **CellsPerML**: Determined number of cells per mL in the original culture.
- **TotalCultureVolume**: Total cultivation volume (µL).
- **ProteinSRM** (Optional): If the enrichment of a subcellular fraction has been calculated using targeted proteomics (SRM). This corresponds to the accession of measured protein in SRM to calculate the enrichment.
- **fmolSRM** (Optional): If the enrichment of a subcellular fraction has been calculated using targeted proteomics (SRM). Fmol of the proteins measured in the targeted proteomics measurements. 
- **Enrichment** `Optional`: Boolean (True or False). Samples that have been enriched should be specified as True.
- **EnrichmentMode** `Optional`: Enrichment or Concentration (see Supplementary Material, [integration protocol](https://www.protocols.io/view/integration-of-sample-preparation-and-subcellular-j8nlk9z46v5r/v1), or the [documentation](https://borfebor.github.io/alpaca_proteomics/) for more details).
- **StdDilution** (Optional): This parameter specifies how many times the stock solution of enrichment standards has been diluted before adding it to the sample. If the standards were not diluted before addition, specify 1. Only used when the enrichment is calculated through the function alpaca.gathers() details of the preparation of the used proteins should be added. 
- **StdVolume** (Optional): Volume of enrichment standards (µL) added to the sample. Only used in case the enrichment is calculated through the function alpaca.gathers() details of the preparation of the used proteins should be added.


## Data Importation & Pre-processing
Functions for data import, cleaning and pre-processing.
> alpaca.**eats**(`File`): this function is meant to offer flexibility on the ProteinGroup file importation as some scientists could have the data in a `.txt`, `.csv` or `.xlsx`.

> alpaca.**spits**(`DataFrame`): this function aims to give coherence to the imported data, as it could be that MaxQuant output organisation is changed by the user or another software like Perseus. It returns our formatted `DataFrame`, and 2 lists: `columns` which contain all df.columns after formatting, and `default` which is a list with all suggested columns for dataframe slicing.

## Protein Quantification

Absolute quantification using Alpaca is optimised for label-free methods, relying on the addition of a set of anchor proteins at a known amount. 

**Table 2.** Format for the file describing the stock solution of anchor proteins.

| Accession  | MW (kDa) | Amount (fmol) |
|------------|---------:|--------------:|
| P02768     |   10.1   |         50    |
| Q9Y6K9     |   65.8   |        100    |
| P05067     |   32.5   |         25    |
| O75475     |   48.2   |         75    |
| Q00653     |   20.9   |         30    |

> alpaca.**census**(): calculates the abundances in moles based on added anchor proteins (Table 2). Anchor proteins should be specified in a dataframe using the format described in Table 2.

### Proteome fraction enrichment (Optional)

In case the study focuses on a fraction of the proteome (e.g., membrane proteome or exoproteome), it is likely that during the sample preparation and enrichment step was performed. This module allows to connect the enrichment step to the data based on how the samples were prepared. 

`Enrichment factors` are calculated based on the fmol quantified in the enriched sample to the raw or non-enriched sample:

$$
ER = \frac{fmol_{enriched}}{fmol_{non-enriched}}
$$

For that purpose, 2 strategies are currently covered under our pipeline:

**1. The quantification of specific proteins of the analysed fraction on both before and after the enrichment step using Targeted MS (SRM).** 

This strategy was described by [Antelo-Varela et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31424929/) and relies on using external protocols (e.g., Skyline) to quantify the enrichment step. Enrichment factors can be added to the parameters table under the column `Enrichment_Factor`. Additionally, the SRM quantified amount for a given protein can be added on the columns `ProteinSRM` (Accession of the quantified protein) and `fmolSRM` (Quantified fmol in the analysed proteome fraction).

**2. The addition of whole proteins at a known concentration before performing the enrichment step.**

This approach was described by [Ferrero-Bordera et al. 2024](https://doi.org/10.1128/spectrum.02616-23) and requires a protein mixture at a known concentration added before the enrichment step. Used standards have to be formatted as specified in the table below:

**Table 3.** Enrichment standards

| Accession | MW (kDa) | StdConcentration (µg/µl) |
|-----------|---------:|-------------------------:|
| P02768    |     10.1 |                     2.5  |
| Q9Y6K9    |     65.8 |                     0.8  |
| P05067    |     32.5 |                     1.2  |
| O75475    |     48.2 |                     3.0  |
| Q00653    |     20.9 |                     2.0  |

> alpaca.**gathers**(): calculates the enrichment factors for each specified fraction based on the sample preparation (Table 1) and the added standards to the sample. Standards should be specified in a dataframe using the format described in Table 3.

## Data Integration

This module connects the protein amounts quantified in the sample and the sample preparation. Thus, allowing the calculation of protein amounts to the original state (e.g. bacterial culture, raw culture supernatant). This step brings deeper insights to the user based on the known experimental parameters, yielding highly valuable data (e.g., molecules per cell, fmol / µmol of protein extract)

> alpaca.**wool**(): This function integrates the sample preparation details with the quantified proteins injected in the MS.
