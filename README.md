![ALPACA](https://github.com/borfebor/alpaca_app/blob/62b6915c377ccc9af4bc85aed6e715ea27c583b3/ALPACA_LOGO2.png)

> An app by [B. Ferrero-Bordera](https://www.linkedin.com/in/borjaferrero/)
### Requirements

Required packages for a proper function of the pipeline:

- pandas 1.3.4
- numpy 1.20.3
- bokeh 2.4.2
- matplotlib.pyplot
- sklearn 1.0.1
- seaborn 0.11.2

## Input Data Requirements

- **Protein Groups** file from MaxQuant as `.csv`, `.txt` or `.xlsx`.
- **Quantification Standards** file as `.csv`, `.txt` or `.xlsx`. It requires 3 columns for a proper execution (Accession, MW (kDa) & StdConcentration (µg/µl). See Quantification for more details.
- **Enrichment Standards** file as `.csv`, `.txt` or `.xlsx`. It requires 3 columns for a proper execution (Accession, MW (kDa) & Mix concentration (µg/µl). See Enrichment for more details.

#### Labwork input

Experimental details (in our example `params.txt`) can be added as txt, csv or xlsx formats. This file can include the columns described in the following table:

**Table 2.** Experimental parameters table. This example covers all possible columns. Nonetheless, not all columns are necessary. For example, Enrichment columns (EnrichmentDirection, StdDilution, StdVolume) are only used if any enrichment step was performed. More information about this is described in the Enrichment section.

| Condition   | SampleVolume | ProteinConcentration | AmountMS | CellsPerML | TotalCultureVolume | ProteinSRM | fmolSRM | Enrichment | EnrichmentDirection | StdDilution | StdVolume |
|-------------|--------------|----------------------|----------|------------|--------------------|------------|---------|------------|---------------------|-------------|-----------|
| Cond1_t0    | 2.31         | 2.99                 | 9.67     | 4.54       | 7.54               | TNAMLN     | 4.44    | False      | Down                | 3.96        | 1.22      |
| Cond2_t1    | 2.50         | 0.20                 | 4.10     | 5.13       | 2.62               | AJFVYC     | 4.85    | True       | Down                | 2.43        | 1.51      |
| Cond3_t2    | 7.38         | 6.56                 | 2.77     | 3.66       | 3.80               | BYEKSC     | 9.71    | True       | Down                | 5.71        | 8.53      |

- **Condition**: Condition in which 
- **SampleVolume**: 
- **ProteinConcentration**: Determined protein concentration (µg/µl) in the sample.
- **AmountMS**: Protein amount (µg) injected in the MS.
- **CellsPerML**: Measured cells per mL of culture.
- **TotalCultureVolume**: Total cultivation volume (µL).
- **ProteinSRM** (Optional): If the enrichment of a subcellular fraction has been calculated using targeted proteomics (SRM). This corresponds to the accession of measured protein in SRM to calculate the enrichment.
- **fmolSRM** (Optional): If the enrichment of a subcellular fraction has been calculated using targeted proteomics (SRM). Fmol of the proteins measured in the targeted proteomics measurements. 
- **Enrichment** (Optional): Boolean (True or False). Samples that have been enriched should be specified as True
- **EnrichmentDirection** (Optional): UP or DOWN. 
- **StdDilution** (Optional): This parameter specifies how many times the stock solution of enrichment standards has been diluted before adding it to the sample. If the standards were not diluted before addition, specify 1. Only used when the enrichment is calculated through the function alpaca.gathers() details of the preparation of the used proteins should be added. 
- **StdVolume** (Optional): Volume of enrichment standards (µL) added to the sample. Only used in case the enrichment is calculated through the function alpaca.gathers() details of the preparation of the used proteins should be added.


## Data Importation & Pre-processing
Functions for data import, cleaning and pre-processing.
> alpaca.**eats**(`File`): this function is meant to offer flexibility on the ProteinGroup file importation as some scientist could have the data in a `.txt`, `.csv` or `.xlsx`.

> alpaca.**spits**(`DataFrame`): formater function aims to give a coherence on the imported data, as it could be that MaxQuant output organisation is changed by the user or another software like Perseus. It returns our formated `DataFrame`, and 2 lists: `columns` which contains all df.columns after formating, and `default` that is a list with all suggested columns for dataframe slicing.

## Protein Quantification

### Proteome fraction enrichment (Optional)

In case the study focuses in a fraction of the proteome (e.g., membrane proteome or exoproteome), it is likely that during the sample preparation and enrichment step was performed. This module allows to translate the enrichment step to the data based on how the samples were prepared. 

`Enrichment factors` are calculated based on the fmol quantified in the enriched sample to the raw or non-enriched sample:

$$
ER = \frac{fmol_{enriched}}{fmol_{non-enriched}}
$$

For that purpose, there are 2 strategies that are currently covered under our pipeline:

**1. The quantification of specific proteins of the analysed fraction on both before and after the enrichment step using Targeted MS (SRM).** 

This strategy was described on [Antelo-Varela et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31424929/) and relies on using external protocols (e.g., Skyline) to quantify the enrichment step. Enrichment factors can be added to the parameters table under the column `Enrichment_Factor`. Additionally, the SRM quantified amount for a given protein can be added on the columns `ProteinSRM` (Accession of the quantified protein) and `fmolSRM` (Quantified fmol in the analysed proteome fraction).

**2. The addition of whole proteins at known concentration before performing the enrichment step.**

This approach was described on [Ferrero-Bordera et al. 2024](https://doi.org/10.1128/spectrum.02616-23) and requires of a protein mixture at known concentration added before the enrichment step. Used standards have to be formatted as specified in the table below:

**Table 3.** Enrichment standards

| Accession | MW (kDa) | StdConcentration (µg/µl) |
|-----------|---------:|-------------------------:|
| P02768    |     10.1 |                     2.5  |
| Q9Y6K9    |     65.8 |                     0.8  |
| P05067    |     32.5 |                     1.2  |
| O75475    |     48.2 |                     3.0  |
| Q00653    |     20.9 |                     2.0  |

> alpaca.**gathers**(): 

## Data Integration

This module connects the protein amounts quantified in the sample and the sample preparation. Thus, allowing to calculate protein amounts to the original state (e.g. bacterial culture, raw culture supernatant). This step brings deeper insights to the user based on the known experimental parameters, yielding high valuable data (e.g., molecules per cell, fmol / µmol of protein extract)

> alpaca.**wool**():
