### Getting Started

1. Install alpaca package on the terminal through 


```
pip install alpaca-proteomics
```


2. Import the package



```python
from alpaca_proteomics import alpaca
```

## Data import and formatting

Alpaca works with unprocessed proteomics datasets from MaxQuant. The package takes the file `proteinGroups.txt`, which can be found in the combined folder of MaxQuant output. This dataset is from the exoproteome data published in [Ferrero-Bordera et al. 2024. Microbiology Spectrum](https://doi.org/10.1128/spectrum.02616-23).


```python
file = 'proteinGroups.txt'

# Data importation
df, id_col, it = alpaca.eats('proteinGroups.txt')
```

    The column Protein IDs was detected to contain your ProteinGroups IDs.
    The following intensity methods were detected in the data: Intensity, iBAQ, LFQ


The function returned:

- **df** is the imported data as a pandas dataframe
- **id_col** corresponds to the column which was detected to contain the Protein IDs
- **it** is a dictionary which groups the columns containing intensity data within each intensity method (e.g. LFQ)

In our example, the data contained 3 intensity methods (Intensity, iBAQ, LFQ)

### Assistance on the analysis (Optional)


```python
standards_file = 'UPS2.xlsx' # Path to the anchor proteins file
st_proteins = alpaca.eats(standards_file) # Importation of the anchor proteins file (more details on these are listed below)
spiked_samples = ['Before_Induction_01', 'Control_01', 'Diamide_01'] # Samples in which anchor proteins were added
values_per_sample = 1/4 # Valid values per condition

suggested = alpaca.Consultant(df,
                         st_proteins,
                         it, 
                         added_samples=spiked_samples,
                         values_per_sample=values_per_sample)
```

    Based on your data, Median-normalized iBAQ is recommended for the quantification.


The function returns an suggested analysis method and a dataframe with the calculated scores for each intensity method with different normalization approaches. These data is visualized below.


```python
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

plt.figure(figsize=(5,5))
sns.set(font_scale=1.5)
sns.heatmap(suggested.pivot(index='Normalization', columns='Intensity method', values='score'),
           annot=True, cmap='viridis', lw=1, cbar_kws={'label': 'Fitting (R$^2$)'})
```




    <Axes: xlabel='Intensity method', ylabel='Normalization'>




    
![png](assets/Tutorial_8_1.png)
    


### Data pre-processing

Based on the suggestions from the function alpaca.Consultant, the analysis can continue with the most suitable parameters.


```python
# Data pre-processing

values_per_sample = 1/4

clean_df = alpaca.spits(df, 
                        lfq_method='iBAQ',
                        formatting=True, 
                        valid_values=values_per_sample,
                        normalization='Median',
                        info_cols=['Accession', 'Gene names'])
clean_df.head(5)
```

    3 experimental conditions were identified in the data: Before_Induction, Control, Diamide
    Items marked on Only identified by site, Reverse, Potential contaminant have been removed from the dataset.
    Dataset formated for further analysis and visualisation.





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Accession</th>
      <th>Protein</th>
      <th>Sample</th>
      <th>iBAQ</th>
      <th>Condition</th>
      <th>Replicate</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>3</th>
      <td>C0SP82</td>
      <td>YoaE</td>
      <td>iBAQ Before_Induction_01</td>
      <td>2.083660</td>
      <td>Before_Induction</td>
      <td>01</td>
    </tr>
    <tr>
      <th>4</th>
      <td>C0SP93</td>
      <td>AccD</td>
      <td>iBAQ Before_Induction_01</td>
      <td>2.312376</td>
      <td>Before_Induction</td>
      <td>01</td>
    </tr>
    <tr>
      <th>5</th>
      <td>C0SP94</td>
      <td>YhfQ</td>
      <td>iBAQ Before_Induction_01</td>
      <td>-2.030633</td>
      <td>Before_Induction</td>
      <td>01</td>
    </tr>
    <tr>
      <th>6</th>
      <td>C0SPA7</td>
      <td>YukB</td>
      <td>iBAQ Before_Induction_01</td>
      <td>-8.627621</td>
      <td>Before_Induction</td>
      <td>01</td>
    </tr>
    <tr>
      <th>7</th>
      <td>C0SPB0</td>
      <td>YtcI</td>
      <td>iBAQ Before_Induction_01</td>
      <td>0.498673</td>
      <td>Before_Induction</td>
      <td>01</td>
    </tr>
  </tbody>
</table>
</div>



## Anchor protein quantification

Absolute quantification using Alpaca is optimised for label-free methods, relying on the addition of a set of anchor proteins at a know amount. 

**Table 1.** Format for the file describing the stock solution of anchor proteins.

| Accession  | MW (kDa) | Amount (fmol) |
|------------|---------:|--------------:|
| P02768     |   10.1   |         50    |
| Q9Y6K9     |   65.8   |        100    |
| P05067     |   32.5   |         25    |
| O75475     |   48.2   |         75    |
| Q00653     |   20.9   |         30    |


```python
# Import the file containing the information about the quantification standards proteins

standards_file = 'UPS2.xlsx'
st_proteins = alpaca.eats(standards_file)

# If applicable, define which samples/replicates contain standards proteins

spiked_samples = ['iBAQ Before_Induction_01', 'iBAQ Control_01', 'iBAQ Diamide_01']

# Quantify the fmol present in the measured samples

quant_df, st_proteins, coef, inter, r2 = alpaca.census(clean_df, st_proteins, lfq_col='iBAQ',
                                                      filter_col = 'Sample', # Defines which column to filter for the spiked samples
                                                      added_samples = spiked_samples) # Adding which samples contain the standars
```


    
![png](assets/Tutorial_12_0.png)
    


`alpaca.census()` adds a column to the processed data with the calculated mol amounts present in the measured samples.

## Experimental details



```python
sample_prep = pd.read_csv('params.csv', sep=',')
sample_prep.sample()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Condition</th>
      <th>SampleVolume</th>
      <th>ProteinConcentration</th>
      <th>AmountMS</th>
      <th>CellsPerML</th>
      <th>TotalCultureVolume</th>
      <th>ProteinSRM</th>
      <th>fmolSRM</th>
      <th>Enrichment</th>
      <th>EnrichmentMode</th>
      <th>StdDilution</th>
      <th>StdVolume</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>Control</td>
      <td>9.310715</td>
      <td>7.766196</td>
      <td>1.060855</td>
      <td>1.658673</td>
      <td>9.223972</td>
      <td>P54576</td>
      <td>6.914357</td>
      <td>True</td>
      <td>Sampling</td>
      <td>9.943855</td>
      <td>9.530689</td>
    </tr>
  </tbody>
</table>
</div>



Experimental details (in our example `params.txt`) can be added as txt, csv or xlsx formats. This file can include the columns described in the following table:

**Table 2.** Experimental parameters table. This example covers all possible columns. Nonetheless, not all columns are necessary. For example, Enrichment columns (EnrichmentMode, StdDilution, StdVolume) are only used if any enrichment step was performed. More information about this is described in the Enrichment section.

| Condition   | SampleVolume | ProteinConcentration | AmountMS | CellsPerML | TotalCultureVolume | ProteinSRM | fmolSRM | Enrichment | EnrichmentMode | StdDilution | StdVolume |
|-------------|--------------|----------------------|----------|------------|--------------------|------------|---------|------------|---------------------|-------------|-----------|
| Cond1_t0    | 2.31         | 2.99                 | 9.67     | 4.54       | 7.54               | TNAMLN     | 4.44    | False      |                 | 3.96        | 1.22      |
| Cond2_t1    | 2.50         | 0.20                 | 4.10     | 5.13       | 2.62               | AJFVYC     | 4.85    | True       | Sampling                | 2.43        | 1.51      |
| Cond3_t2    | 7.38         | 6.56                 | 2.77     | 3.66       | 3.80               | BYEKSC     | 9.71    | True       | Amplification                | 5.71        | 8.53      |

## Proteome fraction enrichment (Optional)

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



```python
# Enrichment Standard importation

enrichment_std = pd.read_excel('enrichment_std.xlsx')
enrichment_std.sample(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Protein</th>
      <th>Accession</th>
      <th>Chain length</th>
      <th>MW (kDa)</th>
      <th>StdConcentration</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>4</th>
      <td>Lisozyme</td>
      <td>P00698</td>
      <td>129 aa</td>
      <td>14.3</td>
      <td>17.780320</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ADH</td>
      <td>P00330</td>
      <td>348 aa</td>
      <td>36.8</td>
      <td>45.756347</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Soybean Trypsin Inhibitor</td>
      <td>P01071</td>
      <td>181 aa</td>
      <td>20.0</td>
      <td>24.867580</td>
    </tr>
  </tbody>
</table>
</div>




```python
enrichment_std, sample_prep_updated = alpaca.gathers(clean_df, enrichment_std, sample_prep)
```

    Enrichment factor for condition Before_Induction: 17.97
    Enrichment factor for condition Control: 34.65


## Data integration

This module connects the protein amounts quantified in the sample and the sample preparation. Thus, allowing to calculate protein amounts to the original state (e.g. bacterial culture, raw culture supernatant). This step brings deeper insights to the user based on the known experimental parameters, yielding high valuable data (e.g., molecules per cell, fmol / µmol of protein extract)


```python
results = alpaca.wool(quant_df, sample_prep_updated)

results.sample(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Accession</th>
      <th>Protein</th>
      <th>Sample</th>
      <th>iBAQ</th>
      <th>Condition</th>
      <th>Replicate</th>
      <th>fmol</th>
      <th>Molecules</th>
      <th>fmolSample</th>
      <th>MoleculesPerCell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>14063</th>
      <td>P29252</td>
      <td>FolK</td>
      <td>iBAQ Diamide_03</td>
      <td>-1.290911</td>
      <td>Diamide</td>
      <td>03</td>
      <td>17.039546</td>
      <td>5.304567e+10</td>
      <td>88.071837</td>
      <td>2.809453e+09</td>
    </tr>
    <tr>
      <th>14186</th>
      <td>P39120</td>
      <td>CitZ</td>
      <td>iBAQ Diamide_03</td>
      <td>3.895811</td>
      <td>Diamide</td>
      <td>03</td>
      <td>533.536478</td>
      <td>1.660948e+12</td>
      <td>2757.675511</td>
      <td>8.796864e+10</td>
    </tr>
    <tr>
      <th>6599</th>
      <td>Q01464</td>
      <td>MinD</td>
      <td>iBAQ Control_01</td>
      <td>-0.565776</td>
      <td>Control</td>
      <td>01</td>
      <td>9.159222</td>
      <td>3.760164e+11</td>
      <td>624.300851</td>
      <td>2.457695e+10</td>
    </tr>
  </tbody>
</table>
</div>




```python

```
