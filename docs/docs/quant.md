# 🔢 Quantification

The quantification module of the **ALPACA (Absolute Protein Quantification)** proteomics pipeline is designed to facilitate the analysis of absolute protein quantification data. This Python-based tool streamlines the processing of mass spectrometry data, enabling researchers to accurately determine protein abundances in complex biological samples. By integrating with various proteomics data formats and employing robust statistical methods, the module ensures precise quantification, thereby enhancing the reliability of downstream analyses.

This module is especially valuable for researchers working with large-scale proteomics datasets, as it simplifies the steps necessary to generate quantitative data that can be used for a wide range of biological insights.

---

## 🧮 **Function: `alpaca.census()`**

The primary function in the quantification module, `alpaca.census()`, is responsible for performing the absolute quantification of protein abundances using various statistical methods.

---

### `alpaca.census(df, standards, concentration=0.5, in_sample=6.0, lfq_col='iBAQ', ratio=1, total_protein=1, filter_col='Replicate', added_samples=None, valid_values=0, plot=True, save=False)`

Performs protein quantification using regression between quantified protein intensities and dynamic standards (e.g., UPS2).

**Parameters:**

- **`df`** (`pandas.DataFrame`):  
  DataFrame containing clean data for quantified proteins.

- **`standards`** (`str` or file-like object):  
  Path to a file containing UPS2 dynamic standards information or the standards file itself (e.g., .csv, .txt, .xlsx).

- **`concentration`** (`float`, optional, default=`0.5`):  
  Stock concentration of the standards in mg/mL.

- **`in_sample`** (`float`, optional, default=`6.0`):  
  Volume (in microliters) of standards added to each sample.

- **`lfq_col`** (`str`, optional, default=`'iBAQ'`):  
  Column name in `df` representing the label-free quantification (LFQ) values.

- **`ratio`** (`float`, optional, default=`1`):  
  A multiplier for adjusting the calculated concentration of each protein.

- **`total_protein`** (`float`, optional, default=`1`):  
  Total protein concentration in the sample.

- **`filter_col`** (`str`, optional, default=`'Replicate'`):  
  Column name used to filter data for specific replicates or conditions.

- **`added_samples`** (`list`, optional, default=`None`):  
  List of samples or conditions where standards were added. If `None`, assumes standards were added to all samples.

- **`valid_values`** (`int`, optional, default=`2`):  
  Minimum number of valid (non-missing) values required for regression.

- **`plot`** (`bool`, optional, default=`True`):  
  Whether to generate and display regression plots.

- **`save`** (`bool`, optional, default=`False`):  
  Whether to save the regression plots.

**Returns:**

- **`df`** (`pandas.DataFrame`):  
  DataFrame with quantified proteins, adjusted for concentration.

- **`ups_red`** (`pandas.DataFrame`):  
  DataFrame containing measured standards in the sample.

- **`coef`** (`float`):  
  Regression slope used for quantification.

- **`inter`** (`float`):  
  Regression intercept used for quantification.

- **`R2`** (`float`):  
  R-squared value representing the goodness-of-fit for the regression.

**Notes:**

- The function performs regression analysis between protein intensities in `df` and spiked-in dynamic standards (e.g., UPS2) to quantify proteins.
- It provides an option to visualize and save regression plots.
- A low R² value (less than 0.8) indicates a poor fit for the regression and should be carefully reviewed.


### 🧾 **Signature**
```python
alpaca.census(data, method="sum", target_protein="total_protein", normalize=False)``

