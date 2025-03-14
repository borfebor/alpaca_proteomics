# 🧼 Data Pre-processing

Data import and pre-processing are foundational steps in proteomics analysis, ensuring that raw experimental data is transformed into a clean, standardized format for downstream analysis. In the context of `alpaca_proteomics`, these steps are critical and applicable for both relative and absolute quantification in quantitative proteomics.  

Proper pre-processing removes contaminants, handles missing values, and normalizes datasets, minimizing technical variability and enhancing the accuracy of protein quantification. This ensures robust and reproducible results, enabling meaningful biological insights. The library streamlines these tasks, making it highly applicable for workflows involving large-scale proteomics data.

---

## 🔁 **Overview of the Workflow**

The pre-processing pipeline in `alpaca` consists of three main steps:

1. **Import** raw data using `alpaca.eats()`  
2. **Consult (Optional)** an automated advisor via `alpaca.Consultant()` to choose optimal parameters
3. **Clean** the dataset with `alpaca.spits()`  


---

## 📥 **1. Importing Data with `alpaca.eats()`**

📚 **Function Overview**

`alpaca.eats()` ingests proteomics data files in various formats and returns a structured `pandas.DataFrame`.

### `alpaca.eats(file, inspect=True)`

Loads proteomics data from a file in `.txt`, `.tsv`, `.csv`, or `.xlsx` format.

**Parameters:**

- **`file`** (`str` or file-like):  
  Path to the data file or an open file object. Supported formats include:
  	- `.txt` or `.tsv`: Tab-separated values  
 	- `.csv`: Comma-separated values  
 	- `.xlsx`: Excel spreadsheet

- **`inspect`** (`bool`, default=`True`):  
  If `True`, the function tries to identify the column containing protein group IDs and the available intensity methods using `detective.alpacaHolmes`.

**Returns:**

- If `inspect=True` and detection is successful:  
  Returns a tuple `(DataFrame, id_col, intensity_dict)`

- If `inspect=True` but detection fails:  
  Returns the `DataFrame` with a warning message.

- If `inspect=False`:  
  Returns the `DataFrame` only.

**Notes:**

- If the file format is not recognized, the function prints a

 

**Example Usage:**

```python
from alpaca_proteomics import alpaca

file = "results.csv"

data = alpaca.eats(file)
```

## 🤖 **2. Quantification assistance with`alpaca.Consultant()`**

### `alpaca.Consultant(df, st_proteins, it, added_samples='all', norm_options=['None', 'Relative', 'Median', 'Quantile'], values_per_sample=0.3)`

Suggests the best normalization and intensity method for protein quantification based on regression performance with spiked-in standards.

**Parameters:**

- **`df`** (`pandas.DataFrame`):  
  DataFrame containing the quantified protein data to be analyzed.

- **`st_proteins`** (`pandas.DataFrame`):  
  DataFrame containing the spiked-in standards used for regression-based quantification.

- **`it`** (`dict`):  
  Dictionary where keys are intensity methods (e.g., 'LFQ') and values are lists of corresponding sample names.

- **`added_samples`** (`list` or `'all'`, default=`'all'`):  
  List of sample names to apply the analysis to. If `'all'`, applies to all samples.

- **`norm_options`** (`list`, default=`['None', 'Relative', 'Median', 'Quantile']`):  
  List of normalization options to test during analysis.

- **`values_per_sample`** (`float`, default=`0.3`):  
  Minimum proportion of valid values (non-missing) required per sample for analysis.

**Returns:**

- `pandas.DataFrame`:  
  DataFrame with tested normalization and intensity methods, and their corresponding R² scores.

**Notes:**

- The function runs regression analysis for each combination of normalization and intensity method, evaluating the quality of fit using R².
- Based on the highest R² score, the function recommends the best combination of intensity method and normalization for protein quantification.



## 🧹 **3. Data pre-processing with `alpaca.spits()`**

### `alpaca.spits(data, lfq_method, id_col='auto', replicate_dict='auto', intensity_dict='auto', info_cols=None, contamination_cols=['identified by site', 'contaminant', 'Reverse'], cleaning=True, formatting='auto', transformation=np.log2, normalization=None, valid_values=0.7, imputation='', **imp_kwargs)`

Processes a proteomics DataFrame for quantitative analysis, performing tasks such as cleaning, transformation, imputation, normalization, and formatting for downstream analysis or visualization.

**Parameters:**

- **`data`** (`pandas.DataFrame`):  
  Input DataFrame containing the raw intensity data.

- **`lfq_method`** (`str`):  
  Column name or label-free quantification (LFQ) method used to extract intensities.

- **`id_col`** (`str`, default=`'auto'`):  
  Column name containing unique protein identifiers. If `'auto'`, it is inferred using `detective.alpacaHolmes`.

- **`replicate_dict`** (`dict` or `'auto'`, default=`'auto'`):  
  Dictionary mapping sample names to `Condition;Replicate` values. If `'auto'`, it is inferred using `detective.alpacaWatson`.

- **`intensity_dict`** (`dict` or `'auto'`, default=`'auto'`):  
  Dictionary mapping LFQ methods to their corresponding columns. Inferred automatically if set to `'auto'`.

- **`info_cols`** (`list`, optional):  
  Columns containing additional metadata to retain. If `None`, selected based on common patterns.

- **`contamination_cols`** (`list`, default=`['identified by site', 'contaminant', 'Reverse']`):  
  Columns used to identify and remove contaminants or reverse hits.

- **`cleaning`** (`bool`, default=`True`):  
  Whether to filter out contaminants and reverse hits.

- **`formatting`** (`str` or `bool`, default=`'auto'`):  
  Controls output format. If `True`, returns long-format data suitable for downstream processing.  
  If `'auto'`, formatting is inferred.

- **`transformation`** (`callable`, default=`np.log2`):  
  Function used to transform intensity values (e.g., log transformation).

- **`normalization`** (`str`, default=`None`):  
  Normalization method. Accepted values: `'None'`, `'Relative'`, `'Median'`, `'Quantile'`.

- **`valid_values`** (`float`, default=`0.7`):  
  Minimum proportion of valid values (non-missing) required to retain a row per condition.

- **`imputation`** (`str`, default=`''`):  
  Imputation method for missing values. Available options:
  `'None'`, `'LOD'`, `'ND'`, `'kNN'`, `'LLS'`, `'SVD'`.  
  If empty, no imputation is applied.

- **`**imp_kwargs`** (`dict`):  
  Additional keyword arguments passed to the imputation function, if specified.

**Returns:**

- `pandas.DataFrame`:  
  A processed DataFrame, either in long format (for analysis and plotting) or wide format (for inspection), depending on the `formatting` parameter.

**Notes:**

- Requires an intensity column and a way to define biological replicates.
- If any key components (`id_col`, `replicate_dict`, or `intensity_dict`) are set to `'auto'`, the function attempts to detect them automatically.
- Designed to integrate with the `alpaca` pipeline and its internal modules like `detective`, `tools`, `Imputation`, and `Normalization`.
