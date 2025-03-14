### `alpaca.gathers(df, enrichment_standards, preparation, lfq_method='iBAQ', plot=False, save_plot=False)`

Calculates enrichment factors for protein standards spiked into samples and optionally plots the results.

**Parameters:**

- **`df`** (`pandas.DataFrame`):  
  DataFrame containing the quantified protein data (e.g., with iBAQ or LFQ values).

- **`enrichment_standards`** (`pandas.DataFrame` or file-like object):  
  DataFrame or file containing the standard information for calculating enrichment.

- **`preparation`** (`pandas.DataFrame`):  
  DataFrame containing information about the experimental preparation (e.g., conditions).

- **`lfq_method`** (`str`, optional, default=`'iBAQ'`):  
  Label-free quantification method used in the analysis (e.g., 'iBAQ', 'LFQ').

- **`plot`** (`bool`, optional, default=`False`):  
  Whether to plot the enrichment factors.

- **`save_plot`** (`bool`, optional, default=`False`):  
  Whether to save the plot as a file.

**Returns:**

- **`e_test`** (`pandas.DataFrame`):  
  DataFrame containing calculated enrichment values for each sample.

- **`preparation`** (`pandas.DataFrame`):  
  Updated preparation DataFrame with enrichment factors added.

**Notes:**

- The function calculates the enrichment factors by comparing the protein levels in the experimental samples to the spiked-in standards.
- It groups the data by experimental conditions and replicates to calculate median enrichment values.
- The function also optionally generates and saves a plot of the enrichment values.

