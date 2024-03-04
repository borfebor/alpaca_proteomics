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

**Other packages**
- math / pi
- random / seed, randint

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

- **Condition**:
- **SampleVolume**:
- **ProteinConcentration**:
- **AmountMS**:
- **CellsPerML**:
- **TotalCultureVolume**:
- **ProteinSRM**:
- **fmolSRM**:
- **Enrichment**:
- **EnrichmentDirection**:
- **StdDilution**:
- **StdVolume**:


## Data Importation & Pre-processing
Functions for data import, cleaning and pre-processing.
> alpaca.**importer**(`File`): this function is meant to offer flexibility on the ProteinGroup file importation as some scientist could have the data in a `.txt`, `.csv` or `.xlsx`.

> alpaca.**formater**(`DataFrame`): formater function aims to give a coherence on the imported data, as it could be that MaxQuant output organisation is changed by the user or another software like Perseus. It returns our formated `DataFrame`, and 2 lists: `columns` which contains all df.columns after formating, and `default` that is a list with all suggested columns for dataframe slicing.

> alpaca.**data_cleaner**(`DataFrame`): removes `Reverse`, `Possible contaminants` & `Identifications by site` from the given data and returns the same `dataframe` droping the empty columns. 

> alpaca.**log_transform**(`DataFrame`): transforms `iBAQ intensities` into a log2 scale for intesity linearization. Returns a `DataFrame` 

> alpaca.**experimenter**(`DataFrame`): seeks for the different `Experimental Conditions` in the data and recognises how many (n),
 which ones (condition) and how many replicates (r). Experimenter function returns `condition` list which the recongnised experimental conditions, and 2 integers: `n`and `r`. `n` is the amount of different experimental conditions and `r`the number of replicates.

> alpaca.**replicator**(`Int`): Gets the amount of replicates from experimenter() - `r` - and creates tags for the columns.
It helps solving the problem of variable amount of replicates in different experiments.
It returns a list of replicate names with numbers, that will be used later by condition_format() function.creates `list` with names for the replicates [`Replicate_1`, `Replicate_2`,...,`Replicate_R`] with R being the total amount of replicates.

> alpaca.**condition_format**(`DataFrame`, `list_1`, `list_2`): transforms the data by chunking it according to their experimental condition. `DataFrame` is the working df, `list_1` is our experimental condition list - could be in the output from experimenter() -, and `list_2` is the list of replicates created by replicator() function. It returns a `formated DataFrame`

> alpaca.**namer**(`integer`): takes an `integer` and returns a list of alphabetical characters for names.

> alpaca.**sample_volumes**(`dict, list`): creates a `dictionary` from the dictionary in which differnt enrichment preparations and experimental conditions are stored - `enrichment_type_dict` - and the list of enrichment protocols - `prep` -. Returns a dictionary with Enrichment type as `key` and sample volume as `value`. This function is integrated into `alpaca.getMolecules_cell()` and `alpaca.alpaca()`.
