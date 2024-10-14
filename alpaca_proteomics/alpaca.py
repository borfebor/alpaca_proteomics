import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import re
import sys, os
from alpaca_proteomics_test.spits import detective, Normalization, tools, clean, Imputation
from alpaca_proteomics_test.census import Quantification, advisor
from alpaca_proteomics_test.gathers import gathering
from alpaca_proteomics_test.wool import yarnScissors

class alpaca:
    
    def eats(file):
        """
        Loads a file (either in CSV, TSV, or XLSX format) into a Pandas DataFrame and attempts 
        to detect protein ID and intensity methods.
    
        Parameters:
        ----------
        file : str or file-like object
            The file can be provided as a string representing the file path (if it's on disk)
            or as a file-like object (e.g., uploaded file in Streamlit).
    
        Returns:
        -------
        tuple : (pandas.DataFrame, str, dict) or pandas.DataFrame
            Returns the DataFrame, the detected ID column, and the intensity methods (if detected).
            If detection fails, returns the DataFrame alone.
        """
        
        if isinstance(file, str):
            file_name = file
        else:
            file_name = file.name

        if 'TXT' in file_name.upper() or 'TSV' in file_name.upper():
            df = pd.read_csv(file, sep='\t')  # Tab-separated values
        elif 'CSV' in file_name.upper():
            df = pd.read_csv(file, sep=',')  # Comma-separated values
        elif 'XLSX' in file_name.upper():
            df = pd.read_excel(file)  # Excel file
        else:
            print('Not compatible format')
            return None  # Return None explicitly if format is unsupported
    
        try:
            id_col, is_pivot, it = detective.alpacaHolmes(df) 
            print(f"The column {id_col} was detected to contain your ProteinGroups IDs.\nThe following intensity methods were detected in the data: {', '.join(i for i in it.keys())}")
            return df, id_col, it
        except:
            #print("This is a bit awkward, I could not identify the columns containing your Protein IDs and intensity methods.\nYou'll have to define them in the spits function parameters to process the data correctly.")
            return df
    
    def spits(data, lfq_method, id_col='auto', replicate_dict='auto', intensity_dict='auto', info_cols=None, 
              contamination_cols=['identified by site', 'contaminant', 'Reverse'],
              cleaning=True, formatting='auto', 
             transformation=np.log2, normalization=None, valid_values=0.7, imputation='', **imp_kwargs):
        """
        Processes a DataFrame for quantitative analysis, performing tasks such as cleaning, 
        transformation, imputation, normalization, and reformatting based on provided parameters.
        
        Parameters:
        ----------
        data : pandas.DataFrame
            Input dataframe containing raw data for processing.
        id_col : str
            Column name in `df` representing unique identifiers for the dataset (e.g., 'Accession').
        lfq_method : str
            Column name or method representing LFQ (Label-Free Quantification) values in the data.
        replicate_dict : dict
            A dictionary mapping sample names to condition and replicate information.
        cleaning : bool, optional, default=True
            Whether to clean the data by removing contaminants and decoys.
        formatting : str or bool, optional, default='auto'
            Controls the format of the output dataframe ('auto', True, or False).
        transformation : callable, optional, default=np.log2
            A transformation function to apply to the data (e.g., log2).
        normalization : str, optional, default='None'
            Method for normalization. Accepted values: 'None', 'Relative', 'Median', 'Quantile'.
        valid_values : float, optional, default=0.7
            Proportion of valid values required for each row to be retained.
        imputation : str, optional, default=''
            Method for data imputation. If empty, no imputation is performed.
        **imp_kwargs : dict
            Additional keyword arguments passed to the imputation function.
        
        Returns:
        -------
        pandas.DataFrame
            Processed DataFrame ready for further analysis or visualization.
        """
        df = data.copy()
        # Check and update id_col and intensity_dict if they are 'auto'
        if id_col == 'auto' or intensity_dict == 'auto':
            # Call alpacaHolmes, but only update the 'auto' fields
            detected_id_col, is_pivot, detected_intensity_dict = detective.alpacaHolmes(df)
            
            # Update only the fields that are still 'auto'
            if id_col == 'auto':
                id_col = detected_id_col
            if intensity_dict == 'auto':
                intensity_dict = detected_intensity_dict
        
        # Raise an error if both id_col and intensity_dict are still 'auto'
        if id_col == 'auto' or intensity_dict == 'auto':
            raise ValueError("'id_col' and 'intensity_dict' must not both be 'auto'.")
        
        # Check and update replicate_dict if it is 'auto'
        if replicate_dict == 'auto':
            replicate_dict, conditions, it = detective.alpacaWatson(df, intensity_dict, id_col, lfq_method)

        df.columns = df.columns.str.replace('.: ', '')
        samples = list(replicate_dict.keys())
        conditions = list(tools.invert_dict(replicate_dict).keys())
        print(f"{len(conditions)} experimental conditions were identified in the data: {', '.join(conditions)}")

        if 'Accession' not in df.columns:
            df.columns = [re.sub(id_col, 'Accession', id_col) 
                          if id_col == col else col for col in df.columns]
  
        # Removal of contaminants, and decoys
        cont_key = [col for col in df.select_dtypes(exclude=np.number).columns for item in contamination_cols if item in col]
    
        # Select columns with additional information 
        all_ids, priority = ['Accession'],  ['name', 'kDa']
        [all_ids.append(col) for col in df.columns for i in priority if i.upper() in col.upper()]
        [all_ids.append(col) for col in df.select_dtypes(exclude=np.number).columns if col not in all_ids]
        if info_cols is None:
            info_cols = all_ids
        ids = [col for col in df.select_dtypes(exclude=np.number).columns if col in info_cols]
        
        if cleaning is True:
            df = df[df[cont_key].isna().all(axis=1)].reset_index(drop=True)
            print(f"Items marked on {', '.join(cont_key)} have been removed from the dataset.")
    
        df = df[ids+samples].replace(0, np.nan)
        log_samples = [i for i in samples if 'LOG' not in i.upper()]
        df[log_samples] = df[log_samples].apply(lambda x: transformation(x)) 
        df[samples] = clean.filter_rows_by_missing_values(df[samples], replicate_dict, valid_values)

        if imputation != '':
            imputation_methods = {
                'None': ('', {''}),
                'LOD': (Imputation.impute_lod, {'lod': 0.01}),
                'ND': (Imputation.impute_nd, {'lod': 0.01}),# 'mean': mean, 'std': std}),
                'kNN': (Imputation.impute_knn, {'n_neighbors': 5}),
                'LLS': (Imputation.impute_lls, {'max_iter': 10}),
                'SVD': (Imputation.impute_svd, {'n_components': 2}),
                }
            if imputation not in list(imputation_methods.keys()):
                raise ValueError(f"Invalid imputation method. Available options: {', '.join(imputation_methods.keys())}")
            imp_method = imputation_methods[imputation][0]
            imp_kwargs = imputation_methods[imputation][1]
            df[samples] = Imputation.by_condition(df[samples], replicate_dict, imp_method, **imp_kwargs)
            
        norm_methods = {'None':None, 'Relative': Normalization.Relative, 'Median': Normalization.Median, 'Quantile': Normalization.Quantile}
        if normalization not in list(norm_methods.keys()):
                raise ValueError(f"Invalid normalization method. Available options: {', '.join(norm_methods.keys())}")
        df = Normalization.normalizer(df, samples, norm_methods[normalization])
        
        if formatting == 'auto':
            formatting = tools.check_format(df, 'Accession')
    
        if formatting is True:
            df = df.melt(id_vars=ids, value_vars=samples, 
                        var_name='Sample', value_name=lfq_method).replace(0, np.nan)
            df[['Condition', 'Replicate']] = df['Sample'].replace(replicate_dict).str.split(';', expand=True)
            df = df.dropna(subset=lfq_method)

            if 'Gene names' in ids:
                df['Gene names'] = df['Gene names'].str[0].str.upper() + df['Gene names'].str[1:]
                df = df.rename(columns={'Gene names':'Protein'})
                
            print('Dataset formated for further analysis and visualisation.')   
            return df
        else:
            print('Data is formated for human vision.\nThat could lead to errors or incompatibilities in further analysis using Alpaca pipeline.\nConsider formating your dataset if you see any anomally.')
            return df
    
    def census(df, standards, concentration=0.5, in_sample=6.0, lfq_col='iBAQ', ratio=1, 
               total_protein= 1, filter_col='Replicate', added_samples=None, valid_values=0,
               plot=True, save=False):
        '''
        Performs protein quantification using regression between quantified protein intensities and 
        dynamic standards (e.g., UPS2).
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing clean data for quantified proteins.
        standards : str or file-like object (.csv, .txt, .xlsx)
            Path to a file containing UPS2 dynamic standards information or the standards file itself.
        concentration : float, optional, default=0.5
            Stock concentration of the standards in mg/mL.
        in_sample : float, optional, default=6.0
            Volume (in microliters) of standards added to each sample.
        lfq_col : str, optional, default='iBAQ'
            Column in `df` representing the label-free quantification (LFQ) values.
        ratio : float, optional, default=1
            A multiplier for adjusting the calculated concentration of each protein.
        total_protein : float, optional, default=1
            Total protein concentration in the sample.
        filter_col : str, optional, default='Replicate'
            Column to filter data for specific replicates or conditions.
        added_samples : list or None, optional, default=None
            List of samples or conditions where standards were added. If None, assume standards were added to all samples.
        valid_values : int, optional, default=2
            Minimum number of valid (non-missing) values required for regression.
        plot : bool, optional, default=True
            Whether to generate and display regression plots.
        save : bool, optional, default=False
            Whether to save regression plots.
    
        Returns
        -------
        df : pandas.DataFrame
            DataFrame with quantified proteins, adjusted for concentration.
        ups_red : pandas.DataFrame
            DataFrame containing measured standards in the sample.
        coef : float
            Regression slope (used for quantification).
        inter : float
            Regression intercept (used for quantification).
        R2 : float
            R-squared value representing the goodness-of-fit for the regression.
        '''
        try:
            ups2 = Quantification.abacus(standards, concentration, in_sample, total_protein=total_protein)
        except Exception as e:
            raise ValueError(f"Error processing standards: {e}")
        if lfq_col not in df.columns:
            raise ValueError(f"The specified lfq_col '{lfq_col}' does not exist in the dataframe.")
        if added_samples is None:
            added_samples = df[filter_col].unique().tolist()  # Apply to all replicates by default
        ups_red, coef, inter, R2 = Quantification.regression(df, ups2, lfq_col=lfq_col, filter_col=filter_col,
                                                     added_samples=added_samples, valid_values=valid_values) # Regression between intensities and standards
        if R2 < 0.8:
            print(f"Low R² value ({R2:.2f}), indicating a poor fit for the regression.")
        if plot == True:
            Quantification.abacusreg(ups_red, lfq_col=lfq_col, R2=R2, save=save) # Plots the Regression
        df = Quantification.moles(df, coef, inter, lfq_col=lfq_col, ratio = ratio)
    
        return df, ups_red, coef, inter, R2
    
    def gathers(df, enrichment_standards, preparation, lfq_method='iBAQ', plot=False, save_plot=False):
        """"
        Calculates enrichment factors for protein standards spiked into samples and optionally plots the results.
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing the quantified protein data (e.g., with iBAQ or LFQ values).
        enrichment_standards : pandas.DataFrame
            DataFrame or file containing the standard information for calculating enrichment.
        preparation : pandas.DataFrame
            DataFrame containing information about the experimental preparation (e.g., conditions).
        lfq_method : str, optional
            Label-free quantification method used in the analysis (e.g., 'iBAQ', 'LFQ'). Default is 'iBAQ'.
        plot : bool, optional
            Whether to plot the enrichment factors. Default is False.
        save_plot : bool, optional
            Whether to save the plot as a file. Default is False.
    
        Returns
        -------
        e_test : pandas.DataFrame
            DataFrame containing calculated enrichment values for each sample.
        preparation : pandas.DataFrame
            Updated preparation DataFrame with enrichment factors added.
        """
        required_columns = ['Condition', 'Replicate', lfq_method]
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns from the input DataFrame: {', '.join(missing_columns)}")
            
        e_test = gathering.enrichment_calculator(df, enrichment_standards, preparation, lfq_method).dropna(subset='Enrichment')
        grouping = ['Condition', 'Replicate']
        col_grouper = [columns for columns in df.columns if columns in grouping]

        enrichments = e_test.groupby(col_grouper)['Enrichment'].median().reset_index()
        col_grouper.remove('Replicate')
        enrichment_factors = enrichments.groupby(col_grouper).apply(lambda x: pd.Series({
                                                                    'EnrichmentFactor':x["Enrichment"].median()})
            ).reset_index()

        preparation = preparation.merge(enrichment_factors, on='Condition', how='left')
            
        for index, row in enrichment_factors.iterrows():
            print(f'Enrichment factor for condition {row["Condition"]}: {np.round(row["EnrichmentFactor"],2)}')
     
        if plot == True:
            gathering.enrichment_plot(enrichments, save_plot)
        
        return e_test, preparation   
     
    def wool(df, preparation):
        """"
        Integrates experimental preparation to the measured molar amounts. Applies enrichment factor corrections to a dataframe, converts amounts to molecules using Avogadro's number,
        and computes sample-specific and cell-specific molecule concentrations.
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing the quantified protein data (e.g., fmol values).
        preparation : pandas.DataFrame
            DataFrame containing experimental setup and enrichment information (e.g., Enrichment factors, sample volumes).
    
        Returns
        -------
        df : pandas.DataFrame
            Updated DataFrame with corrected fmol values and molecule counts.
        """
        enrichment_params = ['Enrichment', 'EnrichmentMode', 'ProteinSRM', 'fmolSRM', 'EnrichmentFactor']
        sample_params = ['SampleVolume', 'ProteinConcentration', 'AmountMS']
        cells_params = ['CellsPerML', 'TotalCultureVolume']
        
        required_columns = ['Condition', 'fmol']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"Required columns {required_columns} are missing in the input DataFrame")
        
        if 'EnrichmentFactor' in preparation.columns:
        
            for condition, values in preparation.set_index('Condition')[enrichment_params].fillna(1).iterrows():
                
                if values['EnrichmentMode'] == 'Amplification':
                    """
                    This calculation is made for samples which correspond to a higher fraction 
                    compared to the original proteome. E.g., Membrane
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                        df.fmol / values['EnrichmentFactor'], 
                                        df.fmol)
                    
                elif values['EnrichmentMode'] == 'Sampling': 
                    """
                    This calculation is made for samples which correspond to a smaller fraction
                    to the original proteome. E.g., Secretome
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                          df.fmol * values['EnrichmentFactor'], 
                                          df.fmol)
                    
            preparation = yarnScissors.correctionSRM(df, preparation)

            if "CorrectionSRM" in preparation.columns:
                
                for condition, values in preparation.set_index('Condition').fillna(1).iterrows():
                        
                            df['fmol'] = np.where(df.Condition == condition, 
                                                df.fmol * values['CorrectionSRM'], 
                                                df.fmol)
        
        df['Molecules'] = df['fmol'] * 6.023e8  # Avogadro's number fixed for fmol (-15)
        
        if all(item in preparation.columns.to_list() for item in sample_params):
            
            df['fmolSample'] = np.nan
            for condition, values in preparation.set_index('Condition')[sample_params].fillna(1).iterrows():
                
                total_protein = values['SampleVolume'] * values['ProteinConcentration'] #calculate the µg of protein in the sample
                MS_to_sample = total_protein / values['AmountMS']
                
                df['fmolSample'] = np.where(df['Condition'] == condition,
                                                    df['fmol'] * MS_to_sample, 
                                                    df['fmolSample'])
                
                df['Molecules'] = df['fmolSample'] * 6.023e8 # Avogadro's number fixed for fmol (-15)
            
        if all(item in preparation.columns.to_list() for item in cells_params):
               
            df['MoleculesPerCell'] = np.nan
            for condition, values in preparation.set_index('Condition')[cells_params].fillna(1).iterrows():
                
                cells = values['CellsPerML'] * values['TotalCultureVolume']
                
                df['MoleculesPerCell'] = np.where(df['Condition'] == condition,
                                                    df['Molecules'] / cells, 
                                                    df['MoleculesPerCell'])
        
        return df
    
    def Consultant(df, st_proteins, it, added_samples='all',
                norm_options = ['None', 'Relative', 'Median', 'Quantile'],
                values_per_sample=0.3, ):
        
        """
        Suggests the best normalization and intensity method for protein quantification based on regression performance
        with spiked-in standards.
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing quantified protein data.
        st_proteins : pandas.DataFrame
            DataFrame containing the spiked-in standards for regression-based quantification.
        it : dict
            Dictionary with intensity methods as keys and corresponding sample lists as values.
        added_samples : list or str, optional
            List of sample names to apply or 'all' to apply to all samples. Default is 'all'.
        norm_options : list of str, optional
            List of normalization options to test. Default includes ['None', 'Relative', 'Median', 'Quantile'].
        values_per_sample : float, optional
            Minimum proportion of valid values required per sample for analysis. Default is 0.3.
    
        Returns
        -------
        output : pandas.DataFrame
            DataFrame with tested normalization, intensity methods, and corresponding R² score.
        """

        export = []
        for k, v in it.items():
            if added_samples != 'all':
                samples = [sample for i in added_samples for sample in v if i in sample]
            else:
                samples = v
            for n in norm_options:
                mock = alpaca.ninja_spits(df, k, formatting=True, valid_values=values_per_sample,
                                   normalization=n)
                anchor = st_proteins.copy()
                quant_df, measured_anchor, coef, inter, r2 = alpaca.census(mock, anchor, lfq_col=k,
                                                                     filter_col = 'Sample', # Defines which column to filter for the spiked samples
                                                                     added_samples=samples,
                                                                    plot=False) # Adding which samples contain the standars
                export.append([n, k, r2])
        output = pd.DataFrame(export, columns=['Normalization', 'Intensity method', 'score'])
        met = output['score'].idxmax()
        lfq_advise, norm_advise = output.iloc[met]['Intensity method'], output.iloc[met]['Normalization']
        if norm_advise != 'None':
            print(f'Based on your data, {norm_advise}-normalized {lfq_advise} is recommended for the quantification.')
        else:
            print(f'Based on your data, {lfq_advise} without normalization is recommended for the quantification.')
        return output
    
    def ninja_spits(data, lfq_method, id_col='auto', replicate_dict='auto', intensity_dict='auto', info_cols=None, 
              contamination_cols=['identified by site', 'contaminant', 'Reverse'],
              cleaning=True, formatting='auto', 
             transformation=np.log2, normalization=None, valid_values=0.7, imputation='', **imp_kwargs):

        df = data.copy()
        # Check and update id_col and intensity_dict if they are 'auto'
        if id_col == 'auto' or intensity_dict == 'auto':
            # Call alpacaHolmes, but only update the 'auto' fields
            detected_id_col, is_pivot, detected_intensity_dict = detective.alpacaHolmes(df)
            
            # Update only the fields that are still 'auto'
            if id_col == 'auto':
                id_col = detected_id_col
            if intensity_dict == 'auto':
                intensity_dict = detected_intensity_dict
        
        # Raise an error if both id_col and intensity_dict are still 'auto'
        if id_col == 'auto' or intensity_dict == 'auto':
            raise ValueError("'id_col' and 'intensity_dict' must not both be 'auto'.")
        
        # Check and update replicate_dict if it is 'auto'
        if replicate_dict == 'auto':
            replicate_dict, conditions, it = detective.alpacaWatson(df, intensity_dict, id_col, lfq_method)

        df.columns = df.columns.str.replace('.: ', '')
        samples = list(replicate_dict.keys())

        if 'Accession' not in df.columns:
            df.columns = [re.sub(id_col, 'Accession', id_col) 
                          if id_col == col else col for col in df.columns]
  
        # Removal of contaminants, and decoys
        cont_key = [col for col in df.select_dtypes(exclude=np.number).columns for item in contamination_cols if item in col]
    
        # Select columns with additional information 
        all_ids, priority = ['Accession'],  ['name', 'kDa']
        [all_ids.append(col) for col in df.columns for i in priority if i.upper() in col.upper()]
        [all_ids.append(col) for col in df.select_dtypes(exclude=np.number).columns if col not in all_ids]
        if info_cols is None:
            info_cols = all_ids
        ids = [col for col in df.select_dtypes(exclude=np.number).columns if col in info_cols]
        
        if cleaning is True:
            df = df[df[cont_key].isna().all(axis=1)].reset_index(drop=True)
            #print(f"Items marked on {', '.join(cont_key)} have been removed from the dataset.")
    
        df = df[ids+samples].replace(0, np.nan)
        log_samples = [i for i in samples if 'LOG' not in i.upper()]
        df[log_samples] = df[log_samples].apply(lambda x: transformation(x)) 
        df[samples] = clean.filter_rows_by_missing_values(df[samples], replicate_dict, valid_values)
        if imputation != '':
            imputation_methods = {
                'None': ('', {''}),
                'LOD': (Imputation.impute_lod, {'lod': 0.01}),
                'ND': (Imputation.impute_nd, {'lod': 0.01}),# 'mean': mean, 'std': std}),
                'kNN': (Imputation.impute_knn, {'n_neighbors': 5}),
                'LLS': (Imputation.impute_lls, {'max_iter': 10}),
                #'Random Forest (RF)': (Imputation.impute_rf, {'max_iter': 10}),
                'SVD': (Imputation.impute_svd, {'n_components': 2}),
                #'Bayesian Principal Component Analysis (BPCA)': (Imputation.impute_bpca, {'n_components': 2, 'max_iter': 100})
                }
            if imputation not in list(imputation_methods.keys()):
                raise ValueError(f"{imputation} not in the accepted methods/input.\nThe accepted methods are {', '.join(list(imputation_methods.keys()))}")
            imp_method = imputation_methods[imputation][0]
            imp_kwargs = imputation_methods[imputation][1]
            df[samples] = Imputation.by_condition(df[samples], replicate_dict, imp_method, **imp_kwargs)
            
        norm_dict = {'None':None, 'Relative': Normalization.Relative, 'Median': Normalization.Median, 'Quantile': Normalization.Quantile}
        if normalization not in list(norm_dict.keys()):
            raise ValueError(f"{normalization} not in the accepted methods/input.\nThe accepted methods are {', '.join(list(norm_dict.keys()))}")
        df = Normalization.normalizer(df, samples, norm_dict[normalization])
        
        if formatting == 'auto':
            formatting = tools.check_format(df, 'Accession')
    
        if formatting is True:
            df = df.melt(id_vars=ids, value_vars=samples, 
                        var_name='Sample', value_name=lfq_method).replace(0, np.nan)
            df[['Condition', 'Replicate']] = df['Sample'].replace(replicate_dict).str.split(';', expand=True)
            df = df.dropna(subset=lfq_method)

            if 'Gene names' in ids:
                df['Gene names'] = df['Gene names'].str[0].str.upper() + df['Gene names'].str[1:]
                df = df.rename(columns={'Gene names':'Protein'})
                
            #print('Dataset formated for further analysis and visualisation.')   
            return df
        else:
            #print('Data is formated for human vision.\nThat could lead to errors or incompatibilities in further analysis using Alpaca pipeline.\nConsider formating your dataset if you see any anomally.')
            return df
   
    def generate_example_params(df):
        """
        Given a dataframe `df` containing columns 'Condition', 'Replicate', and 'Sample',
        generate an example parameter table with random values for the other columns.
        """
    
        # Get unique values for the 'Condition', 'Replicate', and 'Sample' columns
        conditions = df['Condition'].unique()
    
        # Generate random values for the other columns
        n_conditions = len(conditions)
    
        param_names = ['SampleVolume', 'ProteinConcentration', 'AmountMS',
                       'CellsPerML', 'TotalCultureVolume', 
                       'ProteinSRM', 'fmolSRM', 
                       'Enrichment', 'EnrichmentMode', 'StdDilution', 'StdVolume']
        param_types = [float, float, float,
                       float, float, 
                       str, float,
                       bool, str, float, float]
    
        # Generate random values for each parameter
        data = {}
        for name, dtype in zip(param_names, param_types):
            if dtype == bool:
                data[name] = np.random.choice([True, False], size=(n_conditions))
            elif name == 'EnrichmentDirection':
                data[name] = np.random.choice(['Amplification', 'Sampling'], size=(n_conditions))
            elif name == 'ProteinSRM':
                data[name] = np.random.choice(df.Accession, size=(n_conditions))
            else:
                data[name] = np.random.rand(n_conditions) * 10
    
        # Create dataframe
        index = conditions
        df_params = pd.DataFrame(data=data, index=index, columns=param_names).reset_index(names='Condition')
    
        return df_params
        
