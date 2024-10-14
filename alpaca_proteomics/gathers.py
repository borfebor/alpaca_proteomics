#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:19:34 2024

@author: borfebor
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

class gathering:
    
    def standards_preparation_check(preparation):
    
        try:
            vol_col = [col for col in preparation.columns if ("STDV" or "STD_V" or "STD V") in col.upper()][0]
        except:
            print('Added volume is missing in sample parameters')
        
        try:
            splv_col = [col for col in preparation.columns if ("SAMPLEV" or "SAMPLE_V" or "SAMPLE V") in col.upper()][0]
        except:
            print('Sample volume is missing in sample parameters')
        
        try:
            dil_col = [col for col in preparation.columns if ("STDD" or "STD_D" or "STD D") in col.upper()][0]
        except:
            print('Dilution of the standard stock solution is missing in sample parameters')
            
        return vol_col, splv_col, dil_col
    
    def multiplier(enriched_conditions, standards):
    
        arranger = pd.DataFrame()
    
        for condition in enriched_conditions:
    
            standards['Condition'] = condition
            arranger = pd.concat([arranger, standards])
            
        return arranger
    
    def enrichment_calculator(df, standards, preparation, lfq_method='iBAQ', subproteome=None):
        
        enriched_conditions = [item[1]['Condition'] for item in preparation.iterrows() if item[1]['Enrichment'] == True]
        
        if enriched_conditions == []:
            pass
        else:
            vol_col, splv_col, dil_col = gathering.standards_preparation_check(preparation)
    
            standards_mod = gathering.multiplier(enriched_conditions, standards)
    
            standards_mod = standards_mod.merge(preparation, on=['Condition'])
            
            MW = [col for col in standards_mod.columns if 'Da' in col][0]
    
            standards_mod['StdMass'] = standards_mod['StdConcentration'] / standards_mod[dil_col] * standards_mod[vol_col]
            standards_mod['StdFmol'] = standards_mod['StdMass'] / standards_mod[MW] * 1000
    
            standards_mod['StdFmol'] = np.where(standards_mod['StdFmol'] == 0, np.nan, standards_mod['StdFmol'])
    
            standards_mod = standards_mod[['Accession', 'Condition', MW, 'StdMass','StdFmol']]
    
            ID_standards = df[df.Condition.isin(enriched_conditions)
                             ].merge(standards_mod, how='right', on=["Accession", "Condition"])
            
            ID_standards['Enrichment'] = ID_standards['StdFmol'] / ID_standards['fmol']  
    
            e_reduction = ['Accession', MW, 'Sample', lfq_method, 
                                           'Condition', 'Replicate', 'fmol', 'StdFmol', 'Subproteome', 'Enrichment']
    
            ID_standards = ID_standards[[col for col in ID_standards.columns if col in e_reduction]]
                    
            return ID_standards
        
    def enrichment_plot(enrichments, save_plot=False):
        sns.set(font_scale=1.5)
        ggplot_styles = {
            'axes.edgecolor': 'black',
            'axes.facecolor': 'EBEBEB',
            'axes.grid': True,
            'axes.grid.which': 'both',
            'axes.grid.axis': 'y',
            'axes.spines.left': True,
            'axes.spines.right': True,
            'axes.spines.top': True,
            'axes.spines.bottom': True,
            'grid.color': 'black',
            'grid.linewidth': '0.2',
            'xtick.color': 'black',
            'xtick.major.bottom': True,
            'xtick.minor.bottom': False,
            'ytick.color': 'black',
            'ytick.major.left': True,
            'ytick.minor.left': False,
        }
        
        plt.rcParams.update(ggplot_styles)
        sns.catplot(data=enrichments, x='Condition', y='Enrichment', linecolor='k',
                    hue='Condition',
                    kind='box', width=0.5)
        sns.despine(top=False, right=False)
        if save_plot == True:
            plt.savefig('enrichment_std.svg', bbox_inches='tight', pad_inches=0.5)