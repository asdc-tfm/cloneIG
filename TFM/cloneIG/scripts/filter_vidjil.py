"""Vidjil FILTERING

This script allows the user to filter a list of Vidjil patient-specific 
outputs and obtain a single data frame with the rearrangement data for
all given conditions.

    * Just tab-separated value files (.tsv) from Vidjil software are accepted.
    * `pandas` and `numpy` installation in working environment is required.
    * The script is consisted of one function: filter_vidjil.
    * Implemented in a snakemake pipeline.

Author: Andrea Sanchez de la Cruz - 8.6.2022
"""

import pandas as pd
import numpy as np

## SNAKEMAKE I/O ##

vidjil_output = snakemake.input.get("vidjil_output")
vidjil_filtered = snakemake.output.get("vidjil_filtered")
log = snakemake.log_fmt_shell(stdout = True, stderr = True)

## FUNCTIONS ##

def filter_vidjil(vidjil_list):
    """Filter Vidjil patient-specific output.

    Parameters
    ----------
    vidjil_list : list
        List of Vidjil outputs at different conditions. 

    Return
    ------
    pandas DataFrame
        a pandas DataFrame containing Vidjil filtered data.
    """

    df_list = []
    missing_cond = []

    for file in vidjil_list:

        #sample and condition from filename
        sample = file.split('/')[-1].split('_')[0]
        condition = file.split('/')[-1].split('_')[1].split('.')[0]

        #load file content in data frame
        tsv_df = pd.read_csv(file, sep = '\t')

        #filter by productive and KDE rearrangements
        prod = tsv_df[tsv_df.productive == 'T']
        KDE = tsv_df[tsv_df.j_call == 'KDE']
        tsv_df = pd.concat([prod, KDE])

        #continue with next file if no rows left
        if not len(tsv_df): 
            missing_cond.append(condition)
            continue

        #select cols of interest
        cols_to_keep = ['locus','sequence','duplicate_count','v_call',
                        'd_call','j_call','productive','cdr3_aa',
                        'cdr3_sequence_start','cdr3_sequence_end']

        #keep cols of interest
        tsv_df = tsv_df[cols_to_keep]

        #get global reads from vidjil file
        vidjil_file = file.split('.')[0] + '.vidjil'
        with open(vidjil_file, 'r') as vfile:
            for line in vfile:
                if 'junction detected in' in line:
                    global_reads = line.split()[5] 

        #insert tool, condition, sample_id and reads in df
        tsv_df.insert(loc = 0, column = 'Tool', value = 'Vidjil')
        tsv_df.insert(loc = 0, column = 'Condition', value = condition)
        tsv_df.insert(loc = 0, column = 'Sample', value = sample)
        tsv_df.insert(loc = 5, column = 'Reads', value = global_reads)

        #convert variables to integers
        cdr3_coordinates = ['cdr3_sequence_start','cdr3_sequence_end']
        cdr3_df = tsv_df[cdr3_coordinates].fillna(0)
        cdr3_df = cdr3_df.astype({'cdr3_sequence_start':'int'})
        cdr3_df = cdr3_df.astype({'cdr3_sequence_end':'int'})
        tsv_df = tsv_df.astype({'Reads':'int'})

        #add int CDR3 coordinates to df
        tsv_df['start'] = cdr3_df['cdr3_sequence_start']
        tsv_df['end'] = cdr3_df['cdr3_sequence_end']

        #find nCDR3 in sequence using coordinates
        nCDR3 = tsv_df.apply(lambda x: x['sequence'][x['start']-1:x['end']], 1)
        nCDR3 = nCDR3.replace('', np.nan, regex = False)
        tsv_df.insert(loc = 11, column = 'nCDR3', value = nCDR3)

        #calculate frequency and round to third decimal
        frequency = tsv_df['duplicate_count'] / tsv_df['Reads']
        tsv_df.insert(loc = 6, column = 'Frequency', value = frequency)
        tsv_df = tsv_df.round({'Frequency': 3})

        #rename df columns
        tsv_df.rename(inplace = True, 
                      columns = {'locus':'Locus',
                                 'sequence': 'Sequence',
                                 'v_call': 'V',
                                 'd_call': 'D',
                                 'j_call': 'J',
                                 'cdr3_aa': 'aaCDR3'})

        #delete IGH|IGK and *00 allele info in V, D and J cols
        VDJ = ['V','D','J']
        VDJ_df = tsv_df[VDJ].replace('IG[H|K]', '', regex = True)
        VDJ_df = VDJ_df[VDJ].replace('\*[0-9]{2}', '', regex = True)

        #join VDJ columns in one V(D)J column and insert in df
        VDJ_col =VDJ_df['V'] + '_' + VDJ_df['D'].fillna('') + '_' + VDJ_df['J']
        tsv_df.insert(loc = 8, column = 'VDJ', value = VDJ_col)

        #delete unnecessary columns
        del_cols = ['duplicate_count','productive','start','end',
                    'cdr3_sequence_start','cdr3_sequence_end']
        tsv_df.drop(inplace = True, columns = VDJ + del_cols)

        df_list.append(tsv_df)

    #available rearrangement data
    if len(df_list) > 0:
        
        #if no data in a specific condition
        if len(missing_cond) > 0:
            
            #insert a dummy row with zero values
            for c in missing_cond:
                
                dummy_row = df_list[0].head(1).copy()
                dummy_row['Condition'] = c
                dummy_row[['Locus', 'Sequence', 'VDJ', 'nCDR3', 'aaCDR3']] ='-'
                dummy_row[['Reads', 'Frequency']] = 0
                df_list.append(dummy_row)

        #join conditions
        return pd.concat(df_list)

    #no rearrangements after filtering
    else: return pd.DataFrame([])

    
## MAIN CODE ##

vidjil_data = filter_vidjil(vidjil_output)

#save data frame to CSV file
vidjil_data.to_csv(vidjil_filtered, sep = ';')
