"""MiXCR FILTERING

This script allows the user to filter a list of MiXCR patient-specific 
outputs and obtain a single data frame with the rearrangement data for
all given conditions.

    * Just tab-separated value files (.tsv) from MiXCR software are accepted.
    * `pandas` previous installation in working environment is required.
    * The script is consisted of one function: filter_mixcr.
    * Implemented in a snakemake pipeline.

Author: Andrea Sanchez de la Cruz - 8.6.2022
"""

import pandas as pd

## SNAKEMAKE I/O ##

mixcr_output = snakemake.input.get("mixcr_output")
mixcr_filtered = snakemake.output.get("mixcr_filtered")
log = snakemake.log_fmt_shell(stdout = True, stderr = True)

## FUNCTIONS ##

def filter_mixcr(mixcr_list):
    """Filter MiXCR patient-specific output.

    Parameters
    ----------
    vidjil_list : list
        List of MiXCR outputs at different conditions. 

    Return
    ------
    pandas DataFrame
        a pandas DataFrame containing MiXCR filtered data.
    """

    df_list = []
    missing_cond = []

    for file in mixcr_list:

        #sample and condition from filename
        sample = file.split('/')[-1].split('_')[0]
        condition = file.split('/')[-1].split('_')[1]
        
        #load file content in data frame
        fullclones_df = pd.read_csv(file, sep = '\t')

        #filter by min number of reads supporting the clone
        fullclones_df = fullclones_df[fullclones_df.cloneCount > 4.0]

        #remove rows with non-functional or incomplete CDR3
        fullclones_df = fullclones_df[
            ~fullclones_df.aaSeqImputedCDR3.str.contains('\*|_')]

        #continue with next file if no clones left
        if not len(fullclones_df): 
            missing_cond.append(condition)
            continue

        #select cols of interest
        cols_to_keep = ['targetSequences','cloneCount',
                        'cloneFraction','allVHitsWithScore',
                        'allDHitsWithScore','allJHitsWithScore',
                        'nSeqImputedCDR3','aaSeqImputedCDR3']

        #keep cols of interest
        fullclones_df = fullclones_df[cols_to_keep]

        #insert condition and sample_id in df
        fullclones_df.insert(loc = 0, column = 'Tool', value = 'MiXCR')
        fullclones_df.insert(loc = 0, column = 'Condition', value = condition)
        fullclones_df.insert(loc = 0, column = 'Sample', value = sample)

        #calculate global reads, insert in df and convert to int
        global_reads = fullclones_df['cloneCount'] / fullclones_df[
            'cloneFraction']
        fullclones_df.insert(loc = 4, column = 'Reads', value = global_reads)
        fullclones_df = fullclones_df.astype({'Reads':'int'})

        #rename df columns
        fullclones_df.rename(inplace = True, 
                             columns = {'targetSequences':'Sequence',
                                        'cloneFraction': 'Frequency',
                                        'allVHitsWithScore': 'V',
                                        'allDHitsWithScore': 'D',
                                        'allJHitsWithScore': 'J',
                                        'nSeqImputedCDR3': 'nCDR3',
                                        'aaSeqImputedCDR3': 'aaCDR3'})

        #remove first and last codon of CDR3 seq
        fullclones_df['nCDR3'] = fullclones_df.nCDR3.str[3:-3]
        fullclones_df['aaCDR3'] = fullclones_df.aaCDR3.str[1:-1]

        #round frequency values to third decimal
        fullclones_df = fullclones_df.round({'Frequency': 3})

        #extract locus from V gene
        locus_col = fullclones_df['V'].str.split('V', 1, expand = True)[0]
        fullclones_df.insert(loc = 3, column = 'Locus', value = locus_col)

        #delete IGH|IGK and *00(*) allele info in V, D and J cols
        VDJ = ['V','D','J']
        VDJ_df = fullclones_df[VDJ].replace('IG[H|K]', '', regex = True)
        VDJ_df = VDJ_df[VDJ].replace('\*00\([0-9]+(\.[0-9])*\)', '', 
                                     regex = True)

        #join VDJ columns in one V(D)J column and insert in df
        VDJ_col =VDJ_df['V'] + '_' + VDJ_df['D'].fillna('') + '_' + VDJ_df['J']
        fullclones_df.insert(loc = 8, column = 'VDJ', value = VDJ_col)

        #delete unnecessary columns
        fullclones_df.drop(inplace = True, columns = ['cloneCount'] + VDJ)

        df_list.append(fullclones_df)

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

mixcr_data = filter_mixcr(mixcr_output)

#save data frame to CSV file
mixcr_data.to_csv(mixcr_filtered, sep = ';')
