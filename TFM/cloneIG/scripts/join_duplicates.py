"""DUPLICATE MANAGEMENT

This script allows the user to detect rearrangement duplicates (identical
CDR3 or VDJ) for the same condition, selecting the rearrangement with the 
highest frequency and replacing its value with the frequency sum of all rows.

    * Just semicolon-separated value files (.csv) from MiXCR or Vidjil 
    filtering are accepted.
    * `pandas` previous installation in working environment is required.
    * The script is consisted of two functions: sum_freq and join_duplicates.
    * Implemented in a snakemake pipeline.

Author: Andrea Sanchez de la Cruz - 8.6.2022
"""

import pandas as pd

## SNAKEMAKE I/O ##

filtered_tool = snakemake.input.get("filtered")
joined_tool = snakemake.output.get("joined")
log = snakemake.log_fmt_shell(stdout = True, stderr = True)

## FUNCTIONS ##

def sum_freq(df, colname):
    """Select row with the highest Frequency value and
    replace that value with the frequency sum of all rows.

    Parameters
    ----------
    df : pandas DataFrame
        Condition-specific data frame with MiXCR or Vidjil data. 
    colname : str
        Column name of given data frame.

    Return
    ------
    list
        a list of one-row data frames with different colname values.
    """
    
    elem_list = []

    for elem in df[colname].unique():
        
        #subset according to element value
        elem_df = df[df[colname] == elem]

        #if more than 1 row in subset
        if len(elem_df) > 1:
            
            #select row with highest frequency
            high_freq = elem_df[
                elem_df.Frequency == elem_df.Frequency.max()][0:1]

            #sum frequencies of all rows
            sum_freq = elem_df.Frequency.sum()

            #set reference rearrangement
            ref = high_freq.copy()

            #insert freq sum in ref Frequency column
            ref['Frequency'] = sum_freq

            elem_list.append(ref)
        
        #return intact row
        else: elem_list.append(elem_df)

    return elem_list


def join_duplicates(tool_filtered):
    """Join rearrangement duplicates by CDR3 or VDJ column.

    Parameters
    ----------
    tool_filtered : CSV file
        Filtered rearrangement data from MiXCR or Vidjil tool. 

    Return
    ------
    pandas DataFrame
        a pandas DataFrame containing the joined data.
    """
    
    #load data from csv file
    df =  pd.read_csv(tool_filtered, sep = ';')

    if len(df) > 0:

        #remove D suffix in VDJ subgenes
        df['VDJ'] = df['VDJ'].replace(r'D([-|_])', r'\1', regex = True)

        cond_list = []

        #unique condition names in data frame
        cond_array = df['Condition'].unique()
        
        for cond in cond_array:
            #condition subset
            cond_df = df[df.Condition == cond]

            #sum frequencies according to aaCDR3
            CDR3_list = sum_freq(cond_df,'aaCDR3')

            #non-CDR3 data subset
            na_df = cond_df[cond_df.aaCDR3.isna()]
            
            if len(na_df) > 0:

                #sum frequencies according to VDJ
                CDR3_list += sum_freq(na_df,'VDJ')

            #join unique rearrangements in data frame
            CDR3_merge = pd.concat(CDR3_list)

            #join rows with different CDR3 but identical VDJ
            if len(CDR3_merge) > len(CDR3_merge['VDJ'].unique()):
                
                #sum frequencies according to VDJ
                VDJ_list = sum_freq(CDR3_merge,'VDJ')

                cond_list.append(pd.concat(VDJ_list))

            else: cond_list.append(CDR3_merge) 

        #join condition data frames in a single one
        return pd.concat(cond_list)

    else: return pd.DataFrame([])


## MAIN CODE ##

joined_data = join_duplicates(filtered_tool)

#save data frame to CSV file
joined_data.to_csv(joined_tool, sep = ";")
