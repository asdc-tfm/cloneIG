"""MIXCR AND VIDJIL MERGE

This script allows the user to merge the obtained rearrangements of 
MiXCR- and Vidjil-analysed sequencing data.

    * Just semicolon-separated value files (.csv) from MiXCR and Vidjil 
    harmonized outputs are accepted.
    * `pandas` previous installation in working environment is required.
    * The script is consisted of one function: merge_tools.
    * Implemented in a snakemake pipeline.

Author: Andrea Sanchez de la Cruz - 8.6.2022
"""

import pandas as pd

# SNAKEMAKE I/O

m_std = snakemake.input.get("mixcr_df")
v_std = snakemake.input.get("vidjil_df")
merge_df = snakemake.output.get("merge_df")
warnings_file = snakemake.params.get("warnings_file")
log = snakemake.log_fmt_shell(stdout = True, stderr = True)

## FUNCTIONS ##

def merge_tools(mixcr_harmonized, vidjil_harmonized):
    """Merge VDJ-standardized MiXCR and Vidjil data frames.

    Parameters
    ----------
    mixcr_harmonized : CSV file
        Rearrangement VDJ-standardized data for MiXCR tool. 
    vidjil_harmonized : CSV file
        Rearrangement VDJ-standardized data for Vidjil tool.

    Return
    ------
    pandas DataFrame or zero value.
        a pandas DataFrame containing the merged data 
        or zero value when no rearrangements found.
    """

    #load data frames
    mixcr_df = pd.read_csv(mixcr_harmonized, sep = ';')
    vidjil_df = pd.read_csv(vidjil_harmonized, sep = ';')

    #rearrangement data in both tools
    if len(vidjil_df) > 0 and len(mixcr_df) > 0:

        #Vidjil df: remove Tool column and rename Reads and Freq names
        vidjil_df = vidjil_df.drop(columns = 'Tool')
        vidjil_df.rename(inplace = True, columns = {'Reads':'vReads',
                                                    'Frequency':'vFreq'})

        #MiXCR df: remove Tool colum and rename Reads and Freq names
        mixcr_df = mixcr_df.drop(columns = 'Tool')
        mixcr_df.rename(inplace = True, columns = {'Reads':'mReads',
                                                   'Frequency':'mFreq',
                                                   'VDJ':'mVDJ', 
                                                   'Sequence':'mSeq'})

        #merge both data frames according to 'on' list
        all_merge = pd.merge(vidjil_df, mixcr_df, 
                             how = 'outer', 
                             on = ['Sample','Condition',
                                   'Locus','nCDR3','aaCDR3'])

        #replace NaN values with '-'
        all_merge = all_merge.fillna('-')

        merge_list = []
        
        ##subset of available D gene in Vidjil
        if len(all_merge[~all_merge.VDJ.str.contains(r'_Unk_')]) > 0:
            
            complete_df = all_merge[~all_merge.VDJ.str.contains(r'_Unk_')].copy()
        
            #add MiXCR Seq and VDJ name in Sequence and VDJ columns
            complete_df['Sequence'] = complete_df.apply(
                lambda x: x.Sequence.replace('-', x.mSeq), 1)
            complete_df['VDJ'] = complete_df.apply(
                lambda x: x.VDJ.replace('-', x.mVDJ), 1)
            
            #append to list if non-empty data frame
            if len(complete_df) > 0 : merge_list.append(complete_df)
        
        ##subset of unknown D gene in Vidjil
        if len(all_merge[all_merge.VDJ.str.contains(r'_Unk_')]) > 0:
            
            unk_df = all_merge[all_merge.VDJ.str.contains(r'_Unk_')].copy()
            
            #leave just D gene in MiXCR VDJ name when just 1 hit
            unk_df['mVDJ'] = unk_df.mVDJ.replace(
                r'V\d[\d/]*_(D\d{1,2})_J\d[\d/]*', 
                r'\1', regex = True)

            #otherwise, replace MiXCR VDJ call with 'Unk'
            unk_df.loc[unk_df.mVDJ.str.contains(r'/|Unk|-'), 'mVDJ'] = 'Unk'

            #replace Vidjil D gene with MiXCR D gene
            unk_df['VDJ'] = unk_df.apply(lambda x: x.VDJ.replace('Unk', 
                                                                x.mVDJ), 1)

            #append to list if non-empty data frame
            if len(unk_df) > 0: merge_list.append(unk_df)
            
        ##join subsets and select columns of interest
        merged_df = pd.concat(merge_list)
        merged_df = merged_df[['Sample','Condition','Locus','Sequence',
                               'mReads', 'mFreq','vReads','vFreq','VDJ',
                               'nCDR3','aaCDR3']]

        #sort rows by sample, locus and condition
        merged_df = merged_df.sort_values(by = ['Sample','Locus','Condition'])

        #delete no-meaning rows with zeros
        merged_df = merged_df.drop(merged_df[(merged_df['mReads'] == 0.0) & 
                                            (merged_df['vReads'] == '-')].index)
        merged_df = merged_df.drop(merged_df[(merged_df['vReads'] == 0.0) & 
                                            (merged_df['mReads'] == '-')].index)

    #rearrangement data just in MiXCR
    elif len(vidjil_df) == 0 and len(mixcr_df) > 0:
        
        #MiXCR df: remove Tool col and rename reads and freq colnames
        merged_df = mixcr_df.drop(columns = 'Tool')
        merged_df.rename(inplace = True, columns = {'Reads':'mReads',
                                                    'Frequency': 'mFreq'})

        #replace NaN values with '-'
        merged_df = merged_df.fillna('-')
        
        #add vReads, vFreq with zero values
        merged_df.insert(6, 'vReads', 0)
        merged_df.insert(7, 'vFreq', 0.0)

    #rearrangement data just in Vidjil
    elif len(vidjil_df) > 0 and len(mixcr_df) == 0:
        
        #Vidjil df: remove Tool col and rename reads and freq colnames
        merged_df = vidjil_df.drop(columns = 'Tool')
        merged_df.rename(inplace = True, columns = {'Reads':'vReads',
                                                    'Frequency': 'vFreq'})
        #replace NaN values with '-'
        merged_df = merged_df.fillna('-')
        
        #add mReads and mFreq with zero values
        merged_df.insert(4, 'mReads', 0)
        merged_df.insert(5, 'mFreq', 0.0)

    #no rearrangements in any of the tools
    else: merged_df = 0

    return merged_df


## MAIN CODE ##

merged_data = merge_tools(m_std, v_std)

#variable contains a data frame
if type(merged_data) != int:
    #save data frame in CSV file
    merged_data.to_csv(merge_df, sep = ';', index = False)

#variable == 0
else:
    #get sample name
    sample = merge_df.split("/")[-1].split("_")[0]

    #save empty data frame in CSV file
    pd.DataFrame([]).to_csv(merge_df, sep = ";")
    
    #print warning in TXT file
    with open(warnings_file, "a") as f:
        f.write(f"No rearrangements found for {sample}.\n")
