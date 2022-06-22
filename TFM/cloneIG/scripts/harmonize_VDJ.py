"""HARMONIZATION OF V(D)J NOMENCLATURE

This script allows the user to modify V(D)J names in MiXCR and Vidjil to
a standardized nomenclature, only taking information of major VDJ genes.

    * Just semicolon-separated value files (.csv) from MiXCR or Vidjil 
    duplicate-joined output are accepted.
    * `pandas` previous installation in working environment is required.
    * The script is consisted of one function: harmonize_VDJ.
    * Implemented in a snakemake pipeline.

Author: Andrea Sanchez de la Cruz - 8.6.2022
"""

import pandas as pd

## SNAKEMAKE I/O ##

joined_tool = snakemake.input.get("joined")
harmonized_tool = snakemake.output.get("harmonized")
log = snakemake.log_fmt_shell(stdout = True, stderr = True)

## FUNCTIONS ##

def harmonize_VDJ(tool_joined):
    """Standardize MiXCR and Vidjil VDJ nomenclature.

    Parameters
    ----------
    tool_joined : CSV file
        Non-duplicated rearrangement data from MiXCR or Vidjil tool. 

    Return
    ------
    pandas DataFrame
        a pandas DataFrame containing the harmonized data.
    """
    
    #load data from csv file
    df_std = pd.read_csv(tool_joined, sep = ';')

    if len(df_std) > 0:

        ##store subset without rearrangement data aside
        df_missing = df_std[df_std.VDJ == '-']

        ##modify VDJ colummn: remove subgenes
        vdj_std = df_std.VDJ.replace('-([A-Z]{1,2})?[0-9]{1,2}', '', 
                                     regex = True)
        #remove commas (MiXCR)
        vdj_std = vdj_std.replace(',', '', regex = True)
        #remove repeated VDJ genes in the same row (MiXCR)
        vdj_std = vdj_std.replace(r'([V|D|J]{1}\d{1,2})\1+', r'\1', 
                                  regex = True)
        #add '/' when several hits in V, D and J genes in the same row
        vdj_std = vdj_std.replace(r'(\d{1,2})[V|D|J]', r'\1/', regex = True)

        #replace VDJ column
        df_std['VDJ'] = vdj_std

        ##select IGH locus (Vidjil)
        df_stdH = df_std[df_std.Locus.isin(['IGH','IGH+'])].copy()

        #insert 'Unk' between repeated underscore
        vdj_std = df_stdH.VDJ.replace(r'__', r'_Unk_', regex = True)
        #replace several V, D or J with 'Unk' when more than 2 hits (MiXCR)
        vdj_std = vdj_std.replace(r'[V|D|J]\d(/\d){2,}', 'Unk', regex = True)

        #replace VDJ column
        df_stdH['VDJ'] = vdj_std

        ##select IGK locus (Vidjil)
        df_stdK = df_std[df_std.Locus.isin(['IGK','IGK+'])].copy()

        #remove repeated underscore
        vdj_std = df_stdK.VDJ.replace(r'__', r'_', regex = True)

        #replace VDJ column
        df_stdK['VDJ'] = vdj_std

        return pd.concat([df_stdH, df_stdK, df_missing])

    else: return pd.DataFrame([])


## MAIN CODE ##

std_data = harmonize_VDJ(joined_tool)

#save data frame to CSV file
std_data.to_csv(harmonized_tool, sep = ";")
