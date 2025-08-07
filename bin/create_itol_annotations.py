#!/usr/bin/env python

"""
Creates various annotation files for use with ITOL. A few files are required

accessions.txt - single column text file containing list of all accessions in tree

comP_mutants.txt - tab-delimited file (with header) detailing mutations in comP. 
    columns required: Genome Accession

rapP_straints.txt - tab delimited file (with header) detailing rapP status
    columns required: Genome Accession, rapP (boolean: 1), 'N or T': (str: N/T)
"""

import pandas as pd

def get_line(accession, colour):
    """
    Returns output line for background colour files

    Required parameters:
        accession (str): genome accession to highlight
        colour (str):  rgb hex colour code
    """
    return f"{accession},label,node,#FFFFFF,1,normal,{colour}"

def get_comP_rapP_both_header():
    """
    Returns text header for comP/rapP/both status 

    Required parameters:
        None

    Returns:
        header(str): header block
    """

    HEADER = '''
DATASET_BINARY
SEPARATOR COMMA

DATASET_LABEL,Overall_status

COLOR,#000000
FIELD_SHAPES,2,2,2
FIELD_LABELS,comP_mutant,rapP,both
FIELD_COLORS,#6a51a3,#4292c6,#ff7f00

HEIGHT_FACTOR,1.5
SHOW_LABELS,1
LABEL_SHIFT,10

DATA
'''

    return HEADER

def get_comP_header():
    """
    Returns text header for comP status 

    Required parameters:
        None

    Returns:
        header(str): header block
    """

    HEADER = '''
DATASET_STYLE
SEPARATOR COMMA
DATASET_LABEL,comP_mutant
COLOR,#ff0000

DATA
'''

    return HEADER

def get_rapP_header():
    """
    Returns text header for rapP status 

    Required parameters:
        None

    Returns:
        header(str): header block
    """

    HEADER = '''
DATASET_STYLE
SEPARATOR COMMA
DATASET_LABEL,rapP_Status
COLOR,#ff0000

DATA
'''

    return HEADER

def main():

    """
    Main function to create ITOL annotation files.
    """

    overall_df = pd.read_csv('accessions.txt', header=None, names = ['Accession'])
    comP_df = pd.read_csv('comP_mutants.txt', sep="\t")

    comP_df['comP_mutant']=1
    comP_df = comP_df[['Genome Accession', 'comP_mutant']]

    overall_df = overall_df.merge(comP_df, left_on='Accession', right_on='Genome Accession', how = "left")

    rapP_df = pd.read_csv('rapP_strains.txt', sep="\t")
    rapP_df = rapP_df[['Genome Accession', 'rapP', 'N or T']]

    overall_df = overall_df.merge(rapP_df, left_on='Accession', right_on='Genome Accession', how = "left")
    overall_df = overall_df[['Accession','comP_mutant','rapP', 'N or T']]
    overall_df.fillna(value='0', inplace=True)
    overall_df = overall_df.astype({"comP_mutant": int, 'rapP': int})

    overall_df['both']=overall_df.sum(numeric_only=True, min_count=0, axis=1)
    overall_df.loc[overall_df['both']==1, 'both'] = 0
    overall_df.loc[overall_df['both']==2, 'both'] = 1

    # This block for selects comP/rapP/both
    #both = overall_df[overall_df['both']==2]
    # Select subsets for compP/rapP only
    #comP = overall_df.loc[(overall_df['comP_mutant']==1) & (overall_df['rapP']==0)]
    #rapP = overall_df.loc[(overall_df['comP_mutant']==0) & (overall_df['rapP']==1)]

    header = get_comP_rapP_both_header()
    data_block = overall_df[['Accession', 'comP_mutant', 'rapP', 'both']]
    data_block = data_block.to_csv(index=False, header=False)

    with open('overall_status.txt', 'w', encoding='UTF-8') as fh:
        for block in (header, data_block):
            fh.write(block)

    # For rapP N/T status label highlights
    rapP_N = overall_df.loc[overall_df['N or T']=='N']
    rapP_T = overall_df.loc[overall_df['N or T']=='T']

    header = get_rapP_header()
    N_lines = list(rapP_N['Accession'].apply(lambda x: get_line(x, '#6a51a3')))
    T_lines = list(rapP_T['Accession'].apply(lambda x: get_line(x, '#4292c6')))

    with open('rapP_NT_status.txt', 'w', encoding='UTF-8') as fh:
        for block in (header, "\n".join(N_lines), "\n", "\n".join(T_lines)):
            fh.write(block)

    # for comP only
    header = get_comP_header()
    comP = overall_df.loc[overall_df['comP_mutant']==1]
    comP_lines = comP['Accession'].apply(lambda x: get_line(x, '#6a51a3'))
    with open('comP_status.txt', 'w', encoding='UTF-8') as fh:
        for block in (header, "\n".join(comP_lines) ):
            fh.write(block)


if __name__ == "__main__":
    main()
