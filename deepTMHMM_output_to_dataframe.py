# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 15:33:28 2024

@author: cassp
"""

import pandas as pd

def reformat_tmr_df(file):
    df = pd.read_csv(file, skip_blank_lines = False, header = None, skiprows = 1)
    df = df[df[0].str.startswith('#')]
    df = df[0].str.split(' ', expand=True)
    df_collapsed = df.groupby([1]).first().reset_index()
    df_short = df_collapsed[[1, 3, 6]]
    df_short.columns = ['accession', 'aa_length', 'TMRs']
    df_short['accession'] = df_short['accession'].str.split(':').str.get(0)
    return(df_short)

df1 = reformat_tmr_df(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\Bella\Bioinformatics\spo0B_TM_110525\TMRs (1).gff3")

df2 = reformat_tmr_df(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\Bella\Bioinformatics\spo0B_TM_110525\TMRs (2).gff3")

df3 = reformat_tmr_df(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\Bella\Bioinformatics\spo0B_TM_110525\TMRs (3).gff3")

df_combined = pd.concat([df1, df2, df3], axis=0)

df_combined.to_csv(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\Bella\Bioinformatics\spo0B_TM_110525\TMR_df_all_111025.csv")

