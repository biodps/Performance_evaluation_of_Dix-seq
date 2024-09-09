import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype

cr_97_stats_df_all = pd.read_csv(
    'cr_97_stats_df_all.csv')

cr_97_stats_df_all['software'] = cr_97_stats_df_all['software'].astype(CategoricalDtype(categories = ['q1_open_ref_97', 'q1_open_ref_99', 'dix_seq', 'ez_amplicon', 'q2_dada2', 'q2_deblur'], ordered = True))

cr_97_stats_df_all['cr_q1_97_otus_count_pct'] = cr_97_stats_df_all['cr_q1_97_otus_count'] / (cr_97_stats_df_all['cr_q1_97_otus_count'] + cr_97_stats_df_all['de_novo_otus_count'])
cr_97_stats_df_all.drop(columns = ['cr_q1_97_otus_count', 'de_novo_otus_count'], inplace = True)
cr_97_stats_df_all['cr_q1_97_otus_pct'] = cr_97_stats_df_all['cr_97_counts'] / (cr_97_stats_df_all['cr_97_counts'] + cr_97_stats_df_all['de_novo_counts'])
cr_97_stats_df_all.drop(columns = ['cr_97_counts', 'de_novo_counts'], inplace = True)

# wo factoring reads count
cr_97_otus_counts_all = cr_97_stats_df_all.groupby(['depth', 'software'])['cr_q1_97_otus_count_pct'].mean().reset_index()
cr_97_otus_counts_all = cr_97_otus_counts_all.loc[cr_97_otus_counts_all['software'] != 'q1_open_ref_97', :]
# w factoring reads count
cr_97_otus_all = cr_97_stats_df_all.groupby(['depth', 'software'])['cr_q1_97_otus_pct'].mean().reset_index()
cr_97_otus_all = cr_97_otus_all.loc[cr_97_otus_all['software'] != 'q1_open_ref_97', :]
stats_df_all = pd.merge(cr_97_otus_counts_all, cr_97_otus_all, on = ['depth', 'software'])

depth_lst = [cr_97_stats_df_all['depth'].unique()[0]]

depth_lst = depth_lst + \
    np.flip(cr_97_stats_df_all['depth'].unique()[1:]).tolist()

dfs = {}

for depth in cr_97_stats_df_all['depth'].unique():

    df = stats_df_all.loc[stats_df_all['depth'] == depth, [
        'software', 'cr_q1_97_otus_count_pct','cr_q1_97_otus_pct']]

    dfs[depth] = df.melt(id_vars = ['software'], value_vars = ['software', 'cr_q1_97_otus_count_pct','cr_q1_97_otus_pct'], var_name = 'indicator', value_name = 'value')
    dfs[depth]['value'] = dfs[depth]['value'] * 100

