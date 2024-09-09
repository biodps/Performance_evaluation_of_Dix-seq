import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype

alpha_div_df_all = pd.read_csv(
    'alpha_div_df_all.csv')

alpha_div_df_all['software'] = alpha_div_df_all['software'].astype(CategoricalDtype(categories = ['q1_open_ref_97', 'q1_open_ref_99', 'dix_seq', 'ez_amplicon', 'q2_dada2', 'q2_deblur'], ordered = True))


depth_lst = [alpha_div_df_all['depth'].unique()[0]]

depth_lst = depth_lst + \
    np.flip(alpha_div_df_all['depth'].unique()[1:]).tolist()

dfs = {}

for depth in alpha_div_df_all['depth'].unique():

    df = alpha_div_df_all.loc[alpha_div_df_all['depth'] == depth, [
        'software', 'sample_id', 'observed_otus', 'pielou_e',
       'shannon', 'pd_wt']]

    dfs[depth] = df.melt(id_vars = ['software', 'sample_id'], value_vars = ['observed_otus', 'pielou_e', 'shannon', 'pd_wt'], var_name = 'alpha_div', value_name = 'value')
    

