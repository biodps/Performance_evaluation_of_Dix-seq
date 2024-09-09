import pandas as pd
import numpy as np

tax_diff_df_sw_rank = pd.read_csv(
    'tax_diff_df_sw_rank.csv')

depth_lst = [tax_diff_df_sw_rank['depth'].unique()[0]]

depth_lst = depth_lst + \
    np.flip(tax_diff_df_sw_rank['depth'].unique()[1:]).tolist()
    
depth_lst_num = [100, 98, 90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10]

depth_rp_dict = dict(zip(depth_lst, depth_lst_num))

tax_diff_df_sw_rank['depth'] = tax_diff_df_sw_rank['depth'].replace(depth_rp_dict)

dfs = {}

for sw in tax_diff_df_sw_rank['software'].unique():

    df = tax_diff_df_sw_rank.loc[tax_diff_df_sw_rank['software'] == sw, [
        'depth', 'rank', 'tax_diff_pct']]

    dfs[sw] = df
    

    

