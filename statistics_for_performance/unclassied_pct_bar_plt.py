import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype

unclassified_pct_df = pd.read_csv('unclassified_pct_df.csv')
# unclassified_pct_df = pd.read_csv(
#     'unclassified_features_comp_total_reads_counts_df.csv')
unclassified_pct_df.columns = ['software', 'depth', 'rank', 'pct']

depth_lst = [unclassified_pct_df['depth'].unique()[0]]

depth_lst = depth_lst + \
    np.flip(unclassified_pct_df['depth'].unique()[1:]).tolist()

dfs = {}

for rank in unclassified_pct_df['rank'].unique():

    df = unclassified_pct_df.loc[unclassified_pct_df['rank'] == rank, [
        'software', 'depth', 'pct']]

    df['depth'] = df['depth'].astype(
        CategoricalDtype(categories=depth_lst, ordered=True))

    dfs[rank] = df
