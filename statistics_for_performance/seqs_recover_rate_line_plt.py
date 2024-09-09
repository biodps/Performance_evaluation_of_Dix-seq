import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype

seqs_count_long_df = pd.read_csv('seqs_count_long_tab.csv')

sequencing_depth_df = pd.read_csv('sequencing_depth.csv')

recover_rate_df = pd.merge(seqs_count_long_df, sequencing_depth_df, on = ['sample_id', 'depth'], how = 'left', validate = 'many_to_one')

recover_rate_df['valid_data_recover_rate'] = recover_rate_df['reads_mapped_to_features'] / recover_rate_df['sequencing_depth']

recover_rate_df['valid_data_recover_rate'] = round(recover_rate_df['valid_data_recover_rate'] * 100, 2)

depth_lst = [recover_rate_df['depth'].unique()[0]]

depth_lst = depth_lst + np.flip(recover_rate_df['depth'].unique()[1:]).tolist()

software_lst = recover_rate_df['software'].unique()

df = recover_rate_df.loc[:, ['software', 'depth', 'sample_id', 'valid_data_recover_rate']]

df['depth'] = df['depth'].astype(CategoricalDtype(categories = depth_lst, ordered = True))
