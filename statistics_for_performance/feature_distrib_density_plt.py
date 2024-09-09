import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype

asvs_count_long_df = pd.read_csv('asvs_count_long_tab.csv')

depth_lst = [asvs_count_long_df['depth'].unique()[0]]

depth_lst = depth_lst + np.flip(asvs_count_long_df['depth'].unique()[1:]).tolist()

df = asvs_count_long_df.loc[asvs_count_long_df['depth'] == 'full', ['software', 'Feature_ID', 'cumulative_abundance']]


