import pandas as pd
import numpy as np
import os
import biom
from qiime2 import Artifact
import skbio
import hashlib

md = pd.read_excel('diff_software_diff_depth_res.xlsx', sheet_name = 'diff_software_diff_depth_res')

wd = '/home/navi/shared_wd/Dix-Seq'

unclasified_pct_df = pd.DataFrame(columns = ['software', 'depth', 'rank', 'pct'])

for index, row in md[['software', 'depth']].drop_duplicates().iterrows():
    software = row['software']
    software_lst = []
    software_lst.append(software)
    depth = row['depth']
    depth_lst = []
    depth_lst.append(depth)
    print('now working on: ' + software + ' ' + depth + '...')
    fn = software + "_" + depth + ".pkl"
    fp = os.path.join(wd, 'mega_tabs', fn)
    df = pd.read_pickle(fp)
    df.set_index('md5_id', inplace = True)

    rank_unclasified_pct_df = pd.DataFrame(columns = ['rank', 'pct'])
    
    for rank in [x.lower() for x in ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]:
        rank_lst = []
        rank_lst.append(rank)
        print('current rank: ' + rank + '...')
        tmp_se_current_rank = df.loc[:, df.columns.str.lower().isin(rank_lst)].iloc[:, 0]
        unclassified_sum = (tmp_se_current_rank.isin(['unclassified', 'unassigned', '']) | (tmp_se_current_rank.isna())).sum()
                
        unclassfied_pct_current_rank = unclassified_sum / len(tmp_se_current_rank)
        unclassfied_pct_current_rank = unclassfied_pct_current_rank * 100
        unclassfied_pct_current_rank_lst = []
        unclassfied_pct_current_rank_lst.append(unclassfied_pct_current_rank)
        tmp_stats_df = pd.DataFrame({'rank': rank_lst, 'pct': unclassfied_pct_current_rank_lst})
        rank_unclasified_pct_df = pd.concat([rank_unclasified_pct_df, tmp_stats_df], axis = 0, ignore_index = True)
        rank_unclasified_pct_df[['software', 'depth']] = [software, depth]
        rank_unclasified_pct_df = rank_unclasified_pct_df[['software', 'depth', 'rank', 'pct']]
    unclasified_pct_df = pd.concat([unclasified_pct_df, rank_unclasified_pct_df], axis = 0, ignore_index = True)

unclasified_pct_df.to_csv('/home/navi/Data-Visualization-Projects/Dix-Seq/unclassified_pct_df.csv', index = False)
