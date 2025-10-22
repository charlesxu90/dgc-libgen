import pandas as pd
import numpy as np
from loguru import logger

from .utils.seq_mutations import NAT_AAs

def get_aa_pos_weight(df_topn_uniq_pos_top):
    """Get the amino acid position weight matrix for the topn unique positions dataframe.
    Args:
        df_topn_uniq_pos_top (pd.DataFrame): The dataframe containing selected positions.
    Returns:
        list: A list of fixed positions.
        pd.DataFrame: A dataframe containing the amino acid position weight matrix.
    """
    # obtain fixed positions and aa 
    fix_positions = []
    fix_aa_mut = []
    for pos, fix_aa in zip(df_topn_uniq_pos_top.pos, df_topn_uniq_pos_top.fix_aa):
        if fix_aa:
            fix_positions.append(int(pos))
            fix_aa_mut.append(fix_aa)
    
    weight_mat = {aa:[] for aa in NAT_AAs}
    for i, row in df_topn_uniq_pos_top.iterrows():
        pos = int(row.pos)
        if pos in fix_positions: 
            # assign fixed positions frequencies directly
            idx = fix_positions.index(pos)
            for aa in NAT_AAs:
                if aa == fix_aa_mut[idx]:
                    weight_mat[aa].append(1.0)
                else:
                    weight_mat[aa].append(0.0)
        else:
            # associated positions
            for aa in NAT_AAs:
                if aa not in row.sel_mut_aa:
                    weight_mat[aa].append(0.0)
                else:
                    if aa == row['wt_aa']:
                        weight_mat[aa].append(row['wt_count'])
                    else:
                        weight_mat[aa].append(row['mut_aa_count'][row['mut_aa'].index(aa)])

    weight_mat_df = pd.DataFrame(weight_mat)
    weight_mat_df['pos'] = [int(pos) for pos in df_topn_uniq_pos_top.pos]
    weight_mat_df.sort_values(by=['pos'], inplace=True)
    weight_mat_df[NAT_AAs] = weight_mat_df[NAT_AAs].div(weight_mat_df[NAT_AAs].sum(axis=1), axis=0).replace(np.nan, 0)
    
    return fix_positions, weight_mat_df.reset_index()


def get_topn_variants_pos_df(df_variants_topn, 
                             wt_seq, 
                             topn=10000):
    """Get the topn unique positions dataframe from the variants dataframe.
    Args:
        df_variants_topn (pd.DataFrame): The dataframe containing the topn variants.
        wt_seq (str): The reference sequence to compare against. Usually is the starting sequence for
            design or wild-type sequence.
        topn (int, optional): The number of top variants to consider. Defaults to 10000.
    Returns:
        pd.DataFrame: The dataframe containing the topn unique positions.
    """
    topn_variants_uniq_pos = list(set([pos for pos_list in df_variants_topn.pos_list for pos in pos_list]))
    topn_variants_uniq_pos.sort()
    logger.info(f'total unique positions in top {topn}: {len(topn_variants_uniq_pos)}')

    topn_variants_uniq_pos_count = {pos: {} for pos in topn_variants_uniq_pos}
    for mut_list in df_variants_topn.mut_list:
        for mut in mut_list:
            pos = int(mut[1:-1])
            aa = mut[-1]
            if aa not in topn_variants_uniq_pos_count[pos]:
                topn_variants_uniq_pos_count[pos][aa] = 0
            topn_variants_uniq_pos_count[pos][aa] += 1

    # logger.info(topn_variants_uniq_pos_count)

    df_topn_variants_uniq_pos = pd.DataFrame(topn_variants_uniq_pos, columns=['pos'])
    df_topn_variants_uniq_pos['wt_aa'] = [wt_seq[pos-1] for pos in df_topn_variants_uniq_pos.pos]
    df_topn_variants_uniq_pos['mut_aa'] = [list(topn_variants_uniq_pos_count[pos].keys()) for pos in df_topn_variants_uniq_pos.pos]
    df_topn_variants_uniq_pos['mut_aa_count'] = [list(topn_variants_uniq_pos_count[pos].values()) for pos in df_topn_variants_uniq_pos.pos]
    df_topn_variants_uniq_pos['wt_count'] = topn - df_topn_variants_uniq_pos['mut_aa_count'].apply(sum)

    df_topn_variants_uniq_pos.sort_values(by='wt_count', inplace=True)
    return df_topn_variants_uniq_pos