
import pandas as pd
import numpy as np
from loguru import logger

from .codon import get_codons_from_aa
from .weblogo import draw_weblogo
from .position_freqs import get_aa_pos_weight, get_topn_variants_pos_df
from .utils.seq_mutations import get_mutations
from .utils import plot_style_utils

 
def get_topn_variants(df_variants, 
                      ref_seq, 
                      topn=10000, 
                      seq_col='aa_seqs'):
    """Get the topn variants from the dataframe and annotate mutations and number of mutations.
    Args:
        df_variants (pd.DataFrame): The dataframe containing variant data.
        ref_seq (str): The reference sequence to compare against.
        topn (int, optional): The number of top variants to consider. Defaults to 10000.
        seq_col (str, optional): The column name for sequences in the data. Defaults to 'aa_seqs'.
    Returns:
        pd.DataFrame: The dataframe containing the topn variants with mutations and number of mutations annotated.
    """
    df_variants_topn = df_variants.head(topn).copy()
    df_variants_topn['mutations'] = [get_mutations(aa_seq, wt_seq=ref_seq) for aa_seq in df_variants_topn[seq_col]]
    df_variants_topn['n_muts'] = [len(muts.split('/')) if '/' in muts else 0 for muts in df_variants_topn.mutations]

    df_variants_topn['mut_list'] = [muts.split('/') if '/' in muts else muts for muts in df_variants_topn.mutations]
    df_variants_topn['pos_list'] = [list([int(m[1:-1]) for m in muts]) for muts in df_variants_topn.mut_list]
    return df_variants_topn


def load_variants_data(data_path, 
                       pred_columns=['pred_activity'], 
                       seq_col='aa_seqs'):
    """Load the variants data from a CSV file and rank the variants based on the mean of specified prediction columns.
    Args:
        data_path (str): The path to the CSV file containing variant data.
        pred_columns (list, optional): List of prediction columns to consider for ranking. Defaults to ['pred_activity'].
        seq_col (str, optional): The column name for sequences in the data. Defaults to 'aa_seqs'.
    Returns:
        pd.DataFrame: The dataframe containing the variants with an additional 'rank' column.
    """

    df_variants = pd.read_csv(data_path)
    df_variants['rank'] = df_variants[pred_columns].rank(ascending=False).mean(axis=1)
    df_variants = df_variants[[seq_col, 'rank'] + pred_columns].copy()
    df_variants.sort_values(by='rank', inplace=True, ascending=True)
    return df_variants


def expand_ref_mutations(df_top_uniq_pos, 
                         ref_seq, 
                         wt_seq):
    """Expand the reference mutations to ensure all wild-type positions are included in the dataframe.
    Args:
        df_top_uniq_pos (pd.DataFrame): The dataframe containing selected positions
        df_top1k_uniq_pos_top (pd.DataFrame): The dataframe containing selected positions
        ref_seq (str): The reference sequence to compare against. Usually is the starting sequence for design or wild-type sequence.
        wt_seq (str): The wild-type sequence.
    Returns:
        pd.DataFrame: The expanded dataframe containing all wild-type positions.
        str: The reference mutations string.
    """
    ref_mutations = get_mutations(ref_seq, wt_seq=wt_seq) 

    for wt_mut in ref_mutations.split('/'):
        pos = int(wt_mut[1:-1])
        if pos in df_top_uniq_pos.pos.values:
            continue
        
        aa = wt_mut[-1]
        logger.info(f"Position {wt_mut} of wild-type is not in the selected positions.")
        df_new_row = pd.DataFrame({'pos':[pos], 'wt_aa': [aa], 'fix_aa': [aa], 'n_codons': [1], 'n_expanded_codons': [1],})
        df_top_uniq_pos = pd.concat([df_top_uniq_pos, df_new_row], ignore_index=True)

    logger.info(f"total positions: {len(df_top_uniq_pos)}")
    df_top_uniq_pos.sort_values(by='pos',ascending=True)
    return df_top_uniq_pos, ref_mutations


def filter_pos_aa_for_library(df_top_uniq_pos, 
                              total_mutants=10000, 
                              fix_threshold=0.75, 
                              aa_freq_threshold=0.1):
    """Filter the top variants position dataframe to select positions and amino acids for library design.
    Args:
        df_top_uniq_pos (pd.DataFrame): The dataframe containing position and amino acid frequencies.
        total_mutants (int, optional): The total number of top variants considered. Defaults to 10000.
        fix_threshold (float, optional): The frequency threshold to fix an amino acid at a position. Defaults to 0.75.
        aa_freq_threshold (float, optional): The frequency threshold to select amino acids at a position. Defaults to 0.1
    Returns:
        pd.DataFrame: The filtered dataframe containing selected positions and amino acids for library design.
    """

    df_top_uniq_pos_top = df_top_uniq_pos[df_top_uniq_pos['ref_count']/total_mutants <= fix_threshold].copy()
    logger.info(f"Totally {len(df_top_uniq_pos_top)} positions are selected.")

    df_top_uniq_pos_top['fix_aa'] = [df_top_uniq_pos_top['mut_aa'].iloc[i][0] if df_top_uniq_pos_top['mut_aa_count'].iloc[i][0]/total_mutants >= fix_threshold else '' for i in range(len(df_top_uniq_pos_top)) ]
    df_top_uniq_pos_top['sel_mut_aa'] = [list(aa  for aa, count in zip(row['mut_aa'], row['mut_aa_count']) if count/total_mutants >= aa_freq_threshold) + [row['ref_aa']] for i, row in df_top_uniq_pos_top.iterrows()]
    df_top_uniq_pos_top['dg_codons'], df_top_uniq_pos_top['expanded_codons'] = zip(*df_top_uniq_pos_top.apply(get_codons_from_aa, axis=1))
    df_top_uniq_pos_top['n_codons'] = [len(codons) for codons in df_top_uniq_pos_top.dg_codons]
    df_top_uniq_pos_top['n_expanded_codons'] = [len(codons) for codons in df_top_uniq_pos_top.expanded_codons]
    library_size = np.prod([ n for n in df_top_uniq_pos_top.n_expanded_codons if n > 0])

    logger.info(f"Library size: {library_size}")
    return df_top_uniq_pos_top


def generate_library_for_topn_with_weblogo(data_path, 
                          ref_seq, 
                          wt, 
                          pred_cols=['raw_activity'], 
                          seq_col='aa_seqs', 
                          save_path=None, 
                          topn=1000, 
                          fix_threshold=0.75, 
                          aa_freq_threshold=0.1, ):
    """Plot the topn positions weblogo and save the position table for library design.
    Args:
        data_path (str): The path to the CSV file containing variant data.
        ref_seq (str): The reference sequence to compare against. Usually is the starting sequence for design or wild-type sequence.
        wt (str): The wild-type sequence.
        pred_cols (list, optional): List of prediction columns to consider for ranking. Defaults to ['raw_activity'].
        seq_col (str, optional): The column name for sequences in the data. Defaults to 'aa_seqs'.
        save_path (str, optional): The path to save the figure and position table. If None, it will be derived from data_path. Defaults to None.
        topn (int, optional): The number of top variants to consider. Defaults to 1000.
        fix_threshold (float, optional): The frequency threshold to fix an amino acid at a position. Defaults to 0.75.
        aa_freq_threshold (float, optional): The frequency threshold to select amino acids at a position. Defaults to 0.1.
    Returns:
        pd.DataFrame: The dataframe containing selected positions and amino acids for library design.
    """
    df_variants = load_variants_data(data_path, pred_columns=pred_cols, seq_col=seq_col)
    df_variants_topn = get_topn_variants(df_variants, ref_seq=ref_seq, topn=topn, seq_col=seq_col)
    df_topn_uniq_pos = get_topn_variants_pos_df(df_variants_topn, ref_seq=ref_seq, topn=topn)
    df_topn_uniq_pos = filter_pos_aa_for_library(df_topn_uniq_pos, total_mutants=topn, fix_threshold=fix_threshold, aa_freq_threshold=aa_freq_threshold)
    df_topn_uniq_pos, ref_mutations = expand_ref_mutations(df_topn_uniq_pos, ref_seq=ref_seq, wt_seq=wt)
    logger.info(f"Ref mutations: {ref_mutations}")
    
    plot_style_utils.set_pub_plot_context(context="talk")
    fig, ax = plot_style_utils.simple_ax(figsize=(14, 3))
    plot_style_utils.prettify_ax(ax)
    fix_positions, weight_mat_df = get_aa_pos_weight(df_topn_uniq_pos)
    draw_weblogo(weight_mat_df, fix_positions, ax=ax)
    if not save_path:
        save_path = data_path.replace('.csv', f'_library')

    plot_style_utils.save_for_pub(fig, path=save_path)
    df_topn_uniq_pos.to_csv(save_path + '.csv', index=False)
    logger.info(f"Figure and position table are saved to {save_path}.*")
    
    return df_topn_uniq_pos