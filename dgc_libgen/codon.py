
from .dynamcc.dynamcc_0 import get_dg_codon_dict, load_ecoli_codons, get_aa_top1_codon


sorted_dict = load_ecoli_codons()


def get_codons_from_aa(df_row):
    """Get the codons for a given amino acid or list of amino acids.
    Args:
        df_row (pd.Series): A row from the dataframe containing amino acid information.
    
    Returns:
        list: A list of codons for the given amino acid or list of amino acids.
        list: A list of expanded codons for the given amino acid or list of amino acids.
    """
    if df_row.fix_aa: 
        # print(df_row.fix_aa)
        codon, _ = get_aa_top1_codon(df_row.fix_aa, sorted_dict)
        return [codon], [codon]
    else:
        # print(df_row.mut_aa)
        codon_dict, _ = get_dg_codon_dict(df_row.sel_mut_aa)
        return list(codon_dict.keys()), [c for l in  list(codon_dict.values()) for c in l]