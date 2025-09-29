
NAT_AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def get_mutations(aa_seq, wt_seq):
    muts = []
    for i in range(len(wt_seq)):
        if wt_seq[i] != aa_seq[i]:
            muts.append(f"{wt_seq[i]}{i+1}{aa_seq[i]}")
    return '/'.join(muts)