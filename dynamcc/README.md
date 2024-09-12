# DYNAMCC 0


Usage:

## In command line:
```shell
$ python -m dynamcc.dynamcc_0  -i A,C,G
Exploded codons: {'GCC': ['GCC'], 'KGC': ['GGC', 'TGC']}
DG codon: GCC
codon: GCC, rank: ('2', 0.26, 'A')
DG codon: KGC
codon: GGC, rank: ('1', 0.37, 'G')
codon: TGC, rank: ('1', 0.54, 'C')
```

## calling functions
```python
from dynamcc.dynamcc_0 import get_dg_codon_dict

get_dg_codon_dict(['A', 'G', 'C'], 'keep', 'rank', 2)
```
