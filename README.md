# DGC-LIBGEN: Design degenerate codon libraries from Machine learning predictions

## Install

```commandline
git clone git@github.com:charlesxu90/dgc-libgen.git
cd dgc-libgen/
pip install -e .
```

## Examples
See [`examples`](examples) directory for examples.
See [`example.ipynb`](example.ipynb) to see how to execute the example.

Simple usage as follows:
```python

wt_seq = 'AQSVPWGISRVQAPAAHNRGLTGSGVKVAVLDTGISTHPDLNIRGGASFVPGEPSTQDGNGHGTHVAGTIAALNNSIGVLGVAPSAELYAVKVLGASGSGSVSSIAQGLEWAGNNGMHVANLSLGSPSPSATLEQAVNSATSRGVLVVAASGNSGAGSISYPARYANAMAVGATDQNNNRASFSQYGAGLDIVAPGVNVQSTYPGSTYASLNGTSMATPHVAGAAALVKQKNPSWSNVQIRNHLKNTATSLGSTNLYGSGLVNAEAATR'

data_path = 'example/run1_top10k_data_pred.csv'

generate_library_for_topn(data_path, ref_seq=wt_seq, wt=wt_seq, 
                        pred_cols=['raw_activity', 'pred_esm2_t36_3B_UR50D', 'pred_ProtT5-XL', 'pred_ProtAlbert'],
                        seq_col='aa_seqs',
                        topn=1000, fix_threshold=0.75, aa_freq_threshold=0.1)
```
Note: You can adjust `fix_threshold` and `aa_freq_threshold` to fit your own library size.

## License
This code is licensed under [MIT License](./LICENSE.txt).