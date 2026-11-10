# Input files

Put these three files in this folder. The notebook expects exactly these names.

**`query_proteins.faa`**: the proteins the students will search, as a protein FASTA file (for example, predicted proteins from a genome, a MAG, or a set of contigs).

**`core_gene.hmm`**: one profile HMM in HMMER3 format. If your HMM file holds several models, pull out the one you want by its NAME:

```
hmmfetch all_models.hmm PolB > core_gene.hmm
```

**`references.faa`**: reference sequences of the same gene from known organisms, as an unaligned protein FASTA file. The notebook aligns them together with the students' hits.

## Tips for the reference file

- The first word of each header becomes the label in the tree, so make it short, unique, and informative, for example `>Mimiviridae_APMV_PolB`. Anything after the first space is dropped.
- Use only letters, numbers, `_`, `.` and `-` in that first word. The notebook replaces any other character with `_`.
- Make sure no reference shares an ID with a protein in `query_proteins.faa`.

## Size

- Use at least 10 sequences in total (references plus hits). With 9 or fewer, the trimming setting in the notebook (`trimal -gt 0.1`) removes nothing, because a single residue is already more than 10% of a column.
- Phylogeny.fr accepts at most 200 sequences, each at most 2,000 amino acids long. For an uploaded alignment, keep the trimmed alignment under 2,000 columns.
