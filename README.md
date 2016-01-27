# Sequence_Transition_Scorer

REQUIRES BIOPYTHON

Simple command line program for comparing protein sequences and scoring their differences using a custom amino acid scale.

To run, open the commandline and type - `python sts.py PARAMETERS`

Required parameters:
* Input sequence file. Must be Fasta format.
* Output results location. Any directory to save results to.
* Scale file. Must be tab-delimited file with two columns. First column must be one of the 20 standard amino acids (Capital letter). Second column must be the integer score.

Optional parameters:
* Positions file. Single column of position number to limit comparison to. Assumes indexing starts at 1.
* Tree file. Must be in newick format. If provided, only parent-child transitions are looked at.
* Sequence names file. Two-column tab-delimited file, where each column must include sequence names for sequences in the input sequence file. Pairs of names provided here indicate transitions of interest and limits the comparisons to these sequences. If this file and a tree file is provided, only transitions in this file will be used.

If no tree or sequence name file is provided, all sequences will be compared to each other.
If no positions file is provided, all positions will be compared.

Output files:
* Score differences at all specified positions for all specified pairs of protein sequences
* Mean score differences for all specified pairs of protein sequences.

Additional Notes:
* Currently only positions that change in each transition are scored. So if sequence A and B vary by 10 amino acids, than their score will be based on these 10 positions. However if comparing sequence A and sequence C we observed the difference at only 5 locations, the score will be based on these 5 locations.
* Make sure sequence names in all files are identical.
* Any positions that are a gap are ignored for the specific transition they are observed in.