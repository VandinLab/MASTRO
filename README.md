# MASTRO: Discovering Significant Evolutionary Trajectories in Cancer Phylogenies

This repository contains the implementation of MASTRO, our algorithm to discover statistically significant conserved evolutionary trajectories of alterations from phylogenetic tumor trees. All details on MASTRO are in the paper (see the `mastro.pdf` file).

Please address questions and bug reports to Leonardo Pellegrina (leonardo.pellegrina@unipd.it) and Fabio Vandin (fabio.vandin@unipd.it).

Our implementation of MASTRO leverages LCM (version 5.3 by Takeaki Uno http://research.nii.ac.jp/~uno/codes.htm ) to mine frequent itemsets.
The folder `/MASTRO/` contains the source files for the implementation of MASTRO and for reproducing the experiments described in the paper. The folder `/data/` contains the data analyzed by MASTRO (and scripts to preprocess it).

To run MASTRO use the python script `run_MASTRO.py`. Such script accepts various parameters, listed and described by using the `-h` flag:

```
usage: run_MASTRO.py [-h] [-g G] [-p P] [-s S] [-minp MINP]

optional arguments:
  -h, --help  show this help message and exit
  -g G        input file with graphs
  -p P        permutation type: 0 = indipendent, 1 = permutation, 2 = ind. in random topology (def=0)
  -s S        minimum support of trajectories (def=2)
  -minp MINP  path to append minimum pvalue (optional)
```

The `-g` option specifies the path to the input file containing (see below for the input format).
The parameter `-s` allows to set the minimum support of the trajectories to find (the minimum number of tumor trees in which any trajectory must be observed to be frequent).
The flag `-p` allows to specify the type of statistical test to use for evaluating the significance of a frequent trajectory.
`0` is the statistical test in which alterations are assumed to be inserted uniformly and independently at random on the nodes, preserving the set of alterations of each patient and the topology of the tumor trees.
`1` is the statistical test assuming that alterations are randomly permuted (preserving the number of alterations in each node).
`2` is the statistical test that consider uniform and independent assignment of each alteration on a tree with random topology (that is sampled uniformly from the set of topologies of the cohort).

To correct for multiple hypothesis testing, MASTRO uses resampling-based procedures: the Westfall-Young permutation testing procedure to bound the Family-Wise Error Rate (FWER) and the Storey-Tibshirani method to estimate the False Discovery Rate (FDR).
The script `run_wy_correction.py` can be used to generate resampled datasets and running MASTRO on them.
The script `plot_results.py` generates various figures with the results, and also computes the corrected significance thresholds for bounding the FWER.
The script `plot_emp_FDR.py` estimates the FDR from the resampled datasets of the top-k results for various values of k.

### Input format
MASTRO takes in input a file containing, in each line, the edges of the complete expanded tumor graph of the tumor phylogenetic trees, separated by a space.
Alterations can be of any type (e.g., SNV or CNA) but must have a distict id in each tree (e.g., `TP53_SNV` and `TP53_DEL` are allowed in the same tumor tree).
The complete expanded tumor graph is composed by three types of edges:
1. a directed edge `A->-B` denotes that the alteration `A` is an anchestor of the alteration `B`
2. an undirected edge `A-/-B` denotes that `A` and `B` belong to different branches of the tree
3. an undirected edge `A-?-B` denotes that the order between `A` and `B` is not known (they belong to the same node)
Note that the last two edges are undirected, therefore the alterations are meant to be sorted alphabetically (i.e., `A-?-B` is correct, `B-?-A` is not).
(it is not necessary to include the edges incident to the germline node to find trajetories composed by at least 2 alterations; otherwise, use `g` to denote the germline node)

For example, the following file denotes the input of MASTRO of three tumor trees with edges `{g->-A, A->-B, B->-C}`, `{g->-A, A->-B, A->-C}`, `{g->-A, g->-B, A->-C}`:
```
A->-B A->-C B->-C
A->-B A->-C B-/-C
A-/-B A->-C B-/-C
```

### Output format
The output of MASTRO is a file containing a list of frequent trajectories. If the input file is `input.txt`, the trajectories are written in `input_final.csv`.
The trajectories are represented with the same format of the imput tumor trees (always including edges that are incident to the germline node).
