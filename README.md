# RFDistributionR
This code implements an improved version of the RF distribution computation by Bryant et al.
Here, we modified the dynamic programming algorithm introduced by Bryant et al for computing the distribution of RF distance 
for a given tree by leveraging the Number-Theoretic Transform (NTT), and improve the running time from O(l^5) to O(l^3log(l)), 
where l is the number of tips of the tree.
Given an unrooted phylogenetic tree T with l tips, the procedure for computing the RF distribution of this tree is as follows:
Denote the node adjacent to tip l in T by v_0. Remove tip l, and root the resulting tree with v_0 as the root. We use this 
rooted tree as the input to the dynamic programming algorithm.

## Installation

To install the R package, open a terminal and type:<br><br>
```
sudo add-apt-repository -y ppa:cran/poppler
sudo apt-get update
sudo apt-get install -y libpoppler-cpp-dev
sudo apt-get install libcurl4-openssl-dev libxml2-dev
R
install.packages("devtools")
library("devtools")
install_github("WGS-TB/RFDistributionR")
```

## Usage Examples

```
library("rfdistr")
library("ape")
rfdistr::ntt_polynomial(rtree(5),6)
```


```
library("rfdistr")
library("ape")
rfdistr::polynomial(rtree(5),6)
```
If you are interested in the Python version instead, please visit https://github.com/WGS-TB/RFDistribution.
