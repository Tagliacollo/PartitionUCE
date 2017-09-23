# Estimating improved partitioning schemes for UltraConserved Elements (UCEs)

Most phylogenomic studies of UCEs use partitioning methods to account for variation in rates and patterns of molecular evolution among sites. The accuracy of phylogenetic inferences often depends on choosing an appropriate model of molecular evolution, and many studies have demonstrated that accounting for variation in rates and patterns of evolution among sites is of primary importance. Given that the rates of molecular evolution within each UCE are known to vary predictably, with near-zero variation in the conserved core of the UCE increasing to large amounts of variation at the edges, we developed new partitioning methods for phylogenomic studies of UCE to improve model-fit and parameter estimates. Here, we present the scripts of the two methods: **_Sliding-Window Site Characteristics_** and **_UCE Site Position_**. 
 
**Method 1: _Sliding-Window Site Characteristics_ (SWSC)**

 This first partitioning method uses a sliding-window approach and site characteristics such as entropy (SWSC-EN), GC content (SWSC-GC), and multinomial (SWSC-MU) to automatically select partitions that account for heterogeneity in rates and patterns of molecular evolution within UCE alignments. The SWSC methods are graphically summarized below:

![SWSC](/Tables-Figures/Figures/Figure2.png)

The diagram above includes three hypothetical UCE alignments: UCE-1 (green), UCE-2 (yellow), and UCE-3 (blue). The alignments include the same patterns of molecular evolution (i.e. entropies, GC content, and multinomial rates) as seen in the UCE markers. The SWSC algorithm includes three steps. First, it proposes all combinations of three-window models in the alignments. It delimitates windows by locating every conceivable pairs of nucleotide sites in the alignments. Second, it estimates site-wise molecular evolution and nucleotide proportions/counts across the alignments; in case, it uses three alternatives properties of UCEs: (i) Entropy (SWSC-EN), (ii) GC content (SWSC-GC), and (iii) multinomial likelihoods (SWSC-MU). Third, it calculates the Sum of Squared Errors (SSE) statistics across window and sum them up to obtain the sum of SSEs for every three data block model. The SWSC algorithm selects the best models by minimizing the squared residuals among three windows. 

**Method 2: _UCE Site Position_ (UCESP)**
The UCESP partitioning method groups together sites that are found in similar physical locations within the UCEs. This method relies on the assumption the UCE regions (e.g. cores and flacking regions) have evolved under similar evolutionary processes. The UCESP method is graphically summarized below.

![UCESP](/Tables-Figures/Figures/Figure3.png)

 The diagram above includes three hypothetical alignments: UCE-1 (green), UCE-2 (yellow), and UCE-3 (blue). The positions of nucleotide sites in the alignments are indicated by the numbers of each site. The UCESP algorithm includes two steps. First, it defines the central site of the UCE as site 0 (zero) and sites to the left of this site are labelled with negative numbers, and sites to the right of this site are labelled with positive numbers. Second, it creates one data block for each label, by combining all the sites with the same label into a single data block. Note that there are nine partitions, all sharing the same locations of sites among the UCEs. This method is rate-free; i.e. it does not require estimating proxies of molecular evolution

## Input

The input of both methods is a concatenated nexus alignment (.nex) comprised of UCE markers and including charsets with the locations of the UCEs (example inputs are in the `/raw_data` folder of this repository) 

## Output Files

The output of both methods is a PartitionFinder configuration file (.cfg) to be used as the input file for [PartitionFinder 2](https://academic.oup.com/mbe/article/34/3/772/2738784/PartitionFinder-2-New-Methods-for-Selecting). The SWSC methods also outputs a csv file (.csv) with values of entropy (SWSC-EN), GC content (SWSC-GC), and multinomial (SWSC-EN) for each site of the UCEs. (example outputs are in the `/processed_data/UCE-subsets` folder of this repository). 

## Python version

Python 3.6.x or higher

> Required Python modules: 

> - [Bio](https://pypi.python.org/pypi/biopython)
> - [collections](https://docs.python.org/3/library/collections.html)
> - [io](https://docs.python.org/3/library/io.html)
> - [itertools](https://docs.python.org/3/library/itertools.html) 
> - [math](https://docs.python.org/3/library/math.html) 
> - [numpy](https://pypi.python.org/pypi/numpy)
> - [os](https://docs.python.org/3/library/os.html)
> - [pathlib2](https://pypi.python.org/pypi/pathlib2/)
> - [subprocess](https://docs.python.org/3/library/subprocess.html)
> - [time](https://docs.python.org/3/library/time.html) 
> - [tqdm](https://pypi.python.org/pypi/tqdm)


## How to run the two new partitioning methods?

1- Open Terminal (on most Macs, this is found in Applications/Utilities)
2- Type “python“ followed by a space (remember, it needs to be python 3.6.x or higher)
3- Drag and drop the `py_script/analysis.py` file onto the command prompt.
4- Type another space
5- Drag and drop the `example/example_dataset.nex` folder onto the command prompt
6- Hit Enter/Return to run



