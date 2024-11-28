# MultiplexNDD

Identifying mechanisms and modulators of proteotoxicity in neurodegeneration through multiplexed gene-gene interaction screens.

All code necessary to reproduce the results of our manuscript are available as scripts in this repository.

Instructions and code for preprocessing the raw sequencing reads are found in the file preproc_readme.md

Installation
Get the source
git clone https://github.com/alejandrochavezlab/MultiplexNDD.git
The repository is fewer than 100Mb, so installation generally takes a few seconds depending on internet speed.

Environment and Dependencies
Execution is OS agnostic and does not require special hardware.

All code was run with R Version 4.0.2 using the following dependencies (versions are parenthesized):

dplyr
Rsamtools
reshape2
tibble
tidyr
ggplot2
stringr
metap
matrixStats
stringr

Usage
To preprocess the raw sequencing reads, please see the instructions in preproc_readme.md

To run the hit caller, execute the following from within a local version of R studio: python multiplexNDD.R
	- This will require user provided inputs. These inputs have been specified within the script.

License
See the LICENSE file
