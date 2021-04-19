# MS.variant.simulations
A RScript to simulate single variant GWAS in one population, and calculate correlation of the most significant variant in other populations.

# Input files
The script is intended to simulate a region of the genome that contains one causal variant. In order to generate causal variant all distributed along the genome this same script can be use iteratively across the different regions, which will results in independent loci. 

The program takes the path to the following files as input:
* Raw file with dosages, samples of the discovery population in the rows and variants of the region to simulate in the columns. The raw file can be made with Plink's flag --recodeA, click [here](http://zzz.bwh.harvard.edu/plink/dataman.shtml#recode) for more details on the format.
* 

