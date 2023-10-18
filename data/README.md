# Mosquito aging modulates the development, virulence and transmission potential of pathogens.

[https://doi.org/10.5061/dryad.n02v6wx3n](https://doi.org/10.5061/dryad.n02v6wx3n)

5 datasets + 1 R script corresponding to the analysis of the data

## Description of the data and file structure

analysis of the effect of age on oocyst prevalence and intensity in An. coluzzii ####
"fig2\_oocyst\_data.txt" has 4 columns:
\# 	mosquito\_id : a unique code for each dissected mosquito \(n= 605\)
\# 	age\_class : a 3\-level categorical variable corresponding to each of the three age group \(12\, 8 and 4\-day\-old\)
\# 	isolate : a 4\-level categorical variable corresponding to each of the four gametocyte carriers \(replicates\)
\# oocyst : number of developed oocyst in each mosquito midgut

analysis of the effect of age on sporozoite prevalence and intensity in An. coluzzii
"fig2\_sporozoite\_data.txt" has 5 columns:
\# 	mosquito\_id : a unique code for each dissected mosquito \(n= 729\)
\# 	age\_class : a 3\-level categorical variable corresponding to each of the three age group \(12\, 8 and 4\-day\-old\)
\# 	isolate : a 4\-level categorical variable corresponding to each of the four gametocyte carriers \(replicates\)
\# Ct : number of cycle during the qPCR \(this is a proxy of the quantity of parasite DNA in salivary glands\)
\# NB: this variable includes NA values which mean that no amplification occured and hence that the sample is negative
\# positive: a binary variable corresponding to the infection status of mosquitoes \(1: presence of\, 0: absence of sporozoite\)

<br>
analysis of the effect of age and infection on mosquito survival
"fig3\_survival\_data.txt" has 8 columns:
\# 	mosquito\_id : a unique code for each tracked mosquito \(n= 657\)
\# cup: a code for each paper cup containing mosquitoes \(n=48 cups of \~15 mosquitoes \(range = 4\-18\)\)
\# The naming convention for each cup follows a structured pattern: the first two letters represent one of 
\# the three age groups \(G1\, G2\, or G3\)\, the third letter indicates the ID of the parasite isolate \(A to D\)\, 
\# the fourth and fifth letters signify whether the mosquitoes were exposed \(yes\) or not \(no\) to the infectious blood meal\,
\# and the final letter specifies the cup number \(1 or 2\) for each group\, isolate\, exposure status\.
\# 	age\_class : a 3\-level categorical variable corresponding to each of the three age group \(12\, 8 and 4\-day\-old\)
\# 	isolate : a 4\-level categorical variable corresponding to each of the four gametocyte carriers \(replicates\) 
\# exposure : a two\-level categorical variable corresponding to whether mosquitoes received an infectious \(yes\) or uninfectious blood\-meal \(no\)
\# infection\_status: a 3\-level categorical variable corresponding to whether mosquitoes were uninfected controls\, became infected upon exposure or remained uninfected upon exposure
\# daysPI: the day \(time post\-infection\) of mosquito death
\# censor: censoring indicator \(1 indicates that the response is a time at death\, 0 indicates that the individual was alive when last seen\)\. because all mosquitoes were followed until death\, this indicator is 1 for every mosquito

<br>
analysis of the effect of age on the parasite's EIP
fig4\_EIP\_data.txt" has 7 columns:
\# 	mosquito\_id : a unique code for each dissected mosquito \(n= 802\)
\# 	age\_class : a 2\-level categorical variable corresponding to each of the two age group \(12 and 4\-day\-old\)
\# 	isolate : a 3\-level categorical variable corresponding to each of the three gametocyte carriers \(replicates\)
\# dpi : day \(time post\-infection\) at which mosquitoes were dissected
\# oocyst : number of developed oocyst in each mosquito midgut
\# broken\_oocyst: number of ruptured oocysts \(NB: there are some NA values for dpi 6 because broken oocysts are unlikely 
\# to occur at this timepoint and were therefore not recorded\)
\# spz: a binomial variable corresponding to the presence \(1\) or absence \(0\) of sporozoite in mosquito salivary glands
\# \(NB: there are some NA values for dpi 6 because sporozoites are unlikely to occur at this timepoint and were 
\# therefore not recorded\. A few additional NA are observed beyond dpi 6 when mosquito salivary glands were lost during the dissection\)

analysis of the effect of age on vectorial capacity
"fig5\_vectorial\_capacity\_data.txt" has 3 columns:
\# 	scenario : a 3\-level categorical variable corresponding to each scenario
\# 	age\_class : a 2\-level categorical variable corresponding to each of the two age group \(12 and 4\-day\-old\)
\# vectorial\_capacity: the response variable simulated from the mathematical model

<br>
## Sharing/Access information

## Code/Software

1 R script "Rscript\_mosquito\_age\_infection" corresponding to the analysis of the 5 data set presented above