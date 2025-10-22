# ComoUnGdgt

Suite of functions to work, manipulate, and apply GDGTs for paleoclimate applications. The package is optimized to take GDGT data as HPLC-MS peak areas, but can handle also fractional abundances. While in principle it can handle concentration values, this is not recommended and should be avoided.

For this package it is important to have the source data with the appropriate format (samples in rows, variables in columns), and GDGT variables named as follow (order is irrelevant).

1) brGDGTs: IIIa5, IIIa6, IIIa56, IIIa7, IIIb5, IIIb6, IIIc5, IIIc6, IIa5, IIa6, IIa7, IIb5, IIb6, IIb7 IIc5, IIc6, Ia, Ib, Ic
2) isoGDGTs: GDGT0, GDGT1, GDGT2, GDGT3, Cren, Creni (Cren')
3) isoGMGTs: GMGT0, GMGT1, GMGT2, GMGT3, GMGT4
4) OHGDGTs: GDGTX0, GDGTX1, GDGTX2, GDGTX3

(GMGTs and OHGDGTs are recognized but not indices are available yet)

The package has three main functionalities:

1) Calculate directly indices
2) Generate fractional abundances from the data using different approaches
3) Generate environmental reconstructions directly from the peak areas/fractional abundances

## Calculate Indices

A range of the most common GDGT indices have been included in the package (MBT'5Me, CBT', BIT, TEX86, etc.). All calculations by default require the presence of all the structures in the original formulation of the index by default. The parameter 'complete' can be manually set to FALSE to override this, and the formula for the index will be modified to remove all missing compounds. A warning will be issued listing all the compounds removed from the formula. By default the function won't ignore NA values in the data, but this can be overritten by turning 'na.ignore' as TRUE, in this case, any NA values in the data will be turned to 0.

### List of currently supported indices

1) BITx: BIT index
2) CBTp: CBT'
3) ComIdx: Community Index (CI)
4) DegCyc: Degrees of Cylisation (DC)
5) DegCycp: DC'
6) G2G3: GDGT2/GDGT3
7) IBT: Isomerization of Branchet tetraethers
8) IsoRat: Isomer Ratio 6-methyl
9) IsoRatsevn: Isomer Ratio 7-methyl
10) IsoRatSS: Isomer 6 + 7Me Ratio
11) MBT5: MBT'5Me
12) MBT6: MBT'6Me
13) MetIdx: Methane Index
14) pGDGT0: %GDGT0
15) pGDGT2: %GDGT2
16) tG0C: GDGT0/Cren
17) Rib: Ratio isobrenoid/branched
18) RinIdx: Ring Index
19) RinTet: #ringstetra
20) S3a2a: Sum IIIa/ Sum IIa
21) TEX86: TEX86

BayMBT:
1) bayrmbt_forward
2) bayrmbt_model
3) bayrmbt_predict

BaySPAR
1) bayspR_tex
2) bayspR_tex_analog
3) TEX_forward


## Generate fractional abundances

The function FracA will take a dataset of peak areas with the correct format and it can calculate the fractional abundance of one or multiple groups of GDGTs. Supported groups are "br", "iso", "br_extended" (5/6 isomer and 7-methyl), "isoGMGTs", and "OHGDGTs". The function by default calculates the fractional abundance of each group independently, but it can also calculate the fractional abundance relative to the sum of all the selected GDGTs.

The function Calculate fractional abundances for the structural sets defined in Raberg et al. (2021). This function depends on the formatting of the original data using the format_GDGT_sets, then calculating each set with calculate_single_set and presenting all the fractional abundances of all sets in a long format table.

## Generate environmental reconstructions

Three options for environmental reconstructions are integrated in the package (for v.0.1.0): BaySPR, BAYMBT, and linearCalib. BaySPR and BayMBT are just the application in R of the previously published models (https://github.com/jesstierney/BAYSPAR and https://github.com/jesstierney/BayMBT). On the other hand linearCalib has the linear equations to convert MBT'5Me and fractional abundances to temperature using the equations presented in the Suplementary table 1 of [[Review paper]]. Both the bayesian and linear calibrations are modified to take the dataset with peak areas or fractional abundances as input and directly output a reconstruction.

# Example

## Install package

```
library(devtools) #Make sure that the devtools library is loaded

options(timeout=400) #The package is a little heavy currently. Just making sure it doesn't time out.

install_github("Martoxa/ComoUnGDGT")

example #subset of 100 samples from the dataset from Martinez-Sosa, 2023
```
## Fractional abundances
```
FAexample<-FracA(example,group = c("iso","br"),how="each") #calculates the fractional abundance of br and iso GDGTs separately (default)

Setexample<- calculate_FAs_by_set(example,set_type = "Meth") #calculates the fractional abundance of the Meth set according to Raberg et al. (2021)
```

## Indices
```
exampleMBT<-MBT5(example);

exampleCBT<-CBTp(example,complete=FALSE) #missing IIIc6 in example dataset. Warning of the modified formula is issued.

exampleCBTfa<-CBTp(FAexample) #Indices can also be calculated from fractional abundances. Fractional abundance calculation can also import missing variables as 0s. A warning is issued.

```
## Calibrations
```
linearCalib(example,calibration = 1)

bayrmbt_predict(example[example$Type=="L",],10,10,Tmodel = "T0",Type = "lake")

out<-bayspR_tex(example[c(63,71),],example[c(63,71),]$Longitude,example[c(63,71),]$Latitude,6,runname = "SST") #some marine samples as example
out
```
