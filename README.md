# ComoUnGDGT

Suite of functions to work, manipulate, and apply GDGTs for paleoclimate applications. The package is optimized to take GDGT data as HPLC-MS peak areas, but can handle also fractional abundances. While in principle it can handle concentration values, this is not recommended and should be avoided.

For this package it is important to have the source data with the appropriate format (samples in rows, variables in columns), and GDGT variables named as follow (order is irrelevant).

brGDGTs: IIIa5, IIIa6, IIIa56, IIIa7, IIIb5, IIIb6, IIIc5, IIIc6, IIa5, IIa6, IIa7, IIb5, IIb6, IIb7 IIc5, IIc6, Ia, Ib, Ic isoGDGTs: GDGT0, GDGT1, GDGT2, GDGT3, Cren, Creni (Cren') isoGMGTs: GMGT0, GMGT1, GMGT2, GMGT3, GMGT4 OHGDGTs: GDGTX0, GDGTX1, GDGTX2, GDGTX3 (GMGTs and OHGDGTs are recognized but not indices are available yet)

The package has three main functionalities:

    Calculate directly indices
    Generate fractional abundances from the data
    Generate environmental reconstructions directly from the peak areas/fractional abundances

## Calculate Indices

A range of the most common GDGT indices have been included in the package (MBT'5Me, CBT', BIT, TEX86). All calculations by default require the presence of all the structures in the original formulation of the index by default. The parameter 'complete' can be manually set to FALSE to override this, and the formula for the index will be modified to remove all missing compounds. A warning will be issued listing all the compounds removed from the formula. By default the function won't ignore NA values in the data, but this can be overritten by turning 'na.ignore' as TRUE, in this case, any NA values in the data will be turned to 0.


## Generate fractional abundances

The function FracA will take a dataset of peak areas with the correct format and it can calculate the fractional abundance of one or multiple groups of GDGTs. Supported groups are "br", "iso", "br_extended" (5/6 isomer and 7-methyl), "isoGMGTs", and "OHGDGTs". The function by default calculates the fractional abundance of each group independently, but it can also calculate the fractional abundance relative to the sum of all the selected GDGTs.

##Generate environmental reconstructions

Three options for environmental reconstructions are integrated in the package (for v.0.1.0): BaySPR, BAYMBT, and linearCalib. BaySPR and BayMBT are just the application in R of the previously published models (https://github.com/jesstierney/BAYSPAR and https://github.com/jesstierney/BayMBT). On the other hand linearCalib has the linear equations to convert CBT' to pH, and MBT'5Me to temperature as published by DeJonge, 2014; Naafs, 2017; Russell, 2018; Martinez-Sosa, 2022 (as of v.0.1.0). Both the bayesian and linear calibrations are modified to take the dataset with peak areas or fractional abundances as input and directly output a reconstruction.
