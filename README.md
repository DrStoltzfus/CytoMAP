# CytoMAP

Histo-Cytometric Multidimensional Analysis Pipeline (CytoMAP) 

## Introduction 

This program is meant to make advanced data analytic techniques
accessible for Histo-Cytometry data. The goal is to take established
analytical techniques, such as neural network based clustering
algorithms, and package them in a user friendly way, allowing researchers
to use them to explore complex datasets. For an in depth walkthrough, [check out the wiki.]( https://gitlab.com/gernerlab/cytomap/wikis/Home )

Authored by: 

Dr. Caleb Stoltzfus - Lead Author

Yajun Wu (Cherie)- Technical Assistant, September 2019 to June 2020, *University of Washington, Mechanical Engineering*

Jakub Filipek - Technical Assistant, October 2018 to August 2019, *University of Washington, Paul G. Allen School of Computer Science and Engineering*

## Citation

If you use CytoMAP please cite our Cell Reports paper:

**Stoltzfus, C. R., Filipek, J., Gern, B. H., Olin, B. E., Leal, J. M.,  Wu, Y., ... & Gerner, M. Y. (2020). CytoMAP: a spatial analysis  toolbox reveals features of myeloid cell organization in lymphoid  tissues. *Cell reports*, *31*(3), 107523.**

## Install 
You can download the necessary MATLAB files to run [CytoMAP here](https://gitlab.com/gernerlab/cytomap/-/archive/master/cytomap-master.zip?path=CytoMAP) 

To run CytoMAP in MATLAB you will also need to install MATLAB 2018b or later and the following toolboxes:
* Optimization Toolbox
* Deep Learning Toolbox
* Image Processing Toolbox
* Statistics and Machine Learning Toolbox
* Bioinformatics Toolbox

You can download an installer for the compiled version of [CytoMAP here](https://github.com/DrStoltzfus/CytoMAP/raw/main/StandaloneInstaller/CytoMAP_Installer_WindowsV1.4.21.exe)

For the compiled version there is less built in functionality in the figure window. However, you do not need to install a full version of MATLAB.

## Input:

### .csv files:

Any dataset structured as a .csv, with at least x-y positional data can be loaded into this program. The data should be formatted as Columns of parameters and rows of individual cells/objects. Individual .csv files will be treated as separate phenotypes. For best results when loading multiple phenotypes, each csv(phenotype) should have the same markers/channels and the first row of all csv files should have the same channel names in the same order. 

#### Example:
> If you use the 'load csv' function and select four .csv files, each containing columns of channel intensities 
> and rows of individual cell objects, CytoMAP will create one new sample with four cell types. The names of the cells will be the names of the .csv files. The name of the sample will be the common string among the four .csv file names.

### .wsp+.fcs files:

Under limited conditions .wsp files, with associated .fcs files can be directly uploaded.
* There must be at least two positional channels, i.e. X,Y, and/or Z.
* There can not be any ellipse or histogram gates.
* Some of the more advanced components of .wsp, like compensation matrices and derived parameters are supported, though not thoroughly tested and can break things.

### .mat Previous workspace files:
You can save a workspace as a .mat file, which can be loaded later. 
When saving a workspace, all data, neighborhoods, and trained clustering 
networks are saved. Any open figures are not saved.

> Always check your data after uploading to confirm everything worked properly.

## Troubleshooting

If you run into bugs, or have questions about features; post a question on the [Scientific Community Image Forum.](https://forum.image.sc/) Be sure to add the *cytomap* tag to your topic for visibility.

**Have fun exploring your data!**# CytoMAP
