# IRAA: A statistical tool for investigating a protein-protein interaction interface from multiple structures

It is a tool to assess protein interface of a any complex [AB] with components [A] and [B].
Multiple 3D structures of complex [AB] in the bound state, as well as multiple 3D structrues of the components [A] and [B] are combined together to identify and understand the properties of the interface residues.

*A GUI within ipynb*

![A GUI within ipynb](https://github.com/kastritislab/IRAA/blob/main/IRAA_GUI_within_notebook.png)


Github repos is structure as below:

Notebook -> Main file: jupyter notebook 
For a quick sneak peek the notebook file is converted to .html and .pdf format only for viewing purposes. To actually interact and run the analysis please open the Run_IRAA.ipynb and go through it cell-by-cell. At the beginning there is a simple GUI panel created within the notebook for easy interaction.

* iraa_utils -> contains scripts to offer some functionality ot the notebook
* figures -> all output figure will be saved here.
* data -> Contains downloaded mmCIF files data as well as processed SASA data. 
* jpnb_overview_files -> Here are four .txt files, each containing lists of PDB IDs of structures of A, B, AB, and allow skipping any files a file \_drop.txt. 

Dependencies:
- numpy
- scipy
- pandas
- biopython
- tqdm
- pylab
- seaborn
- freesasa

You need to install freesasa and make sure it is added to the path, and available in your python run session.