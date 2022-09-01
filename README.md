# explicit-mouse

**Identify TF regulators for mouse gene modules.**

The EXPLICIT-Mouse package is developed in companion with the [EXPLICIT-Human](https://github.com/MaShisongLab/explicit-human) package. Both packages are modified from the [EXPLICIT](https://github.com/MaShisongLab/explicit) package we published previously ([Geng *et al.* 2021](https://github.com/MaShisongLab/explicit-mouse#References)). EXPLICIT-Mouse uses the expression of 1,593 mouse TF genes to predict the expression of 23,481 mouse non-TF genes. It also predicts the expression of every TF gene using the expression of other 1,592 TF genes. It further infers TF regulators for mouse gene modules. For more details, please refer to [Wang et al, 2022](https://github.com/MaShisongLab/explicit-mouse#References) and the [EXPLICIT-Human](https://github.com/MaShisongLab/explicit-human) package.

## Table of Contents
- [Install](https://github.com/MaShisongLab/explicit-mouse#Install)
- [Usage](https://github.com/MaShisongLab/explicit-mouse#Usage)
   - Identify TF regulators for mouse gene modules
   - Visualize TF-target genes interactions using Chord diagrams
   - Build the EXPLICIT-Mouse predictor model
- [References](https://github.com/MaShisongLab/explicit-mouse#References)

## Install
[Perl](https://www.activestate.com/products/perl/downloads/), [R](https://www.r-project.org/), and the [circlize](https://www.rdocumentation.org/packages/circlize/) package in R are required. [MATLAB](https://www.mathworks.com/products/matlab.html) is required only for building gene expression predictors from gene expression matrices. Once these softwares are installed, just download the package and start using it.

## Usage

### 1. Infer TF regulators for mouse gene modules

#### (1). Prepare the module file
The file `modules_to_analyze.txt` is used to store gene module information. It has two tab-separated columns for gene ids and module names, respectively. Only standard mouse Ensembl gene ids are supported. The `modules_to_analyze.txt` file came with the package can be used as a test case. The file is preloaded with mouse gene modules whose genes are homologues of the human genes within the human gene co-expression modules as described in [Wang *et al*, 2022](https://github.com/MaShisongLab/explicit-mouse#References).  
```shell
Gene_Name   ModuleID
ENSMUSG00000029831	Module0037
ENSMUSG00000041534	Module0037
ENSMUSG00000023982	Module0037
ENSMUSG00000035270	Module0037
ENSMUSG00000074991	Module0037
ENSMUSG00000032343	Module0037
ENSMUSG00000049353	Module0037
ENSMUSG00000031293	Module0037
ENSMUSG00000040554	Module0037
ENSMUSG00000044375	Module0037
ENSMUSG00000023979	Module0037
ENSMUSG00000047034	Module0037
.........   ........
```
#### (2). Identify potential TF regulators for the modules
The Perl script `getMouseRegulatorTFs.pl` takes the modules from `modules_to_analyze.txt`, performs enrichment assays to identify potential TF regulators, and saves the results to a file `results.regulator.tfs.txt`. The file can be opened and previewed in EXCEL, which lists the potential TF regulators for the input modules. 
```shell
perl getMouseRegulatorTFs.pl
``` 

### 2. Visualize TF-target genes interactions using Chord diagrams

The `getChordDiagram` function in R draws Chord diagrams to visualize the interactions between TF regulators and their target genes. The function extracts TF-target gene pairs from `results.regulator.tfs.txt` for a input module and draws a chord Diagram accordingly. <br><br>
`getChordDiagram('module', ratio, tfnum, targetnum)` <br/>
`module` - module name <br>`ratio` - relative size of the target gene area <br>`tfnum` - maximum number of TF genes to be included in the diagram <br>`targetnum` - maximum number of target genes to be included in the diagram<br>

Open an R console, change the working directory to the home directory of the EXPLICIT-Mouse package, and type in the following commands:
```R
# R code
# Load the scripts that define the getChordDiagram function.
source("Rscripts.R")  

# The 'circlize' package is required.
library("circlize")

# Draw a Chord diagram for Module0037
getChordDiagram( module="Module0037", ratio = 1, tfnum = 50, targetnum = 15)

# Change the relative size of the target gene area
getChordDiagram( module="Module0037", ratio = 0.5, tfnum = 50, targetnum = 15)

# Draw Chord diagrams for other modules
getChordDiagram( module="Module0018", ratio = 1, tfnum = 50, targetnum = 15)
getChordDiagram( module="Module0057", ratio = 1, tfnum = 25, targetnum = 20)
getChordDiagram( module="Module0085", ratio = 1, tfnum = 20, targetnum = 15)
```

### 3. Build the EXPLICI-Mouse predictor model

The following steps explain how we built the EXPLICIT-Mouse predictor model. One can use a similar process to build a custom predictor model using custom gene expression matrix. For more details on the process, please refer to the EXPLICIT and EXPLICIT-Human packages. 

We used a mouse gene expression matrix extracted from the ARCHS4 database ([Lachmann *et al*, 2018](https://github.com/MaShisongLab/explicit-mouse#References)) to build the EXPLICIT-Mouse model. A mouse gene expression file `mouse_transcript_v7.h5` was downloaded from [the ARCHS4 database](https://maayanlab.cloud/archs4/download.html) and processed into CPM gene expression values. The expression values are then log-transformed via log<sub>2</sub>(CPM + 1). After quality control, 25,074 mouse genes in 54,634 bulk RNA-seq samples are selected to generate a mouse gene expression matrix, contained within a file `mouse_expression_extracted_from_archs4_v7.h5`, which is available via [Figshare](https://figshare.com/s/ec58e5b149c3060e1a6f). 

```shell
mouse_expression_extracted_from_archs4_v7.h5
 ├─ expression  					(54,634 samples [row] X 25,074 genes [column])
 ├─ gene_name						(25,074 genes)
 ├─ sample_geo_accession			(54,634 samples)
 ├─ idx_tf_gene						(specify TF genes used for model construction)
 └─ idx_non_tf_gene 				(specify non-TF genes used for model construction)
```

Download the file `mouse_expression_extracted_from_archs4_v7.h5` from [Figshare](https://figshare.com/s/ec58e5b149c3060e1a6f) and place it in the home directory of the EXPLICIT-Mouse package. Build the EXPLICIT-Mouse model using the following commands within MATLAB.

```matlab
% MATLAB code
% Navigate to and start within the home directory of the EXPLICIT-Mouse package.

%%%%%%%%%%%%%
% Step A. Build a predictor for non-TF genes
% Obtain the expression matrix and gene names.
mtx = h5read('mouse_expression_extracted_from_archs4_v7.h5','/expression');
gene_name = h5read('mouse_expression_extracted_from_archs4_v7.h5','/gene_name');

% idx_tf specifies which of the 25,074 genes are TFs to be used.  1,593 TFs are selected in total.
% idx_non_tf specifies which of the 25,074 genes are TFs to be used.  23,481 non-TFs are selected in total.
idx_tf = h5read('mouse_expression_extracted_from_archs4_v7.h5','/idx_tf_gene') == 1;
idx_non_tf = h5read('mouse_expression_extracted_from_archs4_v7.h5','/idx_non_tf_gene') == 1;

% Obtain the TF expression matrix and non_TF expression matrix
tf_mtx = mtx(:,idx_tf);
non_tf_mtx = mtx(:,idx_non_tf);

% Obtain the TF gene names and non_TF gene names
tf_name = gene_name(idx_tf);
non_tf_name = gene_name(idx_non_tf);

% Build the predictor model
mdl = explicit( tf_mtx, non_tf_mtx, tf_name, non_tf_name);

% Inspect the model
mdl

% The first 5 significant TF-target gene pairs
mdl.SigEdges(1:5,:)

%%%%%%%%%%%%%
% Step B. Build a predictor for TF genes
% Note that this step might take a long time to run.
mdl_tfxtf = explicit_tfxtf(tf_mtx, tf_name);

% Inspect the model
mdl_tfxtf

%%%%%%%%%%%%%
% Step C. Combine the SigEdges from both models
% Combine the SigEdges from the two models to obtain all SigEdges for the EXPLICIT-Mouse model. 
% They are the SigEdges used in Step 1 to infer TF regulators for gene modules. 
% Note that the file Hs.SigEdges.1e-12.txt are split into two files: Mouse.SigEdges.part1.txt and 
% Mouse.SigEdges.part2.txt so that their file sizes are less than 100M. These two files are 
% contained within the data folder.
i = mdl.SigEdges{:,4} <= 1e-12;
i2 = mdl_tfxtf.SigEdges{:,4} <= 1e-12;
writetable([mdl.SigEdges(i,:) ; mdl_tfxtf.SigEdges(i2,:)],'Mouse.SigEdges.1e-12.txt','Delimiter','tab')
```

## References

Geng H, Wang M, Gong J, Xu Y, and Ma S. 2021. An Arabidopsis expression predictor enables inference of transcriptional regulators for gene modules. Plant J 107:597-612.

Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, and Ma'ayan A. 2018. Massive mining of publicly available RNA-seq data from human and mouse. Nature Communications 9:1366.

Wang Y, Zhang Y, Yu N, Li B, Gong J, Mei Y, Bao J, Ma S. 2022. Decoding transcriptional regulation via a human gene expression predictor. submitted 

