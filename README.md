# newCAR
Immunotherapy relies heavily upon the identification of a tumor-associated antigens (TAAs) which are highly expressed and tumor-specific.  To identify suitable TAAs, newCAR.R compares RNA and protein expression between tumor and control samples. Users must provide data matrices for tumor and control to determine TAAs in their disease of interest. Tumor and control datasets can be from different data cohorts, or even different platforms (e.g. RNA-seq vs microarray), since gene expression value transformations are performed to make the distribution of expression values equivalent.

newCAR allows for reproducibility of the workflow that we used to predict tumor-associated antigens in pediatric AMKL.  A function to perform the Weighted Cumulative Percentage (WCP) data transformation is also available in this repository to unify heterogenous data.

## Installation
Install `newCAR.R` by unpacking the source code to a working directory.

### Prerequisites

* [R](https://www.r-project.org/) ^3.5.2
    * [Bioconductor](https://bioconductor.org/) ~3.8
    * [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) ~1.16.6
    * [preprocessCore](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html) ~1.44.0
    * [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.4) ~1.8.4

## Usage
```bash
Rscript newCAR.R config.txt
```

### Input
#### Configuration File
The configuration file should contain 9 lines in the following order:

| Value                                | Description                                                                          |
| ------------------------------------ | ------------------------------------------------------------------------------------ |
| 1. Disease                           | [Example Disease Matrix](./example_matrices/disease_data_example.txt)        |
| 2. Control                           | [Example Control Matrix](./example_matrices/control_data_example.txt)        |
| 3. Microarray Expression             | [Example Microarray Expression Matrix](./example_matrices/microarray_expression_example.txt)        |
| 4. Microarray Detection (Optional)   | [Example Microarray Detection Matrix](./example_matrices/microarray_detection_example.txt)            |
| 5. [Gene Membrane Association File](./reference_information/membrane_association.txt)    | Information regarding whether or not a gene has been associated with the plasma membrane |
| 6. [Cancer-Testis Gene Reference File](./reference_information/cancer_testis.txt) | Information regarding whether or not a gene is a known cancer-testis antigen    |
| 7. [CD Family Identification File](./reference_information/cd_family.txt)     | Information regarding whether or not a gene is a member of the CD family          |
| 8. [Protein Expression Data](./reference_information/protein.txt)           | Protein expression information (from the [HPA](https://https://www.proteinatlas.org/about/download))        |
| 9. Name of Run                       | example_name                                                                          |

An configuration file [example](example_config.txt) is provided.

### Output
The name of the `newCAR.R` run (line 9 of the configuration file) is used to create an output directory.  The output directory contains *data matrices* and *boxplots* regarding the following TAAs:

| Output                        | Description                                                                                                            |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `Disease.TEVs.txt`               | Disease transformed expression values used for TAA prediction |
| `Control.TEVs.txt`               | Control transformed expression values used for TAA prediction |
| `$name.Detection_Calling_Summary.txt` | If microarray detection information is provided, this file contains the number of samples for which each gene was predicted to be present or absent |
| `$name.AllMetrics.txt`           | Reports all genes considered and their corresponding metrics for TAA prediction                                |
| `$name.MembraneGenes.txt`     | Reports all genes considered and their corresponding metrics for TAA prediction for membrane-associated genes (as defined in Ensembl)                              |
| `$name.MembraneGenes_p05.txt`     | Reports all genes considered and their corresponding metrics for TAA prediction for membrane-associated genes (as defined in Ensembl) with a p-value < 0.05                             |
| `$name.CD_Family.txt`            | TAAs in the CD family                                                                                                      |
| `$name.CTgenes.txt`              | TAAs in the Cancer-testis family                                                                                                  |                                                                          |
| `$name.TopTargets.txt`           | TAAs that have met the criteria described in our workflow that are ready to be subjected to manual review                                             |

## Maintainers
[Patrick Schreiner](https://github.com/pschreiner) and [Yiping Fan](https://www.stjude.org/directory/f/yiping-fan.html)
