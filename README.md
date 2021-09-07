# newCAR
newCAR allows for reproducibility of the workflow that we used to predict tumor-associated antigens in pediatric AMKL.  A function to perform the Weighted Cumulative Percentage (WCP) data transformation is also available in this repository to unify heterogenous data.

Immunotherapy relies heavily upon the identification of an TAAs which are highly expressed and tumor-specific.  To identify suitable TAAs, newCAR.R compares RNA and protein expression between tumor and control samples.

Users must provide data matrices for tumor and control to determine TAAs. Tumor and control datasets can be from different data cohorts, or even different platforms (e.g. RNA-seq vs microarray), since gene expression value transformations are performed to make the distribution of expression values equivalent.

## Installation
Install `newCAR.R` by unpacking the source code to a working directory.

### Prerequisites

* [R](https://www.r-project.org/) ^3.5.2
    * [Bioconductor](https://bioconductor.org/) ~3.8
    * [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) ~1.16.6
    * [preprocessCore](https://bioconductor.org/packages/release/bioc/html/preprocessCore.html) ~1.44.0
    * [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.4) ~1.8.4

**Note**: Tumor and control expression values are necessary, but not provided in this package, in order to identify TAAs.

## Usage
```bash
Rscript newCAR.R config.txt
```

### Input
#### Configuration File
The configuration file should contain 9 lines in the following order:

| Value                                | Description                                                                          |
| ------------------------------------ | ------------------------------------------------------------------------------------ |
| 1. Disease                           | [Expression Value Matrix](./reference_files/Example_Disease-Value_Matrix.xls)        |
| 2. Control                           | [Expression Value Matrix](./reference_files/Example_Control-Value_Matrix.xls)        |
| 3. Microarray                        | [Expression Value Matrix](./reference_files/Example_Control-Value_Matrix.xls)        |
| 4. Microarray                        | [Gene Status Matrix](./reference_files/Example_Control-Status_Matrix.xls)            |
| 5. Gene Membrane Association File    | [Gene Membrane Association File](./reference_files/MembraneAssociated-Reference.xls) |
| 6. CD Family Identification File     | [CD Family Identification File](./reference_files/CD_Family-Reference.xls)           |
| 7. Cancer-Testis Gene Reference File | [Cancer-Testis Gene Reference File](./reference_files/CancerTestis-Reference.xls)    |
| 8. Protein Expression Data           | [Protein Expression Data](./reference_files/ProteinExpression-Reference.xls)         |
| 9. Name of Run                       | example_run                                                                          |

An configuration file [example](config.txt) is provided.  If you do not wish to provide any of the __optional__ files above, put "NA" on the corresponding line in your configuration file. Though the files have the `.xls` extension, they are tab delimited and can be edited as such.

### Output
The name of the `newCAR.R` run (line 9 of the configuration file) is used to create an output directory.  The output directory contains *data matrices* and *boxplots* regarding the following TAAs:

| Output                        | Description                                                                                                            |
| ----------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `$name~AllMetrics`           | a. All metrics observed in the intersection of genes in tumor and control data matrices                                |
| `$name~MembraneGenes`        | b. Metrics observed for genes with any association to plasma membrane localization                                     |
| `$name~Membrane_Ensembl`     | c. Metrics observed for genes with association to plasma membrane localization in Ensembl                              |
| `$name~Membrane_HPRD`        | d. Metrics observed for genes with association to plasma membrane localization in the Human Protein Reference Database |
| `$name~CD_Family`            | e. CD family TAAs                                                                                                      |
| `$name~CTgenes`              | f. Cancer-testis TAAs                                                                                                  |
| `$name~AllCriteriaSatisfied` | g. TAAs satisfying criteria in files `a-f`                                                                             |
| `$name~TopTargets`           | h. TAAs in file `g` after p-value and expression filters have been applied                                             |

## Maintainers

* [Patrick Schreiner](https://github.com/pschreiner)
