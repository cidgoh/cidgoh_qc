# cidgoh_qc
Nextflow workflow for quality check for NGS data

## Introduction

**cidgoh_qc** is a bioinformatics analysis workflow based on nextflow to perform QC analysis.

## How to install

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

> **_NOTE:_**  You can load nextflow on the Cedar cluster like this: 

``` 
$ module load nextflow/21.04.3 
```


2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`Conda`](https://conda.io/miniconda.html) as package manager. 

> **_NOTE:_**  Docker and Conda are not allowed on the Cedar cluster. By default, singularity is in the default tools, so you don't need to install on the Cedar cluster.


3. Download source code from github

``` 
$ git clone https://github.com/cidgoh/cidgoh_qc.git
```

> **_NOTE:_**  We have set up a default version on the Cedar cluster at /project/rrg-whsiao-ab/shared_tools/cidgoh_qc

## How to use


