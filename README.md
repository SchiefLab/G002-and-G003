<h1 align="center">
  <br>
 G00x - Generalizable Germline-Targeting Clinical Trial Pipeline
</h1>

<div class="flex-container" align="center">
    <a href="https://img.shields.io/badge/Python-3.10-blue">
    <img src="https://img.shields.io/badge/Python-3.10-blue"
        alt="Python Version"></a>
    <a href="https://github.com/psf/black">
    <img src="https://img.shields.io/badge/code%20style-black-000000.svg"
        alt="Format Version"></a>
    <a href="https://github.com/pre-commit/pre-commit">
    <img src="https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white"
        alt="pre commit"></a>
</div>

- [About / Summary](#about--summary)
- [Data Access](#data-access)
  - [G002](#g002)
  - [G003](#g003)
- [Pipeline](#pipeline)
  - [Installation pre-requisites](#installation-pre-requisites)
  - [Installation](#installation)
  - [Validation](#validation)
    - [FACS/Sorting validation](#facssorting-validation)
      - [G002](#g002-1)
      - [G003](#g003-1)
    - [Sequencing files validation](#sequencing-files-validation)
      - [G002](#g002-2)
      - [G003](#g003-2)
  - [FACS/Flow](#facsflow)
    - [G002](#g002-3)
    - [G003](#g003-3)
  - [BCR sequence analysis](#bcr-sequence-analysis)
    - [G002](#g002-4)
    - [G003](#g003-4)
  - [Development](#development)
    - [Testing](#testing)
- [Figures and Tables](#figures-and-tables)
  - [Main Figures](#main-figures)
  - [Supplementary Figures](#supplementary-figures)
- [Issues](#issues)
- [License](#license)

# About / Summary

The G00x pipeline is designed to facilitate the analysis of germline-targeting vaccine clinical trials. It is an all-in-one pipeline and analysis that parses, validates, calculates B cell frequencies, runs 10X and combines all analyses into a plottable dataframe. The pipeline supports datasets for both G002 and G003 trials (related to germline targeting to elicit VRC01-class antibodies against HIV) and ensures that all data is processed and validated to maintain the integrity of the clinical trial results. The pipeline is readily modifiable to enable handling of other types of germline-targeting trials, with other types of analyses, and can also be modified to enable data storage/data flow pathways different from that used in G002 and G003.

---

# Data Access

## G002
To run testing, you will need a `g002` directory in the root of the project. You can get this with the sync command:

```bash
# this will sync the entire contents of the files into g002 directory
# Beware, this file is 1.5TB
aws sync --delete s3 s3://iavig002public/g002/ ./g002'
```

However, we have also setup so you don't need to download the entire directory. You can run the following to get a subset of the data.

```bash
# get the sorting directory
aws s3 cp --recursive s3://iavig002public/g002/G002/sorting ./g002/G002/sorting

# get the sequencing directory excluding large bcl and fastq files
aws s3 cp --recursive s3://iavig002public/g002/G002/sequencing ./g002/G002/sequencing --exclude *working_directory/* --exclude *.fastq.gz --exclude *.tif --exclude *.cbcl --exclude *.imf1 --exclude *.filter --exclude *.bin --exclude *Logs/* --exclude *_stdout --exclude *_stderr --exclude *Autofocus/* --exclude *Intensities/*
```

## G003
To run testing, you will need a `g003` directory in the root of the project. You can get this with the sync command

```bash
# this will sync the entire contents of the files into g003 directory
# Beware, this file is 1.5TB
aws sync --delete s3 s3://iavig003public/g003/ ./g003'
```

However, we have also setup so you don't need to download the entire directory. You can run the following to get a subset of the data.

```bash
# get the sorting directory
aws s3 cp --recursive s3://iavig003public/g003/G003/sorting ./g003/G003/sorting

# get the sequencing directory excluding large bcl and fastq files
aws s3 cp --recursive s3://iavig003public/g003/G003/sequencing ./g003/G003/sequencing --exclude *working_directory/* --exclude *.fastq.gz --exclude *.tif --exclude *.cbcl --exclude *.imf1 --exclude *.filter --exclude *.bin --exclude *Logs/* --exclude *_stdout --exclude *_stderr --exclude *Autofocus/* --exclude *Intensities/*
```

# Pipeline

## Installation pre-requisites

While not necessary, we highly recommend using the [conda](https://docs.conda.io/en/latest/) open-source package and environment manager. For the purposes of this repository, only a minimal installer for anaconda is necessary (Miniconda).

[Miniconda command line installers](https://docs.anaconda.com/miniconda/#quick-command-line-install)

## Installation

This installation assumes that `git` and `conda` are in your path.

```
# clone the repository
git clone https://github.com/SchiefLab/G00x.git

# change directory
cd G00x

# this will create a conda environment called g00x and install the package
./install.sh
```

## Validation

The most important part of this clinical trial other than it working, is that everything is validated. That means that the data is validated, the code is validated, and the results are validated. This is a very important part of the process and should not be skipped.

### FACS/Sorting validation

#### G002
These are instructions on how to validate the file structure containing the FACS/Sorting data before uploading to box.

```bash
#Run the validator for flow from the command line
g00x g002 validate flow my_path/to/box/G002/


# If you used the AWS sync command above, you can use the command:
g00x g002 validate flow ./g002/G002/sorting/G002
```

Your folder structure on box should look like this.

```bash
# if you used the above commands to sync the data. The folder structure will look like this

g002/G002/sorting
└── G002
    ├── Prescreens
    │   ├── Prescreen_RunDate220825_UploadDate221021
    │   ├── Prescreen_RunDate220826_UploadDate221021
    ....
    └── Sorts
        ├── Sort_RunDate220927_UploadDate221013
        ├── Sort_RunDate220928_UploadDate221014
        ├── Sort_RunDate220929_UploadDate221014
        ├── Sort_RunDate220930_UploadDate221014
```

#### G003
```bash
#Run the validator for flow from the command line
g00x g003 validate flow my_path/to/G003/


# If you used the AWS sync command above, you can use th e
g00x g003 validate flow ./g003/G003/sorting/G003
```

For G003, all the data is available in S3 bucket. The structure should look as below:
```bash
g003/G003/sorting
└── G003
    ├── Prescreens
    │   ├── Sort_RunDate230807_UploadDate230807
    │   ├── Sort_RunDate230808_UploadDate230808
    ....
    └── Sorts
        ├── Sort_RunDate230807_UploadDate230807
        ├── Sort_RunDate230808_UploadDate230808
        ├── Sort_RunDate230809_UploadDate230809
        ├── Sort_RunDate231003_UploadDate231027
```
### Sequencing files validation

#### G002
In order to validate the sequencing files, it must be merged with the flow data. Thus, you will need both Globus and Box access. Validation can be accomplished by validating the merge command

```bash

g00x g002 sequencing -f /path/to/flow -s /path/to/sequencing -o output_file


# if you have the AWS structure
g00x g002 merge -s ./g002/G002/sequencing/G002 -f ./g002/G002/sorting/G002 -o my_merged_file
```

The sequencing folder structure should be as follows:

```bash
G002
├── run0002
│   ├── 221006_VH00497_31_AAAVKCLHV
│   └── sample_manifest.csv
├── run0003
│   ├── 221019_VH00497_32_AAANGGVM5
│   └── sample_manifest.csv
└── run0004
    ├── 221101_VL00414_3_AACFYLCM5
    ├── 221103_VL00414_4_AAATJGYM5
    └── sample_manifest.csv
```

#### G003
In G003, all the data are stored in the S3 bucket; thus, you only need access.

```bash

g00x g003 sequencing -f /path/to/flow -s /path/to/sequencing -o output_file


# if you have the AWS structure
g00x g003 merge -s ./g003/G003/sequencing/G003 -f ./g003/G003/sorting/G003 -o my_merged_file
```

```bash
G003
├── run0001
│   ├── 230818_NB552490_0059_AHL7GFBGXM
│   └── sequencing_manifest.csv
├── run0002
│   ├── 230821_NB552490_0060_AHL7GGBGXM
│   └── sequencing_manifest.csv
├── run0003
│   ├── 231108_NB552490_0065_AHL7CCBGXM
│   └── sequencing_manifest.csv
```

Using the above command, you will generate a file called `my_merged_file.csv` which will have the following fields.

```bash
ptid
group
weeks
visit_id
probe_set
sample_type
run_date
sort_pool
hashtag
run_dir_path
pool_number
sorted_date
vdj_sequencing_replicate
cso_sequencing_replicate
vdj_lirary_replicate
cso_library_replicate
bio_replicate
vdj_index
feature_index
vdj_run_id
cso_run_id
```

No data fields should have NA or empty values. If that happens, that means there are non-merged files.

## FACS/Flow
### G002
```bash
# run in one command
g00x g002 flow /path/to/box/G002/ -o path/to/flow_output
```
### G003
```bash
# run in one command
g00x g003 flow /path/to/G003/ -o path/to/flow_output
```

This will output two dataframes—one in CSV and the other in [feather](https://arrow.apache.org/docs/python/feather.html). The dataframe looks like the following.

|      | run_purpose | run_date   | sort_id | ptid    | group | weeks | visit_id | probe_set | sample_type | sort_software_dv | sort_file_type | sample_tube | file_subset | extention | file_path                                                                                                                                                                  | gate | gate_parent | easy_name   | verbose_name | value_type |       value | hashtag | sort_pool |
| ---: | :---------- | :--------- | :------ | :------ | ----: | ----: | :------- | :-------- | :---------- | :--------------- | :------------- | :---------- | :---------- | :-------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--- | :---------- | :---------- | :----------- | :--------- | ----------: | ------: | --------: |
|    0 | PreS        | 2022-08-26 | S6C     | G002611 |     4 |     8 | V600     | Cg28v2    | PBMC        | DV               | Summary        | T1          | a           | .csv      | /Users/jordanwillis/Box/G002/Prescreens/Prescreen_RunDate220826_UploadDate220902/PopulationSummaryFilesFromDV/PreS_220826_S6C_G002611_V600_Cg28v2_PBMC_DV_Summary_T1_a.csv | P4   | P3          | all_b_cells | All B Cells  | count      | 3.04304e+06 |     nan |       nan |
|    1 | PreS        | 2022-08-26 | S6C     | G002710 |     4 |    -5 | V091     | Cg28v2    | PBMC        | DV               | Summary        | T1          | a           | .csv      | /Users/jordanwillis/Box/G002/Prescreens/Prescreen_RunDate220826_UploadDate220902/PopulationSummaryFilesFromDV/PreS_220826_S6C_G002710_V091_Cg28v2_PBMC_DV_Summary_T1_a.csv | P4   | P3          | all_b_cells | All B Cells  | count      |      589576 |     nan |       nan |
|    2 | PreS        | 2022-08-26 | S6C     | G002611 |     4 |     4 | V160     | Cg28v2    | PBMC        | DV               | Summary        | T1          | a           | .csv      | /Users/jordanwillis/Box/G002/Prescreens/Prescreen_RunDate220826_UploadDate220902/PopulationSummaryFilesFromDV/PreS_220826_S6C_G002611_V160_Cg28v2_PBMC_DV_Summary_T1_a.csv | P4   | P3          | all_b_cells | All B Cells  | count      |      518265 |     nan |       nan |
|    3 | PreS        | 2022-08-26 | S6C     | G002710 |     4 |     8 | V600     | Cg28v2    | PBMC        | DV               | Summary        | T1          | a           | .csv      | /Users/jordanwillis/Box/G002/Prescreens/Prescreen_RunDate220826_UploadDate220902/PopulationSummaryFilesFromDV/PreS_220826_S6C_G002710_V600_Cg28v2_PBMC_DV_Summary_T1_a.csv | P4   | P3          | all_b_cells | All B Cells  | count      | 2.39997e+06 |     nan |       nan |

Each row is a single flow measurement. For example, the count of how many B cells were in the sort. There are many meta data columns as well. A single donor, time point and flow experiment is repeated for each measurement we are after. This makes it very easy to extract and plot data.

These are the fields from the previous dataframe.

| fields           | description                                                                    |
| :--------------- | :----------------------------------------------------------------------------- |
| run_purpose      | Sort or Preclinical (PreS)                                                     |
| run_date         | The date the sort was completed                                                |
| sort_id          | internal sort id, e.g. S6                                                      |
| pubID            | the patient id, e.g. G002611                                                   |
| group            | the vaccine group, 1,2,3,4                                                     |
| weeks            | the week from the sample                                                       |
| visit_id         | the visit id corresponding with the group                                      |
| probe_set        | the group this flow sample was sorted with                                     |
| sample_type      | PBMC or FNA                                                                    |
| sort_software_dv | Diva software, DV                                                              |
| sort_file_type   | Only Summary at this point. But could be FCS                                   |
| sample_tube      | the labeling of the sample tube                                                |
| file_subset      | if the sorter is started or stopped, this can divide the file, e.g a,b,c..     |
| extension        | the file extension                                                             |
| file_path        | the path of the file                                                           |
| gate             | the gate name, e.g P1                                                          |
| gate_parent      | The parent gate of the gate                                                    |
| easy_name        | an easier more programmatic name for the gate measurement, e.g antigenic_count |
| verbose_name     | A more verbose name for labeling                                               |
| value_type       | if the gate a frequency or a count                                             |
| value            | what is the value of the gate or frequency                                     |
| hashtag          | if this sample has been sorted, what is the hashtag                            |
| sort_pool        | the sorting pool, P01 - P10                                                    |

## BCR sequence analysis

### G002
```bash
g00x g002 pipeline demultiplex -f ./g002/G002/sorting/G002/ -s ./g002/G002/sequencing/G002/ -o ./g002/G002/output/demultiplex

g00x g002 pipeline vdj -d ./g002/G002/output/demultiplex.feather -o ./g002/G002/output/vdj

g00x g002 pipeline cso -d ./g002/G002/output/demultiplex.feather -o ./g002/G002/output/cso

# Merge VDJ and CSO dataframes, run the output through SADIE for
g00x g002 pipeline airr -k 3 -c ./g002/G002/output/cso.feather -v ./g002/G002/output/vdj.feather -o ./g002/G002/output/final_df

# Output merged.feather not used to generate figures, but is a sanity check.
g00x g002 validate merge -s ./g002/G002/sequencing/G002/ -f ./g002/G002/sorting/G002 -o ./g002/G002/output/merged

# Output flow_and_sequencing.feather used to generate figures; will contain counts, frequencies, and metadata from flow and sequencing data.
g00x g002 analysis report -s ./g002/G002/output/final_df.feather -f ./g002/G002/output/flow_output.feather -o ./g002/G002/output/flow_and_sequencing
```
### G003
```bash
g00x g003 pipeline demultiplex -f ./g003/G003/sorting/G003/ -s ./g003/G003/sequencing/G003/ -o ./g003/G003/output/demultiplex

g00x g003 pipeline vdj -d ./g003/G003/output/demultiplex.feather -o ./g003/G003/output/vdj

g00x g003 pipeline cso -d ./g003/G003/output/demultiplex.feather -o ./g003/G003/output/cso

# Merge VDJ and CSO dataframes, run the output through SADIE for
g00x g003 pipeline airr -k 3 -c ./g003/G003/output/cso.feather -v ./g003/G003/output/vdj.feather -o ./g003/G003/output/final_df

# Output merged.feather not used to generate figures, but is a sanity check.
g00x g003 validate merge -s ./g003/G003/sequencing/G003/ -f ./g003/G003/sorting/G003 -o ./g003/G003/output/merged

# Output flow_and_sequencing.feather used to generate figures; will contain counts, frequencies, and metadata from flow and sequencing data.
g00x g003 analysis report -s ./g003/G003/output/final_df.feather -f ./g003/G003/output/flow_output.feather -o ./g003/G003/output/flow_and_sequencing
```


## Development

### Testing

If you'd like to help develop. Please use pre-commit before commiting any files.

```bash
# runs pre-commit for syntax,formatting and static type checking
poetry run pre-commit run --all-files
```

Now you can run the tests via poetry

```bash
poetry run pytest -sv --log-cli-level DEBUG tests
```

# Figures and Tables

## Main Figures

Main figures as they are displayed in the G002 and G003  paper. Figure 1 is not included as it is manually generated.

```bash
g00x plots fig2
g00x plots fig3
g00x plots fig4
g00x plots fig5
g00x plots fig6
g00x plots fig7
```


## Supplementary Figures

Supplementary figures are in no particular order. They are generated as needed.

```bash
g00x plots suppfigures --all #TODO: Update later
```

# Issues

Please submit any issues to the [issues page](https://github.com/SchiefLab/G00x/issues) and we are happy to help.

# License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

Copyright © Jordan R. Willis, Troy Sincomb, Caleb Kibet, VISC, and IAVI
