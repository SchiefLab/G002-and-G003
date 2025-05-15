# Validation

The following discusses data structure and validation.

## Data structure

The data stored in our bucket is structured as follows in the `s3://iavig003sabucket/g003/` bucket.

```bash
├── g001
│   ├── sequencing
│   │   ├── run0001
│   │   └── run0002
│   └── sorting
│   │   └── Sorts
│   └── output
└── g003
│   ├── sequencing
│   │   ├── run0001
│   │   └── run0002
│   │   └── run0003
│   │   └── run0004
│   │   └── ...
│   └── sorting
│   │    └── Sorts
│   └── output
```

### Sequencing

The sequencing folder structure takes the following file structure.

```bash
g001/sequencing/
├── run0001
│   ├── 221118_VH00124_107_AAAW2V3HV
│   ├── sequencing_manifest.csv
│   └── working_directory
└── run0002
    ├── 230201_NB552490_0042_AHKY7KBGXM
    ├── sequencing_manifest.csv
    └── working_directory
```

Each sequencing run is uploaded numerically runXXXX. Each run folder contains a `sequencing_manifest.csv` file that looks like:

|     | ptid      | timpoint | sorted_date | cells | hto   | vdj_index | cso_index | pool_number | run_id                       |
| --: | :-------- | :------- | ----------: | ----: | :---- | :-------- | :-------- | ----------: | :--------------------------- |
|   0 | PubID_023 | V08      |      221107 |  5076 | C0251 | SI-TT-G1  | SI-TN-A1  |           1 | 221118_VH00124_107_AAAW2V3HV |
|   1 | PubID_187 | V08      |      221107 |   290 | C0252 | SI-TT-G2  | SI-TN-A2  |           2 | 221118_VH00124_107_AAAW2V3HV |
|   2 | PubID_153 | V08      |      221107 |  1222 | C0253 | SI-TT-G2  | SI-TN-A2  |           2 | 221118_VH00124_107_AAAW2V3HV |
|   3 | PubID_046 | V08      |      221107 |   288 | C0254 | SI-TT-G2  | SI-TN-A2  |           2 | 221118_VH00124_107_AAAW2V3HV |
|   4 | PubID_023 | V08      |      221107 |  4993 | C0251 | SI-TT-G3  | SI-TN-A3  |           3 | 221118_VH00124_107_AAAW2V3HV |
|   5 | PubID_110 | V08      |      221108 |  3183 | C0255 | SI-TT-G4  | SI-TN-A4  |           1 | 221118_VH00124_107_AAAW2V3HV |
|   6 | PubID_154 | V08      |      221108 |  5000 | C0256 | SI-TT-G5  | SI-TN-A5  |           2 | 221118_VH00124_107_AAAW2V3HV |
|   7 | PubID_154 | V08      |      221108 |  1560 | C0256 | SI-TT-G6  | SI-TN-A6  |           3 | 221118_VH00124_107_AAAW2V3HV |
|   8 | PubID_047 | V08      |      221108 |  1480 | C0257 | SI-TT-G7  | SI-TN-A7  |           4 | 221118_VH00124_107_AAAW2V3HV |
|   9 | PubID_047 | V08      |      221108 |  1811 | C0257 | SI-TT-G8  | SI-TN-A8  |           5 | 221118_VH00124_107_AAAW2V3HV |
|  10 | PubID_079 | V08      |      221108 |  1625 | C0258 | SI-TT-G9  | SI-TN-A9  |           6 | 221118_VH00124_107_AAAW2V3HV |

Where each field is defined as:

| Column      | Definition                                                    |
| :---------- | :------------------------------------------------------------ |
| ptid        | The participant id                                            |
| timpoint    | The timepoint of the sequencing                               |
| sorted_date | The date it was sorted `YYMMDD` format                        |
| cells       | The number of cells in the GEM reaction                       |
| hto         | The hashtag number e.g. `HT01` |
| vdj_index   | The index used for the VDJ library                            |
| cso_index   | The index used for the hashtag library                        |
| pool_number | The pool number of the gem reaction                           |
| vdj_run_id  | The Illumina flow-cell id that the VDJ library was run on     |
| cso_run_id  | The Illumina flow-cell id that the hashtag library was run on |

!!! info run_id

 Often, the VDJ and CSO library will be on the same flow-cell ID.

### Sorting

The sorting directory takes the following file structure.

```bash
g001/sorting/
└── Sorts
    ├── Sort_RunDate221118_UploadDate221222
    │   └── ClinicalSamples
    │       ├── DataFilesFromMelody
    │       │   ├── SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Data_T1_P1_a.fcs
    │       │   ├── SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Data_T1_P3_b.fcs
                ..
    │       ├── DataStats
    │       │   ├── SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Summary_T1_P1_a.xlsx
    │       │   ├── SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Summary_T1_P3_b.xlsx
                ..
    │       ├── ScreenshotsCounts
                ..
    │       └── ScreenshotsMelodyStats
    │           ├── SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Summary_T1_P1_a.png
    │           ├── SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Summary_T1_P3_b.png

```

Each sort data upload can have multiple sorts. Each upload will have the upload date and sort date, e.g. `Sort_RunDate221118_UploadDate221222`.

There are four sub-directories:

1. `ClinicalSamples/DataFilesFromMelody` - The FCS files from the sort
2. `ClinicalSamples/DataStats` - The summary statistics from the sort calculated in FlowJo.
3. `ClinicalSamples/ScreenshotsCounts` - The screenshots of the counts of the sort.
4. `ClinicalSamples/ScreenshotsMelodyStats` - The screenshots of the gating.

Notice each file name has the following format, e.g. `SortHCT019a_221107_M1_PubID_023_V8_eODGT8_PBMC_Chorus_Summary_T1_P1_a.xlsx`

| Field       | Definition                               |
| :---------- | :--------------------------------------- |
| Sort        | The sort name e.g. `SortHCT019a` |
| Date        | The date the sort was performed `YYMMDD` |
| Machine     | The machine number e.g. `M1` |
| ptid        | The participant id                       |
| Visit       | The visit number e.g. `V08` |
| sort_probe  | The sort probe e.g. `eODGT8` |
| sample_type | The sample type e.g. `PBMC` |
| software    | The software used e.g. `Chorus` |
| type        | The type of file e.g. `Summary` |
| tube_number | The tube number e.g. `T1` |
| pool_number | The pool number e.g. `P1` |
| file_subset | The file subset e.g. `a` |
