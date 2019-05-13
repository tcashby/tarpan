<p align="center">
  <img src="/www/LogoBig.jpg" width="25%" height="25%">
</p>

## 
TarPan Viewer is a tool used to visually inspect targeted panel sequencing data.


## Built With
- [R](https://www.r-project.org/)
- [R Shiny](https://shiny.rstudio.com/)
- [Python](https://www.python.org/)

## Dependencies
The easiest and fastest way to run TarPan is to install [R Studio](https://www.rstudio.com/). 

Python 3.x + Pandas

## Installation of TarPan Viewer
1. Clone the repository

```sh
git clone https://github.com/tcashby/tarpan.git
```
2. Create a new R Studio project in the cloned directory and open ```ui.r```, ```server.r```, and ```install.r```.

3. Install dependencies

```sh
source("install.R")
```

4. Assuming the dependencies installed correctly, you should now be able to run the application in R Studio.

## Limitations

Due to limitations of SQLite, databases can not be created on NSF/SMB mounted volumes.

## Scripts
In order to utilize the Python scripts, Python must be installed.

Python can be installed standalone or through Anaconda.

#### Python Setup

1. Download Python 3.7.x from:https://www.python.org/downloads/


#### Anaconda setup

1. Install Anaconda (for Python 3.7) from: https://www.anaconda.com/distribution/
2. Install Pandas via Anaconda Navigator or using the Anaconda Prompt using PIP.  
```sh
  python -m pip install Pandas
```

#### Script Arguments

| Argument | Description |  | Argument | Description |
| -------- | ----------- |  | -------- | ----------- |
| `-refgen`	| Reference Genome	|  | `-muttool`	| Mutation File Tool	|
| `-pipeline`	| Pipeline Name	|  | `-svfile`	| Structural Variant File	|
| `-targetbed`	| Target Bed File	|  | `-svtool`	| Structural Variant File Tool	|
| `-groupbed`	| Group Bed File	|  | `-snp`	| SNP File	|
| `-blacklist`	| Blacklist Bed File	|  | `-met`	| Metrics File	|
| `-bedfile`	| Other Bed File	|  | `-depth`	| Depth File	|
| `-bedtype`	| Other Bed File Type	|  | `-log`	| Log File	|
| `-sampid`	| Sample ID	|  | `-sex`	| Patient Sex	|
| `-normid`	| Normal ID	|  | `-patid`	| Patient Identifier	|
| `-mutfile`	| Mutation File	|




**DB Create**
```sh
required_args = ['db','refgen','targetbed','groupbed']
```
```cs
create_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>"
             -targetbed "<TARGET_BED>" -groupbed "<GROUP_BED>" -blacklist "<BLACKLIST_BED>"
             -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"
```

**DB Update**
```sh
required_args = ['db']
```
```cs
update_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>"
             -targetbed "<TARGET_BED>" -groupbed "<GROUP_BED>" -blacklist "<BLACKLIST_BED>"
             -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"
```

**Entry Create**
```sh
required_args = ['db','sampid']
```
```cs
create_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>"
                   -mutfile "<MUTATION file>" -muttool "<MUTATION tool>"
                   -mutfile "<MUTATION2 file>" -muttool "<MUTATION2 tool>"
                   -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>"
                   -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>"
                   -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>"
```

**Entry Update**
```sh
required_args = ['db','sampid']
```
```cs
update_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>"
                   -mutfile "<MUTATION file>" -muttool "<MUTATION tool>"
                   -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>"
                   -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>"
                   -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>"
```

**Entry Delete**
```sh
required_args = ['db','sampid']
```
```cs
delete_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" 
```
