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

## Requirements

Due to limitations of SQLite, databases can not be created on NSF/SMB mounted volumes.

## Scripts

**DB Create**
```
create_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<TARGET_FILE>" -groupbed "<GROUP_FILE>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"

required_args = ['db','refgen','pipeline','targetbed','groupbed']
```

**DB Update**
```
update_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<TARGET_FILE>" -groupbed "<GROUP_FILE>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"

required_args = ['db']
possible_args = ['refgen','pipeline','targetbed','groupbed']
```

**Entry Create**
```
create_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>" -mutfile "<MUTATION file>" -muttool "<MUTATION tool>" -mutfile "<MUTATION2 file>" -muttool "<MUTATION2 tool>" -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>" -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>" -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>" 

required_args = ['db','sampid','mutfile','muttool'] 
```

**Entry Update**
```
update_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>" -mutfile "<MUTATION file>" -muttool "<MUTATION tool>" -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>" -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>" -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>" 

required_args = ['db','sampid'] 
possible_args = ['normid','mutfile','svfile','snp','met','depth','log','patid','sex']
```

**Entry Delete**
```
delete_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" 

required_args = ['db','sampid']
```
