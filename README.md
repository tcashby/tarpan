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

| Argument | Description |
| -------- | ----------- |
| `-refgen`	| Reference Genome	|
| `-pipeline`	| Pipeline Name	|
| `-targetbed`	| Target Bed File	|
| `-groupbed`	| Group Bed File	|
| `-blacklist`	| Blacklist Bed File	|
| `-bedfile`	| Other Bed File	|
| `-bedtype`	| Other Bed File Type	|
| `-sampid`	| Sample ID	|
| `-normid`	| Normal ID	|
| `-mutfile`	| Mutation File	|
| `-muttool`	| Mutation File Tool	|
| `-svfile`	| Structural Variant File	|
| `-svtool`	| Structural Variant File Tool	|
| `-snp`	| SNP File	|
| `-met`	| Metrics File	|
| `-depth`	| Depth File	|
| `-log`	| Log File	|
| `-sex`	| Patient Sex	|
| `-patid`	| Patient Identifier	|

**DB Create**
```
create_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<TARGET_BED>" -groupbed "<GROUP_BED>" -blacklist "<BLACKLIST_BED>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"

required_args = ['db','refgen','targetbed','groupbed']
```

**DB Update**
```
update_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<TARGET_BED>" -groupbed "<GROUP_BED>" -blacklist "<BLACKLIST_BED>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"

required_args = ['db']
possible_args = ['refgen','pipeline','targetbed','groupbed','blacklist','bedfile']
```

**Entry Create**
```
create_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>" -mutfile "<MUTATION file>" -muttool "<MUTATION tool>" -mutfile "<MUTATION2 file>" -muttool "<MUTATION2 tool>" -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>" -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>" -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>" 

required_args = ['db','sampid']
possible_args = ['mutfile','svfile','snp','met','depth'] 
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
