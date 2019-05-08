# Setup Instructions:
1) Open server.r, ui.r, install.r in Rstudio
2) Run install.r to install needed packages
3) Run app, you should see the Viewer

## Requirements

Lacks support for NSF/SMB mounted volumes

## Scrips

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