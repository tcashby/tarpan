<p align="center">
  <img src="/www/LogoBig.jpg" width="25%" height="25%">
</p>

##
TarPan Viewer is a tool used to visually inspect targeted panel sequencing data.

You can access the demo site [here](https://tarpan.shinyapps.io/tarpan/) (it may take a moment to load example data).

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

3. Install dependencies (this command should work, but you may have to manually intervene if a package gets updated)

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

1. Download Python 3.7.x from: https://www.python.org/downloads/
2. Install Python (If Windows, select "Add PATH" on install)
3. Install Pandas via Command Prompt or Shell prompt.  
  ```sh
    python -m pip install Pandas
  ```
If issues occur, this [Python tutorial](https://docs.python.org/3/tutorial/interpreter.html) might help.

#### Script Arguments

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
```sh
required_args = -db, -refgen, -targetbed, -groupbed
```
```cs
create_db.py -db "Tarpan.db" -refgen "hg19" -pipeline "pipeline01"
             -targetbed "target.bed" -groupbed "group.bed" -blacklist "blacklist.bed"
             -bedfile "other.bed" -bedtype "other"
```

**DB Update**
```sh
required_args = -db
```
```cs
update_db.py -db "Tarpan.db" -refgen "hg19" -pipeline "pipeline01"
             -targetbed "target.bed" -groupbed "group.bed" -blacklist "blacklist.bed"
             -bedfile "other.bed" -bedtype "other"
```

**Entry Create**
```sh
required_args = -db, -sampid
```
```cs
create_db_entry.py -db "Tarpan.db" -sampid "sample01" -normid "sample01Norm"
                   -mutfile "sample01Mut.vcf" -muttool "Strelka2"
                   -mutfile "sample01Mut02.vcf" -muttool "MuTect"
                   -svfile "sample01SV.vcf" -svtool "Manta"
                   -snp "sample01SNP.csv" -met "sample01MET.csv" -depth "sample01DEP.csv"
                   -log "sample01.log" -patid "patientId" -sex "Male"
```

**Entry Update**
```sh
required_args = -db, -sampid
```
```cs
update_db_entry.py -db "Tarpan.db" -sampid "sample01" -normid "sample01Norm"
                   -mutfile "sample01Mut.vcf" -muttool "Strelka2"
                   -svfile "sample01SV.vcf" -svtool "Manta"
                   -snp "sample01SNP.csv" -met "sample01MET.csv" -depth "sample01DEP.csv"
                   -log "sample01.log" -patid "patientId" -sex "Male"
```

**Entry Delete**
```sh
required_args = -db, -sampid
```
```cs
delete_db_entry.py -db "Tarpan.db" -sampid "sample01"
```
