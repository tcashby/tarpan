library(config)

#get genome build from config file
configGenomeBuild <- config::get("genomeBuild")
if (!(configGenomeBuild %in% c("hg19", "hg38"))) {
  stop("Incorrect or no reference genome build defined.")
}

#build file path for p/q boundaries
pqDataFilePath <- file.path("data", configGenomeBuild, "pqboundaries.txt")
if (!file.exists(pqDataFilePath)) {
  stop("p/q boundaries file not found.")
}

#build file path for chromosome sizes
chromSizesFilePath <- file.path("data", configGenomeBuild, "chromsizes.txt")
if (!file.exists(pqDataFilePath)) {
  stop("Chromosome sizes file not found.")
}

#get path of SQLite DB
configSQLiteDB <- config::get("SQLiteDB")
# SQLLiteDBFilePath <-file.path(configSQLiteDB)
# if (!file.exists(SQLLiteDBFilePath)) {
#   stop("SQLite DB file not found.")
# }