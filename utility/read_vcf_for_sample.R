#test vcf reconstruction
library(DBI)
library(RSQLite)
library(VariantAnnotation)

read_vcf_for_sample = function(mydb, sample, tool)
{
  #read in the vcf header for a specific tool's vcf file
  header_query = paste0("SELECT file_header FROM GENOM_FILES where sample_id = '", sample,"' AND file_tool = '",tool,"'")
  header = dbGetQuery(mydb, header_query, stringsAsFactors = F)
  header = header$file_header
  
  #read in the data for that vcf
  vcf_data_query = paste0("SELECT * from GENOM_MUTATIONS where sample_id = '", sample,"' AND mut_tool = '",tool,"'")
  vcf_data = dbGetQuery(mydb, vcf_data_query, stringsAsFactors = F)
  vcf_data = vcf_data[,-1:-3]
  
  #write a temporary file and read it into VariantAnnotator format
  temp_filename = tempfile(fileext = ".vcf")
  write.table(vcf_data,'test.txt')
  writeChar(header, temp_filename)
  colnames(vcf_data)[1] = "#CHROM"
  suppressWarnings(write.table(vcf_data, file=temp_filename, sep="\t", row.names=F, quote=F, col.names=T, append = T))
  vcf = readVcf(temp_filename)
  unlink(temp_filename)
  vcf
}

make_table = function(vcf)
{
  ref = as.character(ref(vcf))
  alt = as.character(unlist(alt(vcf)))
  filter = as.character(filt(vcf))
  ranges = rowRanges(vcf)
  chrom = as.character(ranges@seqnames)
  pos = as.numeric(ranges@ranges@start)
  
  vcf_table = data.frame(chrom = chrom, pos = pos, ref = ref, alt = alt, filter = filter)
  
  #if CSQ is there, load the annotations
  if(any(rownames(info(header(vcf))) %in% "CSQ"))
  {
    #this is hacky, but seems to work to extract the key for what features are available in VEP annotation
    vep_features = unlist(strsplit(gsub("Consequence annotations from Ensembl VEP. Format: ","",as.character(info(header(vcf))[which(rownames(info(header(vcf))) %in% "CSQ"),3])),"\\|"))
    CSQ = info(vcf)$CSQ
    CSQ = unlist(CSQ)
    CSQ = read.table(text = CSQ, sep="|", stringsAsFactors = F)
    colnames(CSQ) = vep_features
    
    vcf_table$consequence = CSQ$Consequence
    vcf_table$symbol = CSQ$SYMBOL
    vcf_table$ccode = CSQ$HGVSc
    vcf_table$pcode = gsub("%3D","=",CSQ$HGVSp)
  }
  vcf_table
}
