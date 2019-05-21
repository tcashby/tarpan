#dependencies
library(shiny)
library(shinyBS)
library(RColorBrewer)
library(DT)
#sqlite libraries
library(DBI)
library(RSQLite)
#CircosDiagram
library(RCircos)

#load utility scripts
source("utility/read_vcf_for_sample.R")
source("utility/readConfig.R")

#get database from config file
databases <- config::get("SQLiteDB")

#create handle to sqlite file
dbhandle <- dbConnect(RSQLite::SQLite(), databases[1])

#established p and q boundaries
pq <- read.table(pqDataFilePath, sep="\t", stringsAsFactors = FALSE, 
                header = TRUE)
#read in the chromosome sizes for a genome of interest
chrom_sizes <- read.table(chromSizesFilePath, sep = "\t", 
                          stringsAsFactors = FALSE, header = FALSE)
chrom_sizes <- c(chrom_sizes$V2, 0)
#cumulative sums of the chromosome sizes
cumulative_chrom_sizes <- cumsum(chrom_sizes)
cumulative_chrom_sizes <- c(0, cumulative_chrom_sizes)
cumulative_chrom_sizes <- cumulative_chrom_sizes[1:25]

#read groups from database
query <- "SELECT chrom, start, end, id from GENOM_BEDFILES where type = 'groups'"
genome_bed <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
#remove chrom prefix
genome_bed$chrom <- gsub("chr", "", genome_bed$chrom)

#genes in the blacklist bed file will be hidden by default
query <- "SELECT chrom, start, end, id from GENOM_BEDFILES "
query <- paste0(query, "where type = 'blacklist'")
blacklist <- dbGetQuery(dbhandle, query, stringsAsFactors = F)
if (nrow(blacklist) == 0) {
  blacklist <- data.frame(V1 = "N", V2 = -1, V3 = -1, stringsAsFactors = FALSE)
}
blacklist$V1[blacklist$V1 %in% "X"] <- "23"

# Define server logic for application
shinyServer(function(input, output, session) {
  #when the database changes redraw
  observeEvent(input$database, {
    dbDisconnect(dbhandle)
    dbhandle <<- dbConnect(RSQLite::SQLite(), input$database)
    res = dbGetQuery(dbhandle, "SELECT sample_id from GENOM_SAMPLE")

    #update sample list
    updateSelectInput(session, "sample", label="Sample:", 
                      choices = res$sample_id, selected = res$sample_id[1])

    #update gene groups
    query <- "SELECT chrom, start, end, id from GENOM_BEDFILES "
    query <- paste0(query, "where type = 'groups'")
    genome_bed <<- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
    #remove chrom prefix
    genome_bed$chrom <<- gsub("chr", "", genome_bed$chrom)

    #reload db
    query <- "SELECT chrom, start, end, id from GENOM_BEDFILES "
    query <- paste0(query, "where type = 'blacklist'")
    blacklist <<- dbGetQuery(dbhandle, query, stringsAsFactors = F)
    if (nrow(blacklist) == 0) {
      blacklist <<- data.frame(V1 = "N", V2 = -1, V3 = -1, stringsAsFactors = FALSE)
    }
    else {
      colnames(blacklist) <<- c("V1","V2","V3","V4")
    }
    blacklist$V1[blacklist$V1 %in% "X"] <- "23"
  })

  #returns all input chromosomes for normalization
  get_chrom_normalization_input <- function() {
    if (!input$chromnorm) {
      1:22
    } else {
      ids <- paste0("chrom", 1:22)
      x <- sapply(ids, function(x) input[[x]])
      x <- ids[which(x == TRUE)]
      x <- as.numeric(gsub("chrom", "", x))
      if (length(x) == 0) {
        1:22
      } else {
        x
      }
    }
  }

  # generate a plot of copy number data. 
  draw_copy_number_plot <- function() {
    #load the data
    query <- "SELECT chrom, start, end, id, normal_mean, tumor_mean, "
    query <- paste0(query, "tumor_depth from GENOM_DEPTH where sample_id = '")
    query <- paste0(query, input$sample,"'")
    genome_depth <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
    #cast to proper data types
    genome_depth$start <- as.numeric(genome_depth$start)
    genome_depth$end <- as.numeric(genome_depth$end)
    genome_depth$normal_mean <- as.numeric(genome_depth$normal_mean)
    genome_depth$tumor_mean <- as.numeric(genome_depth$tumor_mean)

    #if the chr prefix is on chrom, remove it
    genome_depth$chrom <- gsub("chr", "", genome_depth$chrom)

    #remove the sex chromosomes and where the mean == 0
    genome_depth <- genome_depth[genome_depth$tumor_mean != 0, ]
    genome_depth <- genome_depth[genome_depth$chrom != "Y", ]
    genome_depth$chrom[genome_depth$chrom %in% "X"] <- "23"

    #sort by chrom, start
    genome_depth <- genome_depth[order(genome_depth$chrom, genome_depth$start), ]

    #this hides genes on the blacklist
    if (!input$show_blacklist) {
      for (i in 1:nrow(blacklist)) {
        genome_depth <- genome_depth[!(genome_depth$start >= blacklist$V2[i] & 
                genome_depth$end <= blacklist$V3[i] &
                genome_depth$chrom %in% blacklist$V1[i]), ]
      }
    }
  
    #set up x-axis
    xaxis <- genome_depth$start + 
      cumulative_chrom_sizes[as.numeric(genome_depth$chrom)]
    #set up y-axis, if there's no normalization flag, 
    #generate a base normalization

    #if tumor depth is set in the database, use those values initially
    if (!input$chromnorm && 
        length(genome_depth$tumor_depth[!is.na(genome_depth$tumor_depth)]) > 0) {
      genome_depth <- genome_depth[!is.na(genome_depth$tumor_depth), ]
      xaxis <- genome_depth$start + 
        cumulative_chrom_sizes[as.numeric(genome_depth$chrom)]
      yaxis <- genome_depth$tumor_depth
    } else {    
      chromlist <- get_chrom_normalization_input()
      n_fac <- mean(genome_depth$normal_mean[genome_depth$chrom %in% chromlist])
      t_fac <- mean(genome_depth$tumor_mean[genome_depth$chrom %in% chromlist])
      yaxis <- (genome_depth$tumor_mean / t_fac) /
        (genome_depth$normal_mean / n_fac)
    }

    #set min and max limits
    maxlim <- 2
    minlim <- 0
    if (input$nolim)
    {
      temp_y <- yaxis[!is.infinite(yaxis)]

      if (as.numeric(input$maxlim) > 0) {
        maxlim <- as.numeric(input$maxlim)
      } else {
        maxlim <- max(temp_y)
        if (maxlim < 2) {
          maxlim <- 2
        }
      }
    }

    #this controls upper and lower bounds
    upper_bound <- 1.25
    lower_bound <- 0.75
    
    #this colors the control genes gray
    genecolors <- rep("black", length(xaxis))
    genecolors[yaxis >= upper_bound] = "#e41a1c"
    genecolors[yaxis <= lower_bound] = "#377eb8"

    if (input$type == 1) {
      #display pandq boundaries
      pandq <- c()
      if (input$pandq) {
        pandq <- pq$pqbound + cumulative_chrom_sizes[as.numeric(pq$chrom)]
      }

      #count chromosome label positions
      my_breaks <- chrom_sizes / 2 + cumulative_chrom_sizes

      plot(NA, axes = FALSE, ylim = c(minlim, maxlim),
           xlim = c(0, max(cumulative_chrom_sizes[1:24])), xlab = "Chromosome",
           ylab = "Mean Ratio")
      
      #shading on the chromosomes
      for (i in 1:23) {
        if (i %% 2 == 0) {
          polygon(c(cumulative_chrom_sizes[i], cumulative_chrom_sizes[i], 
                    cumulative_chrom_sizes[i + 1], 
                    cumulative_chrom_sizes[i + 1]),
                  c(0, maxlim, maxlim, 0), col = "#d9d9d9", border = FALSE)
        } else {
          polygon(c(cumulative_chrom_sizes[i], cumulative_chrom_sizes[i], 
                    cumulative_chrom_sizes[i + 1], 
                    cumulative_chrom_sizes[i + 1]),
                  c(0, maxlim, maxlim, 0), col = "#f0f0f0", border = FALSE)
        }
      }
      
      points(xaxis, as.numeric(yaxis), pch = 16, cex = 0.7, col = genecolors)
      axis(side = 1, at = my_breaks[1:23], tick = FALSE, labels = c(1:22, "X"),
           las = 2)
      abline(v = pandq, col = "#FFFFFF", lty = 2)
      abline(h = c(lower_bound, upper_bound), col = "#A6ACA7", lty = 2)
      if (input$normline) {
        abline(h = 1, col = "black", lty = 2)
      }
    } else {
      #set up x-axis
      xaxis <- genome_depth$start[genome_depth$chrom == input$chromnum]
      yaxis <- yaxis[genome_depth$chrom == input$chromnum]

      genecolors <- rep("black", length(xaxis))
      pchtype <- rep(21, length(xaxis))

      #set a common boundary so the aspect ratio will remain constant
      boundary <- 100000
      if (input$continuous & input$showlast) {
        plot(NA, axes = FALSE, ylim = c(-0.25, maxlim + 1), 
             xlim = c(0, boundary), xlab = paste("Chromosome", input$chromnum), 
             ylab = "Mean Ratio")
      } else {
        plot(NA, axes = FALSE, ylim = c(-0.25, maxlim), xlim = c(0, boundary),
             xlab = paste("Chromosome", input$chromnum), ylab = "Mean Ratio")
      }
      polygon(c(0, 0, boundary, boundary), c(0, maxlim, maxlim, 0), 
              col = "#FFFFFF", border = FALSE)
      
      last_gene_line <- c()
      last_gene_name <- c()

      group <- genome_bed[genome_bed$chrom == input$chromnum, ]
      pal <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
               "#757780", "#b15928")
      pal <- rep(pal, 1000)
      pal <- pal[1:dim(group)[1]]
      group <- cbind(group, pal)
      for (i in 1:dim(group)[1]) {
        if (length(which(xaxis >= group$start[i] & xaxis<=group$end[i])) > 0) {
          genecolors[which(xaxis >= group$start[i] & xaxis <= group$end[i])] <- 
            as.character(group$pal[i])
          last_gene_name <- c(last_gene_name, group$id[i])
          last_gene_line <- c(last_gene_line, min(which(xaxis >= group$start[i] 
                                                        & xaxis<=group$end[i])))
        } else {
          genecolors[which(xaxis >= group$start[i] & xaxis <= group$end[i])] <- 
            "black"
        }
      }
      
      if (input$chromnum == 8) {
        legend(0, 0, legend = group$id, col = as.character(group$pal), pch = 16,
               horiz = TRUE, bty = "n", cex = 0.75)
      } else {          
        legend(0, 0, legend = group$id, col = as.character(group$pal), pch = 16,
               horiz = TRUE, bty = "n", cex = 1.0)
      }
      
      if (input$continuous) {
        spacer <- floor(chrom_sizes[input$chromnum] / length(xaxis))
        xaxis <- 1:length(xaxis) * spacer
      }

      pandq <- c()
      if (input$pandq & !(input$continuous)) {
        pandq <- pq$pqbound[as.numeric(pq$chrom) == input$chromnum]
        pandq <- (pandq / chrom_sizes[input$chromnum]) * boundary
      }
      if (input$pandq & input$continuous) {
        pandq <- pq$pqbound[as.numeric(pq$chrom) == input$chromnum]
        pandq <- (xaxis[min(which(genome_depth$start[genome_depth$chrom == 
                                                    input$chromnum] > pandq))]
                  / chrom_sizes[input$chromnum]) * boundary          
        spacer <- floor(chrom_sizes[input$chromnum] / length(xaxis))
        pandq <- pandq - ((spacer / 2) / chrom_sizes[input$chromnum] * boundary)
      }
      
      points((xaxis / chrom_sizes[input$chromnum]) * boundary, as.numeric(yaxis),
             pch = pchtype, cex = 1.3, col = "black", bg = genecolors)

      abline(v = pandq, col = "#000000", lty = 2)
      abline(h = c(lower_bound, upper_bound), col = "#00000080", lty = 5, 
             lwd = 1)
      abline(h = c(0.25,1.75), col = "#000000BF", lty = 5, lwd = 1)

      if (maxlim == 2) {
        axis(side = 2, las = 1, at = c(0:4) * 0.5)
        axis(side = 4, at = c(0, 0.5, 1.0, 1.5, 2.0), 
             labels = c("0", "1", "2", "3", "4"), las = 1)
      } else {
        axis(side = 2, at = 0:(floor(maxlim / 0.4) - 1) * 0.5, las = 1)
        axis(side = 4, at = 0:(floor(maxlim / 0.4) - 1) * 0.5, 
             labels = 0:(floor(maxlim / 0.4) - 1), las = 1)
      }

      stagger <- 0
      if (input$continuous & input$showlast) {
        for (i in 1:length(last_gene_line)) {
          if (stagger %% 5 == 0) {
            text(x = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum] *
                 boundary, y = maxlim + 1, adj = c(0.5, 0.5), 
                 labels = last_gene_name[i], col = "black")
            segments(x0 = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum]
                     * boundary, x1 = xaxis[last_gene_line[i]] / 
                       chrom_sizes[input$chromnum] * boundary, y0 = 0.17,
                     y1 = maxlim + 0.9, col = "#00000080", lty = 3, lwd = 0.5)
          } else if (stagger %% 4 == 0) {
            text(x = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum]
                 * boundary, y = maxlim + 0.8, adj = c(0.5, 0.5),
                 labels = last_gene_name[i], col = "black")
            segments(x0 = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum]
                     * boundary, x1 = xaxis[last_gene_line[i]] / 
                       chrom_sizes[input$chromnum] * boundary, y0 = 0.17,
                     y1 = maxlim + 0.7, col = "#00000080", lty = 3, lwd = 0.5)
          } else if (stagger %% 3 == 0) {
            text(x = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum] * 
                   boundary, y = maxlim + 0.6, adj = c(0.5, 0.5),
                 labels = last_gene_name[i], col = "black")
            segments(x0 = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum]
                     * boundary, x1 = xaxis[last_gene_line[i]] / 
                       chrom_sizes[input$chromnum] * boundary, y0=0.17, 
                     y1=maxlim + 0.5, col = "#00000080", lty = 3, lwd = 0.5)
          } else if (stagger %% 2 == 0) {
            text(x = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum] * 
                   boundary, y = maxlim + 0.4, adj = c(0.5, 0.5), 
                 labels = last_gene_name[i], col = "black")
            segments(x0 = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum]
                     * boundary, x1 = xaxis[last_gene_line[i]] / 
                       chrom_sizes[input$chromnum] * boundary, y0 = 0.17, 
                     y1 = maxlim + 0.3, col = "#00000080", lty = 3, lwd = 0.5)
          } else {
            text(x = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum] * 
                   boundary, y = maxlim + 0.2, adj = c(0.5, 0.5),
                 labels = last_gene_name[i], col = "black")
            segments(x0 = xaxis[last_gene_line[i]] / chrom_sizes[input$chromnum]
                     * boundary, x1 = xaxis[last_gene_line[i]] / 
                       chrom_sizes[input$chromnum] * boundary, y0 = 0.17,
                     y1 = maxlim + 0.1, col = "#00000080", lty = 3, lwd = 0.5)
          }
          stagger <- stagger + 1
          if (stagger > 6) {
            stagger <- 1
          }
        }
      }
      
      if (input$normline) {
        abline(h = 1, col = "black", lty = 2)
      }
    }
  }

  draw_snp_plot <- function() {
    if (input$type != 1 & input$showsnps == 1) {
      query <- "SELECT chrom, pos, ncount, nvaf, tcount, tvaf from "
      query <- paste0(query, "GENOM_SNPDIFF where sample_id = '")
      query <- paste0(query, input$sample,"'")
      genome_snpdiff <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
      #cast
      genome_snpdiff$pos <- as.numeric(genome_snpdiff$pos)
      genome_snpdiff$nvaf <- as.numeric(genome_snpdiff$nvaf)
      genome_snpdiff$tvaf <- as.numeric(genome_snpdiff$tvaf)
      genome_snpdiff$chrom <- as.character(genome_snpdiff$chrom)
      #limit to only heterozygous
      genome_snpdiff <- genome_snpdiff[genome_snpdiff$nvaf < 0.6 & 
                                       genome_snpdiff$nvaf > 0.4, ]

      #remove the chr prefix
      genome_snpdiff$chrom <- gsub("chr", "", genome_snpdiff$chrom)
      #remove the sex chromosomes and where the mean == 0
      genome_snpdiff$chrom[genome_snpdiff$chrom %in% "X"] <- "23"
      genome_snpdiff <- genome_snpdiff[genome_snpdiff$chrom %in% 
                                       input$chromnum, ]
      genome_snpdiff <- genome_snpdiff[order(genome_snpdiff$pos), ]

      #hide blacklist genes unless specified otherwise
      if (!input$show_blacklist) {
        for (i in 1:nrow(blacklist)) {
          genome_snpdiff <- 
            genome_snpdiff[!(genome_snpdiff$pos >= blacklist$V2[i] & 
                             genome_snpdiff$pos <= blacklist$V3[i] & 
                             genome_snpdiff$chrom %in% blacklist$V1[i]), ]
        }
      }

      xaxis <- genome_snpdiff$pos
      boundary <- 100000

      plot(NA, axes = FALSE, ylim = c(0, 1), xlim = c(0, boundary), 
           xlab = paste("Chromosome ", input$chromnum, sep = ""), ylab = "VAF")

      abline(h = 0.5, lty = 2)
      axis(side = 2, las = 1)

      #calculate the boundaries
      draw_list_min <- c()
      draw_list_max <- c()
      draw_list_color <- c()
      draw_list_names <- c()
      group <- genome_bed[genome_bed$chrom == input$chromnum, ]

      pal <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#757780",
               "#b15928")
      pal <- rep(pal, 1000)
      pal <- pal[1:dim(group)[1]]
      group <- cbind(group, pal)

      #add padding
      group$start <- group$start - 100
      group$end <- group$end + 100
      mypalette <- rep("red", dim(genome_snpdiff)[1])
      for (i in 1:dim(group)[1]) {
        if (length(which(genome_snpdiff$pos < group$end[i] & 
                         genome_snpdiff$pos > group$start[i]))) {
          mypalette[which(genome_snpdiff$pos < group$end[i] &
                          genome_snpdiff$pos > group$start[i])] <-
            as.character(group$pal[i])
          draw_list_min <- c(draw_list_min, 
                             min(which(genome_snpdiff$pos < group$end[i] & 
                                       genome_snpdiff$pos > group$start[i])))
          draw_list_max <- c(draw_list_max, 
                             max(which(genome_snpdiff$pos < group$end[i] & 
                                         genome_snpdiff$pos > group$start[i])))
          draw_list_color <- c(draw_list_color, as.character(group$pal[i]))
          draw_list_names <- c(draw_list_names, group$id[i])
        }
      }

      z <- genome_snpdiff
      pandq <- c()
      
      if (input$pandq) {
        pandq <- pq$pqbound[as.numeric(input$chromnum)]
        pandq - (pandq / chrom_sizes[input$chromnum]) * boundary
      }
      if (input$continuous) {
        spacer <- floor(chrom_sizes[input$chromnum] / 
                        length(genome_snpdiff$pos))
        genome_snpdiff$pos <- 1:length(genome_snpdiff$pos) * spacer
        genome_snpdiff$pos <- genome_snpdiff$pos - spacer / 2
      }
      if (input$pandq & input$continuous) {
        spacer <- floor(chrom_sizes[input$chromnum] / 
                        length(genome_snpdiff$pos))
        pandq <- pq$pqbound[as.numeric(input$chromnum)]
        pandq <- ((genome_snpdiff$pos[min(which(z$pos > pandq))] - spacer / 2) / 
                    chrom_sizes[input$chromnum]) * boundary
      }

      abline(v = pandq, col = "#FFFFFF", lty = 2)

      segments((genome_snpdiff$pos / chrom_sizes[input$chromnum]) * boundary, 
               genome_snpdiff$nvaf,
               (genome_snpdiff$pos / chrom_sizes[input$chromnum]) * boundary, 
               genome_snpdiff$tvaf)
      points((genome_snpdiff$pos / chrom_sizes[input$chromnum]) * boundary, 
             genome_snpdiff$nvaf, pch = 22, cex = 1.35, col = "black", 
             bg = "black")
      points((genome_snpdiff$pos / chrom_sizes[input$chromnum]) * boundary, 
             genome_snpdiff$tvaf, pch = 21, cex = 1.35, col = "black", 
             bg = mypalette)
    
      if (input$normline) {
        abline(h = c(0, 0.5, 1), col = "gray", lty = 2)
      }
    }
  }

  get_mutations <- function(all_chroms = 0) {
    #load the raw vcf
    if (input$all_mutation_information) {
      query <- "SELECT * from GENOM_MUTATIONS where sample_id = '"
      query <- paste0(query, input$sample,"'")
      genome_muts <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
      genome_muts$mut_id <- NULL
      genome_muts$sample_id <- NULL
    } else { #otherwise, parse out the data
      #load the data
      #get the samples that are in the database
      query <- "SELECT distinct mut_tool from GENOM_MUTATIONS "
      query <- paste0(query, "where sample_id = '",input$sample,"'")
      mut_tools <- dbGetQuery(dbhandle, query, stringsAsFactors=FALSE)
      mut_tools <- mut_tools[,1]

      genome_muts <- c()
      for (i in mut_tools) {
        vcf <- read_vcf_for_sample(dbhandle, input$sample, i)
        vcf <- make_table(vcf)
        genome_muts <- rbind(genome_muts, vcf)
      }

      if (input$hide_blank_mutations) {
        genome_muts <- genome_muts[!(genome_muts$symbol %in% ""), ]
      }


    }

    #remove chrom prefix
    genome_muts$chrom = gsub("chr", "", genome_muts$chrom)

    if(!input$type && all_chroms == 0) {
      genome_muts <- genome_muts[genome_muts$chrom %in% input$chromnum, ]
    }
    if (!input$show_failed_mutations) {
      genome_muts <- genome_muts[genome_muts$filter %in% "PASS", ]
    }
    #this hides genes on the blacklist
    if (!input$show_blacklist) {
      for (i in 1:nrow(blacklist)) {
        genome_muts <- genome_muts[!(genome_muts$pos >= blacklist$V2[i] & 
          genome_muts$pos <= blacklist$V3[i] &
          genome_muts$chrom %in% blacklist$V1[i]),]
      }
    }

    genome_muts
  }

  output$plot <- renderPlot(
  {
    draw_copy_number_plot()
  })

  output$plot2 <- renderPlot(
  {
    draw_snp_plot()
  })

  output$rcircos <- renderPlot(width=750, height=750, {
    #do struct vars exist for this sample?
    query <- "select * from GENOM_STRUCTVAR where sample_id = '"
    query <- paste0(query, input$sample,"'")
    structvar <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
    #do mutations exist for this sample?
    mutations <- get_mutations(all_chroms = 1)
    #do copy number regions exist for this sample?
    query <- "SELECT chrom, start, end, id, normal_mean, tumor_mean, "
    query <- paste0(query, "tumor_depth from GENOM_DEPTH where sample_id = '")
    query <- paste0(query, input$sample,"'")
    copynumber <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)

    circos_cn_data <- data.frame()
    circos_link_data <- data.frame()
    circos_mut_data <- data.frame()

    chr.exclude <- "Y"
    if (!input$type) {
      chr.exclude <- c(as.character(1:22), "X", "Y")
      chr.exclude <- setdiff(chr.exclude, input$chromnum)
    }

    data(UCSC.HG19.Human.CytoBandIdeogram)
    UCSC.HG19.Human.CytoBandIdeogram$Chromosome <- 
      gsub("chr", "" , UCSC.HG19.Human.CytoBandIdeogram$Chromosome)

    if (nrow(structvar) > 0) {
      if (!input$show_failed_structvars) {
        structvar <- structvar[structvar$filter %in% "PASS", ]
      }

      #cast as integer
      structvar$pos <- as.numeric(structvar$pos)

      #if it's not bnd the end should be in info
      nonbnd <- structvar[!(grepl("BND",structvar$id)), ]
      ends <- gsub("END=([0-9]*?);.*", "\\1", nonbnd$info)
      ends <- as.numeric(ends)

      circos_nonbnd <- data.frame(breakpoint1_chrom = nonbnd$chrom,
                                  breakpoint1_start = nonbnd$pos, 
                                  breakpoint1_end = nonbnd$pos, 
                                  breakpoint2_chrom = nonbnd$chrom, 
                                  breakpoint2_start = ends, 
                                  breakpoint2_end = ends, 
                                  stringsAsFactors=FALSE)

      #is this a BND event?
      structvar <- structvar[grepl("BND",structvar$id), ]
      buffer <- gsub("A|C|G|T|\\[|\\]", "", structvar$alt)
      buffer <- strsplit(buffer, ":")
      chrom <- unlist(lapply(buffer, "[", 1))
      pos <- as.numeric(lapply(buffer, "[", 2))

      circos_link_data <- data.frame(breakpoint1_chrom = structvar$chrom,
                                     breakpoint1_start = structvar$pos, 
                                     breakpoint1_end = structvar$pos,
                                     breakpoint2_chrom = chrom,
                                     breakpoint2_start = pos,
                                     breakpoint2_end = pos,
                                     stringsAsFactors=FALSE)

      circos_link_data <- rbind(circos_link_data, circos_nonbnd)

      #don't keep the non-standard chroms
      circos_link_data <- circos_link_data[circos_link_data$breakpoint1_chrom %in% c(as.character(1:22), "X", "Y"),]
      circos_link_data <- circos_link_data[circos_link_data$breakpoint2_chrom %in% c(as.character(1:22), "X", "Y"),]

      if (!input$type) {
        if (input$show_inter_svs) {
          circos_link_data <-
            circos_link_data[circos_link_data$breakpoint1_chrom %in% 
                               input$chromnum | 
                               circos_link_data$breakpoint2_chrom %in% 
                               input$chromnum, ]
          involved <- unique(c(circos_link_data$breakpoint1_chrom, 
                               circos_link_data$breakpoint2_chrom))
          chr.exclude <- setdiff(chr.exclude, involved)
        } else {
          circos_link_data <-
            circos_link_data[circos_link_data$breakpoint1_chrom %in% 
                               input$chromnum & 
                               circos_link_data$breakpoint2_chrom %in% 
                               input$chromnum, ]
        }
      }
    }

    if (nrow(mutations) > 0) {
      circos_mut_data <- data.frame(breakpoint1_chrom = mutations$chrom, 
                                    breakpoint1_start = mutations$pos, 
                                    breakpoint1_end = mutations$pos, 
                                    symbol = mutations$symbol, 
                                    PlotColor="black", stringsAsFactors=FALSE)

      if (!input$type) {
        circos_mut_data <- circos_mut_data[!(circos_mut_data$breakpoint1_chrom
                                             %in% chr.exclude), ]
      }
    }

    if (nrow(copynumber) > 0) {
      #cast to proper data types
      copynumber$start <- as.numeric(copynumber$start)
      copynumber$end <- as.numeric(copynumber$end)
      copynumber$normal_mean <- as.numeric(copynumber$normal_mean)
      copynumber$tumor_mean <- as.numeric(copynumber$tumor_mean)

      #if the chr prefix is on chrom, remove it
      copynumber$chrom <- gsub("chr", "", copynumber$chrom)

      #remove the Y chromosome and where the mean == 0
      copynumber <- copynumber[copynumber$tumor_mean != 0, ]
      copynumber <- copynumber[copynumber$chrom != "Y", ]

      #this hides genes on the blacklist
      if (!input$show_blacklist) {
        for (i in 1:nrow(blacklist)) {
          copynumber <- copynumber[!(copynumber$start >= blacklist$V2[i] &
                                       copynumber$end <= blacklist$V3[i] & 
                                       copynumber$chrom %in% blacklist$V1[i]), ]
        }
      }
    
      #set up x-axis
      xaxis <- copynumber$start + 
        cumulative_chrom_sizes[as.numeric(copynumber$chrom)]

      #if tumor depth is set in the database, use those values initially
      if (!input$chromnorm && 
          length(copynumber$tumor_depth[!is.na(copynumber$tumor_depth)]) > 0) {
        copynumber <- copynumber[!is.na(copynumber$tumor_depth),]
        xaxis <- copynumber$start + 
          cumulative_chrom_sizes[as.numeric(copynumber$chrom)]
        yaxis <- copynumber$tumor_depth
      } else {    
        chromlist <- get_chrom_normalization_input()
        n_fac <- mean(copynumber$normal_mean[copynumber$chrom %in% chromlist])
        t_fac <- mean(copynumber$tumor_mean[copynumber$chrom %in% chromlist])
        yaxis <- (copynumber$tumor_mean / t_fac) / 
          (copynumber$normal_mean / n_fac)
      }

      circos_cn_data <- data.frame(breakpoint1_chrom = copynumber$chrom, 
                                   breakpoint1_start = copynumber$start,
                                   breakpoint1_end = copynumber$end, 
                                   ratio = yaxis, 
                                   PlotColor = "black",
                                   point.type = 16, stringsAsFactors = FALSE)
      circos_cn_data$PlotColor[circos_cn_data$ratio >= 1.25] <- "red"
      circos_cn_data$PlotColor[circos_cn_data$ratio <= 0.75] <- "blue"

      if (!input$type) {
        circos_cn_data <-
          circos_cn_data[!(circos_cn_data$breakpoint1_chrom %in% chr.exclude), ]
      }
    }

    tracks.inside <- 5
    tracks.outside <- 0
    RCircos.Set.Core.Components(cyto.info = UCSC.HG19.Human.CytoBandIdeogram, 
                                chr.exclude, tracks.inside, tracks.outside)
    RCircos.Set.Plot.Area()
    RCircos.Chromosome.Ideogram.Plot()

    rcircos.params <- RCircos.Get.Plot.Parameters()
    rcircos.params$point.type <- 16
    rcircos.params$point.size <- 1
    rcircos.params$text.size <- 1
    RCircos.Reset.Plot.Parameters(rcircos.params)
    RCircos.List.Plot.Parameters()

    if (nrow(circos_link_data) > 0) {
      RCircos.Link.Plot(circos_link_data, track.num = 1, by.chromosome = TRUE,
                        genomic.columns = 3, is.sorted = FALSE)
    }
    if (nrow(circos_cn_data) > 0) {
      RCircos.Scatter.Plot(circos_cn_data, 4, track.num = 2, min.value = 0,
                           max.value = 2)
    }
    if (nrow(circos_mut_data) > 0 && input$hide_circos_muts == 0) {
      RCircos.Gene.Connector.Plot(circos_mut_data, 3, "in")
      RCircos.Gene.Name.Plot(circos_mut_data, 4, 4, "in")
    }
  })

  output$AutoCN <- 
    DT::renderDataTable(DT::datatable(rownames = FALSE, 
                                      options = list(paging = FALSE, 
                                                     searching = TRUE), {
    # load the data
    query <- "select chrom, start, end, id, normal_mean, tumor_mean, "
    query <- paste0(query, "tumor_depth from GENOM_DEPTH where sample_id = '")
    query <- paste0(query, input$sample,"'")
    genome_depth <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
    #cast to proper data types
    genome_depth$start <- as.numeric(genome_depth$start)
    genome_depth$end <- as.numeric(genome_depth$end)
    genome_depth$normal_mean <- as.numeric(genome_depth$normal_mean)
    genome_depth$tumor_mean <- as.numeric(genome_depth$tumor_mean)
    #remove chr prefix if it exists
    genome_depth$chrom <- gsub("chr","",genome_depth$chrom)

    #remove the sex chromosomes and where the mean == 0
    genome_depth <- genome_depth[genome_depth$tumor_mean != 0, ]
    genome_depth <- genome_depth[genome_depth$chrom != "Y", ]
    genome_depth$chrom[genome_depth$chrom %in% "X"] <- "23"

    #hide blacklist genes unless specified otherwise
    if (!input$show_blacklist) {
      for (i in 1:nrow(blacklist)) {
        genome_depth <- 
          genome_depth[!(genome_depth$start >= blacklist$V2[i] &
                         genome_depth$end <= blacklist$V3[i] & 
                         genome_depth$chrom %in% blacklist$V1[i]), ]
      }
    }

    #load normalization, if there's no normalization flag, 
    #generate a base normalization
    if (!input$chromnorm && 
        length(genome_depth$tumor_depth[!is.na(genome_depth$tumor_depth)]) > 0) {
      genome_depth <- genome_depth[!is.na(genome_depth$tumor_depth), ]
      norm_depth <- genome_depth$tumor_depth
    } else {    
      chromlist <- get_chrom_normalization_input()
      n_fac <- mean(genome_depth$normal_mean[genome_depth$chrom %in% chromlist])
      t_fac <- mean(genome_depth$tumor_mean[genome_depth$chrom %in% chromlist])
      norm_depth <- (genome_depth$tumor_mean / t_fac) / 
        (genome_depth$normal_mean / n_fac)
    }
    
    #regions above this are considered a gain
    upper_bound <- 1.25
    #regions below this are considered a loss
    lower_bound <- 0.75

    genome_depth$norm_depth <- norm_depth

    genome_depth$cn_status <- "normal"
    genome_depth$cn_status[norm_depth >= upper_bound] <- "gain"
    genome_depth$cn_status[norm_depth <= lower_bound] <- "loss"

    genome_depth$normal_mean <- NULL
    genome_depth$tumor_mean <- NULL
    genome_depth$tumor_depth <- NULL

    #only show the chromosome we're looking at if 
    #a specific chromosome is indicated
    if (!input$type) {
      genome_depth <- genome_depth[genome_depth$chrom %in% input$chromnum, ]
    }

    genome_depth
  }))

  output$AutoCNGroups <- 
    DT::renderDataTable(DT::datatable(rownames = FALSE, 
                                      options = list(paging = FALSE, 
                                                     searching = TRUE), {
    # load the data
    query <- "select chrom, start, end, id, normal_mean, tumor_mean, "
    query <- paste0(query, "tumor_depth from GENOM_DEPTH where sample_id = '")
    query <- paste0(query, input$sample,"'")
    genome_depth <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)

    #cast to proper data types
    genome_depth$start <- as.numeric(genome_depth$start)
    genome_depth$end <- as.numeric(genome_depth$end)
    genome_depth$normal_mean <- as.numeric(genome_depth$normal_mean)
    genome_depth$tumor_mean <- as.numeric(genome_depth$tumor_mean)

    #remove the sex chromosomes and where the mean == 0
    genome_depth <- genome_depth[genome_depth$tumor_mean != 0, ]
    genome_depth <- genome_depth[genome_depth$chrom != "Y", ]
    genome_depth$chrom[genome_depth$chrom %in% "X"] = "23"

    #hide blacklist genes unless specified otherwise
    if (!input$show_blacklist) {
      for (i in 1:nrow(blacklist)) {
        genome_depth <- 
          genome_depth[!(genome_depth$start >= blacklist$V2[i] & 
                         genome_depth$end <= blacklist$V3[i] & 
                         genome_depth$chrom %in% blacklist$V1[i]), ]
      }
    }

    #load the curator normalization
    curator_norm <- 1:22
    n_fac <- mean(genome_depth$normal_mean[genome_depth$chrom %in% curator_norm]) 
    t_fac <- mean(genome_depth$tumor_mean[genome_depth$chrom %in% curator_norm])

    norm_depth <- 
      (genome_depth$tumor_mean / t_fac) / (genome_depth$normal_mean / n_fac)

    # load snp data
    query <- "SELECT chrom, pos, ncount, nvaf, tcount, tvaf from GENOM_SNPDIFF "
    query <- paste0(query, "where sample_id = '", input$sample,"'")
    snp <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
    #cast to proper data types
    snp$nvaf <- as.numeric(snp$nvaf)
    snp$tvaf <- as.numeric(snp$tvaf)
    #calc diff
    snp$snp_diff <- abs(snp$nvaf - snp$tvaf)
    
    #regions above this are considered a gain
    upper_bound <- 1.25
    #regions below this are considered a loss
    lower_bound <- 0.75

    genome_depth$norm_depth <- norm_depth

    genome_depth$cn_status <- "normal"
    genome_depth$cn_status[norm_depth >= upper_bound] <- "gain"
    genome_depth$cn_status[norm_depth <= lower_bound] <- "loss"

    #define groups
    results <- c()
    for (i in 1:dim(genome_bed)[1]) {
      cn_status <- 
        genome_depth$cn_status[(genome_depth$start >= genome_bed$start[i] & 
                                genome_depth$end <= genome_bed$end[i] &
                                genome_depth$chrom %in% genome_bed$chrom[i])]
      snp_status <- 
        snp$snp_diff[(snp$pos >= genome_bed$start[i] & 
                      snp$pos <= genome_bed$end[i] &
                      snp$chrom %in% genome_bed$chrom[i])]
      drift <- 0
      if (length(snp_status) > 0) {
        drift <- median(snp_status)
      } 

      #we want at least two consecutive regions of loss/gain
      cn_regions_count <- rle(cn_status)[[1]]
      cn_regions_type <- rle(cn_status)[[2]]
      #get the regions that aren't normal
      if (any(cn_regions_count[cn_regions_type != "normal"] >= 2)) {
        cn_status <- cn_status[cn_status != "normal"]
      } else {
        cn_status <- "normal"
      }
      cn_status <- unique(cn_status)
      if(cn_status == "normal" & drift >= 0.2) {
        cn_status = "CNN_LOH"
      }
      if(cn_status == "gain" & drift >= 0.3) {
        cn_status = "GAIN_LOH"
      }
      cn_status <- paste(unique(cn_status),collapse=",")
      results <- rbind(results, c(genome_bed$id[i], cn_status))
    }

    results
  }))

  # Generate mutations information
  output$mutations_table <- 
    DT::renderDataTable(DT::datatable(rownames = FALSE, 
                                      options = list(paging = FALSE, 
                                                    searching = TRUE), {
    get_mutations()
  }))

  # Generate manta information
  output$manta_table <- 
    DT::renderDataTable(DT::datatable(rownames = FALSE, 
                                      options = list(paging = FALSE, 
                                                     searching = FALSE), {
    query <- paste0("SELECT * from GENOM_STRUCTVAR where sample_id = '",
                    input$sample,"'")
    genome_struct <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)
    genome_struct$sv_id <- NULL
    genome_struct$sample_id <- NULL

    if (!input$show_failed_structvars) {
      genome_struct = genome_struct[genome_struct$filter %in% "PASS", ]
    }

    #if the chr prefix is on chrom, remove it
    genome_struct$chrom <- gsub("chr", "", genome_struct$chrom)

    genome_struct$chrom[genome_struct$chrom %in% "X"] <- "23"
    #only show the chromosome we're looking at 
    #if a specific chromosome is indicated
    if (!input$type) {
      genome_struct <- genome_struct[genome_struct$chrom %in% input$chromnum, ]
    }

    genome_struct
  }))

  # Generate debug information
  output$debug <- renderText({
    paste(get_chrom_normalization_input(), collapse = "\n")
  })

  # Generate metrics information
  output$metrics_table <- 
    DT::renderDataTable(DT::datatable(rownames = FALSE, 
                                      options = list(paging = FALSE, 
                                                     searching = FALSE), {
    #load the data
    query <- paste0("SELECT * from genom_depth where sample_id = '",
                    input$sample,"'")
    genome_struct <- dbGetQuery(dbhandle, query, stringsAsFactors = FALSE)

    #hide blacklist by default
    if (!input$show_blacklist) {
      for (i in 1:nrow(blacklist)) {
        genome_struct <- genome_struct[!(genome_struct$start >= blacklist$V2[i] & 
                genome_struct$end <= blacklist$V3[i] &
                genome_struct$chrom %in% blacklist$V1[i]), ]
      }
    }

    res <- data.frame(mean_normal_depth = mean(genome_struct$normal_mean), 
                      mean_tumor_depth = mean(genome_struct$tumor_mean))
    res
  }))

  output$downloadCopyNumberPlot <- downloadHandler(
    filename = function() { 
      if(input$type == 1)
        paste(input$sample, "_copynumber.pdf", sep = "")
      else
        paste(input$sample, "_chrom_", input$chromnum, 
              "_copynumber.pdf", sep = "")
    },
    content = function(file) {
      pdf(h = 5, w = 18, file = file)
      draw_copy_number_plot()
      dev.off()
    }
  )

  output$downloadSNPPlot <- downloadHandler(
    filename = function() 
    { 
      if(input$type == 1)
        paste(input$sample, "_snp.pdf", sep = "")
      else
        paste(input$sample, "_chrom_", input$chromnum, "_snp.pdf", sep = "")
    },
    content = function(file) {
      pdf(h = 5, w = 18, file = file)
      draw_snp_plot()
      dev.off()
    }
  )
})

shadowtext <- function(x, y = NULL, labels, col = 'white', bg = 'black', 
                       theta = seq(0, 2 * pi, length.out = 50), r = 0.1, ... ) {
    xy <- xy.coords(x, y)
    xo <- r * strwidth('A')
    yo <- r * strheight('A')

    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text(xy$x + cos(i) * xo, xy$y + sin(i) * yo, labels, col = bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col = col, ... )
}

# returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

# returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
 