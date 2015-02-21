require(stringr)
require(RCurl)
require(rentrez)

genbankCheckExists <- function(accession){

  if(accession != "") { 
    ## get the whole annotation file
    URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                 accession, "&rettype=gb", sep = "")
    txt <- try(scan(file = URL, what = "", sep = "\n", quiet = TRUE) )
  }
    
  #if it's longer than 2, then it's a real Accession number
  return(length(txt)>2) 
}

## function returns number of sequences
## requires the virus name (as specified in Taxonomy browser) and can take a common acronym, e.g. DENV-1
numSeqByRules <- function(virusName, virusCode=NA, humanOnly=F) {
  
  virusName <- gsub(" ", "+", virusName)
  
  if(humanOnly == T) {
    humanSection <- ""
  } else {
    humanSection <- "+AND+(human+OR+Homo+sapiens)+NOT+homo+sapiens%5BOrganism%5D"  #add arguments for sorting human hosts
  }
  
  if(virusCode == "") {
    virusURL <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D)", humanSection)
  } else {
    virusURL <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D+OR+", virusCode, "%5BAll+Fields%5D)", humanSection) 
  }

  URL <- sprintf("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=%s&RetMax=500000", virusURL)
  
  txt <- scan(file=URL, what="", sep="\n")
  accCodes <- txt[grep("<Id>", txt)]
  accCodes <- sub("<Id>", "", accCodes)
  accCodes <- sub("</Id>", "", accCodes)
  accCodes <- sub("\t", "", accCodes)
  numSeqs <- length(accCodes)
  
  return(numSeqs)
}

accCodesByRules <- function(virusName, virusCode=NA, humanOnly=F) {
  
  virusName <- gsub(" ", "+", virusName)
  
  if(humanOnly == T) {
    humanSection <- ""
  } else {
    humanSection <- "+AND+(human+OR+Homo+sapiens)+NOT+homo+sapiens%5BOrganism%5D"  #add arguments for sorting human hosts
  }
  
  if(virusCode == "") {
    virusURL <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D)", humanSection)
  } else {
    virusURL <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D+OR+", virusCode, "%5BAll+Fields%5D)", humanSection) 
  }
  
  URL <- sprintf("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=%s&RetMax=500000", virusURL)
  
  txt <- scan(file=URL, what="", sep="\n")
  accCodes <- txt[grep("<Id>", txt)]
  accCodes <- sub("<Id>", "", accCodes)
  accCodes <- sub("</Id>", "", accCodes)
  accCodes <- sub("\t", "", accCodes)
  
  return(accCodes)
}

# virusName <- virusNames[i]
# virusCode <- acronym[i]
# segmentName <- segment[i] 
# minLength <- minLength[i] 
# maxLength <- maxLength[i]

# virusName <- virusNames[i]
# virusCode <- acronym[i]
# segmentName <- segment[i]
# humanOnly <- humanOnly[i]
# minLength <- minLength[i]
# maxLength <- maxLength[i]

accCodesByRulesFullLength <- function(virusName, virusCode=NA, segmentName=NA, humanOnly=F, minLength, maxLength) {
  
  virusName <- gsub(" ", "+", virusName)
  segmentName <- gsub(" ", "+", segmentName)
  
  nameSection <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D)")
  sizeSection <- paste0("+", minLength, "%3A", maxLength, "%5Bslen%5D")
  
  if(humanOnly == T) {
    humanSection <- ""
  } else {
    humanSection <- "+AND+(human+OR+Homo+sapiens)+NOT+homo+sapiens%5BOrganism%5D"  #add arguments for sorting human hosts
  }
  
  if(virusCode == "" & segmentName == "") {
    virusURL <- paste0(nameSection, sizeSection, humanSection)
  } else {
    if(segmentName == "" & virusCode != "") { 
      nameAndCode <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D+OR+", virusCode, "%5BAll+Fields%5D)")
      virusURL <- paste0(nameAndCode, sizeSection, humanSection) 
    } else {
      if(segmentName != "" & virusCode == "") {
        segmentSection <- paste0(segmentName, "%5BAll+Fields%5D")
        virusURL <- paste0(nameSection, segmentSection, sizeSection, humanSection) 
      } else {
        if(segmentName != "" & virusCode != "") {
          segmentSection <- paste0(segmentName, "%5BAll+Fields%5D")
          nameAndCode <- paste0("(", virusName, "%5BOrganism%5D+OR+", virusName, "%5BAll+Fields%5D+OR+", virusCode, "%5BAll+Fields%5D)")
          virusURL <- paste0(nameAndCode, segmentSection, sizeSection, humanSection)     
        }
      }
    }
  }
  
  URL <- sprintf("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=%s&RetMax=500000", virusURL)
  
  txt <- scan(file=URL, what="", sep="\n")
  accCodes <- txt[grep("<Id>", txt)]
  accCodes <- sub("<Id>", "", accCodes)
  accCodes <- sub("</Id>", "", accCodes)
  accCodes <- sub("\t", "", accCodes)
  
  return(accCodes)
}

genbankMetadata <- function(accCodes) {
  
  metadata <- c()
  
  if(length(accCodes) > 500) {     #chopUp

    numChunks <- ceiling(length(accCodes)/500)
    maxCodes <- length(accCodes)

   for(k in 1:numChunks) {
     
     if(k*500 > maxCodes) {
       #handle last case specially
       accCodeChunk <- accCodes[(k*500-499):maxCodes]
       accList <- c()
       for(j in accCodeChunk) {
         accList <- paste(accList, j, sep=",")
       }
     } else {
       #send 500
       accCodeChunk <- accCodes[(k*500-499):(k*500)]
       accList <- c()
       for(j in accCodeChunk) {
         accList <- paste(accList, j, sep=",")
       }
     }
     #collect result
     metadata <- rbind(metadata, genbankQueryFull(accList, accCodeChunk))
   }
    
  } else {
    accList <- c()
    for(j in accCodes) {
      accList <- paste(accList, j, sep=",")
    }
  metadata <- genbankQueryFull(accList, accCodes)  
  } 
  
  return(metadata)
}

genbankQueryFull <- function(accList, accCodeChunk){

  ## get the whole annotation file
  URL <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&retmode=text&rettype=gb&id=", accList)
  txt <- getURL(url = URL)
  
  #split into an entry for a sequence
  entries <- str_split(txt, "\n//")  #list of length 1, with number of entriess = num sequences
  
  #results table
  virusName <- c()
  seqLength <- c()
  gene <- c()
  collectionDate <- c()
  location <- c()
  note1 <- c()
  note2 <- c()
  
  for(i in 1:length(accCodeChunk)) {
    
    #split entry into lines
    entry <- str_split(entries[[1]][i], "\n")
    
    ## some general parsing
    entry <- gsub("^[[:blank:]]*","", entry[[1]]) 
    entry <- entry[entry != ""]
    
    ## get virus name field ##
    virus <- entry[grep("SOURCE", entry)]
    virus <- sub(".*SOURCE\\s*","", virus)
    virus <- sub("[^(]*[(]", "", virus)
    virus <- sub("[)]$","", virus)
    
    virusName[i] <- virus
    
    ## get sequence length ##
    seq.length <- sub("[[:alnum:]]*[[:blank:]]*[[:alnum:]_]*[[:blank:]]*","", entry[1]) #seq.length always on first line
    seq.length <- as.integer(sub("[[:blank:]]*bp.*","", seq.length))
    
    seqLength[i] <- seq.length
    
    ## get collection date ##
    dates <- entry[grep("/collection_date", entry , ignore.case=TRUE)]
    dates <- sub(".*/[Cc]ollection_date[[:blank:]]*[=:][[:blank:]]*","", dates)
    dates <- gsub("\"","", dates)
    dates <- sub(":.*","", dates)
    if(length(dates)==0) dates <- NA
    
    collectionDate[i] <- dates
    
    ## get country ##
    country <- entry[grep("/country", entry, ignore.case=TRUE)]
    country <- sub(".*/[Cc]ountry[[:blank:]]*[=:][[:blank:]]*","", country)
    country <- gsub("\"","", country)
    country <- sub(":.*","", country)
    if(length(country)==0)  country <- NA
    
    location[i] <- country
    
    ## get gene ##
    gene <- entry[grep("/gene", entry, ignore.case=TRUE)]
    gene <- sub(".*/[Cc]ountry[[:blank:]]*[=:][[:blank:]]*","", gene)
    gene <- gsub("\"","", gene)
    gene <- sub(":.*","", gene)
    gene <- sub("/gene=", "", gene)
    if(length(gene)==0)  gene <- NA
    
    if(length(gene) > 1) {
      gene <- "multiple"
    }
    gene[i] <- gene
    
    ## get notes ##
    note <- entry[grep("note", entry, ignore.case=TRUE)]
    note <- sub(".*note","", note)
    note <- sub(".*history","", note)
    note <- gsub("[=:\"[:blank:]]","", note)
    if(length(note)==0) {
      note1[i] <- NA
      note2[i] <- NA
    } else {
      if(length(note) == 1) {
        note1[i] <- note
        note2[i] <- NA
      } else {
        if(length(note) > 1) {
          note2[i] <- note[2] 
          note1[i] <- note[1]
        } else {
          if(length(note) > 2) {
            print("Too Many Notes")
          } 
        }
      }
    }
  } 
  result <- cbind(accCodeChunk, virusName, seqLength, collectionDate, location, gene,note1, note2)
  result <- as.data.frame(result)
  return(result)
}
