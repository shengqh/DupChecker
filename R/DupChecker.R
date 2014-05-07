lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

##' arrayExpressDownload
##'
##' The function downloads array express raw data from EBI ftp server based on 
##' datasets user provided. Once the compressed raw data is downloaded, 
##' individual target file will be extracted from compressed raw data. The 
##' dataset/count table will be returned.
##'
##' @param datasets the dataset names, for example: c("E-TABM-43", "E-TABM-158") 
##' @param targetDir the target directory to store the datasets
##' @param filePattern the file pattern of the expected data file extracted 
##' from gzipped file
##' @param tar the path to the command to be used in untar function
##' @param overwrite If TRUE, overwrite existing files, otherwise ignore such 
##' files. The equivalent of unzip -o.
##' @return a data frame containing dataset and how many target files in that 
##' dataset 
##' @importFrom RCurl getURL
##' @export
##' @examples 
##' #download three datasets from ArrayExpress website
##' #datatable<-arrayExpressDownload(datasets = c("E-TABM-43", "E-TABM-158"), 
##' #                                targetDir=getwd())
arrayExpressDownload<-function(datasets, 
                               targetDir = getwd(), 
                               filePattern=".CEL$", 
                               tar="internal", 
                               overwrite=FALSE){
  counts<-data.frame(dataset = datasets, count = rep(0, length(datasets)))
  for(dataset in datasets){
    subdir<-file.path(targetDir, dataset)
    dir.create(subdir, showWarnings = FALSE)
    
    dname<-unlist(strsplit(dataset, '-'))[2]
    if(is.na(dname)){
      stop(paste0("Wrong array express dataset name : ", dataset, "\n"))
    }
    
    curl<-paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/"
                 , dname, "/", dataset, "/")
    filenames<-getURL(curl, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    filenames<-strsplit(filenames, "\r*\n")[[1]]
    links<-paste0(curl, filenames)
    localfiles<-paste0(subdir, "/", filenames)
    lapply(c(1:length(links)), function(x){
      if(overwrite || !file.exists(localfiles[x])){
        cat("downloading file ", links[x], "\n")
        download.file(links[x], localfiles[x], method="auto", mode="wb")
      }
      
      if(grepl(".zip$", localfiles[x])){
        unzip(localfiles[x], overwrite=overwrite, exdir=subdir)
      }
    })
    celfiles<-list.files(subdir, filePattern, ignore.case=TRUE)
    counts$count[counts$dataset==dataset] = length(celfiles)
    cat("downloading file for dataset ", dataset, "finished.\n")
  }
  return(counts)
}

##' geoDownload
##'
##' The function downloads GEO raw data from ncbi ftp server based on datasets 
##' user provided. Once the compressed raw data is downloaded, 
##' individual gzipped target file will be extracted from compressed raw data, 
##' and individual target file will be extracted from corresponding 
##' gzipped file. The dataset/count table will be returned.
##'
##' @param datasets the GEO dataset names, for example: c("GSE14333") 
##' @param targetDir the target directory to store the datasets
##' @param filePattern the file pattern of the expected data file 
##'        extracted from gzipped file
##' @param tar the path to the command to be used which is used in 
##'        untar function
##' @return a data frame containing dataset and how many target files in 
##'         that dataset 
##' @importFrom RCurl getURL
##' @importFrom R.utils gunzip
##' @export
##' @examples 
##' #download three datasets from GEO website
##' #datatable<-geoDownload(datasets = c("GSE14333", "GSE13067", "GSE17538"), 
##' #                       targetDir=getwd())
geoDownload<-function(datasets, 
                      targetDir = getwd(), 
                      filePattern=".CEL$", 
                      tar="internal"){
  counts<-data.frame(dataset = datasets, count = rep(0, length(datasets)))
  for(dataset in datasets){
    subdir<-file.path(targetDir, dataset)
    dir.create(subdir, showWarnings = FALSE)
    curl<-paste0("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/"
                 , dataset, "/")
    filenames<-getURL(curl, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    filenames<-strsplit(filenames, "\r*\n")[[1]]
    links<-paste0(curl, filenames)
    localfiles<-paste0(subdir, "/", filenames)
    lapply(c(1:length(links)), function(x){
      if(!file.exists(localfiles[x])){
        cat("downloading file ", links[x], "\n")
        download.file(links[x], localfiles[x], method="auto", mode="wb")
      }
      
      if(grepl(".tar$", localfiles[x])){
        celfiles<-list.files(subdir, filePattern, ignore.case=TRUE)
        if(length(celfiles) == 0){
          gzfiles<-list.files(subdir, ".gz$")
          if(length(gzfiles) == 0){
            cat("de-compress file ", localfiles[x], "\n")
            untar(localfiles[x], exdir=subdir, tar=tar)
          }
        }
        gzfiles<-list.files(subdir, ".gz$")
        if(length(gzfiles) > 0){
          gzfiles<-file.path(subdir, gzfiles)
          for(file in gzfiles){
            cat("de-compress file ", file, "\n")
            gunzip(file, overwrite = TRUE, remove=TRUE)
          }
        }
      }
    })
    celfiles<-list.files(subdir, filePattern, ignore.case=TRUE)
    counts$count[counts$dataset==dataset] = length(celfiles)
    cat("downloading file for dataset ", dataset, "finished.\n")
  }
  return(counts)
}

##' buildFileTable
##' 
##' The function build file table in the subdirectories under root directories 
##' user provided. The result table contains two columns, dataset and filename
##' 
##' @param rootDir the root of directories whose sub directories contains file 
##'        waiting for validation. It can be vector of directories, or just 
##'        one directory 
##' @param subDirPattern the pattern of sub directory name
##' @param filePattern the pattern of file waiting for validation. 
##'        default is "CEL$"
##' @param ignore.case ignore the case difference when list files from sub 
##'        directory using filePattern
##' @return a data frame containing full file name and its corresponding 
##'         dataset, which will be used at validateFile
##' @export
##' @examples 
##' datafile<-buildFileTable(rootDir=getwd())
##' #or
##' datafile<-buildFileTable(rootDir=c(paste0(getwd(), 
##'                          c("/E-TABM-43", "/E-TABM-158") )))
buildFileTable<-function(rootDir, 
                         subDirPattern = NULL, 
                         filePattern = "CEL$", 
                         ignore.case=TRUE){
  lst<-list()
  
  if(is.character(rootDir)){
    dirs<-c(rootDir)    
  }else{
    dirs<-rootDir
  }
  
  for(dir in dirs){
    subdirs<-list.dirs(dir, recursive=FALSE)
    if(!is.null(subDirPattern)){
      basedirs<-basename(subdirs)
      acceptSubDirs<-grepl(subDirPattern, basedirs, perl=TRUE)
      subdirs<-subset(subdirs, acceptSubDirs)
    }
    
    subdirs<-c(dir, subdirs)
    for(subdir in subdirs){
      basedir<-basename(subdir)
      files<-list.files(subdir, filePattern, full.names=TRUE, no..=TRUE, 
                        ignore.case=ignore.case)
      for(file in files){
        lst<-lappend(lst, c(basedir, file))
      }
    }
  }
  
  result = data.frame(dataset = character(length(lst)), 
                      file = character(length(lst)), stringsAsFactors=FALSE)
  if(length(lst) > 0){
    for(i in c(1:length(lst))){
      result$dataset[i]<-lst[[i]][1]
      result$file[i]<-lst[[i]][2]
    }
  }
  
  return (result)
}

##' validateFile
##' 
##' The function calculate MD5 fingerprint for each file in table and then 
##' check to see if any two files have same MD5 fingerprint. The files with 
##' same fingerprint will be treated as duplication. The function will return 
##' a table contains all duplicated files and datasets.
##' 
##' @param fileTable a table with column name "dataset" and "file", 
##'        here column "file" should contain full name of file.
##' @param saveMd5File if calculated MD5 fingerprint should be save to 
##'        local file
##' @return a list contains two tables. One is the table contains three 
##'         columns: "dataset", "file" and "md5". Another one is the 
##'         duplication table whose row indicates MD5 fingerprint and 
##'         whose column indicates dataset, table cell indicates the 
##'         corresponding filename.
##' @importFrom tools md5sum
##' @export
##' @examples 
##' datafile<-buildFileTable(rootDir=getwd())
##' if(nrow(datafile) > 0){
##'   result<-validateFile(datafile)
##'   if(result$hasdup){
##'     duptable<-result$duptable
##'     write.csv(duptable, file="duptable.csv")
##'   }
##' }
validateFile<-function(fileTable, saveMd5File=TRUE){
  cat("calculate and validate files, make sure the directory is readable\n")
  filemd5<-apply(fileTable, 1, function(x){
    celfile<-as.character(x["file"])
    if(!file.exists(celfile)){
      stop(paste0("File not exists : ", celfile))
    }
    md5file<-paste0(celfile, ".md5")
    if(file.exists(md5file)){
      md5<-readLines(md5file, 1, warn=FALSE)
    }else{
      md5<-""
    }
    
    if(nchar(md5) == 0){
      cat("calculate md5 for", celfile, "...\n")
      md5<-md5sum(celfile)
      if(saveMd5File){
        fileConn<-file(md5file)
        writeLines(c(md5), fileConn)
        close(fileConn)
      }
    }
    
    md5
  })
  
  oldtable<-data.frame(fileTable)
  oldtable$md5<-filemd5
  
  md5table<-table(filemd5)
  dupmd5<-names(md5table[md5table > 1])
  
  if(length(dupmd5) > 0){
    dup<-oldtable[oldtable$md5 %in% dupmd5,]
    warning(paste0(nrow(dup), " entries out of total ", nrow(oldtable), 
                   " entries are duplicated at least once. \n"))
    
    dupdatasets<-unique(as.character(dup$dataset))
    dupdatasets<-dupdatasets[order(dupdatasets)]
    x<-"GSE14333"
    dupdsnames<-unlist(lapply(dupdatasets, function(x){
      totalcount<-table(oldtable$dataset==x)["TRUE"]
      dupcount<-table(dup$dataset==x)["TRUE"]
      paste0(x, "[", dupcount, "/", totalcount, "]")
    }))
    
    restable<-matrix(rep("", length(dupmd5) * length(dupdatasets)), 
                     nrow=length(dupmd5), ncol=length(dupdatasets))
    rownames(restable)<-dupmd5
    colnames(restable)<-dupdsnames
    
    for(i in c(1:length(dupmd5))){
      md5<-dupmd5[i]
      ds<-oldtable[oldtable$md5==md5,]
      for(j in c(1:length(dupdatasets))){
        if(dupdatasets[j] %in% ds$dataset){
          restable[i,j]<-paste0(basename(ds[ds$dataset==dupdatasets[j],]$file),
                                collapse=";")
        }    
      }      
    }
    
    result<-list(filetable=oldtable, duptable=restable, hasdup=TRUE)
  }else{
    result<-list(filetable=oldtable, hasdup=FALSE)
  }
  class(result)<-"FileDup"
  result
}
