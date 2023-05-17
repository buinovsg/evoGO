## This file includes functions for retrieving the annotations from
## external resources and to transform them into the input structures,
## which are required by the evoGO constructor.



##' Retrieve the ontologyIndex object from the obolibrary.
##'
##' @return ontologyIndex object
##' @export
getOBO <- function(goRelease = NULL) {
  url <- ifelse(is.null(goRelease), "http://purl.obolibrary.org/obo/go/go-basic.obo",
                paste0("http://release.geneontology.org/", goRelease, "/ontology/go-basic.obo"))
  filename <- tempfile()
  
  curl::curl_download(url, filename)
  
  go <- ontologyIndex::get_ontology(filename,
                                    extract_tags = "everything",
                                    propagate_relationships = c(
                                      "is_a", "part_of", "regulates",
                                      "positively_regulates", "negatively_regulates"
                                    )
  )
  go
}


##' Create an input data.frame for evoGO constructor function.
##'
##' The data.frame includes "child" - "parent" relationships between the GO ids.
##' @param go an ontologyIndex object returned by getOBO function (also see ontologyIndex package)
##' @param rootName the name of the root elements from which to build
##' the hierarchy.
##'
##' @return a data.frame to be used as graphTable argument in evoGO function.
##' @export
ontologyIndex2graph <- function(go, rootName = "cellular_component") {
  parents <- go$parents
  ids <- go$id
  hasParent <- vapply(go$parents, length, numeric(1)) > 0
  rootIndex <- which(go$name == rootName)
  if (length(rootIndex) == 0) {
    stop(sprintf("no such root element found: %s", rootName))
  }
  if (length(rootIndex) > 1) {
    stop(sprintf("more than 1 element found for the root: %s", rootName))
  }
  root <- go$id[rootIndex]
  with_root <- ontologyIndex::get_descendants(go, root)
  o <- which((go$id %in% with_root) & hasParent)
  
  res <- lapply(o, function(i) {
    data.frame(child = go$id[[i]], parent = go$parents[[i]])
  })
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res
}


##' Get gene-GO mapping from biomaRt
##'
##' @param biomart the name of the biomart database (default "ensembl").
##' @param dataset the name of the dataset (default "hsapiens_gene_ensembl"),
##'   use `biomaRt::listDatasets` to see the available ones.
##'
##' @return an object to use as a geneSets argument for evoGO constructor
##'   function.
##' @export
getGeneSets <- function(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", ensemblRelease = NULL) {
  ensembl <- biomaRt::useEnsembl(biomart = biomart, dataset = dataset, version = ensemblRelease)
  
  # Split the query into chunks to avoid timeouts
  gene_ids <- biomaRt::getBM(
    attributes = c("ensembl_gene_id"),
    mart = ensembl
  )
  gene_ids <- unname(unlist(gene_ids))
  gene_ids <- split(gene_ids, ceiling(seq_along(gene_ids) / 1000))
  
  gene_annot <- lapply(1:length(gene_ids), function(i) {
    message("Retrieving GO to gene mapping... (", i, "/", length(gene_ids), ")")
    mart <- biomaRt::getBM(
      attributes = c("go_id", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = gene_ids[[i]],
      mart = ensembl
    )
  })
  gene_annot <- do.call(rbind, gene_annot)
  
  empty <- gene_annot$go_id == "" | gene_annot$ensembl_gene_id == ""
  gene_annot <- gene_annot[!empty, , drop = FALSE]
  split(gene_annot$ensembl_gene_id, gene_annot$go_id)
}

listGOAnnotations <- function(species = NULL, database = "ensembl", goRelease = NULL,
                              ensemblRelease = NULL, custom = NULL, fileName = NULL,
                              path = NULL, returnLatest = FALSE) {
  # Check inputs-
  #browser()
  assertthat::assert_that(is.null(path) | isValidString(path),
                          msg = "Argument 'path' should be a character vector with length of one"
  )
  if (!is.null(fileName)) {
    assertthat::assert_that(isValidString(fileName),
                            msg = "Argument 'fileName' should be a character vector with length of one"
    )
  } else {
    assertthat::assert_that(is.null(species) | isValidString(species),
                            msg = "Argument 'species' should be a character vector with length of one"
    )
    assertthat::assert_that(isValidString(database),
                            msg = "Argument 'database' should be a character vector with length of one"
    )
    assertthat::assert_that(database == "ensembl",
                            msg = "Only Ensembl database is currently supported"
    )
    assertthat::assert_that(!any(class(tryCatch(as.Date(goRelease), error = function(e) e)) == "error"),
                            msg = "Argument 'goRelease' is expected to be a a character vector with length of one in the format \"2023-04-01\""
    )
    assertthat::assert_that(is.null(ensemblRelease) || (is.numeric(ensemblRelease) && ensemblRelease%%1==0),
                            msg = "Argument 'ensemblRelease' is expected to be an integer"
    )
    assertthat::assert_that(is.null(custom) | isValidString(custom),
                            msg = "Argument 'custom' should be a character vector with length of one"
    )
  }
  # Check path
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
                          msg = paste0("Path ", path, " is not available")
  )
  
  # Get file name
  goRelease <- ifelse(is.null(goRelease), "", paste0("__go", gsub('-','',goRelease)))
  ensemblRelease <- ifelse(is.null(ensemblRelease), "", paste0("__ensembl", ensemblRelease))
  species <- ifelse(is.null(species), "", species)
  custom <- ifelse(is.null(custom), "", custom)
  all_files <- list.files(path)
  
  # Assume newer versions of GO data may only have same or newer versions of annotation and vice versa
  selected_files <- all_files[
    grepl("^evogo__", all_files) &
      grepl("\\.rds$", all_files) &
      grepl(goRelease, all_files) &
      grepl(ensemblRelease, all_files) &
      grepl(species, all_files) &
      grepl(custom, all_files)
  ]
  if (length(selected_files) == 0) {
    message("No matching annotation files found. Run listGOAnnotations() to view available files")
    return(NULL)
  }

    selectedFileInfo <- function(selected_files) {
    selected_files <- gsub(pattern = "\\.rds$", "", selected_files)
    file_info <- data.frame()
    for (file in selected_files) {
      name_split <- strsplit(file, split = "__")[[1]]
      go <- substring(name_split[startsWith(name_split, "go")], 3)
      #address files produced using previous version of evoGO
      if (length(name_split) == 3) {
        ensembl <- NA
        custom <- name_split[3]
      } else {
        ensembl <- name_split[startsWith(name_split, "ensembl")]
        ensembl <- ifelse(length(ensembl) == 0, NA, substring(ensembl, 8))
        custom <- name_split[startsWith(name_split, "custom")]
        custom <- ifelse(length(custom) == 0, NA, substring(custom, 7))
      }
      file_info <- rbind(file_info, c(go, ensembl, custom, file))
    }
    
    names(file_info) <- c("goRelease", "ensemblRelease", "custom", "fileName")
    file_info$goRelease <- as.Date(file_info$goRelease, "%Y%m%d")
    file_info$ensemblRelease <- as.integer(file_info$ensemblRelease)
    
    file_info <- file_info[order(file_info$goRelease, decreasing=TRUE),]
    file_info
    
  }
  
  if (returnLatest == TRUE) {
    file_info <- selectedFileInfo(selected_files)
    if (nrow(file_info) == 1) {
      return(file_info)
    } else {
      latest_go <- max(file_info$goRelease)
      latest_go_files <- file_info[file_info$goRelease %in% latest_go,]
      if (custom != ""){
        return(latest_go_files)
      } else {
        
        latest_ensemlb <- max(na.omit(latest_go_files$ensemblRelease))#, na.rm = TRUE)????????????????
        latest_file <- latest_go_files[latest_go_files$ensemblRelease %in% latest_ensemlb,]
        assertthat::assert_that(nrow(latest_file) == 1,
                                msg = cat(c("Multipe annotation file matches:", "", capture.output(latest_file)
                                            ), sep = "\n")
        )
        
        return(latest_file)
      }
    }
    } else {
      file_info <- selectedFileInfo(selected_files)
    } 
  
  
  file_info
}



##' Retrieve and save the latest annotation for species of interest
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param database database to be used for acquiring of genesets. Only 'ensembl'
##' is currently supported.
##' @param nCores maximum number of cores to be used
##' @param save logical value indicating whether the result should be saved
##' @param path path to a directory where the result will be saved. If NULL, the file
##' will be saved to the package 'extdata' directory
##' @param deletePrevious logical value indicating whether previous annotation for
##' that species should be removed when newer is saved to save disk space. Only
##' used in case results are saved to 'extdata' directory
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest version is retrieved
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl). If not specified,
##' latest version is retrieved
##'
##' @return list of evoGO objects containing the annotation (one per GO domain).
##' @export
getGOAnnotation <- function(species = NULL, database = "ensembl", nCores = 1, save = TRUE,
                            path = NULL, deletePrevious = FALSE, goRelease = NULL,
                            ensemblRelease = NULL, customAnnotation = NULL, customName = NULL) {
  
  # Check inputs
#  browser()
  if (is.null(customAnnotation) & is.null(customName)) {
    assertthat::assert_that(isValidString(database),
                            msg = "'database' should be a character vector with length of one"
    )
    assertthat::assert_that(database == "ensembl",
                            msg = "Only Ensembl database is currently supported"
    )
    assertthat::assert_that(isValidString(species),
                            msg = "'species' should be a character vector with length of one"
    )
    custom = FALSE
  } else if (!is.null(customAnnotation) | !is.null(customName)){
    assertthat::assert_that(!is.null(customAnnotation) & !is.null(customName),
                            msg = "Both 'customAnnotation' and 'customName' must be provided in order to use custom annotation"
    )
    assertthat::assert_that(is.list(customAnnotation) && length(names(customAnnotation)) > 0,
                            msg = paste0("Argument 'customAnnotation' should be a named list where each element is named after a GO term ",
                                         "(eg. \"GO:0000002\") and the element is a character vector containing the names of the annotation terms")
    )
    assertthat::assert_that(isValidString(customName) & !grepl( "__", customName, fixed = TRUE),
                            msg = "'customName' should be a character vector with length of one and should not contain \"__\""
    )
    if (!is.null(species)) {
      message("'species' specified by the user but is ignored when using custom annotation")
    }
    custom = TRUE
  }
  assertthat::assert_that(is.null(path) | isValidString(path),
                          msg = "'path' should be a character vector with length of one"
  )
  assertthat::assert_that(is.null(goRelease) || !any(class(tryCatch(as.Date(goRelease), error = function(e) e)) == "error"),
                          msg = "Argument 'goRelease' is expected to be a a character vector with length of one in the format \"2023-04-01\""
  )
  if (!is.null(goRelease)) {
    assertthat::assert_that(as.Date(goRelease) >= as.Date("2013-04-01"),
                            msg = "GO consortium releases prior to 2013-04-01 are not supported"
    )
  }
  assertthat::assert_that(is.null(ensemblRelease) || (is.numeric(ensemblRelease)),
                          msg = "Argument 'ensemblRelease' is expected to be an integer"
  )
  
  # Get/create save directory
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
    if (!dir.exists(path)) {
      dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    assertthat::assert_that(dir.exists(path), msg = "Provided path does not exist")
  }
  
  if (save) {
    assertthat::assert_that(file.access(path, mode = 2) == 0,
                            msg = paste0("No permission to write to ", path, ". Please specify other path.")
    )
  }
  
  # Get GO data
  message("Looking for the annotation...")
  go <- getOBO(goRelease = goRelease)
  go_annotation <- as.data.frame(go[c("id", "name", "def")])
  goVersion <- strsplit(attr(go, "version")[2], "\\/")[[1]][2]
  out_list <- list()
  go_domains <- c("cellular_component", "biological_process", "molecular_function")
  
  if (database == "ensembl" & !custom) {
    # Get geneset data
    message("Waiting for response from Ensembl...")
    ensembl <- biomaRt::useEnsembl(biomart = database)
    ensembl_datasets <- biomaRt::searchDatasets(ensembl, paste0(species, "_gene_ensembl"))
    assertthat::assert_that(!is.null(ensembl_datasets),
                            msg = "Ensembl dataset is not found. Make sure that 'species' value format is correct (example: 'hsapiens' for human)"
    )
    
    dataset <- ensembl_datasets[1, "dataset"]
    
    if (is.null(ensemblRelease)) {
      ensembl_archive <- listEnsemblArchives()
      database_version <- sub(".*? ", "", ensembl_archive$name[2])
    } else {
      database_version <- ensemblRelease
    }
    
    #check if file exists
    file_name <- paste0(
      "evogo__go", gsub(
        pattern = "-", replacement = "",
        x = goVersion
      ), "__ensembl", database_version, "_", species, ".rds" ################add__species
    )
    full_file_name <- file.path(path, file_name)
    if (file.exists(full_file_name)) {
      message("Previously saved up-to-date annotation found")
      return(invisible(readRDS(full_file_name)))
    }
    # tryCatch(listGOAnnotations(ensemblRelease = 33),
    #          +          error = function(e)
    #            +              message("You can't calculate the log of a character"))
    
    # Create evoGO object with latest annotation
    genesets <- getGeneSets(biomart = database, dataset = dataset, ensemblRelease = ensemblRelease)
    attr(genesets, "Annotation_version") <- database_version
    
  } else if (custom) {
    genesets <- customAnnotation
    attr(genesets, "Annotation_version") <- customName
    
    file_name <- paste0(
      "evogo__go", gsub(
        pattern = "-", replacement = "",
        x = goVersion
      ), "__custom", customName, ".rds"
    )
    full_file_name <- file.path(path, file_name)
    if (file.exists(full_file_name)) {
      message("Previously saved up-to-date annotation found")
      return(invisible(readRDS(full_file_name)))
    }
    
  }
  for (domain in go_domains) {
    message("Retrieving graph for '", domain, "'...")
    
    # Prepare evoGO object
    this_graph <- ontologyIndex2graph(go, rootName = domain)
    attr(this_graph, "GO_version") <- goVersion
    attr(this_graph, "GO_domain") <- domain
    
    out_list[[domain]] <- evoGO(this_graph, genesets, annotation = go_annotation, nCores = nCores)
  }
  
  
  message("New annotation retrieved")
  message("GO release: ", goVersion)
  if (custom) {
    message("Custom version: ", customName)
  } else {
    message("Ensembl version: ", ensembl_datasets[1, "description"])
  }
  
  # Save results
  if (save) {
    saveRDS(out_list, full_file_name)
    message("Results were successfully saved")
    
    # Delete previous annotation file to save space if needed (only look in 'extdata')
    if (deletePrevious & path == file.path(path.package("evoGO"), "extdata")) {
      all_annotation_files <- list.files(path = path, pattern = "^evogo__")
      old_annotation_files <- all_annotation_files[
        grepl("\\.rds$", all_annotation_files) &
          #          grepl(species, all_annotation_files) &
          all_annotation_files != file_name
      ]
      if (length(old_annotation_files) > 0) {
        file.remove(file.path(path, old_annotation_files))
        message(
          "Older annotation files were deleted: ",
          paste0(old_annotation_files, collapse = ", ", ".")
        )
      }
    }
  }
  
  return(invisible(out_list))
}


##' Load previously retrieved GO annotation as evoGO object
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param database database to be used for acquiring of genesets. Only 'ensembl' is
##'  currently supported
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest saved version is loaded
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl). If not specified,
##' latest saved version is loaded
##' @param custom name used to save a custom annotation
##' @param path directory where evoGO annotations are stored. If NULL, the default package
##' directory is used ('extdata')
##'
##' @return list of evoGO objects containing the latest annotation (one per GO domain).
##' @export
loadGOAnnotation <- function(species = NULL, database = "ensembl", goRelease = NULL,
                             ensemblRelease = NULL, customName = NULL, fileName = NULL, path = NULL) {
 # browser()
  # Check inputs
  assertthat::assert_that(is.null(path) | isValidString(path),
                          msg = "Argument 'path' should be a character vector with length of one"
  )
  
  if (!is.null(fileName)) {
    assertthat::assert_that(isValidString(fileName),
                            msg = "Argument 'fileName' should be a character vector with length of one"
    )
  } else {
    assertthat::assert_that(is.null(species) | isValidString(species),
                            msg = "Argument 'species' should be a character vector with length of one"
    )
    assertthat::assert_that(isValidString(database),
                            msg = "Argument 'database' should be a character vector with length of one"
    )
    assertthat::assert_that(database == "ensembl",
                            msg = "Only Ensembl database is currently supported"
    )
    assertthat::assert_that(!any(class(tryCatch(as.Date(goRelease), error = function(e) e)) == "error"),
                            msg = "Argument 'goRelease' is expected to be a a character vector with length of one in the format \"2023-04-01\""
    )
    assertthat::assert_that(is.null(ensemblRelease) || (is.numeric(ensemblRelease) && ensemblRelease%%1==0),
                            msg = "Argument 'ensemblRelease' is expected to be an integer"
    )
  }
  # Check path
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
                          msg = paste0("Path ", path, " is not available")
  )
  
  if (!is.null(fileName)) {
    file_name <- fileName
  } else {
    files <- listGOAnnotations(species = species, database = database, goRelease = goRelease,
                              ensemblRelease = ensemblRelease, custom = customName,  path = NULL)
    if (nrow(files) == 1) {
      file_name <- files$fileName
    } else if (nrow(files) > 1) {
      message(cat(c("Multipe annotation file matches:", "", capture.output(files),
                    "", "File with latest releases will be loaded"), sep = "\n"))
      file_name <- listGOAnnotations(species = species, database = database, goRelease = goRelease,
                                ensemblRelease = ensemblRelease, custom = customName,  path = NULL,
                                returnLatest = TRUE)$fileName
    }
 
    message(paste0(file_name, ".rds", " will be loaded"))
  }
  
  
  full_file_name <- file.path(path, paste0(file_name, ".rds"))
  #  assertthat::assert_that(file.access(full_file_name, mode = 4) == 0,
  #    msg = paste0("File ", full_file_name, " is not available")
  #  )
  readRDS(full_file_name)
}


##' Download example annotation for Homo sapiens
##'
##' This function obtains an example annotation and saves it to the package 'extdata'
##' directory (a default directory for stored annotations). The purpose of the function is to
##' provide the means for quick testing of the package functions. For performing a proper GO
##' enrichment analysis please acquire an up-to-date annotation using \code{\link{getGOAnnotation}}.
##' 
##' @param path directory where example evoGO annotation will be stored. If NULL, the default
##' package directory is used ('extdata').
##'
##' @return list of evoGO objects containing the example annotation (one per GO domain).
##' @export
exampleGOAnnotation <- function(path = NULL) {
  
  # Check the path
  assertthat::assert_that(is.null(path) | isValidString(path),
                          msg = "'path' should be a character vector with length of one"
  )
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
    if (!dir.exists(path)) {
      dir.create(path = path, recursive = TRUE, showWarnings = FALSE)
    }
  } else {
    assertthat::assert_that(dir.exists(path), msg = "Provided path does not exist")
  }
  assertthat::assert_that(file.access(path, mode = 2) == 0,
                          msg = paste0("No permission to write to ", path, ". Please specify other path.")
  )
  
  # Download the annotation
  url <- "https://evogodata.s3.eu-central-1.amazonaws.com/evogo__go20220701__ensembl_hsapiens_GRCh38.p13.rds"
  filename <- file.path(path, basename(url))
  curl::curl_download(url, filename)
  readRDS(filename)
}

