## This file includes functions for retrieving the annotations from
## external resources and to transform them into the input structures,
## which are required by the evoGO constructor.



##' Retrieve the ontologyIndex object from the obolibrary.
##'
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest saved version is loaded
##'
##' @return ontologyIndex object
##' @export
getOBO <- function(goRelease = NULL) {
  assertthat::assert_that(isValidGO(goRelease))

  url <- ifelse(is.null(goRelease), "http://purl.obolibrary.org/obo/go/go-basic.obo",
    paste0("http://release.geneontology.org/", goRelease, "/ontology/go-basic.obo")
  )
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
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl). If not specified,
##' latest saved version is loaded
##'
##' @return an object to use as a geneSets argument for evoGO constructor
##'   function.
##' @export
getGeneSets <- function(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", ensemblRelease = NULL) {
  assertthat::assert_that(isValidEnsembl(ensemblRelease))

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


##' List previously saved annotation files based on specified filter arguments
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param path path to the directory where the annotations files were saved
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31")
##' @param ensemblRelease numeric gene annotation version (e.g., 100 for Ensembl)
##' @param customName a string that will be matched as a substring against the names of
##' available custom annotation files
##' @param returnLatest logical value indicating whether to return information for the 
##' file with latest GO and Ensembl releases out of the ones matching the specified filter arguments
##'
##' @return a data.frame containing information about previously saved annotations which
##' matched specified filter arguments
##' @export
listGOAnnotations <- function(species = NULL, goRelease = NULL, database = "ensembl",
                              ensemblRelease = NULL, customName = NULL, path = NULL, returnLatest = FALSE) {
  # Check inputs
 # browser()
  assertthat::assert_that(is.null(path) | isValidString(path),
    msg = "Argument 'path' should be a character vector with length of one"
  )
  assertthat::assert_that(is.null(species) | isValidString(species),
    msg = "Argument 'species' should be a character vector with length of one"
  )
  assertthat::assert_that(isValidString(database),
    msg = "Argument 'database' should be a character vector with length of one"
  )
  assertthat::assert_that(isValidGO(goRelease))
  assertthat::assert_that(isValidEnsembl(ensemblRelease))
  assertthat::assert_that(database == "ensembl",
    msg = "Only Ensembl database is currently supported"
  )
  assertthat::assert_that(is.null(customName) | isValidString(customName),
    msg = "Argument 'customName' should be a character vector with length of one"
  )

  # Check path
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
    msg = paste0("Path ", path, " is not available")
  )

  all_files <- list.files(path)
  all_rds_files <- all_files[grepl("^evogo__", all_files) & grepl("\\.rds$", all_files)]
  assertthat::assert_that(length(all_rds_files) != 0,
    msg = "No previously saved annotations found"
  )

  # Data.frame with information about all annotation files
  file_names <- gsub(pattern = "\\.rds$", "", all_rds_files)
  file_info <- list()
  
  for (i in 1:length(file_names)) {
    file <- file_names[i]
    name_split <- strsplit(file, split = "__")[[1]]
    go <- substring(name_split[startsWith(name_split, "go")], 3)
    # address files produced using previous version of evoGO
    if (length(name_split) == 3) {
      ensembl <- NA
      if (startsWith(name_split[3], "custom")) {
        custom_name <- substring(name_split[3], 7)
      } else {
        custom_name <- name_split[3]
      }
      sp <- NA
    } else {
      ensembl <- name_split[startsWith(name_split, database)]
      ensembl <- ifelse(length(ensembl) == 0, NA, substring(ensembl, 8))
      custom_name <- NA
      sp <- name_split[startsWith(name_split, "species")]
      sp <- ifelse(length(ensembl) == 0, NA, substring(sp, 8))
    }
    info <- c(go, ensembl, custom_name, sp, file)
    file_info[[i]] <- info
  }
  
  file_info <- as.data.frame(do.call(rbind, file_info))
  file_info[, 1] <- as.Date(file_info[, 1], "%Y%m%d")
  file_info[, 2] <- as.integer(file_info[, 2])

  # Filter available file information based on user-provided arguments
  selected_info <- dplyr::filter(file_info, if (!is.null(goRelease)) grepl(goRelease, file_info[, 1]) else TRUE) %>%
    dplyr::filter(., if (!is.null(ensemblRelease)) grepl(paste0("\\b", ensemblRelease, "\\b"), .[, 2]) else TRUE) %>%
    dplyr::filter(., if (!is.null(customName)) grepl(customName, .[, 3]) else TRUE) %>%
    dplyr::filter(., if (!is.null(species)) grepl(species, .[, 4]) else TRUE)

  if (nrow(selected_info) == 0) {
    message("No matching annotation files found. Run listGOAnnotations() to view available files")
    return(NULL)
  }
  
  names(selected_info) <- c("goRelease", "ensemblRelease", "customName", "species", "fileName")
  selected_info <- selected_info[order(selected_info$goRelease, decreasing = TRUE), ]

  # Select the file with latest releases of GO (and Ensembl)
  if (returnLatest == TRUE) {
    latest_go <- max(selected_info$goRelease)
    latest_go_files <- selected_info[selected_info$goRelease %in% latest_go, ]
    if (!is.null(customName)) {
      assertthat::assert_that(nrow(latest_go_files) == 1,
        msg = paste("Multiple annotation files with the same releases:", "",
          paste(capture.output(latest_go_files), collapse = "\n"), "",
          "Please provide a more explicit 'customName' argument",
          sep = "\n"
        )
      )
      return(latest_go_files)
    } else {
      latest_ensemlb <- max(na.omit(latest_go_files$ensemblRelease))
      latest_file <- latest_go_files[latest_go_files$ensemblRelease %in% latest_ensemlb, ]
      assertthat::assert_that(nrow(latest_file) == 1,
        msg = paste("Multiple annotation files with the same releases:", "",
          paste(capture.output(latest_file), collapse = "\n"), "",
          "Please provide more filter arguments",
          sep = "\n"
        )
      )
      return(latest_file)
    }
  }

  selected_info
}



##' Retrieve and save the latest annotation for species of interest
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param database database to be used for acquiring of genesets. Only 'ensembl'
##' is currently supported.
##' @param nCores maximum number of cores to be used
##' @param save logical value indicating whether the result should be saved (default TRUE)
##' @param path path to a directory where the result will be saved. If NULL, the file
##' will be saved to the package 'extdata' directory
##' @param deletePrevious logical value indicating whether previous annotation for
##' that species should be removed when newer is saved to save disk space.
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' If not specified, latest version is retrieved
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl). If not specified,
##' latest version is retrieved
##' @param customAnnotation a named list where the elements are character vectors containing
##' annotation terms and these list elements are named after GO terms associated with the annotation
##' terms. See output of getGeneSets() as format reference
##' @param customName a string used to save the evoGO annotation object
##'
##' @return list of evoGO objects containing the annotation (one per GO domain).
##' @export
getGOAnnotation <- function(species = NULL, database = "ensembl", nCores = 1, save = TRUE,
                            path = NULL, deletePrevious = FALSE, goRelease = NULL,
                            ensemblRelease = NULL, customAnnotation = NULL, customName = NULL) {
  # Check inputs
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
    custom <- FALSE
  } else if (!is.null(customAnnotation) | !is.null(customName)) {
    assertthat::assert_that(!is.null(customAnnotation) & !is.null(customName),
      msg = "Both 'customAnnotation' and 'customName' must be provided in order to use custom annotation"
    )
    assertthat::assert_that(is.list(customAnnotation) && length(names(customAnnotation)) > 0,
      msg = paste0(
        "Argument 'customAnnotation' should be a named list where each element is named after a GO term ",
        "(eg. \"GO:0000002\") and the element is a character vector containing the names of the annotation terms"
      )
    )
    assertthat::assert_that(isValidString(customName) & nchar(customName) > 0 & !grepl("__", customName, fixed = TRUE),
      msg = "'customName' should be a non-empty character vector with length of one and should not contain \"__\""
    )
    if (!is.null(ensemblRelease)) {
      message("'ensemblRelease' specified by the user but is ignored when using custom annotation")
    }
    if (!is.null(species)) {
      message("'species' specified by the user but is ignored when using custom annotation")
    }
    custom <- TRUE
  }
  assertthat::assert_that(is.null(path) | isValidString(path),
    msg = "'path' should be a character vector with length of one"
  )

  assertthat::assert_that(isValidEnsembl(ensemblRelease))

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

    # check if file exists
    file_name <- paste0(
      "evogo__go", gsub(
        pattern = "-", replacement = "",
        x = goVersion
      ), "__ensembl", database_version, "__species", species, ".rds"
    )
    full_file_name <- file.path(path, file_name)
    if (file.exists(full_file_name)) {
      message(paste0("Previously saved ", file_name, " annotation found and will be loaded"))
      return(invisible(readRDS(full_file_name)))
    }

    # Get gene-GO mapping
    genesets <- getGeneSets(biomart = database, dataset = dataset, ensemblRelease = ensemblRelease)
    attr(genesets, "Annotation") <- paste0("Ensembl ", database_version)
  } else if (custom) {
    genesets <- customAnnotation
    customName <- sub(" ", "_", customName)
    attr(genesets, "Annotation") <- customName

    # check if file exists
    file_name <- paste0(
      "evogo__go", gsub(
        pattern = "-", replacement = "",
        x = goVersion
      ), "__custom", customName, ".rds"
    )
    full_file_name <- file.path(path, file_name)
    if (file.exists(full_file_name)) {
      message(paste0("Previously saved ", file_name, " annotation found and will be loaded"))
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
    message("Custom annotation: ", customName)
  } else {
    message("Ensembl version: ", ensembl_datasets[1, "description"])
    message("Ensembl release: ", database_version)
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


##' Load file with latest release from previously retrieved GO annotation matching the specified
##' filter arguments
##'
##' @param species species identifier (example: 'hsapiens' for human)
##' @param database database used for acquiring of genesets. Only 'ensembl' is
##'  currently supported
##' @param goRelease GO release date in yyyy-mm-dd format (e.g., "2022-01-31").
##' @param ensemblRelease gene annotation version (e.g., 100 for Ensembl).
##' @param custom a string that will be matched as a substring against the names of
##' available custom annotation files
##' @param path directory where evoGO annotations are stored. If NULL, the default package
##' directory is used ('extdata')
##'
##' @return list of evoGO objects containing the latest annotation (one per GO domain).
##' @export
loadGOAnnotation <- function(species = NULL, database = "ensembl", goRelease = NULL,
                             ensemblRelease = NULL, customName = NULL, path = NULL) {
  # Check inputs
  assertthat::assert_that(!is.null(species) || !is.null(customName),
    msg = "Either 'species' or 'customName' argument has to be provided"
  )
  assertthat::assert_that(is.null(species) != is.null(customName),
    msg = "Only one of 'species' and 'customName' arguments can be provided"
  )

  # Check path
  if (is.null(path)) {
    package_dir <- path.package("evoGO")
    path <- file.path(package_dir, "extdata")
  }
  assertthat::assert_that(file.access(path, mode = 4) == 0,
    msg = paste0("Path ", path, " is not available")
  )

  # Find files matching filter arguments
  file_name <- listGOAnnotations(
    species = species, database = database, goRelease = goRelease,
    ensemblRelease = ensemblRelease, custom = customName, path = path, returnLatest = TRUE
  )$fileName


  message(paste0(file_name, ".rds", " will be loaded"))

  full_file_name <- file.path(path, paste0(file_name, ".rds"))
  assertthat::assert_that(file.access(full_file_name, mode = 4) == 0,
    msg = paste0("File ", full_file_name, " is not available")
  )
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
