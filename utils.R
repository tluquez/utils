sink_all <- function(expr, output_file) {
  # Open the output file for writing and messages
  out_file <- file(output_file, open = "wt")

  # Close the sinks and the file connection
  on.exit({
    sink(type = "message")
    sink()
    close(out_file)
    cat("Sink closed. Output redirected to ", output_file, "\n")
  })

  sink(out_file, type = "message")
  sink(out_file, type = "output", append = TRUE)

  # Capture warnings and errors using tryCatch
  withCallingHandlers({
    eval(expr, envir = parent.frame())
  }, warning = function(w) {
    cat("\nWarning: ", w$message, "\n", file = out_file, append = TRUE)
  }, error = function(e) {
    cat("\nError: ", conditionMessage(e), "\n", file = out_file, append = TRUE)
  })
}

Mode <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- x[!is.na(x)]
  }

  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

filter_cell_types <- function(data, grouping_vars, min_n, min_pct, covs_n,
                              donor_col = "donor", cell_type_col = "cell_type",
                              num_cells_col = "num_cells") {
  #' Filter cell types based on specified criteria
  #'
  #' This function filters cell types based on the specified criteria, including
  #'  the number of donors (\code{min_n}), the percentage of non-zero
  #'  \code{num_cells}, and the number of grouping variables (\code{covs_n}).
  #'
  #' @param data A data frame containing the relevant columns.
  #' @param grouping_vars A character vector specifying the grouping variable(s).
  #' @param min_n The minimum number of donors for a cell type to be considered.
  #' @param min_pct The minimum percentage of non-zero \code{num_cells} for a
  #' cell type to be considered.
  #' @param covs_n The minimum number of grouping variables to be considered.
  #' @param donor_col The name of the column containing donor information.
  #' @param cell_type_col The name of the column containing cell type information.
  #' @param num_cells_col The name of the column containing the number of cells information.
  #'
  #' @return A character vector containing the filtered cell types.
  #'
  #' @details This function filters cell types based on the specified criteria. It checks for
  #' the minimum number of donors (\code{min_n}), the minimum percentage of non-zero
  #' \code{num_cells} (\code{min_pct}), and the minimum number of grouping variables (\code{covs_n}).
  #' The resulting character vector contains the cell types that meet the specified criteria.
  #'
  #' @examples
  #' \dontrun{
  #' filtered_types <- filter_cell_types(your_data, c("grouping_var1", "grouping_var2"), 5, 0.1, 3)
  #' }
  #'
  #' @seealso \code{\link{dplyr::group_by}}, \code{\link{dplyr::summarize}}, \code{\link{dplyr::filter}},
  #' \code{\link{tidyr::drop_na}}, \code{\link{tidyr::pull}}
  #'
  #' @importFrom dplyr n_distinct group_by summarize filter pull
  #' @importFrom tidyr drop_na
  #'
  #' @export

  # Safety checks
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }

  if (!is.character(grouping_vars) || length(grouping_vars) == 0) {
    stop("Input 'grouping_vars' must be a non-empty character vector.")
  }

  if (!is.character(donor_col) || !is.character(cell_type_col) || !is.character(num_cells_col)) {
    stop("Columns 'donor_col', 'cell_type_col', and 'num_cells_col' must be character strings.")
  }

  if (!is.numeric(min_n) || !is.numeric(min_pct) || !is.numeric(covs_n)) {
    stop("Inputs 'min_n', 'min_pct', and 'covs_n' must be numeric.")
  }

  # Remove continuous grouping variables
  if (any(purrr::map_lgl(data[,grouping_vars], is.numeric))) {
    grouping_vars <- grouping_vars[!purrr::map_lgl(data[,grouping_vars],
                                                   is.numeric)]
  }

  # Get unique cell types
  cell_types <- unique(data[[cell_type_col]])
  cell_types <- structure(cell_types, names = cell_types)

  # Filter
  cell_types_filtered <- data %>%
    tidyr::drop_na(!!sym(cell_type_col)) %>%
    {
      # Group by multiple categorical grouping variables
      if (length(grouping_vars) > 1) {
        dplyr::group_by(., !!sym(cell_type_col), !!!syms(grouping_vars))
      } else if (length(grouping_vars) == 1) {
        # Group by a single categorical grouping variable
        dplyr::group_by(., !!sym(cell_type_col), !!sym(grouping_vars))
      } else {
        # Group by cell type only
        dplyr::group_by(., !!sym(cell_type_col))
      }
    } %>%
    dplyr::summarize(
      n_donor = dplyr::n_distinct(!!sym(donor_col)),
      absent = mean(!!sym(num_cells_col) != 0) <= min_pct
    ) %>%
    dplyr::filter(n_donor >= min_n, absent == FALSE) %>%
    dplyr::group_by(!!sym(cell_type_col)) %>%
    dplyr::summarize(., n_donor = sum(n_donor)) %>%
    dplyr::filter(n_donor >= covs_n) %>%
    dplyr::pull(!!sym(cell_type_col))

  # Return filtered cell types
  cell_types[cell_types %in% cell_types_filtered]
}

pycat_to_rfactor <- function(x) {
  #' Convert Python-style categorical data to R factor
  #'
  #' This function converts Python-style categorical data with a combination
  #' of `categories` and `codes` into an R factor. \href{https://github.com/scverse/anndata/blob/a96205c58738b7054d7ed663707c2fec6b52e31d/docs/fileformat-prose.rst#id11}{Source}.
  #'
  #' @param x A list containing the categorical data. Expected to have elements
  #'   `categories` and `codes` (numeric indices referring to the categories).
  #' @return An R factor corresponding to the Python-style categorical data.
  #' @examples
  #' # Example input
  #' pycat <- list(categories = c("A", "B", "C"), codes = c(0, 1, 2, 0, 1))
  #' pycat_to_rfactor(pycat)
  #' # [1] A B C A B
  #' # Levels: A B C

  stopifnot(is.list(x))

  if (sum(names(x) %in% c("mask", "values")) == 2) {
    # Category bool to integer vector
    y <- x$values[ifelse(x$mask == T, NA, T)]
  } else if (sum(names(x$categories) %in% c("mask", "values")) == 2) {
    # Category str of int/float to numeric vector
    ## Mask == F if not NA
    x$categories <- as.character(x$categories$values[x$categories$mask == F])
    x$codes[x$codes == -1] <- NA
    y <- factor(x$categories[x$codes + 1])
  } else if (sum(names(x) %in% c("categories", "codes")) == 2) {
    # Category str to character factor
    x$codes[x$codes == -1] <- NA
    y <- factor(x$categories[x$codes + 1])
  } else {
    stop("Could not map python category to R vector")
  }

  return(y)
}

get_h5ad_meta <- function(h5ad, FactorsAsStrings = T) {
  #' Extract Metadata from h5ad File
  #'
  #' This function reads an h5ad file and extracts cell metadata, handling
  #' both SEA-AD and Scanpy formats. Metadata is returned as a `data.frame`,
  #' with cell IDs as row names.
  #'
  #' @param h5ad A character string. The file path to the h5ad file.
  #' @return A `data.frame` containing the cell metadata with cell IDs as
  #'   row names. Columns correspond to metadata attributes from the h5ad file.
  #' @details
  #' The function supports both SEA-AD and Scanpy formats:
  #' - For SEA-AD, categorical variables are converted to factors using the
  #'   `__categories` attribute.
  #' - For Scanpy, Python-style categorical data (`categories`, `codes`) are
  #'   converted to R factors.
  #'
  #' The function validates the file format, ensures consistent cell IDs
  #' and metadata columns, and replaces illegal column names.
  #'
  #' @examples
  #' # Example usage:
  #' metadata <- get_h5ad_meta("path/to/file.h5ad")
  #' head(metadata)
  #'
  #' @seealso [pycat_to_rfactor()] for converting categorical data.
  #' @export

    suppressPackageStartupMessages({
        library(tools)
        library(rhdf5)
        library(purrr)
    })

    stopifnot(file.exists(h5ad))
    stopifnot(tools::file_ext(h5ad) == "h5ad")

    # Import obs and extract cell IDs
    obs <- tryCatch({
      rhdf5::h5read(h5ad, "obs", read.attributes = T)
    }, error = function(e) {
      stop("Error reading 'obs' from h5ad file: ", conditionMessage(e))
    })
    rhdf5::h5closeAll()

    # Extract cell IDs
    if (!"_index" %in% names(attributes(obs))) {
      stop("Missing '_index' attribute in 'obs'.")
    }
    obs_i <- as.character(obs[[attributes(obs)[["_index"]]]])

    # Handle SEA-AD's format
    if ("__categories" %in% names(obs)) {
      cats <- names(obs)[names(obs) %in% names(obs[["__categories"]])]
      cats <- structure(cats, names = cats)
      cats <- purrr::map(cats, ~ factor(obs[[.x]],
                                        labels = obs[["__categories"]][[.x]]))
      num <- names(obs)[!names(obs) %in% names(obs[["__categories"]]) &
                          names(obs) != "__categories"]
      d <- c(cats, obs[num])
    } else {
      # Handle Scanpy's format
      valid_names <- c("categories", "codes", "mask", "values")
      for (i in which(purrr::map_int(obs, purrr::vec_depth) == 3)) {
          if (!sum(names(obs[[i]]) %in% valid_names) == 2) {
              .y <- paste0(names(obs[i]), "/", names(obs[[i]]))
              obs[i] <- unlist(obs[i], F, F)
              names(obs)[i] <- .y
              warning("Illegal column names were replaced")
          }
      }
      # Transform categories (length(.x) == 2) to factors
      d <- purrr::modify_if(obs, .p = ~ length(.x) == 2,
                            .f = ~ pycat_to_rfactor(.x))
    }

    # Validate and enframe metadata
    if (!purrr::every(d, ~ length(.x) == length(obs_i))) {
      stop("Mismatch between cell metadata and cell IDs.")
    }
    if (!"column-order" %in% names(attributes(obs))) {
      stop("Missing 'column-order' attribute in 'obs'.")
    }

    obs <- tryCatch({
      data.frame(d, check.names = F, stringsAsFactors = F)[
        , attributes(obs)[["column-order"]]]
    }, error = function(e) {
      stop("Error constructing data frame from obs: ", conditionMessage(e))
    })
    rownames(obs) <- obs_i

    # Convert factors to string
    if (FactorsAsStrings) {
      i <- sapply(obs, is.factor)
      obs[i] <- lapply(obs[i], as.character)
    }

    return(obs)
}

get_h5ad_expr <- function(h5ad, transpose = T, class = "H5SparseMatrix") {
  #' Extract Expression Matrix from h5ad File
  #'
  #' This function extracts the expression matrix from an h5ad file. It ensures
  #' the row and column indices (e.g., gene names and cell IDs) are correctly
  #' assigned, handling cases where the gene names may be stored as indices
  #' or categorical variables.
  #'
  #' @param h5ad A character string. The file path to the h5ad file.
  #' @param transpose logical. Whether to transpose the cell-by-gene matrix to
  #' return gene-by-cell instead.
  #' @param class character. Either "sparseMatrix" or "H5SparseMatrix". The
  #' latter is useful for large matrices that need to be stored using bit64
  #' (note, spam does not accept dimnames).
  #' @return A `Matrix` object containing the expression data, with genes as rows
  #'   and cells as columns. Row names represent gene names, and column names
  #'   represent cell IDs.
  #' @details
  #' The function ensures robust handling of different storage formats:
  #' - Gene names (`var`) can be stored as `_index`, `feature_name`, or
  #'   `feature_names`.
  #' - Handles categorical storage formats where gene names are represented as
  #'   codes and categories.
  #'
  #' It uses the `rhdf5` package to efficiently read the h5ad file and close all
  #' connections to avoid warnings or resource leaks.
  #'
  #' @examples
  #' # Extract the expression matrix from an h5ad file
  #' expr <- get_h5ad_expr("path/to/file.h5ad")
  #' dim(expr)  # Check dimensions of the expression matrix
  #'
  #' @seealso [get_h5ad_meta()] for extracting metadata from h5ad files.
  #' @export

  suppressPackageStartupMessages({
    library(tools)
    library(rhdf5)
    library(purrr)
  })

  stopifnot(file.exists(h5ad))
  stopifnot(tools::file_ext(h5ad) == "h5ad")

  # Import cell names
  obs_attr <- rhdf5::h5readAttributes(h5ad, "obs")
  obs_i <- tryCatch({
    if (!is.null(obs_attr[["_index"]])) {
      rhdf5::h5read(h5ad, paste0("obs/", obs_attr[["_index"]]))
    } else {
      stop("Critical error: Missing '_index' attribute for obs.")
    }
  }, error = function(e) {
    stop("Unable to load obs names: ", conditionMessage(e))
    NULL
  })

  # Import gene names
  var_attr <- rhdf5::h5readAttributes(h5ad, "var")
  var_i <- tryCatch({
    if (!is.null(var_attr[["_index"]])) {
      var_i_data <- rhdf5::h5read(h5ad, paste0("var/", var_attr[["_index"]]))

      # Handle if the index is categorical
      if (is.list(var_i_data) &&
          all(c("codes", "categories") %in% names(var_i_data))) {
        var_i_data$categories[var_i_data$codes + 1]
      } else {
        var_i_data
      }
    } else {
      stop("Critical error: Missing '_index' attribute for obs.")
    }
  }, error = function(e) {
    stop("Unable to load var names: ", conditionMessage(e))
    NULL
  })
  rhdf5::h5closeAll()

  # Import matrix
  if (class == "sparseMatrix") {
    suppressPackageStartupMessages(library(Matrix))
    X <- rhdf5::h5read(h5ad, "X", read.attributes = T,
                       bit64conversion = "bit64")
    rhdf5::h5closeAll()
    if (purrr::every(X, ~ class(.x) == "array")) {
      X <- purrr::modify_at(X, 1:3, as.numeric)
    }

    # Read in matrix
    if (attributes(X)[["encoding-type"]] == "csr_matrix") {
      X <- Matrix::sparseMatrix(
        j = X$indices,
        p = X$indptr,
        x = X$data,
        dims = attributes(X)$shape,
        dimnames = list(obs_i, var_i),
        index1 = F,
        repr = "R"
      )
    } else if (attributes(X)[["encoding-type"]] == "csc_matrix") {
      X <- Matrix::sparseMatrix(
        i = X$indices,
        p = X$indptr,
        x = X$data,
        dims = attributes(X)$shape,
        dimnames = list(obs_i, var_i),
        index1 = F,
        repr = "C"
      )
    } else if (attributes(X)[["encoding-type"]] == "coo_matrix") {
      X <- Matrix::sparseMatrix(
        i = X$row,
        j = X$col,
        x = X$data,
        dims = attributes(X)$shape,
        dimnames = list(obs_i, var_i),
        index1 = F,
        repr = "T"
      )
    }
    if (transpose) X <- Matrix::t(X)
  } else if (class == "H5SparseMatrix") {
    # Read delayed saprse matrix
    X <- HDF5Array::H5SparseMatrix(h5ad, "X")
    X <- provideDimnames(X, base = list(var_i, obs_i))
    if (!transpose) X <- Matrix::t(X)
  }

  return(X)
}

get_h5ad <- function(h5ad, transpose = T, class = "H5SparseMatrix",
                     FactorsAsStrings = T) {
  # Extracts sparse matrix and cell metadata from h5ad.
  # h5ad character. Path to h5ad file.
  # transpose logical. Whether to transpose the cellsxgenes matrix to return genesxcells.
  # class character. Either "sparseMatrix" or "spam". The latter is useful for large matrices that need to be stored using bit64.
  meta <- get_h5ad_meta(h5ad, FactorsAsStrings = FactorsAsStrings)
  expr <- get_h5ad_expr(h5ad, class = class)
  return(list(expr = expr, meta = meta))
}

merge.sparse <- function(...) {
  # Merge sparse matrices by adding values of unique combinations of indeces
  # ... Sparse matrix.
  # Example:
  # B <- sparseMatrix(c(1,3,5), c(3,6,3), x = c(1,4,1))
  # rownames(B) <- letters[1:nrow(B)]
  # colnames(B) <- LETTERS[1:ncol(B)]
  # merge.sparse(B, B[1:3, 1:3])
  # Reduce(merge.sparse, list(B, B[1:3, 1:3], B[1:5, 1:5], B[3:5, 2:4]))
  #
  # https://www.appsloveworld.com/r/100/71/r-binding-sparse-matrices-of-different-sizes-on-rows

  suppressPackageStartupMessages(library(Matrix))

  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (M in list(...)) {
    cnold <- colnames(M)
    rnold <- rownames(M)

    stopifnot(!is.null(cnold) & !is.null(rnold))

    cnnew <- union(cnnew, cnold)
    rnnew <- union(rnnew, rnold)

    cindnew <- match(cnold, cnnew)
    rindnew <- match(rnold, rnnew)
    ind <- Matrix::summary(M)
    i <- c(i, rindnew[ind[, 1]])
    j <- c(j, cindnew[ind[, 2]])
    x <- c(x, M@x)
  }

  return(Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(length(rnnew), length(cnnew)),
    dimnames = list(rnnew, cnnew)
  ))
}

h5ad_de <- function(fs = fs, replicate_col = "projid", replicate_col_ids = NULL,
                    cell_type_col = "state", label_col = "pmAD",
                    de_family = "pseudobulk", de_method = "limma",
                    de_type = "voom", model = NULL, groups_label = NULL,
                    verbose = F) {
  # Runs limma voom on pseudobulked count matrices from h5ad files
  # fs character. h5ad filenames. If named, the resulting dataframe will use those to refer to the cell type in question.
  # replicate_col character. See Libra.
  # replicate_col_ids character. replicate_col values to be included in the DE analysis. Used for subsetting to specific individuals.
  # label_col character. See Libra.
  # de_family character. See Libra.
  # de_method character. See Libra.
  # de_type character. See Libra.
  # model character. See Libra.
  # groups_label character. h5ad filename label shared by multiple h5ads in fs.
  # verbose character. See Libra.
  #
  # Examples
  # # Get existing files
  # fs <- list.files("data/integrated/p400", "*.h5ad", full.names = T)
  # names(fs) <- fs::path_file(fs::path_ext_remove(fs))
  #
  # # Run limma voom DE at the individual and broad cell type level
  # ict <- h5ad_de(fs = fs, model = ~ age_death + sex + pmi, verbose=T)
  # bct <- h5ad_de(fs = fs, model = ~ age_death + sex + pmi, groups_label = c("Ast", "Mic", "Oli", "Inh", "Exc"))

  suppressPackageStartupMessages({
    library(Matrix)
    library(purrr)
    library(fs)
    library(Libra)
    library(stringr)
    library(dplyr)
    source("~tl3087/utils.R")
  })

  # Define covariates
  if (is.null(model)) {
    covs <- replicate_col
  } else if (!inherits(model, "formula")) {
    stop("Model is not of class formula. Did you forget to add \"~\"?")
  } else {
    covs <- c(replicate_col, labels(terms(model)))
  }

  # Define grouping level for DE
  if (!is.null(groups_label)) {
    if (!all(groups_label %in% stringr::str_extract(fs, groups_label))) {
      stop("At least one element in groups_label is not present in fs.")
    }
  } else {
    groups_label <- fs::path_file(fs::path_ext_remove(fs))
  }

  # Return DE df per group_label
  results <- tryCatch(
    {
      purrr::map_dfr(groups_label, function(g) {
        # Extract files matching group label
        fs <- fs[stringr::str_detect(fs, stringr::fixed(g))]

        # Parse counts and cell-level metadata from h5ads
        cts <- purrr::map(fs, get_h5ad)

        ## Get rid of obs with NAs in any model or ID variable
        if (!is.null(model)) {
          cts <- purrr::map(cts, ~ {
            expr <- .$expr
            meta <- .$meta
            meta <- meta[complete.cases(meta[, covs]), ]
            expr <- expr[, rownames(meta)]
            return(list(expr = expr, meta = meta))
          })
        }

        # Filter to specified replicate_col labels
        if (!is.null(replicate_col_ids)) {
          cts <- purrr::map(cts, ~ {
            expr <- .$expr
            meta <- .$meta
            if (class(meta[[replicate_col]]) != class(replicate_col_ids)) {
              stop("Make sure replicate_col_ids has the same class as replicate_col.")
            }
            meta <- meta[meta[[replicate_col]] %in% replicate_col_ids, ]
            expr <- expr[, rownames(meta)]
            return(list(expr = expr, meta = meta))
          })
        }

        # Pseudobulk each sparse matrix
        if (verbose) {
          message("Pseudobulking...")
        }
        cts <- purrr::map(cts, ~ {
          ps <- tryCatch(
            {
              Libra::to_pseudobulk(
                .$expr,
                meta = .$meta,
                replicate_col = replicate_col,
                cell_type_col = cell_type_col,
                label_col = label_col,
                model = model
              )
            },
            error = function(e) {
              message(e)
              list()
            }
          )
        })

        # Make sure no extra level of list was created
        if (purrr::vec_depth(cts) > 4) {
          cts <- purrr::flatten(cts)
        }

        if (!is.null(groups_label)) {
          # Aggregate pseudocounts accross sub cell types
          expr <- purrr::reduce(purrr::map(cts, "expr"), ~ {
            .x <- as(as.matrix(.x), "sparseMatrix")
            .y <- as(as.matrix(.y), "sparseMatrix")
            merge.sparse(.x, .y)
          })
          meta <- purrr::reduce(purrr::map(cts, "meta"), dplyr::union)
          cts <- list(list(expr = expr, meta = meta)) %>% setNames(g)
        }

        # DE
        cts_de <- tryCatch(
          {
            Libra::run_de(
              cts,
              de_family = de_family,
              de_method = de_method,
              de_type = de_type,
              model = model,
              verbose = verbose
            )
          },
          error = function(e) {
            message(e)
            list()
          }
        )

        return(cts_de)
      })
    },
    error = function(e) {
      message(e)
      data.frame(
        cell_type = NA,
        gene = NA,
        avg_logFC = NA,
        p_val = NA,
        bonf = NA,
        fdr = NA,
        de_family = NA,
        de_method = NA,
        de_type = NA,
        stat = NA,
        lfcse = NA
      )
    }
  )

  return(results)
}

tidy_rma <- function(x, conf.level = 0.95, exponentiate = FALSE,
                     include_studies = FALSE, measure = "GEN", ...) {
  # Parse output from metafor::rma() into a dataframe
  #
  # x metafor object produced after running metafor::rma()

  # tidy summary estimates
  betas <- x$beta
  if (!is.null(nrow(betas)) && nrow(betas) > 1) {
    # get estimate type and fix spelling
    study <- rownames(betas)
    swap <- grepl("intrcpt", study)
    study[swap] <- "intercept"
    betas <- as.double(betas)
  } else {
    study <- "overall"
    betas <- betas[1]
  }

  if (x$level != 1 - conf.level) {
    level <- 1 - conf.level
    if (is.element(x$test, c("knha", "adhoc", "t"))) {
      crit <- if (all(x$ddf > 0)) qt(level / 2, df = x$ddf, lower.tail = FALSE) else NA
    } else {
      crit <- qnorm(level / 2, lower.tail = FALSE)
    }
    conf.low <- c(betas - crit * x$se)
    conf.high <- c(betas + crit * x$se)
  } else {
    conf.low <- x$ci.lb
    conf.high <- x$ci.ub
  }

  results <- tibble::tibble(
    estimate = betas,
    std_error = x$se,
    statistic = x$zval,
    p_value = x$pval,
    conf_low = conf.low,
    conf_high = conf.high,
    tau2 = x$tau2,
    tau2_se = x$se.tau2,
    k = x$k,
    p = x$p,
    m = x$m,
    h2 = x$H2,
    qe = x$QE,
    qep = x$QEp,
    i2 = x$I2,
    r2 = x$R2,
    ll = x$fit.stats$REML[1],
    dev = x$fit.stats$REML[2],
    aic = x$fit.stats$REML[3],
    bic = x$fit.stats$REML[4],
    aicc = x$fit.stats$REML[5],
    reoptimize = if (!is.null(x$reoptimize)) x$reoptimize,
    sigma2 = if (!is.null(x$sigma2)) x$sigma2
  )

  # tidy individual studies
  if (include_studies) {
    # use `metafor::escalc` to standardize estimates and confidence intervals
    estimates <- metafor::escalc(yi = x$yi.f, vi = x$vi.f, measure = measure) %>%
      summary(level = conf.level * 100) %>%
      as.data.frame(stringsAsFactors = FALSE)

    n_studies <- length(x$slab)

    estimates <- dplyr::bind_cols(
      term = as.character(x$slab),
      # dplyr::bind_cols is strict about recycling, so repeat for each study
      type = rep("study", n_studies),
      estimates[, c("yi", "sei", "zi")],
      p.value = rep(NA, n_studies),
      estimates[, c("ci.lb", "ci.ub")]
    )

    names(estimates) <- c(
      "term", "type", "estimate", "std.error", "statistic",
      "p.value", "conf.low", "conf.high"
    )
    estimates <- as_tibble(estimates)
    results <- dplyr::bind_rows(estimates, results)
  }

  if (exponentiate) {
    results <- exponentiate(results)
  }

  return(results)
}

my_tidy_rma <- function(i, lfcs, seis, measure = "ROM", method = "REML",
                        verbose = F, ...) {
  # Run and format radom effects mixed model meta-analysis with metafor
  #
  # i character. Meta-analysis unit. e.g. Gene name or cell type.
  # lfcs numeric. Vector with log fold change for each study.
  # seis numeric. Vector with the standard errors for each study.
  # description
  # It uses tidy_rma(), a modified version of \code{broom.mixed::tidy.rma()} that extracts more model information than the default.

  if (verbose == T & !is.null(i)) {
    message(i)
  }
  x <- tryCatch(
    {
      # Check for NAs
      if (any(is.na(lfcs))) {
        na.pos <- which(is.na(lfcs))
        lfcs <- lfcs[-na.pos]
        seis <- seis[-na.pos]
      } else if (any(is.na(seis))) {
        na.pos <- which(is.na(seis))
        lfcs <- lfcs[-na.pos]
        seis <- seis[-na.pos]
      }

      # Check legths match and there are at least two studies
      stopifnot(
        length(lfcs) == length(seis),
        length(lfcs) > 1
      )

      # Run meta-analysis
      x <- metafor::rma(
        yi = lfcs,
        sei = seis,
        measure = measure,
        method = method,
        verbose = verbose,
        ...
      )

      # Return tidy metafor
      return(tidy_rma(x))
    },
    error = function(e) {
      return(
        tibble::tibble(
          estimate = NA,
          std_error = NA,
          statistic = NA,
          p_value = NA,
          conf_low = NA,
          conf_high = NA,
          tau2 = NA,
          tau2_se = NA,
          h2 = NA,
          qe = NA,
          qep = NA,
          i2 = NA,
          ll = NA,
          dev = NA,
          aic = NA,
          bic = NA,
          aicc = NA
        )
      )
    }
  )
  return(x)
}

influence_to_df <- function(x) {
  # Extract the relevant information from the influence object into a data frame
  # metafor objects can be diagnosed with this function
  data.frame(
    rstudent = x$inf$rstudent,
    dffits = x$inf$dffits,
    cook.d = x$inf$cook.d,
    cov.r = x$inf$cov.r,
    tau2.del = x$inf$tau2.del,
    QE.del = x$inf$QE.del,
    hat = x$inf$hat,
    weight = x$inf$weight,
    dfbs = x$dfbs$intrcpt,
    inf = x$is.infl
  )
}

rma_with_restart <- function(data, type = "RE", yi_col = "log_fc",
                             sei_col = "lfcse", slab = NULL) {
  #' Perform Meta-Analysis with Restart Logic
  #'
  #' This function performs meta-analysis with restart logic based on specified models and tests.
  #'
  #' @param data A data frame containing the relevant variables for the meta-analysis.
  #' @param type A character string specifying the type of meta-analysis.
  #'   Options include "RE" for random effects, "REHSK" for random effects with Hedges' g,
  #'   and "FE" for fixed effects. Default is "RE".
  #' @param yi_col A character string specifying the column name in \code{data} containing
  #'   the effect sizes (yi) for the meta-analysis. Default is "log_fc".
  #' @param sei_col A character string specifying the column name in \code{data} containing
  #'   the standard errors of the effect sizes (sei) for the meta-analysis. Default is "lfcse".
  #' @slab optional vector with labels for the k studies
  #'
  #' @return A data frame containing the meta-analysis results, including effect size estimates,
  #'   standard errors, and additional information about the optimization process.
  #'
  #' @examples
  #' # Example data
  #' data <- tibble::tibble(
  #'   race_eth = c("H", "NHAA", "NHW"),
  #'   estimate = c(-0.0176, -0.162, -0.0537),
  #'   std.error = c(0.415, 0.245, 0.326),
  #'   statistic = c(-0.0424, -0.658, -0.165),
  #'   p.value = c(0.967, 0.515, 0.870)
  #' )
  #' # Example usage:
  #' result <- rma_with_restart(data, type = "RE", yi_col = "log_fc", sei_col = "lfcse", slab = "race_eth")
  #' result$reoptimize
  #'
  #' @references
  #' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package.
  #'   Journal of Statistical Software, 36(3), 1-48.
  #'
  #' @export

  # Define helper function to optimize model
  optimize_model <- function(method, test, control = NULL, slab = slab) {
    args <- list(
      yi = data[[yi_col]],
      sei = data[[sei_col]],
      method = method,
      test = test,
      control = control,
      slab = slab
    )

    # Remove NULL elements from the argument list
    args <- args[!sapply(args, is.null)]

    # Use do.call to invoke the function
    res <- do.call(metafor::rma, args)

    if (!inherits(res, "try-error")) {
      res$reoptimize <- paste(tolower(method), test, sep = "_")
    } else {
      message(paste("Error in optimization for method:", method, "and test:", test))
      res <- NULL
    }

    res
  }

  # Check if yi_col and sei_col exist in data
  if (!yi_col %in% colnames(data) | !sei_col %in% colnames(data)) {
    stop("Columns yi_col and sei_col must exist in the data.")
  }

  # Define optimization settings based on the type
  if (type == "RE") {
    methods <- c("REML", "REML", "ML", "FE")
    tests <- c("knha", "knha", "knha", "t")
    controls <- list(NULL, list(maxiter = 10000, stepadj = 0.5), NULL, NULL)
  } else if (type == "REHSK") {
    methods <- c("HSk", "HS", "HE", "SJ")
    tests <- c("knha", "knha", "knha", "knha")
    controls <- rep(list(NULL), 4)
  } else if (type == "FE") {
    methods <- c("FE", "FE", "FE")
    tests <- c("t", "t", "z")
    controls <- list(NULL, list(maxiter = 10000, stepadj = 0.5), NULL)
  }

  # Iterate over methods and tests to find the first successful optimization
  for (i in seq_along(methods)) {
    current_result <- optimize_model(methods[i], tests[i], controls[[i]],
                                     slab = slab)
    if (!is.null(current_result)) break
  }

  current_result
}

rma_with_diagnostics <- function(data, key = NULL, study_col = NULL,
                                 yi_col = NULL, sei_col = NULL,
                                 type = "RE") {
  #' Perform Meta-Analysis with Diagnostics
  #'
  #' Conducts a meta-analysis using the rma_with_restart function and provides various diagnostics, including leave-one-out results, outlier detection,
  #' and residuals.
  #'
  #' @param data A data frame containing the relevant variables for the meta-analysis.
  #' @param key A character string specifying the key variable in the data.
  #' @param study_col A character string specifying the column name in data containing
  #'   the study identifiers.
  #' @param yi_col A character string specifying the column name in data containing
  #'   the effect sizes (yi) for the meta-analysis.
  #' @param sei_col A character string specifying the column name in data containing
  #'   the standard errors of the effect sizes (sei) for the meta-analysis.
  #' @param type A character string specifying the type of meta-analysis.
  #'   Options include "RE" for random effects, "REHSK" for random effects with Hedges' g, and "FE" for fixed effects. Default is "RE".
  #'
  #' @return A data frame containing the meta-analysis results, including effect size estimates,
  #'   standard errors, and additional information about the optimization process.
  #'
  #' @examples
  #' # Example usage:
  #' # Example data
  #' data <- tibble::tibble(
  #'   region = c("AC"),
  #'   variable = c("braaksc_cat"),
  #'   term = c("Braak3"),
  #'   cell_type = c("Astro"),
  #'   yis_and_seis = list(tibble::tibble(
  #'     race_eth = c("H", "NHAA", "NHW"),
  #'     estimate = c(-0.0176, -0.162, -0.0537),
  #'     std.error = c(0.415, 0.245, 0.326),
  #'     statistic = c(-0.0424, -0.658, -0.165),
  #'     p.value = c(0.967, 0.515, 0.870),
  #'     conf.low = c(-0.872, -0.652, -0.685),
  #'     conf.high = c(0.767, 0.313, 0.598)
  #'   ))
  #' )
  #' result <- rma_with_diagnostics(data, key = "yis_and_seis",
  #'  study_col = "race_eth", yi_col = "estimate", sei_col = "std.error")
  #' print(result)
  #'
  #' @references
  #' Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package.
  #'   Journal of Statistical Software, 36(3), 1-48.
  #'
  #' @export

  if (!key %in% colnames(data)) stop("Column key must exist in the data.")
  key_col_names <- colnames(data[[key]][[1]])
  if (any(!c(study_col, yi_col, sei_col) %in% key_col_names)) {
    stop("Columns study_col, yi_col, and sei_col must exist in the key list column.")
  }

  # Run meta-analysis safely
  safe_rma <- purrr::compose(
    purrr::list_flatten,
    purrr::safely(purrr::quietly(rma_with_restart))
  )

  tictoc::tic.clear()
  tictoc::tic("Total")

  tictoc::tic("Running meta-analysis")
  res <- data %>%
    dplyr::mutate(
      # Run meta-analysis safely
      rma_result = purrr::map(
        !!rlang::sym(key),
        ~ safe_rma(.x, type, yi_col, sei_col, .x[[study_col]])
      ),
      # Extract error
      error = purrr::map_chr(
        rma_result,
        ~ ifelse(!is.null(.x$error),
                 as.character(.x$error$message), NA)
      ),
      warning = purrr::map_chr(
        rma_result,
        ~ ifelse(!is.null(.x$result_warnings),
                 as.character(.x$result_warnings), NA)
      ),
      messages = purrr::map_chr(
        rma_result,
        ~ ifelse(!is.null(.x$result_messages),
                 as.character(.x$result_messages), NA)
      ),
      # Extract results
      rma = purrr::map(rma_result, "result_result"),
      # Wider per-study effects
      purrr::map_dfr(
        !!rlang::sym(key),
        ~ .x %>%
          tidyr::pivot_wider(
            names_from = !!rlang::sym(study_col),
            values_from = -!!rlang::sym(study_col),
            names_glue = paste0(
              "original_{", study_col,
              "}_{.value}"
            )
          )
      )
    ) %>%
    dplyr::select(-rma_result)

  # Set aside conditions with errors
  res_error <- dplyr::filter(res, !is.na(error))

  # Tidy metafor results and per observation information in key column
  res %<>%
    # Remove conditions with errors
    dplyr::filter(is.na(error)) %>%
    dplyr::mutate(
      # Wider meta-analysis results
      purrr::map_dfr(rma, tidy_rma)
    )
  tictoc::toc()

  # Get diagnostics
  tictoc::tic("Running leave one out")
  res %<>%
    dplyr::mutate(
      purrr::map_dfr(
        rma,
        ~ tryCatch(
          {
            metafor::leave1out(.x) %>%
              tibble::as_tibble() %>%
              dplyr::bind_cols("{study_col}" := .x$slab) %>%
              tidyr::pivot_wider(
                names_from = !!rlang::sym(study_col),
                values_from = -!!rlang::sym(study_col),
                names_glue = paste0("l1o_{", study_col, "}_{.value}")
              ) %>%
              dplyr::rowwise() %>%
              dplyr::mutate(
                l1o_mean = mean(
                  dplyr::c_across(dplyr::ends_with("_estimate")),
                  na.rm = T
                ),
                l1o_median = median(
                  dplyr::c_across(dplyr::ends_with("_estimate")),
                  na.rm = T
                ),
                l1o_std = sd(
                  dplyr::c_across(dplyr::ends_with("_estimate")),
                  na.rm = T
                ),
                l1o_cv = l1o_std / l1o_mean
              ) %>%
              dplyr::ungroup()
          },
          error = function(e) {
            message(e)
            data.frame(error = e$message)
          }
        )
      )
    )
  tictoc::toc()

  tictoc::tic("Running outlier detection")
  res %<>%
    dplyr::mutate(
      purrr::map_dfr(
        rma,
        ~ tryCatch(
          {
            stats::influence(.x) %>%
              influence_to_df() %>%
              dplyr::bind_cols("{study_col}" := .x$slab) %>%
              tidyr::pivot_wider(
                names_from = !!rlang::sym(study_col),
                values_from = -!!rlang::sym(study_col),
                names_glue = paste0("influence_{", study_col, "}_{.value}")
              ) %>%
              dplyr::rowwise() %>%
              dplyr::mutate(
                influence_n = sum(
                  dplyr::c_across(dplyr::ends_with("_inf")),
                  na.rm = T
                )
              ) %>%
              dplyr::ungroup()
          },
          error = function(e) {
            message(e)
            data.frame()
          }
        )
      )
    )
  tictoc::toc()

  tictoc::tic("Getting residuals")
  res %<>%
    dplyr::mutate(
      purrr::map_dfr(
        rma,
        ~ tryCatch(
          {
            data.frame(residuals = stats::residuals(.x)) %>%
              dplyr::bind_cols("{study_col}" := .x$slab) %>%
              tidyr::pivot_wider(
                names_from = !!rlang::sym(study_col),
                values_from = -!!rlang::sym(study_col),
                names_glue = paste0("residuals_{", study_col, "}_{.value}")
              )
          },
          error = function(e) {
            message(e)
            data.frame()
          }
        )
      )
    )
  tictoc::toc()
  tictoc::toc()

  # Return all conditions
  dplyr::bind_rows(res, res_error)
}

get_cooksd <- function(M, design = NULL, ndups = 1, spacing = 1,
                       weights = NULL) {
  # 	Fit linear model for each gene to a series of arrays
  # 	Gordon Smyth
  # 	18 Apr 2002. Revised 9 June 2020.
  # Edited by Tain Luquez. 03/30/2023 to output residuals and cook's distance

  library(limma)
  # 	Check expression matrix
  M <- as.matrix(M)
  narrays <- ncol(M)

  # 	Check design
  if (is.null(design)) {
    design <- matrix(1, narrays, 1)
  } else {
    design <- as.matrix(design)
  }
  nbeta <- ncol(design)
  coef.names <- colnames(design)
  if (is.null(coef.names)) coef.names <- paste("x", 1:nbeta, sep = "")

  # 	Check weights
  if (!is.null(weights)) {
    weights <- limma::asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }

  # 	Reform duplicated rows into columns
  if (ndups > 1) {
    M <- limma:::unwrapdups(M, ndups = ndups, spacing = spacing)
    design <- design %x% rep_len(1, ndups)
    if (!is.null(weights)) {
      weights <- limma:::unwrapdups(weights, ndups = ndups,
        spacing = spacing)
    }
  }

  # 	Initialize standard errors
  ngenes <- nrow(M)
  stdev.unscaled <- beta <- matrix(NA, ngenes, nbeta, dimnames = list(rownames(M),
    coef.names))

  # 	Check whether QR-decomposition is constant for all genes
  # 	If so, fit all genes in one sweep
  NoProbeWts <- all(is.finite(M)) && (is.null(weights) || !is.null(attr(weights, "arrayweights")))
  if (NoProbeWts) {
    if (is.null(weights)) {
      fit <- lm.fit(design, t(M))
    } else {
      fit <- lm.wfit(design, t(M), weights[1, ])
      fit$weights <- NULL
    }
    if (fit$df.residual > 0) {
      if (is.matrix(fit$effects)) {
        fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 1):narrays, , drop = FALSE]^2))
      } else {
        fit$sigma <- sqrt(mean(fit$effects[(fit$rank + 1):narrays]^2))
      }
    } else {
      fit$sigma <- rep_len(NA_real_, ngenes)
    }
    fit$fitted.values <- fit$residuals <- fit$effects <- NULL
    fit$coefficients <- t(fit$coefficients)
    fit$cov.coefficients <- chol2inv(fit$qr$qr, size = fit$qr$rank)
    est <- fit$qr$pivot[1:fit$qr$rank]
    dimnames(fit$cov.coefficients) <- list(coef.names[est], coef.names[est])
    stdev.unscaled[, est] <- matrix(sqrt(diag(fit$cov.coefficients)), ngenes,
      fit$qr$rank, byrow = TRUE)
    fit$stdev.unscaled <- stdev.unscaled
    fit$df.residual <- rep_len(fit$df.residual, length.out = ngenes)
    dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
    fit$pivot <- fit$qr$pivot
    return(fit)
  }

  # 	Genewise QR-decompositions are required, so iterate through genes
  beta <- residuals <- stdev.unscaled
  cooksd <- M
  sigma <- rep_len(NA_real_, ngenes)
  df.residual <- rep_len(0, ngenes)
  for (i in 1:ngenes) {
    y <- as.vector(M[i, ])
    obs <- is.finite(y)
    if (sum(obs) > 0) {
      X <- design[obs, , drop = FALSE]
      y <- y[obs]
      if (is.null(weights)) {
        out <- lm.fit(X, y)
      } else {
        w <- as.vector(weights[i, obs])
        out <- lm.wfit(X, y, w)
        class(out) <- "lm"
        out$cooksd <- cooks.distance(out)
      }
      est <- !is.na(out$coefficients)
      beta[i, ] <- out$coefficients
      residuals[i, ] <- out$coefficients
      cooksd[i, ] <- out$cooksd
      stdev.unscaled[i, est] <- sqrt(diag(chol2inv(out$qr$qr, size = out$rank)))
      df.residual[i] <- out$df.residual
      if (df.residual[i] > 0) sigma[i] <- sqrt(mean(out$effects[-(1:out$rank)]^2))
    }
  }

  # 	Correlation matrix of coefficients
  QR <- qr(design)
  cov.coef <- chol2inv(QR$qr, size = QR$rank)
  est <- QR$pivot[1:QR$rank]
  dimnames(cov.coef) <- list(coef.names[est], coef.names[est])

  list(coefficients = beta, residuals = residuals, cooksd = cooksd,
    stdev.unscaled = stdev.unscaled, sigma = sigma, df.residual = df.residual,
    cov.coefficients = cov.coef, pivot = QR$pivot, rank = QR$rank)
}

get_pseudobulk <- function(expr, meta, replicate_col = "replicate",
                           cell_type_col = "cell_type", min_cells = 1,
                           min_replicates = 2) {
  # Barebones of Libra::to_pseudobulk() without the argument "label"
  # expr matrix. Counts matrix with dimensions genes x cell of mode numeric
  # meta data.frame. Must have rownames equal to colnames(expr)
  # replicate character. Unit of aggregation
  # cell_type_col character. Unit of categorization
  # min_cells numeric. Positive integer defining the minimum number of cells per cell type to be pseudobulked.
  # min_replicates numeric. Positive integer defining the minimum number of units of aggregation to pseduobulk on.
  #
  # Example
  # pseudobulk <- get_pseudobulk(expr, meta, replicate_col = "donor",
  # cell_type_col = "cell_type")
  # library(DelayedArray)
  # library(dplyr)
  # library(tibble)
  # library(purrr)

  # Keep only cell types with enough cells
  cell_types <- unique(meta[[cell_type_col]])
  cell_types_filtered <- meta %>%
    tibble::rownames_to_column("cellid") %>%
    dplyr::group_by(!!sym(cell_type_col)) %>%
    dplyr::summarize(n_cell = dplyr::n_distinct(cellid),
                     n_replicate = dplyr::n_distinct(!!sym(replicate_col))) %>%
    dplyr::filter(n_cell > min_cells, n_replicate > min_replicates) %>%
    dplyr::pull(!!sym(cell_type_col)) %>%
    unique() %>%
    as.character()
  if (length(cell_types) > length(cell_types_filtered)) {
    message("Cell types dropped by min_cells or/and min_replicates: ",
            paste(setdiff(cell_types, cell_types_filtered), collapse = ", "))
  }
  cell_types <- cell_types[cell_types %in% cell_types_filtered]
  cell_types <- structure(cell_types, names = cell_types)
  stopifnot(length(cell_types) != 0)

  pseudobulk <- purrr::map_dfr(cell_types, function(cell_type) {
    tryCatch(
      {
        message(cell_type)
        meta <- meta[meta[[cell_type_col]] == cell_type, ]

        # Process data into gene x replicate x cell_type matrices
        model <- formula(paste0("~ 0 +", replicate_col))
        mm <- model.matrix(model, data = meta)
        colnames(mm) <- gsub(replicate_col, "", colnames(mm))
        expr <- expr[, rownames(mm)]
        stopifnot(ncol(expr) == nrow(mm))
        mat_mm <- expr %*% mm

        # filter out cell types with no retained genes
        if (nrow(mat_mm) == 0) {
          return(NA)
        }

        # To long format
        mat_mm <- t(mat_mm)
        if (class(mat_mm) == "DelayedMatrix") {
          mat_mm <- as.matrix(as(mat_mm, "dgCMatrix"))
        } else {
          mat_mm <- as.matrix(mat_mm)
        }
        mat_mm <- cbind(num_cells = colSums(mm), mat_mm) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(replicate_col)

        rm(expr, meta, mm)

        return(mat_mm)
      },
      error = function(e) {
        message(e)
        return(data.frame())
      }
    )
  }, .id = cell_type_col)

  return(pseudobulk)
}

sort_files_by_date <- function(folder_path = '.', search_pattern = NULL,
                               by = 'mtime'){
  # https://stackoverflow.com/a/60860704
  # by character. Must be one of "mtime", "atime", or "ctime".'. See file.info().

  library(purrr, include.only = "map_dfr")
  library(dplyr)

  file_names <- data.frame(file_names = list.files(path = folder_path,
                                                   pattern = search_pattern,
                                                   full.names = TRUE))
  files_by_datetime <- file_names %>%
    purrr::map_dfr(file.info) %>%
    dplyr::select(mtime, atime, ctime) %>%
    dplyr::bind_cols(file_names) %>%
    dplyr::arrange(dplyr::case_when(by == 'mtime' ~ mtime,
                                    by == 'atime' ~ atime,
                                    by == 'ctime' ~ ctime)) %>%
    rownames()

  return(files_by_datetime)
}

my_vroom <- function(...) vroom::vroom(..., show_col_types = F, progress = F)

my_vroom_write <- function(...) vroom::vroom_write(..., progress = F)

my_fread <- function(file_path, select_colnames = NULL, filter_expr = NULL,
                     sep = "\t", ...) {
  #' Custom Fast File Reader with Column Selection and Filtering
  #'
  #' Reads a delimited file with options for selecting specific columns and filtering rows based on an expression.
  #'
  #' @param file_path A character string representing the path to the file to be read. Supports gz compressed files.
  #' @param select_colnames A character vector of column names to be selected. If NULL, all columns are selected. Defaults to NULL.
  #' @param filter_expr A character string representing a logical expression to
  #' filter rows. If NULL, no filtering is applied. Not all filtering
  #' expressions are supported. Only >, < and == and combinations of these with
  #' and or or operators are supported. If filtering using character values,
  #' make sure to quote it (e.g. `"column == 'value'"`). Defaults to NULL.
  #' @param sep A character string representing the column delimiter. Defaults to tab ("\t").
  #' @param ... Additional arguments to be passed to `data.table::fread`.
  #'
  #' @return A data.table containing the selected columns and filtered rows from the input file.
  #'
  #' @details
  #' This function reads a delimited file, optionally selecting specific columns and filtering rows based on a logical expression.
  #' It supports gz compressed files and uses `awk` to efficiently process large files. The function first reads the column names
  #' from the file header, validates the `select_colnames`, and constructs the appropriate `awk` commands for column selection and
  #' row filtering. The final data is read into R using `data.table::fread`.
  #'
  #' @examples
  #' \dontrun{
  #' # Read a file with specific columns and a filter expression
  #' data <- my_fread("data.txt", select_colnames = c("col1", "col3"),
  #'                  filter_expr = "col2 > 10 & col4 == 'value'")
  #'
  #' # Read a gz compressed file
  #' data <- my_fread("data.txt.gz")
  #' }
  #'
  #' @export

  column_from_expr <- function(chr_expr) {
    # Split expression at logical operator other than & and |
    split_chr_expr <- strsplit(
      chr_expr, "\\s*==\\s*|\\s*!=\\s*|\\s*>\\s*|\\s*<\\s*"
    )

    # Trim trailing whitespaces
    split_chr_expr <- trimws(unlist(split_chr_expr))

    # Test if column is present in longer set of columns
    col_in_all_colnames <- split_chr_expr[[1]] %in% all_colnames
    if (col_in_all_colnames) {
      return(split_chr_expr[[1]])
    } else {
      stop("Column ", paste0(split_chr_expr[[1]], collapse = ", "),
           " not found.")
    }
  }

  # Check if file exists and normalize path
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  file_path <- normalizePath(file_path)

  # Allow for gz compressed files
  cmd_base <- if (grepl(".gz$", file_path)) {
    paste0("zcat ", file_path)
  } else {
    paste0("cat ", file_path)
  }

  # Read in all colnames
  header <- system(paste0(cmd_base, " | head -n 1"), intern = T)
  all_colnames <- strsplit(header, sep)[[1]]

  # Check if select_colnames are valid
  if (!is.null(select_colnames) && !all(select_colnames %in% all_colnames)) {
    select_colnames <- paste(select_colnames[!select_colnames %in%
                                               all_colnames],
                             collapse = ", ")
    stop(paste0("Some provided column names do not match the column names ",
                "in ", file_path, ": ", select_colnames))
  }

  # Use all columns if select_colnames not provided
  if (is.null(select_colnames)) {
    select_colnames <- all_colnames
  }

  # Build awk filter condition if filter_expr provided
  awk_filter <- ""
  if (!is.null(filter_expr)) {
    # Split filter_expr by '&' or '|'
    split_expr <- trimws(unlist(strsplit(filter_expr, "\\s*&\\s*|\\s*\\|\\s*")))

    # Get column names from each expression
    filter_cols <- unique(unlist(lapply(split_expr, column_from_expr)))

    if (!is.null(filter_cols)) {
      # Ensure filter columns are included in the selection
      select_colnames <- unique(union(select_colnames, filter_cols))
    }

    # Replace column names in filter_expr with corresponding awk indices
    for (colname in filter_cols) {
      col_index <- which(all_colnames == colname)
      filter_expr <- gsub(colname, paste0("$", col_index), filter_expr,
                          fixed = T)
    }

    # Conform AND and OR operator from R to awk
    filter_expr <- gsub("&", "&&", filter_expr)
    filter_expr <- gsub("\\|", "\\|\\|", filter_expr)

    # Conform quotes
    filter_expr <- gsub("'", "\"", filter_expr)

    # Construct awk filter condition
    awk_filter <- filter_expr
  }

  # Build the awk print statement
  selected_indices <- which(all_colnames %in% select_colnames)
  awk_col_select <- paste0("{print ", paste(paste0("$", selected_indices),
                                            collapse = ", "), "}")

  # Construct full awk command
  cmd <- sprintf("%s | awk -F '%s' '%s %s'",
                 cmd_base, sep, awk_filter, awk_col_select)

  # Read and process the file
  data.table::fread(cmd = cmd, col.names = all_colnames[selected_indices], ...)
}

q <- function() quit(save = "no")

voomByGroup <- function(counts, group = NULL, design = NULL, lib.size = NULL,
                        dynamic = NULL, normalize.method = "none",
                        span = 0.5, save.plot = FALSE, print = TRUE,
                        plot = c("none", "all", "separate", "combine"),
                        col.lines = NULL,
                        pos.legend = c("inside", "outside", "none"),
                        fix.y.axis = FALSE, ...) {
  # 14 June 2017 (Last updated 6 May 2022)
  # Charity Law, Xueyi Dong and Yue You
  # Copied from the github repo https://github.com/YOU-k/voomByGroup/blob/main/voomByGroup.R on May 22 2023.
  # Counts
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(group)) {
      group <- counts$samples$group
    }
    # if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 0)
    #   design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) {
      lib.size <- with(counts$samples, lib.size * norm.factors)
    }
    counts <- counts$counts
  } else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,
                                                         "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) {
        out$genes <- Biobase::fData(counts)
      }
      if (length(Biobase::pData(counts))) {
        out$targets <- Biobase::pData(counts)
      }
      counts <- Biobase::exprs(counts)
    } else {
      counts <- as.matrix(counts)
    }
  }
  if (nrow(counts) < 2L) {
    stop("Need at least two genes to fit a mean-variance trend")
  }
  # Library size
  if (is.null(lib.size)) {
    lib.size <- colSums(counts)
  }
  # Group
  if (is.null(group)) {
    group <- rep("Group1", ncol(counts))
  }
  group <- as.factor(group)
  intgroup <- as.integer(group)
  levgroup <- levels(group)
  ngroups <- length(levgroup)
  # Design matrix
  if (is.null(design)) {
    design <- matrix(1L, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  # Dynamic
  if (is.null(dynamic)) {
    dynamic <- rep(FALSE, ngroups)
  }
  # voom by group
  if (print) {
    cat("Group:\n")
  }
  E <- w <- counts
  xy <- line <- as.list(rep(NA, ngroups))
  names(xy) <- names(line) <- levgroup
  for (lev in 1L:ngroups) {
    if (print) {
      cat(lev, levgroup[lev], "\n")
    }
    i <- intgroup == lev
    countsi <- counts[, i]
    libsizei <- lib.size[i]
    designi <- design[i, , drop = FALSE]
    QR <- qr(designi)
    if (QR$rank < ncol(designi)) {
      designi <- designi[, QR$pivot[1L:QR$rank], drop = FALSE]
    }
    if (ncol(designi) == ncol(countsi)) {
      designi <- matrix(1L, ncol(countsi), 1)
    }
    voomi <- limma::voom(
      counts = countsi, design = designi, lib.size = libsizei,
      normalize.method = normalize.method, span = span, plot = FALSE,
      save.plot = TRUE, ...
    )
    E[, i] <- voomi$E
    w[, i] <- voomi$weights
    xy[[lev]] <- voomi$voom.xy
    line[[lev]] <- voomi$voom.line
  }
  # voom overall
  if (TRUE %in% dynamic) {
    voom_all <- limma::voom(
      counts = counts, design = design, lib.size = lib.size,
      normalize.method = normalize.method, span = span, plot = FALSE,
      save.plot = TRUE, ...
    )
    E_all <- voom_all$E
    w_all <- voom_all$weights
    xy_all <- voom_all$voom.xy
    line_all <- voom_all$voom.line

    dge <- edgeR::DGEList(counts)
    disp <- edgeR::estimateCommonDisp(dge)
    disp_all <- disp$common
  }
  # Plot, can be "both", "none", "separate", or "combine"
  plot <- plot[1]
  if (plot != "none") {
    disp.group <- c()
    for (lev in levgroup) {
      dge.sub <- edgeR::DGEList(counts[, group == lev])
      disp <- edgeR::estimateCommonDisp(dge.sub)
      disp.group[lev] <- disp$common
    }
    if (plot %in% c("all", "separate")) {
      if (fix.y.axis == TRUE) {
        yrange <- sapply(levgroup, function(lev) {
          c(min(xy[[lev]]$y), max(xy[[lev]]$y))
        }, simplify = TRUE)
        yrange <- c(min(yrange[1, ]) - 0.1, max(yrange[2, ]) + 0.1)
      }
      for (lev in 1L:ngroups) {
        if (fix.y.axis == TRUE) {
          plot(xy[[lev]], xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25, ylim = yrange)
        } else {
          plot(xy[[lev]], xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25)
        }
        title(paste("voom: Mean-variance trend,", levgroup[lev]))
        lines(line[[lev]], col = "red")
        legend("topleft", bty = "n", paste("BCV:", round(sqrt(disp.group[lev]),
                                                         3)), text.col = "red")
      }
    }


    if (plot %in% c("all", "combine")) {
      if (is.null(col.lines)) {
        col.lines <- 1L:ngroups
      }
      if (length(col.lines) < ngroups) {
        col.lines <- rep(col.lines, ngroups)
      }
      xrange <- unlist(lapply(line, `[[`, "x"))
      xrange <- c(min(xrange) - 0.3, max(xrange) + 0.3)
      yrange <- unlist(lapply(line, `[[`, "y"))
      yrange <- c(min(yrange) - 0.1, max(yrange) + 0.3)
      plot(1L, 1L, type = "n", ylim = yrange, xlim = xrange,
           xlab = "log2( count size + 0.5 )",
           ylab = "Sqrt( standard deviation )")
      title("voom: Mean-variance trend")
      if (TRUE %in% dynamic) {
        for (dy in which(dynamic)) {
          line[[dy]] <- line_all
          disp.group[dy] <- disp_all
          levgroup[dy] <- paste0(levgroup[dy], " (all)")
        }
      }
      for (lev in 1L:ngroups) {
        lines(line[[lev]], col = col.lines[lev], lwd = 2)
      }
      pos.legend <- pos.legend[1]
      disp.order <- order(disp.group, decreasing = TRUE)
      text.legend <- paste(levgroup, ", BCV: ", round(sqrt(disp.group), 3),
                           sep = "")
      if (pos.legend %in% c("inside", "outside")) {
        if (pos.legend == "outside") {
          plot(1, 1, type = "n", yaxt = "n", xaxt = "n", ylab = "", xlab = "",
               frame.plot = FALSE)
          legend("topleft", text.col = col.lines[disp.order],
                 text.legend[disp.order], bty = "n")
        } else {
          legend("topright", text.col = col.lines[disp.order],
                 text.legend[disp.order], bty = "n")
        }
      }
    }
  }
  # Output
  if (TRUE %in% dynamic) {
    E[, intgroup %in% which(dynamic)] <- E_all[, intgroup %in% which(dynamic)]
    w[, intgroup %in% which(dynamic)] <- w_all[, intgroup %in% which(dynamic)]
  }
  out$E <- E
  out$weights <- w
  out$design <- design
  if (save.plot) {
    out$voom.line <- line
    out$voom.xy <- xy
  }
  new("EList", out)
}

cbind.fill <- function(...) {
  # https://stackoverflow.com/a/7962286
  data_list <- list(...)
  max_length <- max(sapply(data_list, length))
  data_padded <- lapply(data_list, function(x) {
    c(x, rep(NA, max_length - length(x)))
    })
  as.data.frame(data_padded)
}

getR2 <- function(p, mod, mod0 = NULL) {
  # Get R2 and adjusted R2 from limma
  #
  # Modified on 10-25-2023 by Tain Luquez.
  # Source: github.com/LieberInstitute/jaffelab/blob/devel/R/getR2.R
  fit1 <- limma::lmFit(p, mod)
  rss1 <- rowSums((p - fitted(fit1))^2)
  n <- ncol(p)
  k <- ncol(mod) - 1

  if (is.null(mod0)) {
    rss0 <- rowSums((p - rowMeans(p))^2)
  } else {
    fit0 <- limma::lmFit(p, mod0)
    rss0 <- rowSums((p - fitted(fit0))^2)
  }

  r2 <- 1 - (rss1 / rss0)
  r2adj <- 1 - ((1 - r2) * (n - 1)) / (n - k - 1)
  out <- data.frame(R2 = r2, Adjusted_R2 = r2adj)
  return(out)
}

check_nonestimable <- function(design, verbose = F) {
  non_estimable <- limma:::nonEstimable(design)
  if (!is.null(non_estimable)) {
    design <- design[, !colnames(design) %in% non_estimable]
    if (verbose) message("Dropping colinear terms from design matrix: ",
                         paste(non_estimable, collapse = ", "))
  }
  design
}

de_pseudobulk <- function(input, meta, model, contrast, de_method = "limma",
                          de_type = "voombygroup", verbose = F,
                          return_interm = F, contrast_matrix = NULL) {
  #' Perform Differential Expression on Pseudobulk Data
  #'
  #' @description
  #' This function performs differential expression analysis on pseudobulk RNA-seq
  #' data. The input can either be a list of cell types, where each entry is a list
  #' containing an expression matrix (`expr`) and corresponding metadata (`meta`),
  #' or a single count matrix with accompanying metadata. It supports differential
  #' expression methods such as `limma`, `edgeR`, and `DESeq2`, providing flexibility
  #' in analysis approaches and contrast settings.
  #'
  #' @param input A list or matrix. If a list, it should contain pseudobulk data for
  #' different cell types generated by `Libra::to_pseudobulk()` with each entry
  #' containing an expression matrix (`expr`) and metadata (`meta`). Row names of
  #' `meta` must match column names of `expr`. Alternatively, `input` can be a
  #' gene-by-observation count matrix, in which case `meta` must be supplied and
  #' should have row names matching column names of `input`. Counts should be
  #' untransformed.
  #' @param meta A `data.frame` containing per-observation covariates used in the
  #' model. This is required if `input` is a matrix and optional if `input` is a list
  #' (in which case each entry in the list should contain its own metadata). Column
  #' names of `meta` should correspond to terms used in `model`.
  #' @param model A formula specifying the model for differential expression, created
  #' using `as.formula()` or `~`. Terms in the formula should match column names in
  #' `meta` or, if `input` is a list, in each entry’s metadata.
  #' @param contrast A numeric or character specifying the position or name of the
  #' contrast within the `design.matrix()` generated from `model` and `meta`.
  #' @param de_method A character string indicating the differential expression method
  #' to use. Supported options are `"limma"`, `"edgeR"`, and `"DESeq2"`. Default is
  #' `"limma"`.
  #' @param de_type A character string specifying the type of analysis within the
  #' selected method. Options vary based on `de_method`, with `"voombygroup"` as the
  #' default for `limma`.
  #' @param verbose A logical indicating the level of information to display. If
  #' `TRUE`, additional information about progress is printed. Default is `FALSE`.
  #' @param return_interm A logical. If `TRUE`, the function will return intermediate
  #' files (e.g., model matrices, fitted models), which can be useful for debugging or
  #' further examination. Default is `FALSE`.
  #' @param contrast_matrix A matrix specifying multiple contrasts, where each column
  #' represents a contrast. The number of rows must match the number of columns in the
  #' design matrix. This option is only supported when `de_method` is `"limma"`.
  #'
  #' @details
  #' The function is inspired by the differential expression functions in `Libra` and
  #' performs various checks on inputs to ensure consistency. It verifies that
  #' observation names match between expression data and metadata and ensures all model
  #' terms are present in the metadata. Integer count values are required, and only
  #' `limma` supports multiple contrasts via a contrast matrix.
  #'
  #' The `de_method` and `de_type` parameters allow for flexible analysis, with the
  #' option for multiple contrasts in `limma` that could return an F-statistic in the
  #' top table (future implementation).
  #'
  #' @return A `data.frame` of differential expression results for each cell type, or,
  #' if `return_interm` is `TRUE`, a list containing intermediate files for further
  #' examination.
  #'
  #' @examples
  #' \dontrun{
  #' # Example usage with a list of pseudobulk data per cell type
  #' pseudobulk_data <- Libra::to_pseudobulk(sce_data)
  #' de_results <- de_pseudobulk(input = pseudobulk_data, meta = NULL,
  #'                             model = ~ condition, contrast = "condition1")
  #'
  #' # Example with a matrix and metadata
  #' count_matrix <- matrix(sample(0:100, 1000, replace = TRUE), nrow = 100)
  #' rownames(count_matrix) <- paste0("Gene", 1:100)
  #' colnames(count_matrix) <- paste0("Sample", 1:10)
  #' meta <- data.frame(condition = rep(c("control", "treatment"), each = 5),
  #'                    row.names = colnames(count_matrix))
  #' de_results <- de_pseudobulk(input = count_matrix, meta = meta,
  #'                             model = ~ condition, contrast = "conditiontreatment")
  #' }
  #'
  #' @seealso \code{\link[Libra]{to_pseudobulk}}, \code{\link[limma]{lmFit}},
  #' \code{\link[DESeq2]{DESeq}}, \code{\link[edgeR]{DGEList}}
  #'
  #' @export


  # Helper functions
  check_inputs <- function(expr, meta, model, contrast, de_method,
                           contrast_matrix) {
    # Check observation names match between expression and metadata
    if (!identical(sort(colnames(expr)), sort(rownames(meta)))) {
      stop("Observation names don't match between expression and metadata")
    }

    # Check if model terms are in metadata
    covs <- all.vars(model)
    if (!all(covs %in% colnames(meta))) {
      stop(paste0("Model terms absent from meta: ",
                  colnames(meta)[!covs %in% colnames(meta)]))
    }

    # Check if counts are integer
    if (!all(expr %% 1 == 0)) {
      stop("Counts are not integer.")
    }

    # Ensure design matrix is full rank
    design_full <- model.matrix(model, data = meta)
    stopifnot(identical(rownames(design_full), colnames(expr)))
    design_full <- check_nonestimable(design_full)

    # Check for contrast matrix
    if (!is.null(contrast_matrix)) {
      if (de_method != "limma") {
        stop("Only limma supports contrast_matrix.")
      }
      stopifnot(nrow(contrast_matrix) == ncol(design_full))
    }

    # Contrast checks
    if (!is.null(contrast)) {
      if (!contrast %in% colnames(design_full) && !contrast %in% covs) {
        stop(paste0("Contrast \"", contrast, "\" was not present in the model ",
                    "nor the full design with terms: ",
                    paste(colnames(design_full), collapse = ", "), "\n"))
      } else if (contrast %in% colnames(design_full)) {
        is.categorical <- length(unique(design_full[, contrast])) == 2
        group_vec <- design_full[, contrast]
      } else if (contrast %in% covs) {
        is.categorical <- any(class(meta[, contrast]) %in%
                                c("factor", "character"))
        group_vec <- meta[, contrast]

        if (is.null(contrast_matrix)) {
          contrast_in_design <- grep(contrast, colnames(design_full), value = T)
          if (length(contrast_in_design) != 1) {
            stop(length(contrast_in_design), " design matrix columns (",
                 paste(contrast_in_design, collapse = ", "), ") matched the ",
                 "contrast ", contrast, ". Select one.")
          }
        }
      }
    } else {
      stop("Contrast cannot be NULL.")
    }

    # Checks of argument congruency
    if (de_method == "limma" && !is.null(contrast_matrix)) {
      if (nrow(contrast_matrix) != ncol(design_full)) {
        stop("Dimensions between the contrast and design matrices do not match.")
      }
    }
    if (!is.categorical && de_type == "voombygroup") {
      stop("voombygroup does not support continuous contrasts. Use voom.")
    }

    return(list(is.categorical = is.categorical, group_vec = group_vec,
                design_full = design_full))
  }

  run_edgeR <- function(expr, meta, design_full, contrast, de_type, group_vec) {
    tryCatch({
      y <- edgeR::DGEList(expr = expr, group = group_vec) %>%
        edgeR::calcNormFactors(method = "TMM") %>%
        edgeR::estimateDisp(design_full)

      test <- switch(
        de_type,
        QLF = {
          fit <- edgeR::glmQLFit(y, design = design_full)
          edgeR::glmQLFTest(fit, coef = contrast)
        },
        LRT = {
          fit <- edgeR::glmFit(y, design = design_full)
          edgeR::glmLRT(fit, coef = contrast)
        }
      )

      # Format output
      res <- edgeR::topTags(test, n = Inf) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(bonf = p.adjust(PValue, "bonferroni"),
                      contrast = contrast) %>%
        dplyr::rename(dplyr::any_of(c(
          fdr = "FDR", lfc = "logFC", ave_expr = "logCPM",
          p.value = "PValue", stat = ifelse(de_type == "LRT", "LR", "F")
        ))) %>%
        dplyr::arrange(p.value)

      # Keep intermeadite files
      list(res = res, interm = list(y = y, test = test, fit = fit))
    }, error = function(e) {
      message(e)
      list()
    })
  }

  run_DESeq2 <- function(expr, meta, design_full, contrast, de_type) {
    tryCatch({
      dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
        countData = expr, colData = meta, design = design_full
      ))
      dds <- switch(
        de_type,
        Wald = {
          try(DESeq2::DESeq(dds, test = "Wald", quiet = !verbose))
        },
        LRT = {
          try(DESeq2::DESeq(dds, test = "LRT", quiet = !verbose,
                            reduced = design_reduced))
        }
      )

      # Format output
      res <- DESeq2::results(dds,
                             name = contrast, cooksCutoff = F,
                             format = "DataFrame"
      ) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        dplyr::mutate(
          bonf = p.adjust(pvalue, "bonferroni"),
          contrast = contrast
        ) %>%
        dplyr::rename(dplyr::any_of(c(
          fdr = "padj", lfc = "log2FoldChange", lfcse = "lfcSE",
          ave_expr = "baseMean", p.value = "pvalue"
        ))) %>%
        dplyr::arrange(p.value) %>%
        tidyr::drop_na()

      # Keep intermeadite files
      list(res = res, interm = dds)
    }, error = function(e) {
      message(e)
      list()
    })
  }

  run_limma <- function(expr, meta, design_full, contrast, de_type,
                        contrast_matrix, is.categorical, group_vec) {
    tryCatch({
      # Compute normalized library sizes
      norm.factors <- edgeR::calcNormFactors(expr)
      lib.size <- norm.factors * colSums(expr)

      # voom, trend or voombygroup
      v <- switch(
        de_type,
        trend = {
          trend_bool <- T
          dge <- edgeR::DGEList(expr, lib.size = lib.size,
                                norm.factors = norm.factors, group = group_vec)
          v <- methods::new("EList")
          v$E <- edgeR::cpm(dge, log = T)
          v
        },
        voom = {
          trend_bool <- F
          v <- limma::voom(expr, design_full, save.plot = T,
                           lib.size = lib.size)
          v$target <- lib.size
          v$norm.factors <- norm.factors
          v
        },
        voombygroup = {
          trend_bool <- F
          v <- voomByGroup(expr, design_full, save.plot = T, plot = "combine",
                           print = verbose, group = group_vec,
                           lib.size = lib.size)
          v$target <- lib.size
          v$norm.factors <- norm.factors
          v
        }
      )

      # Fit and moderated statistics
      fit <- limma::lmFit(v, design_full) %>%
        limma::eBayes(trend = trend_bool, robust = trend_bool)

      if (!is.null(contrast_matrix)) {
        fit2 <- limma::contrasts.fit(fit, contrast_matrix) %>%
          limma::eBayes(trend = trend_bool, robust = trend_bool)
      }

      # Format output
      res <- limma::topTable(
        fit = if (is.null(contrast_matrix)) fit else fit2,
        coef = if (is.null(contrast_matrix)) contrast else NULL,
        number = Inf, confint = T
      )

      # Add standard errors for the coefficients
      if (!is.null(contrast_matrix)) {
        SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
        if (nrow(res) == nrow(SE)) {
          colnames(SE) <- paste0("lfcse.Coef", seq_len(ncol(SE)))
          res %<>% dplyr::bind_cols(SE[rownames(res), , drop = F])
        } else {
          warning("Standard errors for the LFC were not computed due to ",
                  "mismatch gene numbers between toptable and fit.")
        }
      } else {
        SE <- sqrt(fit$s2.post) * fit$stdev.unscaled
        if (nrow(res) == nrow(SE)) {
          res %<>% dplyr::bind_cols(lfcse = unname(SE[rownames(res), contrast]))
        } else {
          warning("Standard errors for the LFC were not computed due to ",
                  "mismatch gene numbers between toptable and fit.")
        }
      }

      # Add voom plot
      if (purrr::pluck_depth(v$voom.xy) == 3) {
        voom_plot <- purrr::map2(
          purrr::imap(v$voom.xy, ~ {
            xy <- cbind(as.data.frame(.x$x), as.data.frame(.x$y))
            colnames(xy) <- c(paste0("voom_plot_", .y, "_x"),
                              paste0("voom_plot_", .y, "_y"))
            xy
          }),
          purrr::imap(v$voom.line, ~ {
            xy_line <- cbind(as.data.frame(.x$x), as.data.frame(.x$y))
            colnames(xy_line) <- c(paste0("voom_line_", .y, "_x"),
                                   paste0("voom_line_", .y, "_y"))
            xy_line
          }),
          ~tibble::rownames_to_column(cbind(.x, .y), "gene")
        ) %>%
          purrr::reduce(dplyr::full_join, by = "gene")
      } else if (purrr::pluck_depth(v$voom.xy) == 2) {
        voom_plot <- dplyr::tibble(voom_plot_x = unname(v$voom.xy[["x"]]),
                                   voom_plot_y = unname(v$voom.xy[["y"]]),
                                   voom_line_x = unname(v$voom.line[["x"]]),
                                   voom_line_y = unname(v$voom.line[["y"]])) %>%
          dplyr::mutate(gene = names(v$voom.xy[["x"]]))
      }
      res %<>%
        tibble::rownames_to_column("gene") %>%
        dplyr::left_join(voom_plot, by = "gene")

      # Add number of observations
      if (is.categorical) {
        res %<>%
          dplyr::bind_cols(
            table(group_vec) %>%
              as.data.frame() %>%
              pivot_wider(names_from = 1, values_from = 2,
                          names_prefix = "n")
          )
      } else {
        res %<>%
          dplyr::mutate(n = sum(!is.na(group_vec)))
      }

      # Adjust p-values
      res %<>%
        dplyr::mutate(bonf = p.adjust(P.Value, "bonferroni"),
                      contrast = contrast) %>%
        dplyr::rename(dplyr::any_of(c(
          fdr = "adj.P.Val", lfc = "logFC", conf.low = "CI.L", b = "B",
          conf.high = "CI.R", ave_expr = "AveExpr", p.value = "P.Value",
          stat = "t", f = "F"
        ))) %>%
        dplyr::arrange(p.value)

      # Keep intermeadite files
      list(
        res = res,
        interm = list(
          v = v,
          fit = if (is.null(contrast_matrix)) fit else list(fit = fit,
                                                            fit2 = fit2)
        )
      )
    }, error = function(e) {
      message(e)
      list()
    })
  }

  if (is.list(input)) {
    cell_types_expr <- purrr::map(input, "expr")
    cell_types_meta <- purrr::map(input, "meta")
  } else if (is.matrix(input) & is.data.frame(meta)) {
    cell_types_expr <- setNames(list(input), "")
    cell_types_meta <- setNames(list(meta), "")
  } else {
    stop("\"pseudobulks\" must be a list of counts per cell type OR a pair of count matrix and metadata data frame.")
  }
  if (verbose) {
    message("Running differential expression...")
  }

  des <- purrr::pmap(
    list(cell_types_expr, cell_types_meta, names(cell_types_expr)),
    function(expr, meta, cell_type) {
      if (verbose) {
        message(cell_type)
      }
      # Check inputs and define set up variables
      l <- check_inputs(expr, meta, model, contrast, de_method, contrast_matrix)
      is.categorical <- l$is.categorical
      group_vec <- l$group_vec
      design_full <- l$design_full

      # Define reduced design for LRT de_type in DESeq2
      if (de_method == "DESeq2" && de_type == "LRT") {
        if (is.character(contrast)) {
          design_reduced <- design_full[, colnames(design_full) != contrast]
        } else if (is.numeric(contrast)) {
          design_reduced <- design_full[, -contrast]
        }
      }

      if (verbose) print(data.frame(design_full))

      # Run DE
      de <- switch(
        de_method,
        edgeR = run_edgeR(expr, meta, design_full, contrast, de_type,
                          group_vec),
        DESeq2 = run_DESeq2(expr, meta, design_full, contrast, de_type),
        limma = run_limma(expr, meta, design_full, contrast, de_type,
                          contrast_matrix, is.categorical, group_vec),
      )

      de[["design_full"]] <- design_full
      if (de_method == "DESeq2" && de_type == "LRT") {
        de[["design_reduced"]] <- design_reduced
      }

      return(de)
    }
  )

  # Return
  if (return_interm) {
    des
  } else {
    purrr::map_dfr(des, "res", .id = "cell_type")
  }
}

diagnose_ruv <- function(input, meta, file_name = NULL, model, contrast,
                         ruv_type = "RUVr", cIdx = NULL, scIdx = NULL,
                         resids = NULL, min_k = 5, max_k = 5,
                         de_method = "limma", de_type = "voombygroup",
                         plotpcs = F, verbose = 1) {
  #' Diagnose RUVSeq
  #'
  #' @description
  #' Returns diagnostic plots and metrics to evaluate the impact of addign RUVseq factors to a differential expression analysis.on the number of DEGs.
  #'
  #' @param input list or matrix. List of cell types with list of expr and meta generated by [Libra::to_pseudobulk()].`rownames(pseudobulks$meta)` must match `colnames(pseudobulks$expr)`. It can also be a count matrix of gene x observation  satisfying `colnames(pseudobulks) == rownames(meta)` and raw counts.
  #' @param meta data.frame. Table containing the per observation covariates included in model. If input is a matrix, meta is required. Is input is a list, meta is expected per each entry.
  #' @param file_name character. Path to file to write PC analysis PDFs if `plotpcs = T`.
  #' @param model formula. Constructed with [as.formula()] or `~` and containing terms in `colnames(meta)` or `colnames(pseudobulks$meta)`.
  #' @param contrast numeric or character. Position or name of the contrast in `colnames(design.matrix(model, data = meta))`.
  #' @param ruv_type character. One of `c("RUVr", "RUVg", "RUVs")`.
  #' @param cIdx,scIdx,resids As in [RUVSeq::RUVs()].
  #' @param min_k,max_k numeric. Number of factors to use for [RUVSeq]. If min_k != max_k it will run all the factors. If min_k == max_k it will only run one iteration
  #' @param de_method,de_type As in [de_pseudobulk].
  #' @param plotpcs logical. Whether to plot PCs for factor 1 and 2 for all covariates in model.
  #' @param verbose logical. Level fo information displayed. 0 to silence, 1 to display progress information and cell type names, 2 to display per k iteration ifnormation and 3 to print matrix of DEG by k.
  #'
  #' @returns
  #' List per cell type with a list containing de_method x de_type data frame of DEGs (`degs`), a data frame of all genes (`de`), a data frame of k vs number of DEGs (`k_vs_ndegs`), a data frame of k vs total explained variance by covariates (`varparts`), nominal p-value distributions (`pval_quantiles`), pearson correlation of the differential expression statistic without and with RUVSeq factors in the model across all genes (`stat_cors`) and the k that maximizes the variance explained for the term of the contrast (`best_k`).

  # Intiialize output
  if (is.null(file_name)) {
    file_name <- paste0(getwd(), "/diagnose_", ruv_type)
  }
  base_dir <- dirname(file_name)

  if (verbose >= 1) cat("Baseline differential expression")
  des <- de_pseudobulk(input = input, meta = meta, model = model,
                       contrast = contrast, de_method = de_method,
                       de_type = de_type, return_interm = T)

  # Loop over each cell type
  furrr::future_imap(des, function(res, cell_type) {
    if (verbose >= 1) cat(paste0("Diagnosis RUVSeq for cell type: ", cell_type,
                                 "\n"))

    # Get normalized and log transformed counts
    expr <- switch(
      de_method,
      edgeR = edgeR::cpm(res$interm$y, log = T),
      DESeq2 = DESeq2::rlog(DESeq2::counts(res$interm$dds, normalized = T)),
      limma = res$interm$v$E
    )

    # Get residuals if not supplied and ruv_type == "RUVr"
    if (is.null(resids) & ruv_type == "RUVr") {
      resids <- switch(
        de_method,
        edgeR = residuals(res$interm$fit),
        DESeq2 = SummarizedExperiment::assays(res$interm$dds)[["mu"]],
        limma = residuals(res$interm$fit, res$interm$v)
      )
    }

    # Store design matrix in case verbose >= 3 for gene count plots
    res_design_full <- res$design_full

    # Get the same cell_type name as de_pseudobulk
    res <- dplyr::bind_rows(res$res, .id = "cell_type")
    qs <- quantile(res$p.value)

    # Get variance partition
    varpart <- variancePartition::fitExtractVarPartModel(expr, model, meta)
    colnames(varpart) <- make.names(colnames(varpart), unique = T)
    varpartmeds <- Rfast::colMedians(varpart)

    # Store k = 0 results
    k0 <- list(
      degs = res[res$fdr < .05, ],
      de = res,
      varparts = c(varpartmeds,
                   total = sum(varpartmeds[names(varpartmeds) != "Residuals"])),
      k_vs_ndegs = nrow(res[res$fdr < .05, ]),
      pval_quantiles = data.frame(
        q0 = qs[1], q25 = qs[2], q50 = qs[3],
        q75 = qs[4], q100 = qs[5], row.names = NULL
      )
    )

    # Run iterations
    if (verbose >= 1) cat(paste0("RUVseq with max_k: ", max_k, "\n"))
    ks_res <- furrr::future_map(min_k:max_k, function(k) {
      if (verbose >= 3) cat(paste0("\n    k: ", k, "\n"))

      # Define RUV type
      ks <- switch(
        ruv_type,
        RUVg = RUVSeq::RUVg(x = expr, cIdx = cIdx, k = k, round = F, isLog = T),
        RUVs = RUVSeq::RUVs(x = expr, cIdx = rownames(expr), scIdx = scIdx,
                            k = k, round = F, isLog = T),
        RUVr = RUVSeq::RUVr(x = expr, cIdx = rownames(expr), k = k, round = F,
                            residuals = resids, isLog = T)
      )

      if (plotpcs) {
        if (verbose >= 3) cat("    Plotting PCs\n")
        pc <- prcomp(t(ks$normalizedCounts), center = T, scale. = T)
        summary(pc)$importance[, 1:3]
        sort(pc$rotation[, 1], decreasing = T)[1:20]
        suppressMessages(
          purrr::imap(colnames(meta), ~ {
            ggplot(cbind(pc$x[, 1:2], meta),
                   aes(PC1, PC2, color = .data[[.x]])) +
              geom_point() +
              labs(title = paste("Number of ks: ", k)) +
              theme_classic()
            ggsave(paste0(base_dir, "/temp_", .y, "_", k, ".pdf"))
          })
        )
        l <- list.files(path = base_dir, pattern = "temp_", full.names = T)
        qpdf::pdf_combine(l, paste0(file_name, "_pcs", k, "_2.pdf"))
        qpdf::pdf_compress(
          paste0(file_name, "_pcs", k, "_2.pdf"),
          paste0(file_name, "_pcs", k, ".pdf")
        )
        invisible(file.remove(c(l, paste0(file_name, "_pcs", k, "_2.pdf"))))
        rm(l)
      }

      # Include RUV factors in DE
      model2 <- as.formula(
        paste(
          "~",
          paste(
            as.character(model)[2],
            paste(paste0("W_", 1:k), collapse = " + "),
            sep = " + "
          )
        )
      )
      if (is.list(input)) {
        counts <- input[[cell_type]][["expr"]]
        meta2 <- cbind(input[[cell_type]][["meta"]], ks$W)
      } else if (is.matrix(input) & is.data.frame(meta)) {
        counts <- input
        meta2 <- cbind(meta, ks$W)
      }
      res2 <- de_pseudobulk(input = counts, meta = meta2, model = model2,
                            contrast = contrast, de_method = de_method,
                            de_type = de_type, return_interm = F)

      # Get DE diagnostics
      qs <- quantile(res2$p.value)
      stat_cor <- cor(abs(res$stat), abs(res2$stat), method = "spearman")

      # Variance partition
      if (verbose >= 3) cat("    Variance partition\n")
      varpart <- variancePartition::fitExtractVarPartModel(expr, model2, meta2)
      colnames(varpart) <- make.names(colnames(varpart), unique = T)
      varpart_covs <- colnames(varpart)
      if (verbose >= 3) {
        for (cov in varpart_covs) {
          assign(cov, varpart[, cov])
        }
        eval(parse(text = paste0("txtplot::txtboxplot(",
                                 paste(varpart_covs, collapse = ", "), ")")))
        rm(list = varpart_covs)
      }
      varpartmeds <- Rfast::colMedians(varpart)

      # Display information
      if (verbose >= 3) {
        cat(paste0("    DEGs (n): ", nrow(res2[res2$fdr < .05, ]), "\n"))
        cat(paste0("    Median p-value: ", median(res2$p.value), "\n"))
        cat("    p-value distribution: \n")
        print(txtplot::txtdensity(res2$p.value))
        cat(paste0(
          "    Correlation between the uncorrected and corrected t: ",
          stat_cor, "\n"
        ))
        print(txtplot::txtplot(abs(res$stat), abs(res2$stat)))
      }

      # Return
      list(
        degs = res2[res2$fdr < .05, ],
        de = res2,
        varparts = c(varpartmeds, total = sum(varpartmeds[names(varpartmeds)
                                                          != "Residuals"])),
        k_vs_ndegs = nrow(res2[res2$fdr < .05, ]),
        pval_quantiles = data.frame(
          q0 = qs[1], q25 = qs[2], q50 = qs[3],
          q75 = qs[4], q100 = qs[5], row.names = NULL
        ),
        stat_cors = stat_cor
      )
    },
    .options = furrr::furrr_options(seed = T)) %>%
      setNames(min_k:max_k)
    ks_res <- c(list("0" = k0), ks_res)

    # Convert output lists into dataframes
    degs <- purrr::map_dfr(ks_res, "degs", .id = "k") %>%
      dplyr::mutate(k = as.integer(k))
    de <- purrr::map_dfr(ks_res, "de", .id = "k") %>%
      dplyr::mutate(k = as.integer(k))
    varparts <- purrr::map_dfr(ks_res, "varparts", .id = "k") %>%
      dplyr::mutate(k = as.integer(k))
    k_vs_ndegs <- purrr::map_dfr(ks_res, ~{
      tibble::enframe(unlist(.x[["k_vs_ndegs"]]), name = "k", value = "n_degs")
    }, .id = "k") %>%
      dplyr::mutate(k = as.integer(k))
    pval_quantiles <- purrr::map_dfr(ks_res, "pval_quantiles", .id = "k") %>%
      dplyr::mutate(k = as.integer(k))
    stat_cors <- purrr::map_dfr(ks_res, ~{
      tibble::enframe(unlist(.x$stat_cors), name = "k", value = "cor")
    }, .id = "k") %>%
      dplyr::mutate(k = as.integer(k))

    # Get k that maximizes the variance explained for the term of interest
    contrast_term <- all.vars(model)[grep(paste(all.vars(model),
                                                collapse = "|"), contrast)]
    best_k <- varparts$k[which.max(varparts[[contrast_term]])]

    if (plotpcs) {
      # Combine PC PDFs into a single pdf
      l <- list.files(
        path = base_dir, pattern = paste0(basename(file_name), "_pcs"),
        full.names = T
      )
      qpdf::pdf_combine(l, paste0(file_name, "_pcs_2.pdf"))
      qpdf::pdf_compress(
        paste0(file_name, "_pcs_2.pdf"),
        paste0(file_name, "_pcs.pdf")
      )
      invisible(file.remove(c(l, paste0(file_name, "_pcs_2.pdf"))))
      rm(l)
    }

    # Print cross k information
    if (verbose >= 1) {
      cat("RUV k vs n_degs\n")
      print(txtplot::txtplot(k_vs_ndegs$k, k_vs_ndegs$n_degs))

      cat("Total variance explained\n")
      print(txtplot::txtplot(varparts$k, varparts$total))
    }

    if (verbose >= 2) {
      cat("Number of ks supporting a DEG\n")
      k_v_deg <- table(degs$gene, degs$k)
      cbind(k_v_deg, total = rowSums(k_v_deg)) %>%
        as.data.frame() %>%
        dplyr::arrange(total) %>%
        knitr::kable() %>%
        message()
    }

    if (verbose >= 3 & nrow(degs) != 0) {
      is.categorical <- length(unique(res_design_full[, contrast])) == 2
      n_degs <- length(unique(degs$gene))
      max_n_degs <- ifelse(n_degs > 10, 10, n_degs)
      if (is.categorical) {
        # Get sample names per category
        expr_0 <- rownames(res_design_full)[res_design_full[, contrast] == 0]
        expr_1 <- rownames(res_design_full)[res_design_full[, contrast] == 1]

        # Plot each gene
        cat(paste0("Top ", max_n_degs, " DEGs for contrast ", contrast, " with ",
                length(expr_0), " controls and ", length(expr_1), " cases.\n"))
        for (i in 1:max_n_degs) {
          gene <- degs %>%
            dplyr::count(gene) %>%
            dplyr::arrange(dplyr::desc(n)) %>%
            magrittr::extract(i, "gene")
          print(gene)
          print(txtplot::txtboxplot(
            expr[rownames(expr) == gene, expr_0],
            expr[rownames(expr) == gene, expr_1]
          ))
        }
      } else {
        for (i in 1:max_n_degs) {
          gene <- degs %>%
            dplyr::count(gene) %>%
            dplyr::arrange(dplyr::desc(n)) %>%
            magrittr::extract(i, "gene")
          print(gene)
          print(txtplot::txtboxplot(expr[rownames(expr) == gene, ]))
        }
      }
    }

    # Return
    list(
      degs = degs, de = de, k_vs_ndegs = k_vs_ndegs, varparts = varparts,
      pval_quantiles = pval_quantiles, stat_cors = stat_cors, best_k = best_k
    )
  }, .options = furrr::furrr_options(seed = T))
}

wFisher <- function(p, weight = NULL, is.onetail = TRUE, eff.sign) {
  #' wFisher
  #'
  #' Sample size-weighted Fisher's method
  #'
  #' @param p A numeric vector of p-values
  #' @param weight A numeric vector of weight or sample size for each experiment
  #' @param is.onetail Logical. If set TRUE, p-values are combined without considering the direction of effects, and vice versa. Default: TRUE.
  #' @param eff.sign A vector of signs of effect sizes. It works when is.onetail = FALSE
  #' @return A list with combined p-value and, if applicable, the direction of combined effects.
  #' @examples
  #' wFisher(p = c(0.01, 0.2, 0.8), weight = c(50, 60, 100), is.onetail = FALSE, eff.sign = c(1, 1, 1))
  #' @importFrom "stats" "pgamma" "qgamma"
  #' @references
  #' Becker BJ (1994). "Combining significance levels." In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215–230. Russell Sage, New York.
  #' Fisher RA (1925). Statistical methods for research workers. Oliver and Boyd, Edinburgh.
  #'
  #' @source https://github.com/unistbig/metapro/blob/master/R/metapro.R
  #' @copiedon 05/30/2023
  #'
  #' @export

  if (is.null(weight)) {
    weight <- rep(1, length(p))
  }
  idx.na <- which(is.na(p))
  if (length(idx.na) > 0) {
    p <- p[-idx.na]
    weight <- weight[-idx.na]
    if (!is.onetail) {
      eff.sign <- eff.sign[-idx.na]
    }
  }
  NP <- length(p)
  NS <- length(weight)
  if (NP != NS) {
    stop("The length of p and weight vector must be identical.")
  }
  N <- NS
  Ntotal <- sum(weight)
  ratio <- weight / Ntotal
  Ns <- N * ratio
  G <- c()

  if (is.onetail) {
    for (i in 1:length(p))
    {
      G <- append(G, qgamma(p = p[i], shape = Ns[i], scale = 2, lower.tail = F))
    }
    Gsum <- sum(G)
    resultP <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = F)
  } else {
    p1 <- p2 <- p
    idx_pos <- which(eff.sign > 0)
    idx_neg <- which(eff.sign < 0)
    # positive direction
    G <- c()
    p1[idx_pos] <- p[idx_pos] / 2
    p1[idx_neg] <- 1 - p[idx_neg] / 2
    for (i in 1:length(p1))
    {
      G <- append(G, qgamma(p = p1[i], shape = Ns[i], scale = 2, lower.tail = F))
    }
    Gsum <- sum(G)
    resultP1 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = F)
    # negative direction
    G <- c()
    p2[idx_pos] <- 1 - p[idx_pos] / 2
    p2[idx_neg] <- p[idx_neg] / 2
    for (i in 1:length(p2))
    {
      G <- append(G, qgamma(p = p2[i], shape = Ns[i], scale = 2, lower.tail = F))
    }
    Gsum <- sum(G)
    resultP2 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = F)
    resultP <- 2 * min(resultP1, resultP2)
    if (resultP > 1.0) {
      resultP <- 1.0
    }
    overall.eff.direction <- if (resultP1 <= resultP2) {
      "+"
    } else {
      "-"
    }
  }
  RES <- if (is.onetail) {
    list(p = min(1, resultP))
  } else {
    list(p = min(1, resultP), overall.eff.direction = overall.eff.direction)
  }
  return(RES)
}

get_props_tidyverse <- function(data, id, cluster, supercluster = NULL,
                                add_supercluster_prop = T) {
  #' Compute proportions per id across cluster
  #'
  #' This function computes the proportions per donor per broad from the input data.
  #'
  #' @param data A data frame containing at least the columns id, cluster, and supercluster.
  #' @param id The column name for the sample information.
  #' @param cluster The column name for the cluster information.
  #' @param supercluster The column name for the supercluster information.
  #' @param add_supercluster_prop Logical indicating whether to include proportions per donor across broad.
  #'
  #' @return A data frame with computed proportions per donor per broad.
  #'
  #' @examples
  #' # Example usage
  #' df <- data.frame(donor = rep(1:3, times = 4),
  #'                  cell_type = rep(letters[1:4], each = 3, times = 4),
  #'                  broad = rep(c("AB", "CD"), each = 6, times = 2))
  #' #' get_props(df, donor, cell_type, broad)
  #'
  #' @importFrom dplyr group_by summarize n mutate ungroup bind_rows filter distinct
  #' @importFrom tidyr complete fill
  #' @importFrom tibble remove_rownames
  #' @importFrom rlang as_string ensym
  #' @importFrom purrr map
  #' @importFrom magrittr %>%
  #' @importFrom stats sum
  #' @importFrom base unique
  #'
  #' @export

  # Convert column names to symbols if they are strings
  id <- rlang::ensym(id)
  cluster <- rlang::ensym(cluster)
  supercluster <- rlang::ensym(supercluster)

  # Check if input data is a data frame
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  # Check if required columns exist in the input data frame
  required_cols <- c(rlang::as_string(supercluster),
                     rlang::as_string(cluster),
                     rlang::as_string(id))
  if (!all(required_cols %in% names(data))) {
    missing_cols <- required_cols[!required_cols %in% names(data)]
    stop(paste("Input data is missing columns:",
               paste(missing_cols, collapse = ", ")))
  }

  if (!is.null(supercluster)) {
    # Compute proportions per id x cluster x supercluster
    result <- data %>%
      dplyr::group_by(!!supercluster, !!cluster, !!id) %>%
      dplyr::summarize(num_cluster = dplyr::n(), .groups = "drop_last") %>%
      dplyr::group_by(!!supercluster, !!id) %>%
      dplyr::mutate(num_supercluster = sum(num_cluster, na.rm = T),
                    prop = num_cluster / num_supercluster) %>%
      dplyr::ungroup()

    # Add proportions per id across supercluster if requested
    if (add_supercluster_prop) {
      supercluster_prop <- result %>%
        dplyr::group_by(!!supercluster, !!id) %>%
        dplyr::summarize(num_cluster = unique(num_supercluster),
                         .groups = "drop_last") %>%
        dplyr::group_by(!!id) %>%
        dplyr::mutate(num_supercluster = sum(num_cluster),
                      prop = num_cluster / num_supercluster,
                      !!cluster := !!supercluster,
                      !!supercluster := rlang::as_string(supercluster)) %>%
        dplyr::ungroup()

      result <- dplyr::bind_rows(result, supercluster_prop)
    }
  } else {
    # Compute proportions per id x cluster x supercluster
    result <- data %>%
      dplyr::group_by(!!cluster, !!id) %>%
      dplyr::summarize(num_cluster = dplyr::n(), .groups = "drop_last") %>%
      dplyr::group_by(!!id) %>%
      dplyr::mutate(num_supercluster = sum(num_cluster, na.rm = T),
                    prop = num_cluster / num_supercluster) %>%
      dplyr::ungroup()
  }

  # Add props=0 for the ids without a cluster
  result <- result %>%
    dplyr::group_by(!!supercluster) %>%
    tidyr::complete(!!cluster, !!id,
                    fill = list(num_cluster = 0, prop = 0)) %>%
    dplyr::group_by(!!supercluster, !!id) %>%
    tidyr::fill(num_supercluster, .direction = "downup") %>%
    dplyr::ungroup()

  # Safety tests
  # Ensure all levels of required arguments are present in the output
  # same_unique_elements <- purrr::map_lgl(required_cols, function(i) {
  #   length(unique(data[[i]])) == length(unique(results[[i]]))
  #   })
  # if (any(!same_unique_elements)) stop("")

  # Ensure no new combinations were created except id x cluster
  original_combs <- data %>%
    dplyr::distinct(!!supercluster, !!cluster, !!id) %>%
    tibble::remove_rownames()
  result_combs <- result %>%
    dplyr::filter(!!supercluster != rlang::as_string(supercluster))
  new_combs <- dplyr::anti_join(result_combs, original_combs,
                                by = required_cols) %>%
    dplyr::summarize(num_cluster = sum(num_cluster), .groups = "drop_last") %>%
    dplyr::pull(num_cluster) %>%
    unique()
  if (new_combs != 0) stop("New id x cluster combinations were created.")

  # Ensure proportions add up to 1
  props_sum <- result %>%
    dplyr::group_by(!!supercluster, !!id) %>%
    dplyr::summarize(n = round(sum(prop), 2), .groups = "drop_last") %>%
    dplyr::pull(n) %>%
    unique()
  if (props_sum != 1) stop("Porportions don't add up to 1.")

  # Clean up
  rm(original_combs, result_combs, new_combs, props_sum)

  return(result)
}

summarize_by_group <- function(data, group, columns = NULL) {
  #' Summarize data by group
  #'
  #' This function aggregates data by group and summarizes specified columns.
  #'
  #' @param data A data frame.
  #' @param group A character vector specifying the grouping column(s).
  #' @param columns Optional character vector specifying columns to summarize.
  #' @return A data frame with summarized data.
  #' @examples
  #' mixed_df <- data.frame(
  #'   Group = rep(letters[1:3], each = 3),
  #'   Numeric_Value = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  #'   Factor_Value = factor(c("low", "medium", "high"),
  #'    levels = c("low", "medium", "high")),
  #'   Character_Value = c("apple", "banana", "apple", "banana", "apple",
  #'    "banana", "apple", "banana", "apple"),
  #'   Logical_Value = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE,
  #'    TRUE),
  #'   Date_Value = as.Date("2022-01-01") + 0:8,
  #'   Complex_Value = as.complex(1:9)
  #' )
  #'
  #' # Summarize all columns (default) by Group
  #' summarize_by_group(mixed_df, "Group")
  #'
  #' # Summarize column Numeric_Value by Group
  #' summarize_by_group(mixed_df, "Group", "Numeric_Value")
  #'
  #' # Summarize two columns by Group
  #' summarize_by_group(mixed_df, "Group", c("Numeric_Value", "Factor_Value"))
  #'
  #' @export

  # Validate input
  if (!is.data.frame(data)) {
    stop("data must be a dataframe")
  }
  if (!is.character(group) || length(group) == 0) {
    stop("group must be a non-empty character vector")
  }
  if (!is.null(columns) && (!is.character(columns) || length(columns) == 0)) {
    stop("columns must be a non-empty character vector")
  }

  # Get the columns to summarize
  if (is.null(columns)) {
    columns <- setdiff(names(data), group)
  }

  # If no columns to summarize, return unique combinations of group columns
  if (length(columns) == 0) {
    df <- dplyr::distinct(data, dplyr::across(tidyselect::all_of(group)))
  }

  # If columns need to be summarized, summarize columns by group
  df <- data %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(group))) %>%
    dplyr::summarize(
      dplyr::across(
        dplyr::all_of(columns),
        ~ {
          if (all(is.na(.))) NA
          else if (is.numeric(.)) mean(., na.rm = T)
          else if (is.factor(.)) levels(.)[1]
          else if (is.character(.)) unique(.)[which.max(table(.))]
          else if (is.logical(.)) unique(.)[which.max(table(.))]
          else if (inherits(., c("Date", "POSIXt"))) as.character(unique(.))[1]
          else if (is.list(.)) NA
          else if (is.complex(.)) NA
          else NA
        }
      ),
      .groups = "drop_last"
    ) %>%
    dplyr::ungroup()

  return(df)
}

get_props <- function(data, id, cluster, supercluster = NULL,
                      add_supercluster_prop = F, add_other_cols = T) {
  #' Compute Proportions
  #'
  #' This function computes proportions per id across cluster. It can group by
  #'  supercluster and compute proportions per group. It can also add
  #'  supercluster proportions across all its levels. Finally, it summarizes
  #'  other columns of data per id.
  #'
  #' @param data Data frame containing \code{id} and \code{cluster}.
  #' @param id The column name representing the observations identifier.
  #' @param cluster The column name representing the cluster identifier.
  #' @param supercluster Optional. The column name representing the supercluster
  #'  identifier.
  #' @param add_supercluster_prop Logical indicating whether to
  #'  include proportions per id across superclusters.
  #' @param add_other_cols Logical indicating whether to summarize other columns of \code{data} by \code{id}. Numeric columns are summarized to their mean, factor to the reference level and character to the mode.
  #'
  #' @return A long data frame with computed proportions per group and summary
  #'  statistics.
  #' - \code{"num_cluster"}: The number of occurrences of each group combination.
  #' - \code{"props"}: The proportions of each group combination within its
  #'  supercluster.
  #' - \code{"supercluster"}: The supercluster associated with each group.
  #' - \code{"num_supercluster"}: The total number of occurrences of each
  #'  \code{supercluster}-\code{id} combination.
  #'
  #' @details This function calculates proportions per group based on the counts of another
  #' column, and summarizes other columns by group, calculating means for numeric columns,
  #' reference levels for factor columns, and the mode for character columns.
  #'
  #' @examples
  #' df <- data.frame(id = rep(1:3, each = 8),
  #'                  cluster = rep(c("C1", "C2"), each = 4, times = 3),
  #'                  supercluster = rep(c("S1", "S2"), each = 4, times = 3),
  #'                  value1 = runif(24),
  #'                  value2 = as.factor(LETTERS[1:12]),
  #'                  value3 = LETTERS[1:12],
  #'                  value4 = NA)
  #' df
  #' get_props(df, "id", "cluster", "supercluster", add_supercluster_prop = T)
  #' get_props(df, "id", "cluster", "supercluster")
  #' get_props(df, "id", "cluster")
  #'
  #' @import dplyr
  #' @import tidyr
  #' @importFrom magrittr %>%
  #'
  #' @export
  #'
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data frame.")
  }

  df <- data.frame(id = rep(1:3, each = 8),
                   cluster = rep(c("C1", "C2"), each = 4, times = 3),
                   supercluster = rep(c("S1", "S2"), each = 4, times = 3),
                   value1 = runif(24),
                   value2 = as.factor(LETTERS[1:12]),
                   value3 = LETTERS[1:12],
                   value4 = NA)

  # Input Validation
  required_cols <- c(id, cluster)
  if (!is.null(supercluster)) {
    required_cols <- c(required_cols, supercluster)
  }
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("The following columns are missing in the data: ",
         paste(missing_cols, collapse = ", "))
  }
  invalid_args <- c(
    "id" = !is.character(id),
    "cluster" = !is.character(cluster),
    "supercluster" = !is.null(supercluster) && !is.character(supercluster)
  )
  invalid_args <- names(invalid_args)[invalid_args]
  if (length(invalid_args) > 0) {
    stop("Argument(s) ", paste(invalid_args, collapse = ", "),
         " must be character vector(s).")
  }

  # Proportion Calculation
  if (!is.null(supercluster)) {
    # Proportions by supercluster x cluster x id
    data_groups <- split(data, data[[supercluster]])
    props <- lapply(data_groups, function(data_group) {
      # Proportion per supercluster by cluster x id
      nums_cluster <- table(data_group[[id]], data_group[[cluster]])
      prop <- as.data.frame(prop.table(nums_cluster, margin = 1))
      nums_cluster <- as.data.frame(nums_cluster)

      # Count observations per supercluster x id
      nums_supercluster <- data.frame(table(data_group[[id]],
                                            data_group[[supercluster]]))

      # Left join and rename
      prop <- merge(nums_cluster, prop, by = c("Var1", "Var2"))
      prop <- merge(prop, nums_supercluster, by = "Var1")
      colnames(prop) <- c(id, cluster, paste0("num_", cluster), "prop",
                           supercluster, paste0("num_", supercluster))

      return(prop)
    })
    props <- do.call(rbind, props)
    rownames(props) <- NULL

    # Add proportions per supercluster column if requested
    if (add_supercluster_prop) {
      # Proportion by supercluster x id
      nums_supercluster <- table(data[[id]], data[[supercluster]])
      prop <- as.data.frame(prop.table(nums_supercluster, margin = 1))
      nums_supercluster <- as.data.frame(nums_supercluster)

      # Count observations per id across superclusters
      nums_id <- data.frame(table(data[[id]]))
      nums_id$Var2 <- supercluster
      nums_id <- nums_id[, c("Var1", "Var2", "Freq")]

      # Left join and rename
      prop <- merge(nums_supercluster, prop, by = c("Var1", "Var2"))
      prop <- merge(prop, nums_id, by = "Var1")
      colnames(prop) <- c(id, cluster, paste0("num_", cluster), "prop",
                          supercluster, paste0("num_", supercluster))

      props <- rbind(props, prop)
    }
  } else {
    # Proportion by cluster x id
    nums_cluster <- table(data[[id]], data[[cluster]])
    prop <- as.data.frame(prop.table(nums_cluster, margin = 1))
    nums_cluster <- as.data.frame(nums_cluster)

    # Merge and clean column names
    props <- merge(nums_cluster, prop, by = c("Var1", "Var2"))
    colnames(props) <- c(id, cluster, paste0("num_", cluster), "prop")

    props
  }

  # Ensure column types conform to original data
  props[[id]] <- as(props[[id]], class(data[[id]]))
  props[[cluster]] <- as(props[[cluster]], class(data[[cluster]]))
  if (!is.null(supercluster)) {
    props[[supercluster]] <- as(props[[supercluster]],
                                class(data[[supercluster]]))
  }

  # Summarize other columns by id
  other_columns <- setdiff(names(data), required_cols)
  if (add_other_cols && length(other_columns) > 0) {
    df <- summarize_by_group(data, id, other_columns)
    props <- dplyr::left_join(props, df, by = id)
  }

  return(props)
}

get_response <- function(formula) {
  if (class(formula) != "formula") {
    stop("Not formula object")
  }
  all.vars(formula[[length(formula) - 1]])
}

get_covs <- function(formula) {
  if (class(formula) != "formula") {
    stop("Not formula object")
  }
  all.vars(formula[[length(formula)]])
}

robust_glm <- function(formula, data, subset = NULL, family = "quasibinomial",
                       robust_weights = T, sandwich = T, add_ci = T, p = NULL,
                       ...) {
  #' Fit a Generalized Linear Model with Robust Weights and Errors
  #'
  #' This function fits a generalized linear model (GLM) with robust weights per
  #'  observation, to down-weight outliers, and robust sandwich standard errors,
  #'   to account for the lack of independance or heteroscedasticity.
  #'
  #' @param formula A formula specifying the model.
  #' @param data A data frame containing the variables in the model.
  #' @param subset An optional character vector specifying a subset of observations
  #'               to be used in the model.
  #' @param family A character string or function (see \code{lm()}) specifying
  #' the distribution family in the GLM (default is "binomial").
  #' @param robust_weights Logical indicating whether to compute robust model weights
  #'                       (default is TRUE).
  #' @param sandwich Logical indicating whether to compute sandwich standard errors
  #'                 (default is TRUE).
  #' @param add_ci Logical indicating whether to add confidence intervals to the model
  #'               coefficients (default is TRUE).
  #' @param p An optional progressor object to monitor progress (default is NULL).
  #' @param ... Other arguments passed on to \code{stats::glm}.
  #'
  #' @return An object of class \code{"glm"} with additional attributes such as
  #'         confidence intervals, sandwich standard errors, and collinear terms.
  #'
  #' @examples
  #' df <- data.frame(y = runif(50, 0, 1),
  #'                  x1 = rep(1:2, 50),
  #'                  x2 = runif(50, .5, 1),
  #'                  x3 = runif(50, 0, .5))
  #' fit <- robust_glm(y ~ x1 + x2, df, family = "quasibinomial")
  #' summary(fit)
  #'
  #' @importFrom MASS rlm
  #' @importFrom bestNormalize bestNormalize
  #' @importFrom sandwich vcovHC
  #' @importFrom lmtest coeftest coefci
  #' @export

  # Input validation
  if (missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' arguments must be provided.")
  }
  if (!inherits(formula, "formula") && !is.character(formula)) {
    stop("Argument 'formula' must be either a character formula or a formula object.")
  }
  if (!inherits(formula, "formula") && is.character(formula)) {
    formula <- as.formula(formula)
  }
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }
  if (!is.null(subset) && !is.character(subset)) {
    stop("Argument 'subset' must be a character vector.")
  }

  # Subset data
  if (!is.null(subset)) {
    data <- subset(data, eval(parse(text = subset)))
  }

  # Extract model terms
  y <- get_response(formula)
  covs <- get_covs(formula)

  # Check if terms are present in the data
  missing_terms <- c(y, covs)[!(c(y, covs) %in% names(data))]
  if (length(missing_terms) > 0) {
    stop("Variable(s) '", paste(missing_terms, collapse = "', '"),
         "' not found in the data.")
  }

  # Clean data
  data <- data[, c(y, covs)]
  data <- na.omit(data)

  # Extract the response variable values
  prop <- data[[y]]

  # Set up the progressor progress bar
  if (!is.null(p) && inherits(p, "progressor")) {
    p()
  }

  # Check for collinear terms
  x <- model.matrix(formula, data)
  ncovs <- ncol(x)
  QR <- qr(x)
  if (QR$rank < ncovs) {
    collinear_terms <- colnames(x)[QR$pivot[(QR$rank + 1):ncovs]]
    if (all(collinear_terms %in% covs)) {
      formula <- reformulate(covs[!covs %in% collinear_terms], response = y)
    }
    collinear_terms <- paste(collinear_terms, collapse = ", ")

    warning("Dropped collinear terms: ", collinear_terms)
  } else {
    collinear_terms <- NULL
  }

  # Compute robust model weights if requested
  if (robust_weights && is.numeric(data[[y]]) &&
      length(unique(data[[y]])) > 2) {
    # From [0, 1] to (-Inf, +Inf)
    data$prop_norm <- bestNormalize::bestNormalize(prop, loo = T, quiet = T)$x.t
    robust_formula <- update.formula(formula, prop_norm ~ .)
    rweights <- MASS::rlm(robust_formula, data = data)$w

    if (length(rweights) != length(prop)) {
      stop("Weight and response have different lengths. Any NA maybe?")
    }

    data$rweights <- rweights

    # Fit a generalized linear model with robust weights
    fit <- stats::glm(formula, data = data, family = family, weights = rweights,
                      ...)
  } else {
    fit <- stats::glm(formula, data = data, family = family, ...)
  }

  if (add_ci) {
    fit$ci <- suppressMessages(confint(fit))
    colnames(fit$ci) <- c("conf.low", "conf.high")
  }

  # Add sandwich errors
  if (sandwich) {
    vcov_mat <- sandwich::vcovHC(fit, type = "HC3")
    fit$sandwich <- lmtest::coeftest(fit, vcov. = vcov_mat)
    sandwich_ci <- lmtest::coefci(fit, vcov. = vcov_mat)
    fit$sandwich <- cbind(fit$sandwich, sandwich_ci)
    colnames(fit$sandwich) <- c("estimate", "std.error", "statistic", "p.value",
                                "conf.low", "conf.high")
  }

  # Add collinear terms
  fit$collinear_terms <- collinear_terms

  # Add class robust_glm
  class(fit) <- c(class(fit), "robust_glm")

  return(fit)
}

tidy_terms <- function(model, id = NULL, exponentiate = F) {
  #' Tidy Model Terms
  #'
  #' This function tidies model terms, including coefficients, confidence intervals,
  #' and sandwich standard errors, from a list of model objects or a single model object.
  #'
  #' @param model A list of or a single robust_glm object.
  #' @param id A character string specifying the identifier column name
  #'  (default is "id").
  #' @param exponentiate Logical. Exponentiate estimate and confidence
  #'  intervals.
  #'
  #' @return A data frame containing tidied model terms with columns:
  #' \describe{
  #'   \item{\code{id}}{Identifier column name.}
  #'   \item{\code{term}}{Model term names.}
  #'   \item{\code{estimate}}{Estimate of coefficients.}
  #'   \item{\code{std.error}}{Standard error of coefficients.}
  #'   \item{\code{statistic}}{Value of test statistics.}
  #'   \item{\code{p.value}}{p-value of test statistics.}
  #'   \item{\code{conf.low}}{Lower bound of confidence interval.}
  #'   \item{\code{conf.high}}{Upper bound of confidence interval.}
  #'   \item{\code{std.error_hc}}{Sandwich standard error of coefficients.}
  #'   \item{\code{statistic_hc}}{Value of test statistics with sandwich
  #'    standard errors.}
  #'   \item{\code{p.value_hc}}{p-value of test statistics with sandwich
  #'    standard errors.}
  #'   \item{\code{conf.low_hc}}{Lower bound of confidence interval with
  #'    sandwich standard errors.}
  #'   \item{\code{conf.high_hc}}{Upper bound of confidence interval with
  #'    sandwich standard errors.}
  #'    \item{\code{estimate_exp}}{Exponentiated estimate.}
  #'   \item{\code{conf.low_exp}}{Lower bound of exponentiated confidence
  #'    interval.}
  #'   \item{\code{conf.high_exp}}{Upper bound of exponentiated confidence
  #'    interval.}
  #'   \item{\code{std.error_exp}}{Exponentiated standard error.}
  #'   \item{\code{std.error_hc_exp}}{Exponentiated sandwich standard error.}
  #' }
  #'
  #' @importFrom broom tidy
  #' @importFrom tibble rownames_to_column
  #' @importFrom dplyr left_join select rename_with
  #' @importFrom purrr map_dfr reduce
  #' @export

  # If lm or glm supplied, convert to a list
  if (inherits(model, c("lm", "glm"))) {
    model <- list(model)
  }

  # Check if model contains at least one model object
  if (length(model) == 0) {
    stop("model must contain at least one model object.")
  }

  # Tidy terms for each model object and combine results
  terms_df <- purrr::map_dfr(model, function(mod) {
    # Check if mod is a valid model object
    if (!inherits(mod, c("lm", "glm"))) {
      stop("model must be a valid 'lm' or 'glm' object.")
    }
    if (!inherits(mod, "robust_glm")) {
      stop("model must be of class robust_glm.")
    }

    # Check if mod contains ci and sandwich, and they are not NA
    ci_df <- if (!is.null(mod$ci) && !all(is.na(mod$ci))) {
      mod$ci %>%
        as.data.frame() %>%
        tibble::rownames_to_column("term")
    } else {
      data.frame(term = names(coef(mod)), conf.low = NA, conf.high = NA)
    }

    sandwich_df <- if (!is.null(mod$sandwich) && !all(is.na(mod$sandwich))) {
      mod$sandwich %>%
        as.data.frame() %>%
        dplyr::select(-estimate) %>%
        dplyr::rename_with(~ paste0(.x, "_hc")) %>%
        tibble::rownames_to_column("term")
    } else {
      data.frame(term = names(coef(mod)), std.error_hc = NA, statistic_hc = NA,
                 p.value_hc = NA, conf.low_hc = NA, conf.high_hc = NA)
    }

    # Tidy model terms and combine with ci and sandwich
    df <- purrr::reduce(list(broom::tidy(mod), ci_df, sandwich_df),
                        dplyr::left_join, by = "term")

    # Exponentiate estimate and confidence intervals
    if (exponentiate) {
      df <- df %>%
        dplyr::mutate(dplyr::across(dplyr::matches("conf|estimate"), exp,
                                    .names = "{.col}_exp"))
    }

    return(df)
  }, .id = id)

  return(terms_df)
}

tidy_model <- function(model, id = NULL, ...) {
  #' Tidy Model Summary Statistics
  #'
  #' Tidies summary statistics of model objects, such as convergence status,
  #' dispersion, and performance metrics.
  #'
  #' @param model A list of or a single robust_glm object.
  #' @param id A character string specifying the identifier column name
  #'  (default is "id").
  #' @param ... Arguments passed to \code{broom::tidy()}.
  #'
  #' @return A data frame containing tidied model summary statistics with
  #'  columns:
  #' \describe{
  #'   \item{\code{id}}{Identifier column name.}
  #'   \item{\code{converged}}{Convergence status of the model.}
  #'   \item{\code{dispersion}}{Dispersion value of the model.}
  #'   \item{\code{glance}}{Summary statistics of the model.}
  #'   \item{\code{rmse}}{Root mean square error.}
  #'   \item{\code{mse}}{Mean square error.}
  #'   \item{\code{hosmer_chisq}}{Hosmer goodness-of-fit test statistic.}
  #'   \item{\code{hosmer_df}}{Degrees of freedom for Hosmer goodness-of-fit
  #'    test.}
  #'   \item{\code{hosmer_p.value}}{p-value of Hosmer goodness-of-fit test.}
  #'   \item{\code{R2_Tjur}}{Tjur's R-squared.}
  #'   \item{\code{Log_loss}}{Logarithmic loss.}
  #' }
  #'
  #' @importFrom purrr map_dfr
  #' @importFrom broom glance
  #' @importFrom performance performance_rmse performance_mse performance_hosmer
  #' @importFrom performance r2_tjur performance_logloss
  #' @importFrom tibble as_tibble
  #' @importFrom dplyr bind_cols
  #' @export

  # If lm or glm supplied, convert to a list
  if (inherits(model, c("lm", "glm"))) {
    model <- list(model)
  }

  # Check if model contains at least one model object
  if (length(model) == 0) {
    stop("Model must contain at least one model object.")
  }

  # Tidy summary statistics for each model object and combine results
  mods_df <- purrr::map_dfr(model, function(mod) {
    # Check if mod is a valid model object
    if (!inherits(mod, "lm") && !inherits(mod, "glm")) {
      stop("Model must be a valid 'lm' or 'glm' object.")
    }
    if (!inherits(mod, "robust_glm")) {
      stop("model must be of class robust_glm.")
    }

    # Extract summary statistics from the model object
    summary_stats <- list(
      converged = tryCatch(mod$converged, error = function(e) NA),
      dispersion = tryCatch(
        if ("dispersion" %in% names(summary(mod))) summary(mod)$dispersion else NA,
        error = function(e) NA
      ),
      rmse = tryCatch(performance::performance_rmse(mod),
                      error = function(e) NA),
      mse = tryCatch(performance::performance_mse(mod),
                     error = function(e) NA),
      hosmer_chisq = tryCatch(performance::performance_hosmer(mod)$chisq,
                              error = function(e) NA),
      hosmer_df = tryCatch(performance::performance_hosmer(mod)$df,
                           error = function(e) NA),
      hosmer_p.value = tryCatch(performance::performance_hosmer(mod)$p.value,
                                error = function(e) NA),
      R2_Tjur = tryCatch(performance::r2_tjur(mod)[[1]],
                         error = function(e) NA),
      Log_loss = tryCatch(performance::performance_logloss(mod)[[1]],
                          error = function(e) NA)
    )
    glance_df <- tryCatch(broom::glance(mod, ...), error = function(e) NA)
    if (!is.null(glance_df)) summary_stats <- c(as.list(glance_df),
                                                summary_stats)

    # Combine summary statistics into a data frame
    tibble::as_tibble(summary_stats)
  }, .id = id)

  return(mods_df)
}

fcs_to_longmat <- function(fcs) {
  library(flowCore)
  # Check flowSet
  if (class(fcs) != "flowSet") {
    stop("This is not a flowSet object.")
  }
  # flatten matrix
  mat <- matrix(
    flowCore::fsApply(fcs, flowCore::exprs),
    byrow = FALSE,
    ncol = length(flowCore::colnames(fcs)),
    dimnames = list(NULL, flowCore::colnames(fcs))
  )
  # Get metadata
  meta <- flowCore::pData(fcs)
  meta$ncells <- as.numeric(fsApply(fcs, nrow))
  # Select variables to expand
  c <- setdiff(names(meta), "ncells")
  # repeat every value x of column c ncells times
  rd <- data.frame(lapply(meta[c], function(x) {
    v <- rep(x, meta$ncells)
  }), row.names = NULL)
  rd$name <- as.factor(rd$name)

  # Check rownames
  if (nrow(mat) != nrow(rd)) {
    stop("Matrix and metadata dimen")
  }
  return(list(mat = mat, rd = rd))
}

mat_to_flowset <- function(mat, markers, id, pheno=NULL){
  # Load libraries
  library(flowCore)
  library(purrr)
  library(dplyr)
  library(tibble)

  # Extract sample ID
  ids <- as.character(unique(mat[[id]]))

  # Select cells from each sample and return only the markers
  flowframes <- purrr::map(ids, ~ flowCore::flowFrame(exprs = as.matrix(mat[mat[[id]] == .x, colnames(mat) %in% markers]))) %>%
    setNames(ids)

  # Create flowset
  if (is.null(pheno)) {
    flowset <- flowCore::flowSet(flowframes)
  } else if (all(pheno %in% colnames(mat))) {
    # Construct df of ids x pheno
    suppressMessages(
      phenodf <- mat %>%
        dplyr::group_by(id) %>%
        dplyr::slice(1L) %>%
        dplyr::select(dplyr::all_of(covs)) %>%
        dplyr::ungroup() %>%
        tibble::column_to_rownames(id)
    )
    flowset <- flowCore::flowSet(flowframes)
    flowCore::pData(flowset) <- phenodf
  } else{
    stop("Not all pheno are columns of mat")
  }
  return(flowset)
}

logit2prob <- function(logit){
  # Source: https://sebastiansauer.github.io/Rcode/logit2prob.R
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

list2tsv <- function(data, out_file, chunk_size = 2000, overwrite = F,
                     compress = T, idcol = "id", sep = "\t", verbose = F) {
  #' Save a List to a Tab-Delimited File
  #'
  #' This function saves a list of data frames to a tab-delimited file (TSV) in
  #' chunks, with options to compress the output file, overwrite existing files,
  #' and specify additional formatting details.
  #'
  #' @param data A list of data frames to be saved. Each element of the list must
  #'   be a data frame or an object coercible to a data frame.
  #' @param out_file A character string specifying the path of the output file
  #'   (without the `.gz` extension if compression is enabled).
  #' @param chunk_size An integer specifying the number of list elements to write
  #'   in each chunk. Default is 2000.
  #' @param overwrite A logical value indicating whether to overwrite an existing
  #'   output file. Default is `FALSE`.
  #' @param compress A logical value indicating whether to compress the output
  #'   file using `pigz`. Default is `TRUE`.
  #' @param idcol A character string specifying the column name to use for
  #'   identifying the source list element in the output. Default is `"id"`.
  #' @param sep A character string specifying the delimiter for the output file.
  #'   Default is `"\t"`.
  #' @param verbose A logical value indicating whether to enable verbose mode for
  #'   `data.table::fwrite`. Default is `FALSE`.
  #'
  #' @return Returns `NULL` invisibly. The function writes the output file to disk.
  #'
  #' @details
  #' The function splits the input list into chunks of size `chunk_size` and writes
  #' each chunk sequentially to the specified file. Optionally, the resulting file
  #' is compressed using the `pigz` command-line tool for efficient storage.
  #'
  #' If the output file already exists and `overwrite` is set to `FALSE`, the
  #' function will stop with an error. If `overwrite` is `TRUE`, the existing file
  #' will be removed.
  #'
  #' Compression with `pigz` requires that the `pigz` tool is installed and
  #' available on the system path.
  #'
  #' @examples
  #' # Example usage:
  #' my_list <- list(
  #'   data.frame(a = 1:5, b = 6:10),
  #'   data.frame(a = 11:15, b = 16:20)
  #' )
  #' list2tsv(
  #'   data = my_list,
  #'   out_file = "output.tsv",
  #'   chunk_size = 1,
  #'   overwrite = TRUE,
  #'   compress = TRUE
  #' )
  #'
  #' @importFrom data.table rbindlist fwrite
  #' @importFrom tictoc tic toc
  #' @export list2tsv

  # Ensure data is a list
  if (!is.list(data)) stop("Data must be of class list.")

  # Overwrite existing file
  file_to_remove <- if (file.exists(paste0(out_file, ".gz"))) {
    paste0(out_file, ".gz")
  } else if (file.exists(out_file)) {
    out_file
  } else {
    NULL
  }
  if (!is.null(file_to_remove)) {
    if (!overwrite) stop("Output file already exists. Set overwrite = T.")
    if (!file.remove(file_to_remove)) {
      stop("Failed to overwrite existing output file.")
    }
  }

  # Ensure output directory exists
  output_dir <- dirname(out_file)
  if (!dir.exists(output_dir)) stop("Output directory does not exist.")

  # Split data into chunks
  chunks <- split(1:length(data), ceiling(seq_along(data) / chunk_size))
  message("Saving list in ", length(chunks), " chunks to ", out_file, "...")

  # Write each chunk to the output file
  for (i in seq_along(chunks)) {
    indices <- chunks[[i]]
    tryCatch({
      tictoc::tic(paste("Chunk", i, "written successfully"))
      # chunk_df <- collapse::unlist2d(data[indices], idcols = idcol, DT = T)
      chunk_df <- data.table::rbindlist(data[indices], idcol = idcol, fill = T)
      data.table::fwrite(chunk_df, out_file, sep = sep, append = T,
                         col.names = (i == 1), verbose = verbose)
      tictoc::toc()
    }, error = function(e) {
      message("Error in chunk ", i, ": ", e$message)
    })
  }

  # Compress the output file with pigz if requested
  if (compress) {
    compressed_file <- paste0(out_file, ".gz")
    message("Beginning compression as ", compressed_file)
    try({
      tictoc::tic(paste("Done compressing", compressed_file))
      system(paste("pigz -9", out_file))
      tictoc::toc()
    })
  }
  return(invisible(NULL))
}

#EOF