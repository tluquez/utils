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

filter_cell_types <- function(data, grouping_vars, min_n, min_pct, covs_n,
                              donor_col = "donor", cell_type_col = "cell_type",
                              num_cells_col = "num_cells") {
  # Filter cell types with at least min_n donors and min_pct of non-zero
  # num_cells per combinations of grouping variables

  # Get all cell types
  cell_types <- unique(data[[cell_type_col]])
  cell_types <- structure(cell_types, names = cell_types)

  # Filter
  cell_types_filtered <- data %>%
    tidyr::drop_na(!!sym(cell_type_col)) %>%
    dplyr::group_by(!!sym(cell_type_col), !!!syms(grouping_vars)) %>%
    dplyr::summarize(
      n_donor = dplyr::n_distinct(!!sym(donor_col)),
      absent = mean(!!sym(num_cells_col) != 0) <= min_pct
    ) %>%
    dplyr::filter(n_donor >= min_n, absent == FALSE) %>%
    dplyr::group_by(!!sym(cell_type_col)) %>%
    dplyr::summarize(
      n_donor = sum(n_donor),
      n_grouping = dplyr::n_distinct(!!!syms(grouping_vars))
    ) %>%
    dplyr::filter(n_donor >= covs_n, n_grouping >= 2) %>%
    dplyr::pull(!!sym(cell_type_col))

  cell_types[cell_types %in% cell_types_filtered]
}

pycat_to_rfactor <- function(x) {
  # Converts categorical variables saved in a python anndata to R vectors.
  # x list. List of at least two arrays.
  # https://github.com/scverse/anndata/blob/a96205c58738b7054d7ed663707c2fec6b52e31d/docs/fileformat-prose.rst#id11

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

get_h5ad_meta <- function(h5ad) {
    # Extracts cell metadata from h5ad.
    # h5ad character. Path to h5ad file.
    suppressPackageStartupMessages({
        library(tools)
        library(rhdf5)
        library(purrr)
    })

    stopifnot(file.exists(h5ad))
    stopifnot(tools::file_ext(h5ad) == "h5ad")

    # Import obs and extract cell IDs
    obs <- rhdf5::h5read(h5ad, "obs", read.attributes = T)
    rhdf5::h5closeAll()
    obs_i <- as.character(obs[[attributes(obs)[["_index"]]]])

    # SEA-AD's format
    if ("__categories" %in% names(obs)) {
        cats <- names(obs)[names(obs) %in% names(obs[["__categories"]])]
        cats <- structure(cats, names = cats)
        cats <- purrr::map(cats, ~ factor(obs[[.x]],
                                          labels = obs[["__categories"]][[.x]]))
        num <- names(obs)[!names(obs) %in% names(obs[["__categories"]]) & names(obs) != "__categories"]
        d <- c(cats, obs[num])
    } else {
        # Scanpy's format
        valid_names <- c("categories", "codes", "mask", "values")
        for (i in which(purrr::map_int(obs, purrr::vec_depth) == 3)) {
            if (!sum(names(obs[[i]]) %in% valid_names) == 2) {
                .y <- paste0(names(obs[i]), "/", names(obs[[i]]))
                obs[i] <- unlist(obs[i], F, F)
                names(obs)[i] <- .y
                message("Illegal column names were replaced")
            }
        }
        # Transform categories (length(.x) == 2) to factors
        d <- purrr::modify_if(obs, .p = ~ length(.x) == 2, .f = ~ pycat_to_rfactor(.x))
    }

    # Enframe obs if all columns have length equal to ncells
    if (purrr::every(d, ~ length(.x) == length(obs_i))) {
        obs <- data.frame(d, check.names = F)[, attributes(obs)[["column-order"]]]
        rownames(obs) <- obs_i
    }
    return(obs)
}

get_h5ad_expr <- function(h5ad, transpose = T, class = "H5SparseMatrix") {
  # Extracts sparse matrix and cell metadata from h5ad.
  # h5ad character. Path to h5ad file.
  # transpose logical. Whether to transpose the cellsxgenes matrix to return genesxcells.
  # class character. Either "sparseMatrix" or "H5SparseMatrix". The latter is useful for large matrices that need to be stored using bit64 (note, spam does not accept dimnames).
  suppressPackageStartupMessages({
    library(tools)
    library(rhdf5)
    library(purrr)
  })

  stopifnot(file.exists(h5ad))
  stopifnot(tools::file_ext(h5ad) == "h5ad")

  # Import row and column names
  obs_attr <- rhdf5::h5readAttributes(h5ad, "obs")
  obs_i <- rhdf5::h5read(h5ad, paste0("obs/", obs_attr[["_index"]]))
  var_attr <- rhdf5::h5readAttributes(h5ad, "var")
  # Keeping support for old versions of AnnData
  if (length(var_attr[["_index"]]) > 1) {
    var_i <- rhdf5::h5read(h5ad, paste0("var/", var_attr[["_index"]]))
  } else {
    var_i <- rhdf5::h5read(h5ad, "var/feature_name")
    var_i <- var_i$categories[var_i$codes + 1]
  }
  rhdf5::h5closeAll()

  # Import matrix
  if (class == "sparseMatrix") {
    suppressPackageStartupMessages(library(Matrix))
    X <- rhdf5::h5read(h5ad, "X",
      read.attributes = T,
      bit64conversion = "bit64"
    )
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

get_h5ad <- function(h5ad, transpose = T, ...) {
  # Extracts sparse matrix and cell metadata from h5ad.
  # h5ad character. Path to h5ad file.
  # transpose logical. Whether to transpose the cellsxgenes matrix to return genesxcells.
  # class character. Either "sparseMatrix" or "spam". The latter is useful for large matrices that need to be stored using bit64.
  meta <- get_h5ad_meta(h5ad)
  expr <- get_h5ad_expr(h5ad, ...)
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
    if(is.null(group))
      group <- counts$samples$group
    # if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 0)
    #   design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size))
      lib.size <- with(counts$samples, lib.size * norm.factors)
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts)))
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts)))
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  if (nrow(counts) < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  # Library size
  if(is.null(lib.size))
    lib.size <- colSums(counts)
  # Group
  if(is.null(group))
    group <- rep("Group1", ncol(counts))
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
  if(print)
    cat("Group:\n")
  E <- w <- counts
  xy <- line <- as.list(rep(NA, ngroups))
  names(xy) <- names(line) <- levgroup
  for (lev in 1L:ngroups) {
    if(print)
      cat(lev, levgroup[lev], "\n")
    i <- intgroup == lev
    countsi <- counts[, i]
    libsizei <- lib.size[i]
    designi <- design[i, , drop = FALSE]
    QR <- qr(designi)
    if(QR$rank<ncol(designi))
      designi <- designi[,QR$pivot[1L:QR$rank], drop = FALSE]
    if(ncol(designi)==ncol(countsi))
      designi <- matrix(1L, ncol(countsi), 1)
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
  #voom overall
  if (TRUE %in% dynamic){
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
  if(plot!="none"){
    disp.group <- c()
    for (lev in levgroup) {
      dge.sub <- edgeR::DGEList(counts[,group == lev])
      disp <- edgeR::estimateCommonDisp(dge.sub)
      disp.group[lev] <- disp$common
    }
    if(plot %in% c("all", "separate")){
      if (fix.y.axis == TRUE) {
        yrange <- sapply(levgroup, function(lev){
          c(min(xy[[lev]]$y), max(xy[[lev]]$y))
        }, simplify = TRUE)
        yrange <- c(min(yrange[1,]) - 0.1, max(yrange[2,]) + 0.1)
      }
      for (lev in 1L:ngroups) {
        if (fix.y.axis == TRUE){
          plot(xy[[lev]], xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25, ylim = yrange)
        } else {
          plot(xy[[lev]], xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25)
        }
        title(paste("voom: Mean-variance trend,", levgroup[lev]))
        lines(line[[lev]], col = "red")
        legend("topleft", bty="n", paste("BCV:", round(sqrt(disp.group[lev]), 3)), text.col="red")
      }
    }


    if(plot %in% c("all", "combine")){
      if(is.null(col.lines))
        col.lines <- 1L:ngroups
      if(length(col.lines)<ngroups)
        col.lines <- rep(col.lines, ngroups)
      xrange <- unlist(lapply(line, `[[`, "x"))
      xrange <- c(min(xrange)-0.3, max(xrange)+0.3)
      yrange <- unlist(lapply(line, `[[`, "y"))
      yrange <- c(min(yrange)-0.1, max(yrange)+0.3)
      plot(1L,1L, type="n", ylim=yrange, xlim=xrange, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )")
      title("voom: Mean-variance trend")
      if (TRUE %in% dynamic){
        for (dy in which(dynamic)){
          line[[dy]] <- line_all
          disp.group[dy] <- disp_all
          levgroup[dy] <- paste0(levgroup[dy]," (all)")
        }

      }
      for (lev in 1L:ngroups)
        lines(line[[lev]], col=col.lines[lev], lwd=2)
      pos.legend <- pos.legend[1]
      disp.order <- order(disp.group, decreasing = TRUE)
      text.legend <- paste(levgroup, ", BCV: ", round(sqrt(disp.group), 3), sep="")
      if(pos.legend %in% c("inside", "outside")){
        if(pos.legend=="outside"){
          plot(1,1, type="n", yaxt="n", xaxt="n", ylab="", xlab="", frame.plot=FALSE)
          legend("topleft", text.col=col.lines[disp.order], text.legend[disp.order], bty="n")
        } else {
          legend("topright", text.col=col.lines[disp.order], text.legend[disp.order], bty="n")
        }
      }
    }
  }
  # Output
  if (TRUE %in% dynamic){
    E[,intgroup %in% which(dynamic)] <- E_all[,intgroup %in% which(dynamic)]
    w[,intgroup %in% which(dynamic)] <- w_all[,intgroup %in% which(dynamic)]
  }
  out$E <- E
  out$weights <- w
  out$design <- design
  if(save.plot){
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

de_pseudobulk <- function(input, meta, model, contrast,
                          de_method = "limma", de_type = "voombygroup",
                          verbose = F, return_interm = F) {
  # input list or matrix. List of cell types with list of expr and meta generated by Libra::to_pseudobulk().rownames(pseudobulks$meta) must match colnames(pseudobulks$expr).
  # It can also be a count matrix of gene x observation  satisfying colnames(pseudobulks) == rownames(meta) and untransformed counts.
  # meta data.frame. Table containing the per observation covariates included in model. If input is a matrix, meta is required. Is input is a list, meta is expected per each entry.
  # model formula. Constructed with as.formula() or ~ and containing terms in colnames(meta) or colnames(pseudobulks$meta).
  # contrast numeric or character. Position or name of the contrast in colnames(design.matrix(model, data = meta)).
  # verbose logical. Level fo information displayed.
  # return_interm logical. If TRUE, returns intermediate files.
  # Inpired by https://github.com/neurorestore/Libra/blob/main/R/pseudobulk_de.R
  # TODO:
  # - Allow for multiple contrasts in limma. In which case an F statistic will be returned by toptalble(). rdrr.io/bioc/limma/man/toptable.html

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
      # Check inputs
      if (!identical(sort(colnames(expr)), sort(rownames(meta)))) {
        stop("Observation names between expression and metadata do not match.")
      }
      covs <- all.vars(model)
      if (!all(covs %in% colnames(meta))) {
        stop(paste0("Model terms absent from meta: ",
                    colnames(meta)[!covs %in% colnames(meta)]))
      }
      if (!all(expr %% 1 == 0)) {
        stop("Counts are not integer.")
      }
      # Ensure design matrix is full rank
      design_full <- model.matrix(model, data = meta)
      stopifnot(identical(rownames(design_full), colnames(expr)))
      design_full <- check_nonestimable(design_full)

      # Check categorical or continuous contrast
      if (!contrast %in% colnames(design_full)) {
        stop(paste0(contrast,
                    " contrast was not present in the full design with terms: ",
                    paste(colnames(design_full), collapse = ", "), "\n"))
      }
      is.categorical <- length(unique(design_full[, contrast])) == 2
      if (!is.categorical & de_type == "voombygroup") {
        message("voombygroup does not support continuous contrasts. Running voom instead.")
        de_type <- "voom"
      }
      # Define reduced design for LRT de_type in DESeq2
      if (is.character(contrast)) {
        design_reduced <- design_full[, colnames(design_full) != contrast]
      } else if (is.numeric(contrast)) {
        design_reduced <- design_full[, -contrast]
      }

      if (verbose) print(data.frame(design_full))

      # Run DE
      de <- switch(
        de_method,
        edgeR = {
          tryCatch(
            {
              y <- edgeR::DGEList(
                expr = expr,
                group = design_full[, contrast]
              ) %>%
                edgeR::calcNormFactors(method = "TMM") %>%
                edgeR::estimateDisp(design_full)
              test <- switch(de_type,
                QLF = {
                  fit <- edgeR::glmQLFit(y, design = design_full)
                  test <- edgeR::glmQLFTest(fit, coef = contrast)
                },
                LRT = {
                  fit <- edgeR::glmFit(y, design = design_full)
                  test <- edgeR::glmLRT(fit, coef = contrast)
                }
              )
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
              list(res = res,
                   interm = list(y = y, test = test, fit = fit))
            },
            error = function(e) {
              message(e)
              list()
            }
          )
        },
        DESeq2 = {
          tryCatch(
            {
              dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(
                countData = expr, colData = meta, design = design_full
              ))
              dds <- switch(de_type,
                Wald = {
                  dds <- try(DESeq2::DESeq(dds,
                    test = "Wald",
                    quiet = !verbose
                  ))
                },
                LRT = {
                  dds <- try(DESeq2::DESeq(dds,
                    test = "LRT", quiet = !verbose,
                    reduced = design_reduced
                  ))
                }
              )
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
              list(res = res, interm = dds)
            },
            error = function(e) {
              message(e)
              list()
            }
          )
        },
        limma = {
          tryCatch(
            {
              # Scale library size
              norm.factors <- edgeR::calcNormFactors(expr)
              lib.size <- norm.factors * colSums(expr)
              # voom or trend
              switch(de_type,
                trend = {
                  trend_bool <- T
                  dge <- edgeR::DGEList(expr, lib.size = lib.size,
                                        norm.factors = norm.factors,
                                        group = design_full[, contrast])
                  v <- methods::new("EList")
                  v$E <- edgeR::cpm(dge, log = T)
                  v
                },
                voom = {
                  trend_bool <- F
                  v <- limma::voom(expr, design_full, save.plot = T,
                                   lib.size = lib.size)
                  voom_plot <- cbind.fill(
                    voom_plot_x = v$voom.xy$x,
                    voom_plot_y = v$voom.xy$y,
                    voom_line_x = v$voom.line$x,
                    voom_line_y = v$voom.line$y
                  ) %>%
                    tibble::rownames_to_column("gene")
                },
                voombygroup = {
                  trend_bool <- F
                  v <- voomByGroup(expr, design_full,save.plot = T,
                                   plot = "all", print = F,
                                   group = design_full[, contrast],
                                   lib.size = lib.size)
                  voom_plot <- cbind.fill(
                    voom_plot_1_x = v$voom.xy[[1]]$x,
                    voom_plot_1_y = v$voom.xy[[1]]$y,
                    voom_plot_2_x = v$voom.xy[[2]]$x,
                    voom_plot_2_y = v$voom.xy[[2]]$y,
                    voom_line_1_x = v$voom.line[[1]]$x,
                    voom_line_1_y = v$voom.line[[1]]$y,
                    voom_line_2_x = v$voom.line[[2]]$x,
                    voom_line_2_y = v$voom.line[[2]]$y
                  ) %>%
                    tibble::rownames_to_column("gene")
                }
              )
              # get fit
              fit <- limma::lmFit(v, design_full) %>%
                limma::eBayes(trend = trend_bool, robust = trend_bool)
              # calculate standard errors for the coefficients
              SE <- sqrt(fit$s2.post) * fit$stdev.unscaled
              # format the results
              res <- fit %>%
                limma::topTable(number = Inf, coef = contrast, confint = T) %>%
                {
                  if (nrow(.) == nrow(SE)) {
                    dplyr::bind_cols(., lfcse = SE[, contrast])
                  } else {
                    .
                  }
                } %>%
                tibble::rownames_to_column("gene") %>%
                dplyr::left_join(voom_plot, by = "gene") %>%
                {
                  if (is.categorical) {
                    dplyr::mutate(
                      .,
                      n0 = sum(design_full[, contrast, drop = T] == 0),
                      n1 = sum(design_full[, contrast, drop = T] == 1)
                    )
                  } else {
                    .
                  }
                } %>%
                dplyr::mutate(
                  bonf = p.adjust(P.Value, "bonferroni"),
                  contrast = contrast
                ) %>%
                dplyr::rename(dplyr::any_of(c(
                  fdr = "adj.P.Val", lfc = "logFC", conf.low = "CI.L", b = "B",
                  conf.high = "CI.R", ave_expr = "AveExpr", p.value = "P.Value",
                  stat = "t"
                ))) %>%
                dplyr::arrange(p.value)
              list(res = res, interm = list(v = v, fit = fit))
            },
            error = function(e) {
              message(e)
              list()
            }
          )
        }
      )
      de[["design_full"]] <- design_full
      de[["design_reduced"]] <- design_reduced

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

#EOF