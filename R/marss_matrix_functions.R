#' Build a MARSS loading matrix with DFA identification constraints
#'
#' Returns an n x m list matrix suitable for use in a MARSS model list (e.g. as
#' a block of Z). The upper triangle is fixed at zero; the lower triangle and
#' the diagonal are free parameters. This gives exactly m(m-1)/2 zero
#' restrictions, which together with Q = I for the factors supplies the m^2
#' constraints needed to identify an m-factor model -- no more, no less.
#'
#' The diagonal is deliberately left FREE. Fixing it to 1 as well would add m
#' further restrictions and force the first m sites' factor-component variances
#' to equal 1 / (1 - rho^2), which in an AR(1) factor model drags the
#' persistence parameter toward whatever those particular sites' variances
#' happen to be.
#'
#' @param n Number of observation rows (e.g. sites).
#' @param m Number of factors. Must satisfy m <= n.
#' @param prefix Character scalar used to name the estimated parameters, e.g.
#'   "gamma" or "lambda". Allows several independent DFA blocks in one model.
#' @param sep Separator between prefix and indices. Default "_".
#' @param fixed_diag If TRUE, fixes the diagonal at 1 (the over-constrained
#'   parameterisation). Default FALSE. Provided only for comparison.
#' @param zero Value placed in the upper triangle. Numeric 0 by default, which
#'   is what MARSS expects for a fixed element.
#'
#' @return An n x m matrix of mode "list": numeric entries are fixed, character
#'   entries are estimated parameter names.
#'
#' @examples
#' G <- dfa_loadings(6, 2, "gamma")
#' L <- dfa_loadings(6, 3, "lambda")
#' print_loadings(G)
dfa_loadings <- function(n, m, prefix, sep = "_",
                         fixed_diag = FALSE, zero = 0) {
  
  # ---- validation -----------------------------------------------------
  if (!is.numeric(n) || length(n) != 1L || n < 1 || n != as.integer(n))
    stop("`n` must be a single positive integer.", call. = FALSE)
  if (!is.numeric(m) || length(m) != 1L || m < 1 || m != as.integer(m))
    stop("`m` must be a single positive integer.", call. = FALSE)
  if (m > n)
    stop("`m` (", m, ") cannot exceed `n` (", n, "): an m-factor model needs ",
         "at least m rows to place the triangular constraints.", call. = FALSE)
  if (!is.character(prefix) || length(prefix) != 1L || !nzchar(prefix))
    stop("`prefix` must be a single non-empty character string.", call. = FALSE)
  if (m == n)
    warning("m == n: the loading matrix is square, so the factors are simply ",
            "a rotation of the sites and nothing is reduced.", call. = FALSE)
  
  n <- as.integer(n); m <- as.integer(m)
  
  # ---- build ----------------------------------------------------------
  M <- matrix(list(zero), nrow = n, ncol = m)
  
  for (i in seq_len(n)) {
    for (j in seq_len(min(i, m))) {
      if (fixed_diag && i == j) {
        M[[i, j]] <- 1
      } else {
        M[[i, j]] <- paste(prefix, i, j, sep = sep)
      }
    }
  }
  
  dimnames(M) <- list(paste0("row", seq_len(n)), paste0("f", seq_len(m)))
  M
}


#' Preview a list matrix as plain characters
#'
#' MARSS list matrices print awkwardly. This renders one as a character matrix
#' so the structure can be checked at a glance.
#'
#' @param M A list matrix, e.g. from dfa_loadings().
#' @param quote Passed to print(). Default FALSE.
print_loadings <- function(M, quote = FALSE) {
  out <- matrix(vapply(M, function(x) as.character(x[[1]]), character(1)),
                nrow = nrow(M), dimnames = dimnames(M))
  print(out, quote = quote)
  invisible(out)
}


#' Count free parameters in a list matrix
#'
#' @param M A list matrix.
#' @return Integer count of distinct estimated parameter names.
n_free <- function(M) {
  length(unique(Filter(is.character, as.vector(M))))
}


#' Rearrange rows named "param_i_j" into a (lower) triangular matrix
#'
#' @param x       A matrix or data.frame whose row names look like "param_i_j",
#'                or a named vector.
#' @param value   Which column of `x` holds the values (name or index).
#'                Ignored when `x` is a vector.
#' @param param   Which parameter to extract, when the names contain more than
#'                one prefix. If NULL and several are present, a named list of
#'                matrices is returned.
#' @param lower   TRUE (default) puts the larger of (i, j) on the row, so you
#'                get a lower triangular result whether your names are stored
#'                with i <= j or i >= j. FALSE gives upper triangular. Set
#'                `fold = FALSE` if you want the indices used verbatim.
#' @param fold    If FALSE, values go at [i, j] exactly as named and `lower` is
#'                ignored. Use this when both (i,j) and (j,i) can appear.
#' @param fill    Value for the empty triangle. Default 0; NA is useful for
#'                spotting gaps.
#' @param nrows, ncols  Force the output dimensions. Default: inferred from the
#'                largest index seen, so rectangular blocks (e.g. a DFA Z that
#'                is n x m with m < n) come out the right shape.
#' @param pattern Regex with three capture groups: param, i, j.
#'
#' @return A matrix, or a named list of matrices if several params are present.
tri_from_names <- function(x,
                           value   = 1,
                           param   = NULL,
                           lower   = TRUE,
                           fold    = TRUE,
                           fill    = 0,
                           nrows   = NULL,
                           ncols   = NULL,
                           pattern = "^(.+)_([0-9]+)_([0-9]+)$") {
  
  ## --- pull out values + names ----------------------------------------
  if (is.null(dim(x))) {
    nms <- names(x)
    v   <- unname(x)
  } else {
    nms <- rownames(x)
    v   <- x[, value]
  }
  if (is.null(nms)) stop("`x` has no names / row names to parse.")
  
  ## --- parse the names ------------------------------------------------
  m  <- regmatches(nms, regexec(pattern, nms))
  ok <- lengths(m) == 4L
  if (!any(ok)) {
    stop("No names matched `pattern`. First few names: ",
         paste(utils::head(nms, 3), collapse = ", "))
  }
  if (any(!ok)) {
    warning(sum(!ok), " name(s) did not match and were dropped: ",
            paste(utils::head(nms[!ok], 3), collapse = ", "))
  }
  
  d <- data.frame(
    param = vapply(m[ok], `[`, character(1), 2L),
    i     = as.integer(vapply(m[ok], `[`, character(1), 3L)),
    j     = as.integer(vapply(m[ok], `[`, character(1), 4L)),
    value = as.numeric(v[ok]),
    stringsAsFactors = FALSE
  )
  
  ## --- one matrix, or one per parameter --------------------------------
  if (!is.null(param)) {
    d <- d[d$param == param, , drop = FALSE]
    if (nrow(d) == 0L) stop("No rows found for param '", param, "'.")
  } else if (length(unique(d$param)) > 1L) {
    ps  <- unique(d$param)
    out <- lapply(ps, function(p)
      .build_tri(d[d$param == p, , drop = FALSE], lower, fold, fill, nrows, ncols))
    names(out) <- ps
    return(out)
  }
  
  .build_tri(d, lower, fold, fill, nrows, ncols)
}

.build_tri <- function(d, lower, fold, fill, nr, nc) {
  if (fold) {
    r <- if (lower) pmax(d$i, d$j) else pmin(d$i, d$j)
    k <- if (lower) pmin(d$i, d$j) else pmax(d$i, d$j)
  } else {
    r <- d$i
    k <- d$j
  }
  
  if (anyDuplicated(cbind(r, k))) {
    stop("Two entries map to the same cell. Your names likely cover a full ",
         "matrix rather than one triangle -- call with fold = FALSE.")
  }
  
  if (is.null(nr)) nr <- max(r)
  if (is.null(nc)) nc <- max(k)
  if (max(r) > nr || max(k) > nc) stop("An index exceeds the requested dimensions.")
  
  out <- matrix(fill, nr, nc)
  out[cbind(r, k)] <- d$value
  out
}