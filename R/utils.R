##### Formula stuff #####
is.formula <- function(x) inherits(x, "formula")

merge.formula <- function(x, y, ...){
  if(!is.formula(x) || length(x) != 3)
    stop("First argument is invalid")
  if(!is.formula(y)) stop("Second argument is invalid")
  if(length(list(...))) warning("extraneous arguments discarded")
  is.gEnv <- function(e) identical(e, .GlobalEnv)

  str <- paste(c(deparse(x[[2]]), "~",
                 deparse(x[[3]]), "+",
                 deparse(y[[length(y)]])), collapse = "")
  f <- stats::as.formula(str)
  ## MM: try to keep environment (where reasonable)
  ex <- environment(x)
  ey <- environment(y)
  if(!is.gEnv(ex)) {
    environment(f) <- ex
    if(!is.gEnv(ey) && !identical(ex,ey)) {
      warning("`x' and `y' have different environments; x's is used")
    }
  } else if(!is.gEnv(ey))
    environment(f) <- ey
  f
}

## Data helper functions
list2numeric <- function(x){
  as.numeric(as.matrix(x))
}


## parallel functions
plapply <- function (X, FUN, ..., .parallel = 1, .seed = NULL, .verbose = TRUE) {
  if (!(useCluster <- inherits(.parallel, "cluster"))) {
    stopifnot(length(.parallel) == 1L, is.vector(.parallel,
                                                 "numeric"), .parallel >= 1)
    .parallel <- as.vector(.parallel, mode = "integer")
    if (.Platform$OS.type == "windows" && .parallel > 1L) {
      useCluster <- TRUE
      .parallel <- parallel::makeCluster(.parallel)
      on.exit(parallel::stopCluster(.parallel))
    }
  }
  FUN <- match.fun(FUN)
  .FUN <- if (useCluster || is.primitive(FUN)) {
    FUN
  }
  else {
    verboseExpr <- if (isTRUE(.verbose)) {
      if (.parallel == 1L && interactive()) {
        env <- new.env(hash = FALSE, parent = environment(FUN))
        environment(FUN) <- env
        env$pb <- txtProgressBar(min = 0, max = length(X),
                                 initial = 0, style = 3)
        on.exit(close(env$pb), add = TRUE)
        quote(setTxtProgressBar(pb, pb$getVal() + 1L))
      }
      else {
        on.exit(cat("\n"), add = TRUE)
        quote(cat("."))
      }
    }
    else if (is.call(.verbose) || is.expression(.verbose)) {
      .verbose
    }
    else if (is.character(.verbose)) {
      on.exit(cat("\n"), add = TRUE)
      substitute(cat(.verbose))
    }
    do.call(add.on.exit, list(FUN, verboseExpr))
  }
  if (!is.null(.seed)) {
    if (useCluster) {
      parallel::clusterSetRNGStream(cl = .parallel, iseed = .seed)
    }
    else {
      if (!exists(".Random.seed", envir = .GlobalEnv,
                  inherits = FALSE)) {
        set.seed(NULL)
      }
      .orig.seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", .orig.seed, envir = .GlobalEnv),
              add = TRUE)
      if (.parallel == 1L) {
        set.seed(seed = .seed)
      }
      else {
        stopifnot(requireNamespace("parallel", quietly = TRUE))
        set.seed(seed = .seed, kind = "L'Ecuyer-CMRG")
        parallel::mc.reset.stream()
      }
    }
  }
  if (useCluster) {
    parallel::parLapply(cl = .parallel, X = X, fun = .FUN,
                        ...)
  }
  else if (.parallel == 1L) {
    lapply(X = X, FUN = .FUN, ...)
  }
  else {
    parallel::mclapply(X = X, FUN = .FUN, ..., mc.preschedule = TRUE,
                       mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = .parallel)
  }
}

add.on.exit <- function (FUN, expr){
  FUN <- match.fun(FUN)
  if (is.null(expr <- substitute(expr))) {
    return(FUN)
  }
  if (is.primitive(FUN)) {
    stop("not implemented for primitive functions")
  }
  onexitexpr <- substitute(on.exit(expr))
  obody <- body(FUN)
  body(FUN) <- if (is.call(obody) && identical(as.name("{"),
                                               obody[[1L]])) {
    as.call(append(x = as.list(obody), values = onexitexpr,
                   after = 1L))
  }
  else {
    as.call(c(as.name("{"), onexitexpr, obody))
  }
  FUN
}


