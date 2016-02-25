#' @importFrom ncdf4 nc_open nc_close
.varnames <- function(x) {
  names(.ndims(x))
}
.ndims <- function(x) {
  nc <- nc_open(x)
  dims <- sapply(nc$var, "[[", "ndims")
  nc_close(nc)
  dims
}

.dimnames <- function(x, varname) {
  nc <- nc_open(x)
  names(nc$dim[nc$var[[varname]]$dimids])
}

ncatts <- function(x) {
  on.exit(nc_close(ncf))
  ncf <- nc_open(x)
  global <- as_data_frame(ncatt_get(ncf, 0))
  var <- lapply(names(ncf$var), function(vname) as_data_frame(ncatt_get(ncf, vname)))
  names(var) <- names(ncf$var)
  list(global = global, var = var)
}


#' Information about a NetCDF file, in convenient form.
#'
#' @param x path to NetCDF file
#' @export
#' @importFrom ncdf4 nc_open
#' @importFrom dplyr as_data_frame bind_rows data_frame 
NetCDF <- function(x) {
  nc <- ncdf4::nc_open(x)
  dims <- do.call(dplyr::bind_rows, lapply(nc$dim, function(x) dplyr::as_data_frame(x[!names(x) %in% c("dimvarid", "vals", "units", "calendar")])))
  unlimdims <- NULL
  if (any(dims$unlim)) unlimdims <- do.call(dplyr::bind_rows, lapply( nc$dim[dims$unlim], function(x) as_data_frame(x[names(x) %in% c("id", "units", "calendar")])))
  ## do we care that some dims are degenerate 1D?
  ##lapply(nc$dim, function(x) dim(x$vals))
  dimvals <- do.call(bind_rows, lapply(nc$dim, function(x) data_frame(id = rep(x$id, length(x$vals)), vals = x$vals)))
  ## the dimids are in the dims table above
  groups <- do.call(bind_rows, lapply(nc$groups, function(x) as_data_frame(x[!names(x) %in% "dimid"]))) #as_data_frame[x[!names(x) %in% "dimid"]]))
  ## leave the fqgn2Rindex for now
  file <- as_data_frame(nc[!names(nc) %in% c("dim", "var", "groups", "fqgn2Rindex")])
  var <- do.call(bind_rows, lapply(nc$var, function(x) as_data_frame(x[!names(x) %in% c("id", "dims", "dim", "varsize", "size", "dimids")])))
  var$id <- sapply(nc$var, function(x) x$id$id)
  vardim <- do.call(bind_rows, lapply(nc$var, function(x) data_frame(id = rep(x$id$id, length(x$dimids)), dimids = x$dimids)))
  ## read attributes, should be made optional (?) to avoid long read time
  atts <- ncatts(x)
  class(atts) <- c("NetCDF_attributes", "list")
  nc_close(nc)
  x <- list(dims = dims, unlimdims = unlimdims, dimvals = dimvals, groups = groups, file = file, var = var, vardim = vardim, atts = atts)
  class(x) <- c("NetCDF", "list")
  x
}

longlistformat <- function(x, n = 8) {
  if (length(x) <= n) return(x)
  paste(paste(head(x, n), collapse = ", "),  "...",  length(x) - n, "more ...")
}
#' @export
print.NetCDF_attributes <- function(x, ...) {
  print("NetCDF attributes:")
  print("Global")
  print("\n")
  print(x$global)
  print("\n")
  print("Variable attributes:")
  print(sprintf("variable attributes: %s", longlistformat(names(x$var))))
}
#' Return the names of variables in the file
names.NetCDF <- function(x) {
  x$vars$name
}


"[[.NetCDF" <- function(x,i,j,...,drop=TRUE) {
  var <-  x$vars %>% filter(name == i)
  class(var) <- c("NetCDFVariable", class(var))
  var
}

print.NetCDFVariable <- function(x, ...) {
  print(t(as.matrix(x)))
}

#library(lazyeval)
"[.NetCDFVariable" <- function(x, i, j, ..., drop = TRUE) {
  # il <- lazy(i)
  # jl <- lazy(j)
  # dl <- lazy(...)
  #  print(dl)
  #  print( format(dl$expr))
  dots <- list(...)
  #  print(dots)
  ## this is ok, but also need array[i] type indexing, as well as array[matrix]
  if (missing(i)) stop("argument i must be provided")
  
  if (missing(j) & x$ndims > 1L) stop("argument j must be provided")
  #browser()
  nindex <- length(dots) + as.integer(!missing(i)) + as.integer(!missing(j))
  #print(nindex)
  if (!nindex == x$ndims) stop(sprintf("number of index elements must match dimensions of variable: %i", x$ndims))
  #print(i)
  ## now the hard work, see nchelper
  args <- c(list(i), if (missing(j)) list() else list(j), dots)
  # largs <- format(il$expr)
  #return(largs)
  # print(format(il$expr))
  #if (!missing(j)) largs <- sprintf("%s,%s", largs, format(jl$expr))
  
  #if (!missing(...)) sprintf(largs, format(dl$expr))
  # print('after')
  args
  # sprintf("%s[")
}


# nc <- NetCDF("data/mer_his_1992_01.nc")
# Cs_w <- nc[["Cs_w"]]
# lon_u <- nc[["lon_u"]]
# Cs_w[2]
# lon_u[2,3]
#
#