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


#' @export
#' @importFrom nabor WKNNF
rastermesh <- function(x, varname) {
  if (!file.exists(x)) return("no file")
  if (missing(varname)) {
    dimnums <- .ndims(x)
    ## first one with the highest dimension
    varname <- names(dimnums)[which.max(dimnums)]
  }
  d <- brick(x, varname = varname, lvar = 4)
  gl <- new("GeolocationCurvilinear", x = raster(x, varname = "lon_u"), y = raster(x, varname = "lat_u"))
  ## why does this give a different answer??
##  print(head(cbind(values(gl@x), range(values(gl@y)))))
  new("RasterMesh", brick(d[[1]]), geolocation = gl, knnQuery = nabor:::WKNNF(coordinates(gl)))
}


boundary <- function(x) {
  left <- cellFromCol(x, 1)
  bottom <- cellFromRow(x, nrow(x))
  right <- rev(cellFromCol(x, ncol(x)))
  top <- rev(cellFromRow(x, 1))
  ## need XYFromCell method
  coordinates(r)[unique(c(left, bottom, right, top)), ]
}



setMethod("extent", "RasterMesh",
          function(x) {
            warning("bounding box extent from geocation values")
            extent(coordinates(x))
          })
setMethod("coordinates", "GeolocationCurvilinear",
          function(obj, ...) {
            cbind(values(obj@x), values(obj@y))
          })
setMethod("coordinates", "RasterMesh",
          function(obj, ...) {
            coordinates(obj@geolocation)
          }
          )

if (!isGeneric("cellFromXY")) {
  setGeneric("cellFromXY", function(object, xy, ...)
    standardGeneric("cellFromXY"))
}
#setMethod("cellFromXY", signature(object='RasterLayer', xy='matrix'), raster::cellFromXY)
#' @export
setMethod("cellFromXY", signature(object='RasterMesh', xy='matrix'),
          function(object, xy, ...) {
            kn <- object@knnQuery$query(xy, k = 1, eps = 0)
            cell <- as.vector(kn$nn.idx)
            bdy <- boundary(object)
            ## TODO test for distance from the edge
            ## could use distance from knn object
            ## need som pre-analysis of the coordinates and their spacing
            outside <- sp::point.in.polygon(xy[,1], xy[, 2], bdy[,1], bdy[,2], mode.checked = TRUE) < 1
            if (any(outside)) {
              #dist <- rep(0, nrow(xy))
              #dist[outside] <- geosphere::dist2Line(xy[outside, , drop = FALSE], bdy)
              #cell[dist > 1e3] <- NA_integer_
              cell[outside] <- NA_integer_
            }
            cell
          }
)

# r <- rastermesh:::rastermesh()
# y <-  cbind(c(146, 147, 145), c(-65, -64, -60))
setMethod("extract", signature(x='RasterMesh', y='matrix'),
          function(x, y, ...){
            cells <- cellFromXY(x, y)
            bad <- is.na(cells)
            vals <- rep(NA_real_, length(cells))
            if (any(!bad)) {
              vals[!bad] <- extract(x, cells[!bad])
            }
            vals
          }
          )


#
# scl <- function(x) (x - na.omit(min(x)))/diff(range(na.omit(x)))
# setMethod("plot", signature(x = 'RasterMesh', y = "ANY"),
#           function(x, y, ...) {
#             #x <- rasterize(coordinates(x), raster(extent(x), nrow = nrow(x) - 30, ncol = ncol(x) - 30), field = values(x[[1]]),
#             #               fun = mean, na.rm = TRUE)
#            # prj <- "+proj=omerc +lonc=147 +lat_0=-65 +alpha=58 +gamma=58"
#           #  plot(project(coordinates(x), prj), col = rainbow(256)[scl(values(x[[1]])) * 255 + 1], pch = 16, cex = 0.3)
#             g <- graticule(seq(xmin(x), xmax(x), length = 55), seq(ymin(x), ymax(x), length = 35))
#             #g1 <- gris(crop(countriesHigh, extent(x)))
#             g1 <- gris(subset(wrld_simpl, NAME == "Antarctica"))
#             g1$v$x1 <- g1$v$x
#             g1$v$y1 <- g1$v$y
#
#
#             qu <- x@knnQuery$query(g1$v %>% dplyr::select(x, y) %>% as.matrix, k = 1, eps = 0)
#             xx <- coordinates(x@geolocation@x)[qu$nn.idx[,1], ]
#             g1$v$x <- xx[,1]
#             g1$v$y <- xx[,2]
#
#
#             #plot(x)
#           })
# plot(r)



