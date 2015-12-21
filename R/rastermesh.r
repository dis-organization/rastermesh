rastermesh <- function(x = "data/mer_his_1992_01.nc", varname = "u") {
  d <- brick(x, varname = varname, lvar = 4)
  gl <- new(".GeolocationCurvilinear", x = raster(x, varname = "lon_u"), y = raster(x, varname = "lat_u"))
  new("RasterMesh", d, geolocation = gl, knnQuery = WKNNF(cbind(values(gl@x), range(values(gl@y)))))
}



setMethod("extent", "RasterMesh",
          function(x) {
            warning("bounding box extent from geocation values")
            extent(coordinates(x))
          })
setMethod("coordinates", "RasterMesh",
          function(obj, ...) {
            cbind(values(obj@geolocation@x), values(obj@geolocation@y))
          }
          )


setMethod("extract", signature(x='RasterMesh', y='matrix'),
          function(x, y, ...){
            x@knnQuery$query(y)

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



