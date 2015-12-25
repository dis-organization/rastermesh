setClass("GeolocationCurvilinear",
         representation(
           x = "Raster",
           y = "Raster"
         ))
#' @importFrom nabor WKNNF
setClass ('RasterMesh',
          contains = 'RasterBrick',
          representation (
           geolocation = "GeolocationCurvilinear",
           knnQuery = "WKNNF"
          )
)



