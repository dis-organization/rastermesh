setClass("GeolocationCurvilinear",
         representation(
           x = "Raster",
           y = "Raster"
         ))
#' @importClassesFrom nabor WKNNF
#' @importFrom nabor WKNNF
setClass ('RasterMesh',
          contains = 'RasterBrick',
          representation (
           geolocation = "GeolocationCurvilinear",
           knnQuery = "WKNNF"
          )
)



