setClass(".GeolocationCurvilinear",
         representation(
           x = "Raster",
           y = "Raster"
         ))
setClass ('RasterMesh',
          contains = 'RasterBrick',
          representation (
           geolocation = .GeolocationCurvilinear
          )
)


