r <- rastermesh("data/mer_his_1992_01.nc")

library(raadtools)
library(dplyr)
sst <- sstfiles()

r <- NetCDF(sst$fullname[1])
r <- NetCDF("data/mer_his_1992_01.nc")
var0 <- r$var$name
## dimensions and name of the this variable
for (i in seq_along(var0))
  print(r$var  %>% filter(name == var0[i])  %>% inner_join(r$vardim)     %>% select(dimids)  %>% inner_join(r$dims, c("dimids" = "id")))

