

#load("C:/Users/a-pierov/Dropbox/Visconti_Joppa/bestpoint/maxnet_0.1/maxnet/data/bradypus.rda")
setwd("C:/Users/a-pierov/Dropbox/Visconti_Joppa/bestpoint/data/invasives")
read.csv()
library(glmnet)
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,":")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('C:/Users/a-pierov/Dropbox/Visconti_Joppa/bestpoint/maxnet_0.1/maxnet/R')
sourceDir('C:/Users/a-pierov/Dropbox/Visconti_Joppa/bestpoint/scripts')

p <- bradypus$presence
data <- bradypus[,-1]
mod <- maxnet(p, data)
plot(mod, type="cloglog")
mod <- maxnet(p, data, maxnet.formula(p, data, classes="lq"))
plot(mod, "tmp6190_ann")