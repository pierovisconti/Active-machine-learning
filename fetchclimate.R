fetchclimate<-function(centroids,variable){
  
# example temperature=fetchclimate(cbind(c(14,14,14),c(40,43,44)),"airt")  fetches mean annual air temperature at lat 40,43,44 and long 14. for mean monthly prec. the variable is prate

## JSON request object. 
## The text may be downloaded from http://fetchclimate2.cloudapp.net/form
## ('Download request' button)
## Alternatively, one may use R list constructor 
## following with a call to 'toJSON' function
  
library(RCurl)
library(rjson)  

body = paste('{"EnvironmentVariableName": "',variable,'", "ParticularDataSources": ["GFDLAvTemp","GHCNv2"],\
"Domain": {\
"Lats": [',paste(centroids[,2],collapse=','),'], \
"Lons": [',paste(centroids[,1],collapse=','),'], "TimeRegion": \
{"Years": [2000,2012],\
"Days": [1,365],\
"Hours": [0,24], "IsIntervalsGridYears": true, "IsIntervalsGridDays": true, "IsIntervalsGridHours": true}, "SpatialRegionType": "Points"}}',sep='')


## perform request
h = basicTextGatherer()
h$reset()
curlPerform(url = "http://fetchclimate2.cloudapp.net/api/compute",
            httpheader=c(Accept="text/plain", 'Content-Type' = "application/json"),
            postfields=body,
            writefunction = h$update)
reply = h$value()
print(reply)

## wait while processing in progress
while (substr(reply,1,7)=='pending' || substr(reply,1,8)=='progress') {
  Sys.sleep(1)
  hash=strsplit(reply,"hash=")[[1]][2]
  h$reset()
  curlPerform(url = paste("http://fetchclimate2.cloudapp.net/api/status?hash=",curlEscape(hash),sep=""),
              httpheader=c(Accept="text/plain"),
              writefunction = h$update)
  reply = h$value()
  print(reply)
}

## get result data
if (substr(reply,1,9)=='completed') {
  msds = substr(reply,11,nchar(reply))
  h$reset()
  curlPerform(url = paste("http://fetchclimate2.cloudapp.net/jsproxy/data?uri=",curlEscape(msds),"&variables=values",sep=""),
              httpheader=c(Accept="application/json"),
              writefunction = h$update)
  result=fromJSON(h$value())
}
output=as.matrix(lapply(result$values,function(x) x[1]))
return(output)
}