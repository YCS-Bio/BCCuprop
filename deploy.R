library(rsconnect)
rsconnect::setAccountInfo(name='cuproptosis-breastcancer',
                          token='64E3DCB3B25A481B2B9A640D66DC947A',
                          secret='Xk2Lw147pzqlaQ+SCeoN5npagq/vgmUBY5QPjAz8')
library(rsconnect)
rsconnect::deployApp(
  appDir = "D:/cuproptosisdeploy", 
  appPrimaryDoc = "appcup.R", 
  appName = "appcup",
  contentCategory = "application"
)
