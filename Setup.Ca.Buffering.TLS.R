pkg.list <- c("minpack.lm", "mosaic", "doParallel", "foreach")
success <- TRUE
for(pkg.name in pkg.list){
    install.packages(pkg.name)
    success <- success && require(package = pkg.name, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
}
if(success){
print("Packages were successfully installed and loaded.")  
}else{
  for(pkg.name in pkg.list){
    warning(paste("Could not install or load package:",pkg.name))
  }  
}

