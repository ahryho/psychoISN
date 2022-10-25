# Function to load libraries

LoadPackages <- function(pkg.list){
  new.pkg <- pkg.list[!(pkg.list %in% installed.packages()[, "Package"])]
  if (length(new.pkg)){
    tryCatch(
      install.packages(new.pkg, dependencies = TRUE),
      
      error = function(cond){
        tryCatch(
          BiocManager::install(new.pkg, dependencies = TRUE),
          
          error = function(cond){
            message(cond)
          }
        )
      }
    )
  }
  suppressMessages(sapply(pkg.list, library, character.only = TRUE))
}