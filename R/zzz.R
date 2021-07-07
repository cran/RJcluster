# # .onLoad <- function(libname, pkgname) 
# # {
# #   library.dynam("RJcluster", pkgname, libname)
# # }
# 
RJclusterStartupMessage <- function()
{
  # Startup message obtained as
  # > figlet -f slant RJcluster
  msg <- c(paste("version", utils::packageVersion("RJcluster")), "\nType 'citation(\"RJcluster\")' for citing this R package in publications")
  return(msg)
}

.onAttach <- function(lib, pkg)
{
  # unlock .RJcluster variable allowing its modification
  # unlockBinding(".RJcluster", asNamespace("RJcluster"))
  # startup message
  msg <- RJclusterStartupMessage()
  if (!interactive())
    msg[1] <- paste("Package 'RJcluster' version", utils::packageVersion("RJcluster"))
  packageStartupMessage(msg)
  invisible()
}

# .onLoad = function(libname, pkgname)
# {
#   msg = RJclusterStartupMessage()
#   packageStartupMessage(msg)      
#   invisible()
# }

