.First.lib <- function(lib, pkg)
{
   require(coda)
   require(MASS)
   library.dynam("MCMCpack", pkg, lib)
   
}
