.First.lib <- function(lib, pkg)
{
   require(coda)
   require(ts)
   require(MASS)
   library.dynam("MCMCpack", pkg, lib)

   invisible()
}
