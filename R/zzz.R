.First.lib <- function(lib, pkg)
{
   cat("##\n## Markov chain Monte Carlo Package (MCMCpack)\n")
   cat("## Copyright (C) 2003 Andrew D. Martin and Kevin M. Quinn\n##\n")
   require(coda)
   require(MASS)
   library.dynam("MCMCpack", pkg, lib)
   
}

