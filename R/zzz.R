.onAttach <- function(...) {
   cat("##\n## Markov Chain Monte Carlo Package (MCMCpack)\n")
   cat("## Copyright (C) 2003, 2004, Andrew D. Martin and Kevin M. Quinn\n")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0350646 and SES-0350613)\n##\n")
   require(coda, quietly=TRUE)
   require(MASS, quietly=TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("MCMCpack", libpath)
}
