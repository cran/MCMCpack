.onAttach <- function(...) {
 
   # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
   cat("##\n## Markov Chain Monte Carlo Package (MCMCpack)\n")
   cat("## Copyright (C) 2003-", this.year,
      " Andrew D. Martin and Kevin M. Quinn\n", sep="")
   cat("##\n## Support provided by the U.S. National Science Foundation\n")
   cat("## (Grants SES-0350646 and SES-0350613)\n##\n")
   require(coda, quietly=TRUE)
   require(MASS, quietly=TRUE)
   require(stats, quietly=TRUE)
}

.onUnload <- function(libpath) {
    library.dynam.unload("MCMCpack", libpath)
}

