# write a matrix out to an ASCII file that can be read by Scythe
# this puts the number of rows and columns in the first row
# followed by the data
#
# 1/29/2003 (ADM)

write.Scythe <- function(outmatrix, outfile = NA, overwrite=FALSE) {
   outmatrix <- as.matrix(outmatrix)
   
   if(is.na(outfile)) {
      stop("Please specify a file name in the write.Scythe() call.\n")
   }
   if(overwrite==FALSE & file.exists(outfile)) {
      cat("File already exists in the write.Scythe() call.\n")
      stop("Either delete the file, or flip the overwrite switch.\n")
      
   }
   
   outfile <- file(outfile, "w")
   cat(dim(outmatrix), "\n", file=outfile)
   write.table(outmatrix, file=outfile,
      row.names=FALSE, col.names=FALSE, quote=FALSE)
   close(outfile)
   0
} 

