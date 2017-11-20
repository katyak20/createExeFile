BUMPERpp <- function(inputFilename, sampleSite) { 
# 08/20/2016. 
#  Written by Grace M. Hwang.  
#  Computes likelihood function for each species found at each sample site
#  This function allows you to visualize files produced by BUMPER. 
#  Likelihood functions are used to derive the posteriors, from the posteriors 
#  predictions can be derived which can be used to generate model metrics. 
#  Example usage: 
#    llhood50 <- BUMPERpp("ts_likelihood_app.data", 50)
#   Args:
#     fileName: Any one of these files produced from BUMPER.f  
#               ts_likelihood_app.data   
#               ts_likelihood_jack.data  
#               core_likelihood.data               
#               The first row is the environmental variable,
#               followed by likelihood functions. 
#               The first two columns identify sample and location.
#               The third and fourth column identifies the
#               species associated with each likelihood function.
#     sampleSite: Refers to location in which samples were gathered
#      
#   Returns:
#     llhood: Matrix with likelihoods for each species at specified sample site 
#
# Mydir should point to data files produced by Bumper and your speciesFile
  
  Mydir <- "."
  setwd(Mydir) 
  fileName <- file.path(Mydir,inputFilename)
  input0 <- read.table(fileName, quote=NULL, comment.char="#", as.is=TRUE)
   
  xlabel <- as.numeric(input0[1,5:dim(input0)[2]])
  siteVec <- unique(input0[2:dim(input0)[1],1])
  species <- input0[2:dim(input0)[1],4]

  input0 <- input0[2:dim(input0)[1],]
  input <- input0[, 5:dim(input0)[2]]

  rowProd <- function(dat) apply(dat, 2, prod)
  llhood <- matrix(data=0, nrow=length(siteVec), ncol=dim(input)[2])

# Compute likelihood for each site by multiplying all species-specific 
# likelihood functions together
# Normalize based on sum of likelihood of species at each site or sample number.

    i=sampleSite 
      dat <- input[input0$V1==siteVec[i],]
      good <- complete.cases(dat)
      # Define colors for each species found at each site
      col = rainbow(length(good))
      # Compute product of likelihoods normalized by sum
      llhood[i,] <- rowProd(dat[good,])  / sum(rowProd(dat[good,]))
      par(mfrow = c(1, 1))
      plot(xlabel,llhood[i,], main=paste("Sample Number:", i, "-- Filename:", inputFilename, sep =" "), col = "black", 
           type ="l",xlab="Reconstructed variable", ylab = "Normalized Likelihood", 
           cex.lab = 1.25)
        for (j in seq(1,length(good))) {
                lines(xlabel,dat[j,], lwd = 2, col = col[j])  
          }
          lines(xlabel,llhood[i,], col="black", lwd = 4)
          # Extract species index for each site
            speciesName <- input0[input0$V1==siteVec[i], 4]
              speciesName[j+1]={"POSTERIOR"}
              len=length(speciesName)
              # Create Posterior in black
                col[len]="black"
      # Draw legend
        legend("topleft", legend = speciesName, col = col[1:len],
             title = "Species", pch = .1, cex=.8, ncol=1, pt.cex = .8) 

return(llhood[sampleSite,])
}
 
