library(ade4)
library(spdep)
library(sp)
library(geoR)
library(gstat)
library(data.table)
library(adespatial)
require(dplyr)

setwd("path")


culture_list = c("corn1", "corn2", "soy")
for (culture in culture_list) {
  for (year in 2010:2019) {
    if (culture == 'corn2') {
      arq = paste(culture, '_', year, '.csv', sep = "")
      arq_save = substr(arq, 1, 10)
    } else arq = paste(culture, '_', year, '_', year + 1, '.csv', sep = "")
    
    dataset <- read.csv(arq, sep = ';', header = TRUE, encoding = 'UTF-8')
    dataset = dataset %>% relocate (colnames(dataset[6]), .after=last_col())
    
    bound <- read.csv(paste(path, "bound_5880.csv"), sep = ';', header = TRUE)
    bound <- bound[,2:3]
    bound$x <- bound$x/1000
    bound$y <- bound$y/1000
    num_atrib <- ncol(dataset) - 5
    
    coord <- coordinates(dataset[,3:4])
    max_dist <- max(dist(coord))
    
    gri <- dnearneigh(coord, 0, max_dist/2)
    lw2 <- nb2listw(gri, style = "W")
    c.moran <- lapply(dataset[,6:(num_atrib+5)], moran.mc, lw2, 999)   
    
    i=moran.randtest(dataset[,6:(num_atrib+5)], lw2, nrepet = 999)
    normalized <- scale(dataset[,6:(num_atrib+5)])
    pca2 <- dudi.pca(normalized, center=T, scannf = FALSE, nf = num_atrib)
    
    mc.pca <- lapply(pca2$li[, pca2$eig > 1], moran.mc, lw2, 999)
    ms2=multispati(pca2, lw2, scannf = F, nfposi = num_atrib)
    
    sum.ms <- summary(ms2)	
    autovet=ms2$c1[, ms2$eig > 1]
    
    if (culture == 'corn1') { arq_save = substr(arq, 1, 15) } 
    if (culture == 'soy') { arq_save = substr(arq, 1, 13) }
    
    write.table(autovet, paste('autovet_', arq_save, ".txt", sep = ""))
    
    mc.mpca <- lapply(ms2$li[, ms2$eig > 1], moran.mc, lw2, 999)
    CP_finaly <- pca2$li[, pca2$eig > 1]
    
    result_pca <- data.table(cbind(dataset[,1:5], CP_finaly)) 
    CP_finaly <- ms2$li[, ms2$eig > 1]
    
    result_mpca <- data.table(cbind(dataset[,1:5], CP_finaly))
    
    write.table(result_pca, paste('pca_scores_', arq_save, ".txt", sep = ""))
    write.table(result_mpca, paste('pca_scores_', arq_save, ".txt", sep = ""))
    
    output <- matrix(ncol=3, nrow=length(mc.mpca))
    for (cs in (1:length(mc.mpca))) {
      output[cs,] = cbind(paste("CS", cs, sep = ""), mc.mpca$CS1$method,
                          mc.mpca[[paste("CS", cs, sep = "")]]$p.value)
    }
    
    output <- data.frame(output)
    colnames(output) = c('CPE', 'Method', 'Monte-Carlo')
    
    write.table(output, paste('monte_carlo_', arq_save, ".txt", sep = ""))
  }
}
