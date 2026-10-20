miniMap_out <- read_tsv('alignment.sam')

# set the column names
colnames(miniMap_out) <- c('queryName', 'queryLength', 'queryStart', 'queryEnd', 'char', 'targetName', 'targetLength', 'targetStart', 'targetEnd', 'matchingBases', 'matchLength', 'matchQuality')

# check the table
head(miniMap_out)

# First define x-axis coordinates. This will help us later define the genomic window (xlim) we want to zoom in on.
start = 0
end = 3000000

# Plot an empty plot (so we can fill it with rectangles (~genome) and polygons (~alignments)).
# For the y-axis (ylim) coordinates I arbitrarily use c(0,10).
plot(NULL, xlim = c(start, end), ylim = c(0,10), axes=F, ylab = '', xlab = '')

# Now we can draw two simple rectangles, one will define the genomic interval/sequence of H. melpomene (deepskyblue), the other H. erato (mediumseagreen).
rect(0,8,2000000,9, col = 'deepskyblue', border = NA)
rect(0,1,2000000,2, col = 'mediumseagreen', border = NA)

# We also know the positions of the optix gene and we can add these with the same rect() function trick.
rect(705604,9.2,706407,10, col = 'red', border = 'red') # position optix melpomene (only has one exon)

rect(1239943,0,1239972,0.8, col = 'red', border = 'red') # first exon position optix erato
rect(1250591,0,1251211,0.8, col = 'red', border = 'red') # second exon position optix erato
rect(1239943,0.4,1251211,0.4, col = 'red', border = 'red') # A little line between the two and we have a gene model!

# With the text function, we can add the gene and species names at the appropriate coordinates.
text(706407, 9.6, substitute(paste(italic('optix'))), pos = 4)
text(1251211, 0.4, substitute(paste(italic('optix'))), pos = 4)

text(start, 7.5, substitute(paste(italic('H. melpomene'))), pos = 4)
text(start, 2.5, substitute(paste(italic('H. erato'))), pos = 4)

polygon(x = c(miniMap_out$targetStart[1], miniMap_out$targetEnd[1], miniMap_out$queryEnd[1], miniMap_out$queryStart[1]), 
        y = c(8,8,2,2),
        col = adjustcolor('black', alpha.f = 0.1), border = FALSE)
# Assigning [1] behind the vectors extracts only the first element from it. When we plot all the matches, we will loop through these vectors using the variable [e].

for(e in 1:nrow(miniMap_out)){
  
  polygon(x = c(miniMap_out$targetStart[e], miniMap_out$targetEnd[e], miniMap_out$queryEnd[e], miniMap_out$queryStart[e]), 
          y = c(8,8,2,2),
          col = adjustcolor('black', alpha.f = 0.1), border = FALSE)
}


