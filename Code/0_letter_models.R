# Function for generating a data set to represent the letters of the alphabet
libs <- list("raster", "jpeg", "mixtools")
lapply(libs,library,character.only = T)
# Generate a library of images for the letters of the alphabet
# Storing images in a separate folder
setwd("Letters/")
# Generate a string specifying upper case, lower case and special characters
all.letters <- c(letters,toupper(letters), "?", "!")
for (i in 1:length(all.letters)){
  # Open the graphics device to save the images as their_name.jpeg; compression isn't an issue
  jpeg(filename = paste(all.letters[i], ".jpg", sep=""))
  # Plot each letter individually
  plot(0,0,pch=all.letters[i],xaxt="n", yaxt="n", ann=F, bty="n", cex=30)
  # Close the graphics device
  dev.off()
}


# Create a list to store our output
out <- list()
# Get a list of the files in the working directory (the letter images)
files <- list.files()
# max.dim is the max number of rows used to specify the image of each letter
max.dim <- 32
# n.mix is the number of mixtures in each character, tuned manually
n.mix <- c(2, 4, 5, 3, 4, 6, 4, 4, 4, 4, 4, 4, 3, 3, 5, 
           5, 3, 3, 2, 2, 3, 3, 3, 3, 2, 2, 5, 4, 3, 3,
           4, 5, 4, 3, 4, 5, 2, 5, 5, 5, 3, 2, 3, 3, 2,
           2, 4, 4, 4, 4, 3, 3, 3, 3)

# Fitting a multivariate Gaussian mixture model to each character
for (i in 1:54) {
  # Read in the image and convert it into a matrix where each entry represents
  # a pixel and gives the colour in that pixel
  img <- readJPEG(files[i])
  img.matrix <- as.raster(img, max = 255)
  seqn <- floor(seq(1, nrow(img.matrix), length=max.dim))
  # Represent the colours as 1 (black) or 0 (white)
  img.matrix[img.matrix=="#000000"] <- 1
  img.matrix[img.matrix=="#010101"] <- 0
  # Convert the matrix to numeric
  img.matrix <- matrix(as.numeric(unlist(img.matrix)),nrow=nrow(img.matrix))
  
  # Compressing the matrix to save on computation time by taking the average of
  # surrounding "pixels" and rounding so we only have 1's and 0's
  compressed <- matrix(0, max.dim, max.dim)
  for(j in 1:(max.dim - 1)){
    for(k in 1:(max.dim - 1)){
      compressed[j,k] <- round(mean(img.matrix[seqn[j]:seqn[j+1], seqn[k]:seqn[k+1]]))
    }
  }

  # Converting the matrix to an (x,y) co-ordinates matrix
  len.compressed <- nrow(compressed)
  coords <- matrix(0, nrow=len.compressed^2, ncol=3)
  for (j in 1:len.compressed) {
    for (k in 1:len.compressed) {
      coords[k+len.compressed * (j-1),] <- c(j, max.dim-k, compressed[j,k])
    }
  }
  # Fit a mixture of Gaussians to the "black pixels" with # mixtures specified
  gauss.mix <- mvnormalmixEM(coords[coords[,3]==1, 1:2], k = n.mix[i])
  
  # Store the name of the letter, the mu, sigma, and mixing proportions in our 
  # output
  out[[i]] <- c(letter = gsub(".jpg","",files[i]),
                lambda = gauss.mix$lambda,
                mu = gauss.mix$mu,
                sigma = gauss.mix$sigma)
}

# After running store the output in the parent directory
dput(out,file="../output")




