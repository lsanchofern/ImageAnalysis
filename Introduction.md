# Introduction to image analysis in Python 

## Basic workflow
1. Import necessary packages and set working and destination directories
2. If analyzing all images from folder, then read in each image in a loop
3. Run through image analysis (whatever that may be)
4. Save analysis- binary images and any relevant tables
5. Conduct some statistical tests or plotting of the whole dataset

## Basic workflow for colocalization of puncta or two different markers 
1. Import necessary packages and set working and destination directories
2. If analyzing all images from folder, then read in each image in a loop
3. Threshold first channel to make binary image
5. Threshold second channel to make binary image
6. Calculate overlap of binary images
7. Make into objects by connected components labeling
8. Call regionpropstable to get measures of the overlap (size, intensity, etc.)
