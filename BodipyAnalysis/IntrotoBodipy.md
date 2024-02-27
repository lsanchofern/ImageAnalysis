# Bodipy analysis

These are astrocyte iPSCs labeled with DAPI and s100b and bodipy dye. 
This code segments out the astrocytes in each image and quantifies bodipy. 

## Basic pipeline: 

**Draw astrocyte masks:**
- Use cellpose previously trained model to segment astrocytes
- save masks for future plotting and double-checking
- make as connected components
- regionproptable to get coordinates of astrocytes 

**Segment DAPI channel:**
- Segment DAPI channel using voronoi_otsu method from pyclesperanto-prototype
- will only count astrocytes that contain DAPI
- make as connected components
- regionproptable to get coordinates of DAPI nuclei 

**Segment BODIPY channel:**
- Filter and threshold BODIPY channel
- make as objects using connected components labeling
- regionproptable to get coordinates of BODIPY puncta 

**Pandas implementation and dictionary comprehension to select for BODIPY inside DAPI+ astrocytes:**
- Use the magic of pandas to extract coordinates of DAPI, astrocytes, and BODIPY
- Select astrocytes that have DAPI+ nuclei inside
- Then count BODIPY puncta within each astrocyte
- save measures for each image in a destination folder 

**More pandas to compile data and plot**
- Read in all tables for all images into one dataframe to work with
- plot and run statistics 
