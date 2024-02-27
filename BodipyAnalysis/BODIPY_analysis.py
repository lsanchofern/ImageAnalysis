# -*- coding: utf-8 -*-

#Created on Fri Jan 12 11:24:00 2024

#@author: sanchofernandez

#%% This is code written for analysis of BODIPY puncta within iPSC astrocytes 
#The code thresholds and binarizes the astrocyte using an astrocyte marker like
# s100b or GFAP and counts the number of BODIPY puncta within this binary mask 

#%% some set up 
#import some necessary packages

from tifffile import imsave 
import os 
import time
import torch
import czifile
import pandas as pd
import numpy as np 
import pyclesperanto_prototype as cle
from skimage import io, data, segmentation
from skimage.transform import rescale 
from skimage.filters import gaussian, threshold_otsu, threshold_niblack, threshold_sauvola
import cv2 
from aicsimageio import AICSImage

import matplotlib.pyplot as plt
import matplotlib as mpl
#%% some set-up 
#Use the GPU 
devices=cle.available_device_names()

#set file directory- 
os.chdir("J:\\Laura\\Jillytest\\toanalyze\\images\\take3") #change this to be your directory 
os.listdir()
rootdir=os.getcwd()
print (os.getcwd()) # if you would like to check the working directory


#%%some other set-up 
# import some other stuff 
    
import tifffile as tf 
from aicsimageio.writers import OmeTiffWriter
from scipy.ndimage import binary_fill_holes
from skimage.measure import regionprops_table, regionprops 

#set a save directory 
savedir_images="J:\\Laura\\Jillytest\\analysis\\images\\take3" #change this to be where analyzed images live 
savedir_tables="J:\\Laura\\Jillytest\\analysis\\tables\\take3" #change this to be where analysis tables live 
cle.set_wait_for_kernel_finish(True) #for dask processing if wanted 

## SOME IMPORTANT METADATA: 
#voxel_size_x=0.161 #microns 
#voxel_size_y=0.1611 #microns

#channel 0: DAPI
#channel 1: BODIPY
#channel 2: Astrocyte marker     

#%% Cellpose segmentation of astrocytes 
from cellpose import io
from cellpose import models, core, utils, metrics  
from cellpose.io import imread 
import tifffile 
from cellpose import plot


#TO TRAIN MODEL: 
    #lab computer: Anaconda/ Miniconda prompt, type "conda activate napari-env"
    # type "cellpose"
    #make sure all images used for training purposes are in one folder and are tiffs
    #you can run a model from model zoo to see if it can segment cells more or less and then manually draw ROIs and add
    #make sure seg.npy file is saved 
    #draw masks for all files in folder 
    #train new model with files/ masks 
    #make sure you give it a normal name you can identify and know where it's saved 

#cellpose parameters
#location of model
model= models.CellposeModel(pretrained_model="J:\\Laura\\Jillytest\\toanalyze\\training_data\\models\\Jilly_iPSCs_astros")
channels= [[3,1]]

# define CHANNELS to run segementation on
# grayscale=0, R=1, G=2, B=3
# channels = [cytoplasm, nucleus]
# if NUCLEUS channel does not exist, set the second channel to 0
# cytoplasm/ astrocyte= channel 3
# nucleus/ dapi = channel 1 

# if diameter is set to None, the size of the cells is estimated on a per image basis
# you can set the average cell `diameter` in pixels yourself (recommended)
# diameter can be a list or a single number for all images


#%%  RUN ANALYSIS FOR ALL IMAGES 
from skimage.io import imread, imshow

#For all the files in the directory  
for filename in os.listdir(rootdir):
    #select the GPU 
    cle.select_device(devices[1])
    
    #read in the image (can be czi file)      
    image=AICSImage(os.path.join(rootdir, filename))
    
    #separate by channel 
    print(filename)
    image_dapi=image.get_image_data('ZYX', C=0) #get DAPI channel
    image_bodipy=image.get_image_data('ZYX', C=1) #get BODIPY channel
    image_astro=image.get_image_data('ZYX', C=2) #get astrocyte channel   
    
    #images are arranged in Z, Y, then X in 1st, 2nd, and 3rd dimensions of arrays
    #Y is the number of rows, X is number of columns 
    #Z will be 1 because it's a z-projection 
    
    #make astrocyte masks
    masks, flows, diams = model.eval(image_astro, diameter=None, channels=channels)
   
    #make an array 
    masks_array=np.asarray(masks)
    
    #make into objects that can be counted 
    astro_labels_filtered=cle.connected_components_labeling_box(masks_array)
    
    #segment nuclei 
    nuclei_dapi=cle.voronoi_otsu_labeling(image_dapi, spot_sigma=50)
    
    #threshold bodipy channel 
    filtered_bodipy=gaussian(image_bodipy, sigma=1.5) #filter first 
    binary_bodipy=cle.threshold_otsu(filtered_bodipy)
    
    #make into objects 
    labeled_bodipy=cle.connected_components_labeling_box(binary_bodipy) #these are your "objects" 
    
    #filter out big bodipy aka the nucleus or random clumps
    size_min_threshold_bodipy=1 #pixels 
    size_max_threshold_bodipy=200 #pixels 
    
    labeled_bodipy_filtered=cle.exclude_labels_outside_size_range(labeled_bodipy,minimum_size=size_min_threshold_bodipy, maximum_size=size_max_threshold_bodipy)
    
    #save images 
    fname_new=filename.replace('.','_', 1)
    tf.imwrite(os.path.join(savedir_images,str('bodipy_binary_')+fname_new),binary_bodipy)
    
    #save if there are masks, if not then don't 
    #the "try" lets the function run but doesn't error and stop code if it doesn't work 
    try: 
        io.save_masks(image_astro, 
              masks, 
              flows, 
              fname_new, 
              channels=[3, 1],
              tif=True, # save masks as TIFFs
              save_outlines=True,
              savedir= savedir_images, # save flows as TIFFs
              in_folders=True
              )
    except:
        print('no masks')
   
    #set scale of pixels (y,x)
    scale1=(0.161,0.161)
    #bodipy data
    props=['label', 'mean_intensity','area','bbox','centroid','coords'] #properties we want to measure
    
    bodipy_labels=cle.pull(labeled_bodipy_filtered) #convert to numpy array 
    
    #I use np.squeeze function because image arrays are (1, whatever, whatever) dimensions 
    #and the first dimension is unnecessary and throws off the scaling 
    #second dimension is y (so number of rows in array)
    #third dimension is x (so number of columns in array)
    
    data=regionprops_table(
        np.squeeze(bodipy_labels), 
        np.squeeze(image_bodipy),
        spacing=scale1,
        properties=props
        )
    
    data_df=pd.DataFrame.from_dict(data) #convert to dataframe 
    
    #astrocyte data 
    labeled_astro_img=cle.pull(astro_labels_filtered) #convert to numpy array 
    data_astro=regionprops_table(
        labeled_astro_img,
        np.squeeze(image_astro),
        spacing=scale1,
        properties=props
        )
    
    data_astro_df=pd.DataFrame.from_dict(data_astro)
    
    #dapi data 
    dapi_labels=cle.pull(nuclei_dapi) #convert to numpy array 
    props2=['label', 'coords']
    
    data_dapi=regionprops_table(
        np.squeeze(dapi_labels), 
        spacing=scale1,
        properties=props2
        )
    
    data_dapi_df=pd.DataFrame.from_dict(data_dapi)
    
    #count bodipy objects, mean intensity, and area per astrocyte label
    #need to use the coordinates to make sure they are within astrocyte object 
    #for each bodipy puncta, extract x coordinates and y coordinates 
    #for each astrocyte, extract x coordinates and y coordinates 
    
    #----------ASTROCYTES-------------------------------------------------------------------
    #initialize empty list for astroctye coordinates 
    num_astros = len(data_astro_df)
    astro_coords=[]
    
    for i in range (num_astros):
        num_coords=len(data_astro_df.coords[i])
        temp_label=data_astro_df.label[i]
        temp=data_astro_df.coords[i]
        temp_x=[]
        temp_y=[]
        temp_lbl=[]
        for j in range (num_coords):
            temp_x.append(temp[j][1])
            temp_y.append(temp[j][0])
    
        astro_coords.append({'label': temp_label, 'ycoords': temp_y, 'xcoords':temp_x}) #make list of dictionaries 
        
    #convert to dataframe   
    astro_coords_df=pd.DataFrame(astro_coords)      
    
    #----------BODIPY-------------------------------------------------------------------
    #initialize empty list for BODIPY coordinates 
    num_bodipy = len(data_df)
    bodipy_coords=[]
    
    for i in range (num_bodipy):
        num_coords=len(data_df.coords[i])
        temp=data_df.coords[i]
        temp_x=[]
        temp_y=[]
        for j in range (num_coords):
            temp_x.append(temp[j][1])
            temp_y.append(temp[j][0])
        bodipy_coords.append({'ycoords': temp_y, 'xcoords':temp_x}) #make list of dictionaries 
        
    #convert to dataframe   
    bodipy_coords_df=pd.DataFrame(bodipy_coords)      
        
    #-------------------DAPI---------------------------------------------------------------
    num_dapi = len(data_dapi_df)
    dapi_coords=[]
    
    for i in range (num_dapi):
        num_coords=len(data_dapi_df.coords[i])
        temp=data_dapi_df.coords[i]
        temp_x=[]
        temp_y=[]
        for j in range (num_coords):
            temp_x.append(temp[j][1])
            temp_y.append(temp[j][0])
        dapi_coords.append({'ycoords': temp_y, 'xcoords':temp_x}) #make list of dictionaries 
        
    #convert to dataframe   
    dapi_coords_df=pd.DataFrame(dapi_coords)      
    
    
    #--------------ONLY ASTROCYTES CONTAINING DAPI-----------------------------------------
    #select for astrocytes containing DAPI 
    dapi_astro=[]
    #find overlapping puncta with each astrocyte "object"
    for i in range (len(dapi_coords_df)): #for each nuclei 
        y=[]
        x=[]
        y=np.array(dapi_coords_df.ycoords[i])
        x=np.array(dapi_coords_df.xcoords[i])
        array=np.vstack((x,y)).T       
        temp=[]
        for j in range (len(astro_coords_df)):
            y_astro=[]
            x_astro=[]
            y_astro=np.array(astro_coords_df.ycoords[j])
            x_astro=np.array(astro_coords_df.xcoords[j])
            array_astro=np.vstack((x_astro,y_astro)).T
            n=np.concatenate((array,array_astro), axis=0)
            n_unique=np.unique(n, axis=0)
            if len(n_unique)<(len(array)+len(array_astro)):
                temp.append(data_astro_df.label[j])
        dapi_astro.append(temp)         
        
    try: 
        dapi_astro_df=pd.DataFrame(dapi_astro,columns=['label'])         
    
        astros_withdapi=dapi_astro_df['label'].tolist()
        mask = data_astro_df['label'].isin(astros_withdapi)
        filtered_astros=data_astro_df[mask]          #select only astrocytes with dapi 
        filtered_astros.reset_index(drop=True, inplace=True) # reset index 
    
        mask2= astro_coords_df['label'].isin(astros_withdapi)
        filtered_astro_coords_df=astro_coords_df[mask2]
        filtered_astro_coords_df.reset_index(drop=True, inplace=True)
                
    #-------------COUNT BODIPY PUNCTA WITHIN EACH ASTROCYTE--------------------------------
    
        bodipy_avg_per_astro=[]
        #find overlapping puncta with each astrocyte "object"
        for i in range (len(filtered_astro_coords_df)): #for each astrocyte
            y=[]
            x=[]
            y=np.array(filtered_astro_coords_df.ycoords[i])
            x=np.array(filtered_astro_coords_df.xcoords[i])
            array=np.vstack((x,y)).T
            temp_bodipy_intensity=[]
            temp_bodipy_area=[]
            for j in range (len(data_df)):
                y_bodipy=[]
                x_bodipy=[]
                y_bodipy=np.array(bodipy_coords_df.ycoords[j])
                x_bodipy=np.array(bodipy_coords_df.xcoords[j])
                array_bodipy=np.vstack((x_bodipy,y_bodipy)).T
                n=np.concatenate((array_bodipy,array), axis=0)
                n_unique=np.unique(n,axis=0)
                if len(n_unique)<(len(array)+len(array_bodipy)):
                    temp_bodipy_intensity.append(data_df.mean_intensity[j])
                    temp_bodipy_area.append(data_df.area[j])
            bodipy_avg_per_astro.append({'intensity': temp_bodipy_intensity, 'area': temp_bodipy_area})
        
        bodipy_per_astro_df=pd.DataFrame(bodipy_avg_per_astro)         
    except: 
        print('no astros with dapi!')
        
    try: #if there actually is bodipy 
        #store average and all bodipy per astrocyte          
        filtered_astros['mean_bodipy_intensity']=bodipy_per_astro_df['intensity'].apply(np.mean)         
        filtered_astros['mean_bodipy_area']=bodipy_per_astro_df['area'].apply(np.mean)
        filtered_astros['bodipy_count']=bodipy_per_astro_df['area'].apply(len)
        filtered_astros['sum_bodipy_area']=bodipy_per_astro_df['area'].apply(np.sum) 
        filtered_astros['all_bodipy_intensity']=bodipy_per_astro_df['intensity']
        filtered_astros['all_bodipy_area']=bodipy_per_astro_df['area']
    
        #divide total bodipy area by astrocyte area
        area_bodipy=filtered_astros['sum_bodipy_area']/filtered_astros['area']
        filtered_astros['percent_bodipy_area']=area_bodipy *100
        #save dataframe 
        filtered_astros.to_csv(os.path.join(savedir_tables,str(filename)+'_astrocyte_bodipy_stats.csv'))
    except:
        print('No bodipy or dapi+ astrocytes, sorry!')
        
    
    #clean up the GPU space 
    cle.pull(image_dapi)
    cle.pull(image_bodipy)
    cle.pull(nuclei_dapi)
    cle.pull(image_astro)

    cle.pull(astro_labels_filtered)
    cle.pull(binary_bodipy)
    cle.pull(labeled_bodipy)
    torch.cuda.empty_cache() 
    
    del image_dapi
    del image_bodipy 
    del image_astro 
    del nuclei_dapi
    del astro_labels_filtered 
    del binary_bodipy
    del labeled_bodipy

    #clear gpu memory 
    torch.cuda.empty_cache()
    
