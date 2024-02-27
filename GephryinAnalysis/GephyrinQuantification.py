# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 10:47:16 2024

@author: sanchofernandez
"""
#%% 
#import some necessary packages

from tifffile import imsave 
import os 
import time
import torch
import czifile
import pandas as pd
import numpy as np 
import pyclesperanto_prototype as cle
from pyclesperanto_prototype import imshow
from skimage import io, data, segmentation
from skimage.transform import rescale 
from skimage.filters import gaussian
from skimage.io import imread
from aicsimageio import AICSImage
import cv2
import zarr
import dask.array as da
from numcodecs import Blosc

#%%
#Use the GPU 
devices=cle.available_device_names()
print(cle.available_device_names())
#set file directory- BPHO server
os.chdir("J:\\Laura\\Gephyrin_imaging\\to_analyze")
os.listdir()
rootdir=os.getcwd()
print (os.getcwd())

#%%files directories 

files_list=[]
for root, subdirs, files in os.walk(rootdir):
    for subdir in subdirs:
        files_list.append(os.path.join(root, subdir))

files=[]
for file_name in files_list:
    print(file_name)


#%%
import tifffile as tf 
from aicsimageio.writers import OmeTiffWriter
from dask_image.ndmeasure import label as daskimage_label 
from skimage.measure import regionprops_table
import dask.array as da  
from skimage.measure import label
from skimage.util import map_array

savedir="J:\\Laura\\Gephryin_imaging\\analysis"
cle.set_wait_for_kernel_finish(True)
substring='Fill' #fill image of traced dendrite 
substring2='Airyscan'

length=len(files_list)    
for i in range (1,length):
    for filename in os.listdir(files_list[i]):
        if substring in filename: 
            print('dendrite')
            cle.select_device(devices[0]) #select gpu 
            f=imread(os.path.join(files_list[i],filename))

            #filter image 
            filtered=gaussian (f,sigma=1.5)
            
            #no need for background subtraction for this one 
            #segmentation
            cle.push (filtered)
            binary_dendrite=cle.threshold_otsu(filtered)
            cle.pull(filtered)
            cle.pull(binary_dendrite)
            del f
            del filtered 
            #binary fill holes
            #save binary 
            tf.imwrite(os.path.join(savedir,str('binary_')+filename),binary_dendrite)
            
            #clear gpu memory 
            torch.cuda.empty_cache()

        if substring2 in filename: 
            cle.select_device(devices[1]) #select gpu 
            print ('gephryin')
            g=AICSImage(os.path.join(files_list[i], filename)) #if czifile 
            g_czyx=g.get_image_data("ZYX", C=1, S=0,T=0) #get gephryin channel

            binary_puncta=cle.threshold_otsu(g_czyx)
            labeled_puncta=label(binary_puncta)
            scale=(0.0441,0.0441,.110) #set scale of pixels (x,y,z)
            props=['label','area']
            data=regionprops_table(
                labeled_puncta, 
                spacing=scale,
                properties=props
                )
            data_df=pd.DataFrame(data)
            area_labels=data_df['label']*(data_df['area']<1)
            new_labeled_puncta=map_array(labeled_puncta, np.asarray(data_df['label']),np.asarray(area_labels))
            
            del area_labels
            del data_df
            del binary_puncta
            del labeled_puncta 
            
        binary_dendrite=np.asarray(binary_dendrite)
        binary_puncta=np.asarray(binary_puncta)
        overlap_binary=cv2.bitwise_and(binary_dendrite,binary_puncta)
        
        overlap_binary_labels=label(overlap_binary)
        scale=(0.0441,0.0441,.110) #set scale of pixels (x,y,z)
        props=['label','area','mean_intensity','centroid','equivalent_diameter']
        data=regionprops_table(
           new_labeled_puncta,
           overlap_binary,
           spacing=scale,
           properties=props
           )

    data_df=pd.DataFrame.from_dict(data)
    del data
    mask=data_df['mean_intensity']>0
    data_filtered=data_df[mask]
    cle.imshow(overlap_binary)