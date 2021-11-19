#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to plot bias field corrected hists of scans using
 combined, augmented, and thoracic masks
 
Inputs:
    -3 bias field corrected vent niftis:
        -combined
        -augmented
        -thoracic
Output:
    -plot of 3 histograms of 3 vent images

Aryil Bechtel 2021
"""
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt

#read in .nii file and convert to numpy array
def nii2Array(inFile):
    
    #read in .nii
    myNii = nib.load(inFile)
    
    #convert .nii to np array
    myArray = np.asarray(myNii.dataobj)
    
    return(myArray)

#subject + file variables
patientPath = "../data/"
patient = "ven_Sub007-004B_highBW" + ".nii"

#load .nii files as matrices
corrImg = nii2Array(patientPath + patient)
plt.imshow(corrImg[:,:,65], cmap='gray')

#create and plot hists

#save figure
