#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to add central airways to thoracic cavity mask  

Inputs:
    -real-valued reconstructed Xe gas img (.nii)
    -thoracic cavity mask
Outouts:
    -new mask with thoracic cavity and central airways (.nii)
    -.png montage comparing several thoracic mask slices to new augmented mask slices

Aryil Bechtel 2021
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#read in .nii file and convert to numpy array
def nii2Array(inFile):
    
    #read in .nii
    myNii = nib.load(inFile)
    
    #convert .nii to np array
    myArray = np.asarray(myNii.dataobj)
    
    return(myArray)

#find xth percentile intensity in thoracic img
def histPercentile(img,mask,perc):
    
    #get img vals within mask
    mskImg = img[np.where(mask==1.0)]
    thresh = np.percentile(mskImg,perc)
    
    return(thresh)
    
#create augmented mask array 
def findAugMask(img,mask,thresh):
    
    #augMask: augmented section + thoracic mask
    augMask = np.zeros(img.shape)
    augMask[np.where(img>thresh)] = 1.0 #include voxels>thresh
    augMask[np.where(mask==1.0)] = 1.0 #also include voxels in og mask
    
    #airMask: augmented section only (sp to be the "airways")
    airMask = np.zeros(img.shape)
    airMask[np.where(img>thresh)] = 1.0
    airMask[np.where(mask==1.0)] = 0 #exclude voxels in og mask
    
    return(augMask,airMask)

#create matrix of only image voxels in augmented mask
def findAugVoxels(img,mask,thresh):
    
    augVoxels = img[np.where(img>=thresh)] #exclude voxels<thresh
    
    return(augVoxels)
    
#create 2 panel fig of a slice of orig mask and aug mask
def createShortSegMontage(img,mask1,mask2,slice):
    
    img = img[:,:,slice-1]
    mask1 = mask1[:,:,slice-1]
    mask2 = mask2[:,:,slice-1]
    opacity = 0.3
    
    fig = plt.figure()
    fig.subplots_adjust(wspace=0)
    p1 = plt.subplot(1,2,1)
    p1.axes.xaxis.set_visible(False)
    p1.axes.yaxis.set_visible(False)
    p1.text(5, 8, "slice "+str(slice), fontsize = 10, \
        fontweight='bold',color='white')
    p1.text(8, 120, "thoracic segmentation", fontsize = 10, \
        fontweight='bold',color='white')
    p1.imshow(img, cmap='gray')
    p1.imshow(mask1, cmap='jet', alpha=opacity)
    
    p2 = plt.subplot(1,2,2)
    p2.axes.xaxis.set_visible(False)
    p2.axes.yaxis.set_visible(False)
    p2.text(8, 120, "augmented segmentation", fontsize = 10, \
               fontweight='bold', color='white')
    p2.imshow(img, cmap='gray')
    p2.imshow(mask2, cmap='jet', alpha=opacity)
    
    return(fig)
    
#create 2xnumSlice panel fig of numSlice slices of orig mask and aug mask
def createLongSegMontage(sbjct,imgIn,mask1In,mask2In):
    
    sliceArray = np.array([4,30,50,52,63,65,98,124]) #array of slices to display
    numSlice = len(sliceArray) #number of slices to show
    opacity = 0.5 #opacity of mask overlay
    
    #loop through slices and plot them
    nrow = 2
    ncol = numSlice
    fig = plt.figure(figsize=(ncol+1, nrow+1)) 
    gs = gridspec.GridSpec(nrow, ncol,\
                           wspace=0.0, hspace=0.0, 
                           top=1.-0.5/(nrow+1), bottom=0.5/(nrow+1), 
                           left=0.5/(ncol+1), right=1-0.5/(ncol+1))
    
    for i in range(numSlice):
        
        slice = sliceArray[i]
        img = imgIn[:,:,slice-1]
        mask1 = mask1In[:,:,slice-1]
        mask2 = mask2In[:,:,slice-1]
        
        p1 = plt.subplot(gs[0,i])
        p1.axes.xaxis.set_visible(False)
        p1.axes.yaxis.set_visible(False)
        p1.text(8, 15, str(slice), fontsize = 7, \
            fontweight='bold',color='white')
        p1.imshow(img, cmap='gray')
        p1.imshow(mask1, cmap='jet', alpha=opacity)
        
        p2 = plt.subplot(gs[1,i])
        p2.axes.xaxis.set_visible(False)
        p2.axes.yaxis.set_visible(False)
        p2.imshow(img, cmap='gray')
        p2.imshow(mask2, cmap='jet', alpha=opacity)
        
        if i==0:
            prctDelta = prctChange(mask1In,mask2In)
            prctDeltaStr = "{:.1f}".format(prctDelta)
            p1.text(8, 120, "thoracic", fontsize = 5, \
                    fontweight='bold',color='white')
            p2.text(8, 120, "augmented, inc " + prctDeltaStr + "%",\
                    fontsize = 5, fontweight='bold', color='white')
            
        if i==numSlice-1:
            p2.text(8, 120, "patient " + sbjct, \
                    fontsize = 5, fontweight='bold', color='white')
            
    return(fig)
    
#calculate percent change from mask1 to mask2
def prctChange(mask1,mask2):
    
    numEl1 = np.count_nonzero(mask1 == 1.0)
    numEl2 = np.count_nonzero(mask2 == 1.0)
    diff = numEl2 - numEl1
    prct = (diff/numEl1) * 100
    
    return(prct)

#subject + file variables
patientNum = "002158"
scanType = "highBW"
maskType = "_uncor"; #if auto leave blank, if manual write '_manual
patientPath = "../data/"

thresholds = np.array([70]) #np.array([30,40,60,70])

#loop over thresholds, make aug mask, create output files + montage image
for i in thresholds:

    #set threshold intensity percentage
    pThresh = i; #in percent
    
    
    fullPath = patientPath + patientNum + "/" + patientNum + "_" + scanType\
            + "/Gas_Exchange/"
    imgPath = fullPath + "GAS_" + patientNum +"_" + scanType + ".nii"
    maskPath = fullPath + "mask_" + patientNum + "_" + scanType + maskType + ".nii"
    
    outFile1 = fullPath + "mask" + maskType + "_augmented_" + patientNum\
            + "_" + scanType + ".nii"
    outFile2 = fullPath + "seg_montage" + maskType + "_augmented_" + patientNum\
            + "_" + scanType + ".png"
    
    
    fullPath = patientPath + patientNum + "_" + scanType + "/"
    imgPath = fullPath + "gx_vent_uncor.nii"
    maskPath = fullPath + "ute_mask.nii"
    
    #outFile1: .nii of augmented (airways + thoracic) and airway (augmented section only)
    outFile1_1 = fullPath + "mask_augmented_" + \
            str(pThresh) + "prct_" + patientNum+ "_" + scanType + ".nii"
    outFile1_2 = fullPath + "mask_airways_" + \
            str(pThresh) + "prct_" + patientNum+ "_" + scanType + ".nii"
            
    outFile2 = fullPath + "seg_montage" + maskType + "_augmented_" + \
            str(pThresh) + "prct_" + patientNum+ "_" + scanType + ".png"
    outFile3 = fullPath + "hists" + maskType + "_" + \
            str(pThresh) + "prct_" + patientNum+ "_" + scanType + ".png"
            
    #get image and mask arrays from .nii files
    img = nii2Array(imgPath)
    mask = nii2Array(maskPath)
    
    #get threshold intensity above P
    thresh = histPercentile(img,mask,pThresh)
    
    #create augmented mask with pixels>thresh
    augMask, airMask = findAugMask(img,mask,thresh)
    
    #plot hist of voxels in each mask and combined mask
    augVoxels = findAugVoxels(img,mask,thresh).flatten();
    ogMaskedVoxels = img[np.where(mask==1.0)].flatten();
    combinedVoxels = np.concatenate((augVoxels,ogMaskedVoxels),axis=None)
    binwidth = 1e-11
    numBins = round((max(combinedVoxels) - min(combinedVoxels))/binwidth) + 1
    binArray = np.linspace(min(combinedVoxels), max(combinedVoxels), numBins)
    plt.figure()
    #plt.xlim((0,1e-9))
    combinedHist = plt.hist(combinedVoxels,bins = binArray,histtype=u'step', color = 'blue',\
                            linewidth = 3, label = 'combined')
    augHist = plt.hist(augVoxels,bins = binArray,histtype=u'step', color = 'green',\
                       linewidth = 3, label = 'augmented')
    ogMaskedHist = plt.hist(ogMaskedVoxels,bins = binArray,histtype=u'step', color = 'red', \
             linewidth = 3, linestyle = 'dotted', label = 'thoracic')
    plt.xlabel('unscaled intensity')
    plt.ylabel('counts')
    plt.legend()
    plt.savefig(outFile3, dpi = 300)
    
    
    #save augmented and airway masks as .nii
    myNii = nib.Nifti1Image(augMask, affine = np.eye(4))
    nib.save(myNii, outFile1_1)
    myNii = nib.Nifti1Image(airMask, affine = np.eye(4))
    nib.save(myNii, outFile1_2)
    
    
    
    slice = 65
    montageFig = createShortSegMontage(img,mask,augMask,slice);
    montageFig = createLongSegMontage(patientNum, img, mask, augMask)
    #montageFig.savefig(outFile2, dpi=300, bbox_inches='tight', pad_inches = 0)