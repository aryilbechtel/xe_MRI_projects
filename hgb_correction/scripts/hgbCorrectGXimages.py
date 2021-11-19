#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to rescale bar/gas and rbc/gas images (ie for hgb correction)

Workflow:
    -input .mat output file from GX pipeline
    -extract mask, bar2gas, rbc2gas and UTE matrices
    -divide bar2gas and rbc2gas matrices by scaling factor from 
        MOXE hgb dependence on rbc and bar signals
    -define 8 and 6 bins (7 and 5 thresholds) for bar and rbc values
    -sort bar and rbc voxels into color bins
    -combine sorted matrices with UTE matrix
    -display color-coded montage of bar and rbc slices with colormap
    -create histograms of rescaled bar and rbc values
    -subtract uncorrected defect, low, and high % (% for hgb = 15 g/dL) from corrected %
    -plot distribution of delta % values for subject cohort
    
    -save hgb and delta % values in one .txt files
    -export montages as .png and .nii files
    -export hists as .png files
    -export distribution (boxplots) as .png files

Aryil Bechtel 2021
"""

import numpy as np
from scipy.io import loadmat
import matplotlib
import matplotlib.pyplot as plt
from hgbCorrectUtils import makeMontage, decideStartInterval, \
    binning, makeHists, binStats, getDeltaPrct
from hgbCorrectUtils import barmap, long_index2color, mmap, short_index2color,\
    hgbKey
#%%

#create delta prct arrays to store values from each subject
deltaDefectRbc = np.zeros(len(hgbKey))
deltaLowRbc = np.zeros(len(hgbKey))
deltaHighRbc = np.zeros(len(hgbKey))
deltaDefectBar = np.zeros(len(hgbKey))
deltaLowBar = np.zeros(len(hgbKey))
deltaHighBar = np.zeros(len(hgbKey))

#index counter for delta prct arrays
i = 0

for patientName in hgbKey:
    #patientName = "004014_highBW"
    patientPath = "../data/" + patientName + "/"
    matPath = patientPath + patientName + ".mat"
    hgbRef = 15 #healthy reference hgb g/dL (such that hgb correction factor = 1)
    #healthy cohort hgb data in g/dL
    
    hgb = hgbKey[patientName] #get patient hgb from key
    #hgb=6
    
    #load .mat file and extract matrices
    matFile = loadmat(matPath)
    
    mask = matFile['mask_reg']
    bar = np.absolute(matFile['bar2gas'])
    rbc = np.absolute(matFile['rbc2gas'])
    ute = matFile['ute_reg']
    
    #define hgb correction fns and scale bar and rbc matrices
    barPoly3 = [-2.925e-5, 0.0026, -0.1163, 2.2509]
    rbcPoly3 = [2.286e-5, -0.0021, 0.0909, 0.0226]
    
    barScale = barPoly3[0]*hgb**3 + barPoly3[1]*hgb**2 + barPoly3[2]*hgb + barPoly3[3]
    rbcScale = rbcPoly3[0]*hgb**3 + rbcPoly3[1]*hgb**2 + rbcPoly3[2]*hgb + rbcPoly3[3]
    
    barCorrect = np.divide(bar,barScale)
    rbcCorrect = np.divide(rbc,rbcScale)
    
    #define color bin thresholds
    
        #from 'GX_defineColormaps.py' in pipeline ver2
    mean_bar = 0.736
    std_bar = 0.278
    barThresh = [mean_bar-2*std_bar, mean_bar-std_bar, mean_bar, mean_bar+std_bar,
                mean_bar+2*std_bar, mean_bar+3*std_bar, mean_bar+4*std_bar]
    rbcThresh = [0.066, 0.250, 0.453, 0.675, 0.956]
    
    #sort matrices into bins
    barSort = binning(bar,barThresh)
    rbcSort = binning(rbc,rbcThresh)
    barSortCorrect = binning(barCorrect,barThresh)
    rbcSortCorrect = binning(rbcCorrect,rbcThresh)
    
    
    
    #ind_start, ind_inter = decideStartInterval(mask = mask)
    ind_start = 60
    ind_inter = 8
    
    
    barMontageUncorr = makeMontage(barSort, ute, long_index2color, ind_start, \
                             ind_inter, patientName, patientPath, "bar_montage_uncorr_short.png")
    rbcMontageUncorr = makeMontage(rbcSort, ute, short_index2color, ind_start,\
                             ind_inter, patientName, patientPath, "rbc_montage_uncorr_short.png")
    
    barMontageCorrect = makeMontage(barSortCorrect, ute, long_index2color, ind_start, \
                             ind_inter, patientName, patientPath, "bar_montage_hgbCorrect_short.png")
    rbcMontageCorrect = makeMontage(rbcSortCorrect, ute, short_index2color, ind_start,\
                             ind_inter, patientName, patientPath, "rbc_montage_hgbCorrect_short.png")
    
    
    
    barUncorrStats = binStats(bar, barSort, mask, 'bar')
    rbcUncorrStats = binStats(rbc, rbcSort, mask, 'rbc')
    
    barCorrStats = binStats(barCorrect, barSortCorrect, mask, 'bar')
    rbcCorrStats = binStats(rbcCorrect, rbcSortCorrect, mask, 'rbc')
    
    deltaPrctBar = getDeltaPrct("bar", barCorrStats, barUncorrStats)
    deltaPrctRbc = getDeltaPrct("rbc", rbcCorrStats, rbcUncorrStats)
    
    deltaDefectBar[i] = deltaPrctBar["defect"]
    deltaLowBar[i] = deltaPrctBar["low"]
    deltaHighBar[i] = deltaPrctBar["high"]
    
    deltaDefectRbc[i] = deltaPrctRbc["defect"]
    deltaLowRbc[i] = deltaPrctRbc["low"]
    deltaHighRbc[i] = deltaPrctRbc["high"]
    
    i+=1
    
    
    barHistOutFile = "hgbBarHistCompare_croppedAx2.png"
    rbcHistOutFile = "rbcBarHistCompare_croppedAx2.png"
    
    matplotlib.rcParams.update({'font.size': 14})
    
    barHistCorr, barHistUncorr = \
        makeHists(barCorrect[np.where(mask == 1.0)], bar[np.where(mask == 1.0)], "bar", \
                             barCorrStats, barUncorrStats, patientPath + barHistOutFile)
        
    rbcHistCorr, rbcHistUncorr = \
        makeHists(rbcCorrect[np.where(mask == 1.0)], rbc[np.where(mask == 1.0)], "rbc", \
                             rbcCorrStats, rbcUncorrStats, patientPath + rbcHistOutFile)


#%% plot delta% vs hgb for real patient data

i = 0
hgbArray = np.zeros(len(hgbKey))
for patientName in hgbKey:
    hgbArray[i] = hgbKey[patientName]
    i+=1
#create array of hgb and delta% values for export
deltaPrctRbcFile = '../data/deltaPrctRbc'
deltaPrctBarFile = '../data/deltaPrctBar'
outArrayRbc = np.vstack((hgbArray,deltaDefectRbc,deltaLowRbc,deltaHighRbc))
outArrayRbc = np.transpose(outArrayRbc)
np.savetxt(deltaPrctRbcFile,outArrayRbc)
outArrayBar = np.vstack((hgbArray,deltaDefectBar,deltaLowBar,deltaHighBar))
outArrayBar = np.transpose(outArrayBar)
np.savetxt(deltaPrctBarFile,outArrayBar)

plt.figure()
plt.scatter(hgbArray,deltaDefectRbc,color = mmap[1,:], label = r'$\Delta$ defect %')
plt.scatter(hgbArray,deltaLowRbc,color = mmap[2,:], label = r'$\Delta$ low %')
plt.scatter(hgbArray,deltaHighRbc, color = mmap[5,:], label = r'$\Delta$ high %')
plt.xlabel('hgb (g/dL)')
plt.ylabel('(corrected - uncorrected) (%)')
plt.title('RBC')
plt.legend()

plt.figure()
plt.scatter(hgbArray,deltaDefectBar,color = barmap[1,:], label = r'$\Delta$ defect %')
plt.scatter(hgbArray,deltaLowBar,color = barmap[2,:], label = r'$\Delta$ low %')
plt.scatter(hgbArray,deltaHighBar, color = barmap[7,:], label = r'$\Delta$ high %')
plt.scatter('hgb (g/dL)')
plt.ylabel('(corrected - uncorrected) (%)')
plt.title('barrier')
plt.legend()

#%% plot distrs of delta prct for entire cohort of subjects

matplotlib.rcParams.update({'font.size': 14})

#%matplotlib inline
deltaBar = [deltaDefectBar, deltaLowBar, deltaHighBar]
figBar = plt.figure(figsize =(4, 8))
axBar = figBar.add_subplot(111)
bpBar = axBar.boxplot(deltaBar, patch_artist = True)
boxColor = [barmap[1,:], barmap[2,:], barmap[7,:]]
axBar.set_xticklabels([r'$\Delta$ defect %', r'$\Delta$ low %', r'$\Delta$ high %'])
axBar.set_ylim([-15, 15])
axBar.set_ylabel('(corrected - uncorrected) (%)')
axBar.set_title('barrier/gas')

for patch, flier, myColor in zip(bpBar['boxes'], bpBar['fliers'], boxColor):
    patch.set_facecolor(myColor)
    flier.set(marker = 'D', markeredgecolor = myColor, markeredgewidth = 2)

for whisker, median in zip(bpBar['whiskers'], bpBar['medians']):
    median.set(color = 'black', linewidth = 1.5)

plt.show()
figBar.savefig('../plots/chang_eta/lit_constants/tr20_flip20/deltaPrctDistrBar.png', dpi = 300, bbox_inches='tight')

deltaRbc = [deltaDefectRbc, deltaLowRbc, deltaHighRbc]
figRbc = plt.figure(figsize =(4, 8))
axRbc = figRbc.add_subplot(111)
bpRbc = axRbc.boxplot(deltaRbc, patch_artist = True)
boxColor = [mmap[1,:], mmap[2,:], mmap[5,:]]
axRbc.set_xticklabels([r'$\Delta$ defect %', r'$\Delta$ low %', r'$\Delta$ high %'])
axRbc.set_ylim([-15, 15])
axRbc.set_ylabel('(corrected - uncorrected) (%)')
axRbc.set_title('RBC/gas')

for patch, flier, myColor in zip(bpRbc['boxes'], bpRbc['fliers'], boxColor):
    patch.set_facecolor(myColor)
    flier.set(marker = 'D', markeredgecolor = myColor, markeredgewidth = 2)

for whisker, median in zip(bpRbc['whiskers'], bpRbc['medians']):
    median.set(color = 'black', linewidth = 1.5)

plt.show()
figRbc.savefig('../plots/chang_eta/lit_constants/tr20_flip20/deltaPrctDistrRbc.png', dpi = 300, bbox_inches='tight')

#%%
plt.figure()
plt.hist(deltaLowBar,color = 'red')
plt.hist(deltaLowRbc, color = 'blue')
plt.show()