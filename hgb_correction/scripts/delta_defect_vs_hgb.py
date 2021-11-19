#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to plot change in defect, low, and high % vs hgb (for a specific subject)

Workflow:
    -input .mat output file from GX pipeline
    -extract mask, bar2gas, rbc2gas and UTE matrices
    -for a set range of hgb values:    
        -divide bar2gas and rbc2gas matrices by scaling factor from 
            MOXE hgb dependence on rbc and bar signals
        -using binStats method, store bar and rbc defect, low, and high %
    -subtract uncorrected % (% for hgb = 15 g/dL) from corrected %

Aryil Bechtel 2021
"""
import numpy as np
from scipy.io import loadmat
import matplotlib 
import matplotlib.pyplot as plt
from hgbCorrectUtils import binning, binStats, getDeltaPrct
from hgbCorrectUtils import barmap, mmap, hgbKey

#set hgb values to evaluate
hgbStart = 10
hgbEnd = 20
numHgb = 70
hgbRange = np.linspace(hgbStart,hgbEnd,numHgb)

patientName = '002101_highBW'
patientPath = "../data/" + patientName + "/"
matPath = patientPath + patientName + ".mat"
hgbRef = 15 #healthy reference hgb g/dL (such that hgb correction factor = 1)

#create delta % arrays to store % values
deltaDefectBar = np.zeros((hgbRange.size))
deltaLowBar = np.zeros((hgbRange.size))
deltaHighBar = np.zeros((hgbRange.size))
deltaDefectRbc = np.zeros((hgbRange.size))
deltaLowRbc = np.zeros((hgbRange.size))
deltaHighRbc = np.zeros((hgbRange.size))

#load .mat file and extract matrices
matFile = loadmat(matPath)

mask = matFile['mask_reg']
bar = np.absolute(matFile['bar2gas'])
rbc = np.absolute(matFile['rbc2gas'])
ute = matFile['ute_reg']

#define hgb correction fns and scale bar and rbc matrices
barPoly3 = [-2.925e-5, 0.0026, -0.1163, 2.2509]
rbcPoly3 = [2.286e-5, -0.0021, 0.0909, 0.0226]          

#index counter for delta % arrays
i=0

for hgb in hgbRange:
    
    #scale bar and rbc matrices
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
#%% 
#%matplotlib inline

matplotlib.rcParams.update({'font.size': 13})

#load real delta% vs hgb data
inArrayRbc = np.loadtxt('../data/deltaPrctRbc')
inArrayBar = np.loadtxt('../data/deltaPrctBar')

plt.figure()
plt.plot(hgbRange,deltaDefectRbc,color = mmap[1,:], label = r'$\Delta$ defect %')
plt.plot(hgbRange,deltaLowRbc,color = mmap[2,:], label = r'$\Delta$ low %')
plt.plot(hgbRange,deltaHighRbc, color = mmap[5,:], label = r'$\Delta$ high %')
plt.scatter(inArrayRbc[:,0], inArrayRbc[:,1], color = mmap[1,:])
plt.scatter(inArrayRbc[:,0], inArrayRbc[:,2], color = mmap[2,:])
plt.scatter(inArrayRbc[:,0], inArrayRbc[:,3], color = mmap[5,:])
plt.xlabel('hgb (g/dL)')
plt.ylabel('(corrected - uncorrected) (%)')
plt.title('RBC/gas')
plt.legend()
outFileRbc = patientPath + 'deltaPrctRbc_vs_hgb_realData.png'
plt.savefig(outFileRbc, dpi = 300)


plt.figure()
plt.plot(hgbRange,deltaDefectBar,color = barmap[1,:], label = r'$\Delta$ defect %')
plt.plot(hgbRange,deltaLowBar,color = barmap[2,:], label = r'$\Delta$ low %')
plt.plot(hgbRange,deltaHighBar, color = barmap[7,:], label = r'$\Delta$ high %')
plt.scatter(inArrayBar[:,0], inArrayBar[:,1], color = barmap[1,:])
plt.scatter(inArrayBar[:,0], inArrayBar[:,2], color = barmap[2,:])
plt.scatter(inArrayBar[:,0], inArrayBar[:,3], color = barmap[7,:])
plt.xlabel('hgb (g/dL)')
plt.ylabel('(corrected - uncorrected) (%)')
plt.title('barrier/gas')
plt.legend(loc = 2, bbox_to_anchor=(0.4,1))
outFileBar = patientPath + 'deltaPrctBar_vs_hgb_realData.png'
plt.savefig(outFileBar, dpi = 300)
