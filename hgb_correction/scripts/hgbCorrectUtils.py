#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Useful methods for hgb correction of GX images. Some methods taken from "Gx_Map_utils.py"
    from GX pipeline ver2, Driehuys Lab, Duke University


Modified/written by Aryil Bechtel 2021
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from skimage import color
import nibabel as nib

def makeHists(imgCorr, imgUncorr, title, myStatsCorr, myStatsUncorr, outFile):
    
    numVox = imgCorr.size
    binwidth = 0.02
    
    imgCorr = imgCorr.flatten()
    imgUncorr = imgUncorr.flatten()
    minVal = np.min(np.concatenate([imgCorr,imgUncorr]))
    maxVal = np.max(np.concatenate([imgCorr,imgUncorr]))   

    numBins = round((maxVal - minVal)/binwidth) + 1
    binArray = np.linspace(minVal, maxVal, numBins)
    plt.figure()
    histUncorr = plt.hist(imgUncorr, bins = binArray, weights=(1/numVox)*np.ones_like(imgUncorr), \
                             histtype=u'step', label='uncorrected', color='red')
    histCorr = plt.hist(imgCorr, bins = binArray, weights=(1/numVox)*np.ones_like(imgCorr), \
                             histtype=u'step', label='corrected', color='blue')
    
    plt.xlim((0,1.5))
    plt.ylim((0,.1))
    textXpos1 = 0.02
    textYpos = 0.078
    textXpos2 = textXpos1 + 0.3
    textXpos3 = textXpos1 + 0.5
    
    if(title == 'rbc'):
        plt.text(textXpos2, textYpos, str(round(myStatsCorr[title + "_defect"]*100)) + \
                 '%\n' + str(round(myStatsUncorr[title + "_low"]*100)) + \
                 '%\n' + str(round(myStatsUncorr[title + "_high"]*100)) + '%', color='red')
        plt.text(textXpos1, textYpos, "defect:\nlow:\nhigh:")
        plt.text(textXpos3, textYpos, str(round(myStatsCorr[title + "_defect"]*100)) + \
                 '%\n' + str(round(myStatsCorr[title + "_low"]*100)) + \
                 '%\n' + str(round(myStatsCorr[title + "_high"]*100)) + '%', color='blue')
    elif(title == 'bar'):
        plt.text(textXpos2, textYpos, str(round(myStatsCorr[title + "_defect"]*100)) + \
                 '%\n' + str(round(myStatsUncorr[title + "_low"]*100)) + \
                 '%\n' + str(round(myStatsUncorr[title + "_high"]*100)) + '%', color='red')
        plt.text(textXpos1, textYpos, "defect:\nlow:\nhigh:")
        plt.text(textXpos3, textYpos, str(round(myStatsCorr[title + "_defect"]*100)) + \
                 '%\n' + str(round(myStatsCorr[title + "_low"]*100)) + \
                 '%\n' + str(round(myStatsCorr[title + "_high"]*100)) + '%', color='blue')

    if any(map(title.__contains__,{"bar"})):
        plt.title("barrier/gas")
    else:
        plt.title("RBC/gas")
    plt.ylabel('fraction of total pixels')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outFile, dpi = 300)
    
    return(histCorr,histUncorr) 

#get delta % from corrected and uncorrected stats keys
def getDeltaPrct(title, myStatsCorr, myStatsUncorr):
    
    delta = {}
    delta["defect"] = (myStatsCorr[title + "_defect"]*100) - \
        (myStatsUncorr[title + "_defect"]*100)
    delta["low"] = round(myStatsCorr[title + "_low"]*100) - \
        round(myStatsUncorr[title + "_low"]*100)
    delta["high"] = round(myStatsCorr[title + "_high"]*100) - \
        round(myStatsUncorr[title + "_high"]*100)
    
    return delta

#define colormaps (from 'GX_defineColormaps.py')
    
    #bar colormap
barmap= np.zeros((9,3))
barmap[0,:] = [0, 0, 0] # background
barmap[1,:] = [1, 0, 0] # lowest value,red
barmap[2,:] = [1, 0.7143, 0]  # yellow 2
barmap[3,:] = [0.4, 0.7, 0.4] # yellow - GREEN 3
barmap[4,:] = [0, 1, 0] # GREEN 4
barmap[5,:] = [184.0/255.0, 226.0/255.0, 145.0/255.0]
barmap[6,:] = [243.0/255.0, 205.0/255.0, 213.0/255.0]
barmap[7,:] = [225.0/255.0, 129.0/255.0, 162.0/255.0]
barmap[8,:] = [197.0/255.0, 27.0/255.0, 125.0/255.0] # highest value

long_index2color = {
        1: [0, 0, 0],
        2: [1, 0, 0],
        3: [1, 0.7143, 0],
        4: [0.4, 0.7, 0.4],
        5: [0, 1, 0],
        6: [184.0/255.0, 226.0/255.0, 145.0/255.0],
        7: [243.0/255.0, 205.0/255.0, 213.0/255.0],
        8: [225.0/255.0, 129.0/255.0, 162.0/255.0],
        9: [197.0/255.0, 27.0/255.0, 125.0/255.0]
    }

    #rbc colormap
mmap= np.zeros((7,3))
mmap[0,:] =[0, 0, 0] # background
mmap[1,:] = [1, 0, 0] # defect 1
mmap[2,:] = [1, 0.7143, 0]# yellow 2
mmap[3,:] = [0.4, 0.7, 0.4]# yellow - GREEN 3
mmap[4,:] = [0, 1, 0]#  GREEN 4
mmap[5,:] = [0, 0.57, 0.71] # green 5
mmap[6,:] = [0, 0, 1] # high-intense 6

short_index2color = {
        1: [0, 0, 0],
        2: [1, 0, 0],
        3: [1, 0.7143, 0],
        4: [0.4, 0.7, 0.4],
        5: [0, 1, 0],
        6: [0, 0.57, 0.71],
        7: [0, 0, 1]
    }

hgbKey = {'002101_highBW': 12.9,\
          '002105_highBW': 13.7,\
          '002111_highBW': 15.4,\
          '002116_highBW': 13.4,\
          '002132_highBW': 14.1,\
          '002133_highBW': 13.2,\
          '002138_highBW': 16.5,\
          '003020_highBW': 13.5,\
          '004007_highBW': 15,\
          '004008_highBW': 14,\
          '004011_highBW': 14.5,\
          '004012_highBW': 13.5,\
          '004013_highBW': 15.1,\
          '004014_highBW': 10.2}

'''following methods are from "GX_Map_utils.py" '''
    
def binning(volume, thresholds):

    # volume: mask_vented, thresholded 3D volume
    # thresholds: just the middle thresholds

    bvolume = np.ones(np.shape(volume))

    bvolume[(volume > 0) & (volume <= thresholds[0])] = 2

    for k in range(len(thresholds)-1):
        bvolume[(volume > thresholds[k]) & (volume <= thresholds[k+1])] = k+3

    bvolume[volume > thresholds[-1]] = len(thresholds)+2
    return bvolume
    
def mergeRGBandGray(ute_slice, binning_slice):
    # function combine the gray scale UTE with the RGB binning via HSV

    # construct RGB version of gray-level ute
    ute_slice_color = np.dstack((ute_slice, ute_slice, ute_slice))

    # Convert the input image and color mask to HSV
    ute_slice_hsv = color.rgb2hsv(ute_slice_color)
    binning_slice_hsv = color.rgb2hsv(binning_slice)

    # Replace the hue and saturation of the original image
    # with that of the color mask
    ute_slice_hsv[..., 0] = binning_slice_hsv[..., 0]
    ute_slice_hsv[..., 1] = binning_slice_hsv[..., 1]

    mask = ((binning_slice[:, :, 0] == 0) & (
        binning_slice[:, :, 1] == 0) & (binning_slice[:, :, 2] == 0))
    mask = ~mask

    ute_slice_hsv[mask, :] = binning_slice_hsv[mask, :]

    colormap = color.hsv2rgb(ute_slice_hsv)

    return colormap

def save3DRGB2nii(volume, file_name):
    # save 4D volume to a RGB nii
    # the input should be of size a*b*c*3(RGB)
    # There is some order difference between python and nifti that requires pre-process

    # nibabel: 2.2.1
    # numpy: 1.14
    # seems only work for 3D volume a * a * a

    color = (np.copy(volume)*255).astype('uint8')  # need uint8 to save to RGB

    # some fancy and tricky re-arrange
    color = np.transpose(color, [2, 3, 0, 1])
    cline = np.reshape(color, (1, np.size(color)))
    color = np.reshape(cline, np.shape(volume), order='A')
    color = np.transpose(color, [2, 1, 0, 3])

    # stake the RGB channels
    shape_3d = volume.shape[0:3]
    rgb_dtype = np.dtype([('R', 'u1'), ('G', 'u1'), ('B', 'u1')])
    nii_data = color.copy().view(dtype=rgb_dtype).reshape(
        shape_3d)  # copy used to force fresh internal structure
    ni_img = nib.Nifti1Image(nii_data, np.eye(4))

    nib.save(ni_img, file_name)
    
def montage(Img, slices=16):
    # plot montage(2*8) of Img
    # Img has to have 16 slices
    img_w = Img.shape[0]
    img_h = Img.shape[1]
    count = slices
    n_row = 2
    n_col = int(count/2)

    img_montage = np.zeros((n_row * img_h, n_col * img_w, 3))

    image_id = 0
    for j in range(n_row):
        for k in range(n_col):
            if image_id >= count:
                break
            sliceN, sliceM = j * img_h, k * img_w
            img_montage[sliceN:sliceN + img_h, sliceM:sliceM +
                        img_w, :] = Img[:, :, image_id, :]
            image_id += 1

    return img_montage

def getIndexMaxOnes(arr):
    # returns the starting index and ending index of the max consecutive ones
    # intitialize count
    cur_count = 0
    cur_sta = 0

    max_count = 0
    pre_state = 0

    index_sta = 0
    index_end = 0

    for i in range(0, np.size(arr)):

        if (arr[i] == 0):
            cur_count = 0
            if((pre_state == 1) & (cur_sta == index_sta)):
                index_end = i-1
            pre_state = 0

        else:
            if(pre_state == 0):
                cur_sta = i
                pre_state = 1
            cur_count += 1
            if(cur_count > max_count):
                max_count = cur_count
                index_sta = cur_sta

    return index_sta, index_end

def decideStartInterval(mask, scan_type="3D Radial"):

    # determine the starting slice and the interval for the montage
    num_slice = 3.0 #in pipeline this is 15.0

    sum_line = np.sum(np.sum(mask, axis=0), axis=0)
    if(scan_type == "2D GRE"):
        binary_arr = sum_line> -1
    else:
        binary_arr = sum_line > 300


    ind_start, ind_end = getIndexMaxOnes(binary_arr)

    flt_inter = (ind_end-ind_start)/num_slice

    # use 0.4 as a threshold to decide interval number
    if np.modf(flt_inter)[0] > 0.4:
        ind_inter = np.ceil(flt_inter).astype(int)
    else:
        ind_inter = np.floor(flt_inter).astype(int)

    return ind_start, ind_inter
    
def makeMontage(bin_index, ute_reg, index2color, ind_start, ind_inter, Subject_ID, patientPath, \
                mon_name, slices=3):
    # make montage (2*8) from binning map and ute image
    # the montage will pick the image from ind_start
    # normalize ute
    ute_thre = np.percentile(ute_reg, 99)
    ute_reg_m = np.divide(ute_reg, ute_thre)
    ute_reg_m[ute_reg_m > 1] = 1

    img_w, img_h, img_d = np.shape(bin_index)

    num_slice_all = np.shape(ute_reg)[0]
    num_slice_mon = slices
    ind_end = ind_start + ind_inter*num_slice_mon

    colormap = np.zeros((img_w, img_h, num_slice_all, 3))

    # convert each slice from index to RGB, then combine ute_reg and bin_index_RGB to HSV
    # for k in range(ind_start, ind_end, ind_inter):
    for k in range(num_slice_all):

        # convert bin_index to bin_rgb
        bin_rgb = [index2color[x] for x in bin_index[:, :, k].flatten()]
        bin_rgb = np.asarray(bin_rgb)
        bin_rgb = np.reshape(bin_rgb, (img_w, img_h, 3))

        # merge bin_rgb with ute_reg through hsv colorspacew
        # colorslice = bin_rgb # without BHUTE as background
        colorslice = mergeRGBandGray(ute_reg_m[:, :, k], bin_rgb)

        colormap[:, :, k, :] = colorslice

    # save RGB stack to nii
    if(mon_name == 'bar_montage_hgbCorrect.png'):
        nii_name = 'bar_hgbCorrect_' +Subject_ID+'.nii'
    elif(mon_name == 'rbc_montage_hgbCorrect.png'):
        nii_name = 'rbc_hgbCorrect_' + Subject_ID+'.nii'
    elif(mon_name == 'bar_montage_uncorr.png'):
        nii_name = 'bar_uncorr_' + Subject_ID+'.nii'
    elif(mon_name == 'rbc_montage_uncorr.png'):
        nii_name = 'rbc_uncorr_' + Subject_ID+'.nii'    
    else:
        nii_name = ''

    if nii_name != '':
        save3DRGB2nii(colormap, patientPath + nii_name)

    colormap_mon = colormap[:, :, ind_start:ind_end:ind_inter, :]
    # make montage from the image stack

    img_montage = montage(colormap_mon, slices=slices)

    # plot and save the montage
    plt.figure()
    
    # dedicated ventilation removes padding
    if slices == 14:
        plt.gca().set_axis_off()
        plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
                    hspace = 0, wspace = 0)
        plt.margins(0,0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())
        pad_inches = -0.0
    else:
        pad_inches = -0.05
    
    #create colorbar (following section not in GX pipeline)
    keyWord = {"bar"}
    if any(map(mon_name.__contains__,keyWord)):
        myMap = mpl.colors.ListedColormap(barmap[1:,:])
    else:
        myMap = mpl.colors.ListedColormap(mmap[1:,:])
    
    plt.imshow(img_montage, interpolation='none', cmap = myMap)
    
    if any(map(mon_name.__contains__,{"Correct"})):
        plt.colorbar(pad = 0).set_ticks([])
    
    plt.axis('off')
    plt.savefig(patientPath + mon_name, transparent=True,
                bbox_inches='tight', pad_inches=pad_inches, dpi=300)
    plt.clf()
    plt.close()
    return img_montage


def binStats(rawdata, bindata, mask, key):

    statsbox = {}
    maskall = np.sum(mask).astype('float')

    statsbox[key+'_defect'] = np.divide(np.sum((bindata == 2)), maskall)
    statsbox[key+'_low'] = np.divide(np.sum((bindata == 3)), maskall)
    
    statsbox[key+'_mean'] = np.average(rawdata[np.where(mask==1.0)])
    statsbox[key+'_median'] = np.median(rawdata[np.where(mask==1.0)])
    statsbox[key+'_SD'] = np.std(rawdata[np.where(mask==1.0)])
    statsbox[key+'_CV'] = (np.std(abs(rawdata[np.where(mask==1.0)])))\
        /(np.average(abs(rawdata[np.where(mask==1.0)])))
    if ((key == 'rbc') | (key == 'ven') | (key == 'r2b')):
        statsbox[key +
                 '_high'] = np.divide(np.sum((bindata == 6) | (bindata == 7)), maskall)
    else:
        # statsbox[key+'_high'] = np.divide(np.sum((bindata == 8)|(bindata == 9)),maskall)
        # calculate barrier high to be top 3 bins
        statsbox[key+'_high'] = np.divide(np.sum((bindata == 7)
                                                 | (bindata == 8) | (bindata == 9)), maskall)

    return statsbox

