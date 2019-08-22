#Create combined ground truth and prediction data set for each ROI for each organelle

import numpy as np
import zarr
import h5py
import z5py
import json

#read meta data from parse_meta.py

meta_path = 'meta/cosem_metadata.json'
with open(meta_path,'r') as f:
    meta_data = json.load(f)

PATH = meta_data['PATH']
EXTENT = meta_data['EXTENT']
ORIGIN = meta_data['ORIGIN']
#LABELS = meta_data['LABELS']
ROI_NAME = meta_data['ROI_NAME']
COMB_GT = meta_data['COMB_GT_FOUND']
COMB_GT_LABEL = meta_data['COMB_GT_LABEL']
COMB_YHAT = meta_data['COMB_YHAT']
ORG_KEYS = meta_data['ORG_KEYS']

PRED_PATH = '/groups/funke/home/aziza/Documents/COSEM/predictions/cell2_257000_89_207.n5'

#OUTPUT_PATH = '/groups/funke/home/aziza/Documents/cosem_scripts/comp_volumes.h5'
OUTPUT_PATH = '/groups/funke/home/aziza/Documents/cosem_scripts/comp_volumes.zarr'

RAW_PATH = '/groups/funke/home/aziza/Documents/COSEM/raw_data/cell2/test2.n5'

def crop_middle(array,roi,channel=None):

    offset = np.floor((np.array(array.shape)-np.array(roi))/2)
    offset = offset.astype(int)
    return array[offset[0]:offset[0]+roi[0],
           offset[1]:offset[1]+roi[1],
           offset[2]:offset[2]+roi[2]]

#with h5py.File(OUTPUT_PATH,'w') as output:

output = zarr.open(OUTPUT_PATH,mode='w')
for i in range(0,len(ROI_NAME)):

    #load annotated dataset and raw dataset
    f1 = h5py.File(PATH[i])
    ann_data = f1['volumes/labels/gt'][0:-1:2,0:-1:2,0:-1:2]
    raw_data = f1['volumes/raw'][:,:,:]
    #raw_data = np.reshape(raw_data,ann_data.shape)
    raw_data = crop_middle(raw_data,[EXTENT[0][2],EXTENT[0][1],EXTENT[0][0]])
    f1.close()

    #output.create_dataset(ROI_NAME[i] + "/" + "raw", data=raw_data)
    output[ROI_NAME[i] + "/" + "raw"] = raw_data
    output[ROI_NAME[i] + "/" + "raw"].attrs['resolution']=(4,4,4)


    ##CREATE GT DATA
    for org in list(COMB_GT_LABEL.keys()):
        #use combined dict to select required labels for given organelle
        #loc_labels,loc_name = findLabels(COMB_GT[org],LABELS[i])
        ann_data_bin = np.zeros(ann_data.shape)
        #set loc_labels to 1
        for x,lab in enumerate(COMB_GT_LABEL[org]):
            #ann_data_bin[ann_data==lab]=1
            #ann_data_bin = np.zeros(ann_data.shape)
            ann_data_bin[ann_data==lab]=1
            #output.create_dataset(ROI_NAME[i]+"/"+org+"/gt/"+COMB_GT[org][x],data=ann_data_bin)
            output[ROI_NAME[i]+"/"+org+"/gt/"+COMB_GT[org][x]] = ann_data_bin

    ##CREATE PREDICTION DATA
    f2 = z5py.File(PRED_PATH)
    for org in ORG_KEYS:
        for channel in COMB_YHAT[org]:
            pred_data = f2[channel][ORIGIN[i][2]:ORIGIN[i][2]+EXTENT[i][2],
                                ORIGIN[i][1]:ORIGIN[i][1]+EXTENT[i][1],
                                ORIGIN[i][0]:ORIGIN[i][0]+EXTENT[i][0]]
            #output.create_dataset(ROI_NAME[i] + "/" + org + "/yhat/" +channel,data=pred_data)
            output[ROI_NAME[i] + "/" + org + "/yhat/" +channel] = pred_data

#output.close()




















