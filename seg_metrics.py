import numpy as np
import h5py
from scipy import ndimage
from scipy import optimize
from funlib.evaluate import rand_voi
import matplotlib.pyplot as plt
import json
import os
from skimage.morphology import watershed

#load labeled data
ROI = h5py.File('/groups/funke/home/aziza/Documents/cosem_scripts/comp_volumes.h5','r')

with open('meta/cosem_metadata.json','r') as f:
    META_DATA = json.load(f)

ROI_NAME = META_DATA['ROI_NAME']
ORG_KEYS = META_DATA['ORG_KEYS']
COMB_GT = META_DATA['COMB_GT_FOUND']
COMB_YHAT = META_DATA['COMB_YHAT']

gt_bin = ROI[ROI_NAME[0]+"/"+ORG_KEYS[0]+"/gt/"+COMB_GT[ORG_KEYS[0]][-1]][:,:,:]
yhat = ROI[ROI_NAME[0]+"/"+ORG_KEYS[0]+'/yhat/'+COMB_YHAT[ORG_KEYS[0]][-1]][:,:,:]


thresh=127
yhat_bin = np.copy(yhat)
yhat_bin[yhat<thresh]=0
yhat_bin[yhat>=thresh]=1

#lumen only
gt_lum_bin = ROI[ROI_NAME[0]+"/"+ORG_KEYS[0]+"/gt/"+COMB_GT[ORG_KEYS[0]][-2]][:,:,:]
yhat_mem_bin = ROI[ROI_NAME[0]+"/"+ORG_KEYS[0]+'/yhat/'+COMB_YHAT[ORG_KEYS[0]][-2]][:,:,:]
yhat_mem_bin[yhat_mem_bin < thresh] = 0
yhat_mem_bin[yhat_mem_bin >= thresh] = 1
yhat_lum_bin = np.logical_and(yhat_bin,np.logical_not(yhat_mem_bin))
ROI.close()

print('compute connected components')
#connected components
gt_comp, gt_num_comp = ndimage.measurements.label(gt_bin)
yhat_comp, yhat_num_comp = ndimage.measurements.label(yhat_bin)
gt_lum_comp, gt_lum_num_comp = ndimage.measurements.label(gt_lum_bin)
yhat_lum_comp, yhat_lum_num_comp = ndimage.measurements.label(yhat_lum_bin)

if(gt_comp.dtype != np.dtype('uint64')):
    gt_comp = gt_comp.astype('uint64')

if(yhat_comp.dtype != np.dtype('uint64')):
    yhat_comp = yhat_comp.astype('uint64')

print('compute vi')
#vi
rand_voi_arrays = rand_voi(gt_comp,yhat_comp,False)
rand_voi_control = rand_voi(gt_comp,gt_comp,False)

print('compute iou over all')
#IOU over all
iou_score = np.sum(np.logical_and(gt_bin,yhat_bin))/np.sum(np.logical_or(gt_bin,yhat_bin)).item()
iou_control = np.sum(np.logical_and(gt_bin,gt_bin))/np.sum(np.logical_or(gt_bin,gt_bin)).item()

print('compute cleft matching')
#cleft matching
#signed distance transformation on gt_bin to create gt_sdt
sdt = ndimage.morphology.distance_transform_edt(gt_bin,sampling=(4,4,4)) -\
        ndimage.morphology.distance_transform_edt(np.logical_not(gt_bin),sampling=(4,4,4))
stdt = np.tanh(sdt/50.)
stdt_uint8 = stdt+np.abs(np.min(stdt))
stdt_uint8 = stdt_uint8*(255/np.max(stdt_uint8))
stdt_uint8 = stdt_uint8.astype(np.uint8)
#find where gt_bin and yhat_bin do not intersect
false = np.logical_xor(gt_bin,yhat_bin)
#compute sum of distance btw. distances of those voxels
cleft = np.sum(np.abs(yhat[false]-stdt_uint8[false])).item()
false_control = np.logical_xor(stdt_uint8,stdt_uint8)
cleft_control = np.sum(np.abs(stdt_uint8[false_control]-stdt_uint8[false_control])).item()

#image pairs + IOU
#results_volume = h5py.File('results/label_volumes.h5','r')
with open('metrics/match_metadata.json','r') as f:
    meta_data = json.load(f)
row_ind = meta_data['row_ind']
col_ind = meta_data['col_ind']
#tot_cost = meta_data['tot_cost']

# #matching pairs
# gt_ind = row_ind+1
# yhat_ind = col_ind+1
# print(gt_ind)
# print(yhat_ind)

# #counts of voxels per component
# yhat_comp_count = []
# for i in range(1, len(col_ind) + 1):
#     yhat_comp_count.append(len(yhat_comp[yhat_comp == i]))
#
# gt_comp_count = []
# for j in range(1, len(row_ind) + 1):
#     gt_comp_count.append(len(gt_comp[gt_comp == j]))
#
# print('gt')
# print(gt_comp_count)
# print('yhat')
# print(yhat_comp_count)

print('match image pairs')
#IOU per comp (+ control)
iou_per_case = []
mismatch_per_case = []
mismatch = []
for j in range(0,len(row_ind)):
    iou_per_comp = []
    for i in range(0,len(row_ind[j])):
        gt_i = np.zeros(gt_comp.shape)
        if j<len(row_ind)-3:
            gt_i[gt_comp == row_ind[j][i]+1] = 1
        else:
            gt_i[gt_lum_comp == row_ind[j][i]+1] = 1
        yhat_i = np.zeros(yhat_comp.shape)
        if j<len(row_ind)-3:
            yhat_i[yhat_comp == col_ind[j][i]+1] = 1
        else:
            yhat_i[yhat_lum_comp == col_ind[j][i]+1] = 1
        iou_per_comp.append(np.sum(np.logical_and(gt_i,yhat_i))/np.sum(np.logical_or(gt_i,yhat_i)).item())
        #TODO: saving mismatch pixels fucks memory usage
        # if(iou_per_comp[i]<0.8):
        #     mismatch.append(yhat_i[:].tolist())
    iou_per_case.append(iou_per_comp)
    mismatch_per_case.append(mismatch)

#WATERSHED

print('save to dict')
metrics_dict = {
    'vi':[rand_voi_arrays,rand_voi_control],
    'iou':[iou_control,iou_score],
    'iou_per_comp':iou_per_case,
    'cleft':[cleft,cleft_control],
    #'mismatch':mismatch_per_case
}


print('dumping in json...')
output_path = 'metrics/'+ORG_KEYS[0]+'_'+ROI_NAME[0]+'.json'
if os.path.isfile(output_path):
    os.remove(output_path)

with open(output_path,'w') as f:
    json.dump(metrics_dict,f)
