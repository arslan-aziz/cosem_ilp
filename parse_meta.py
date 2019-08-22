'''
Regex basics
Identifiers
\d any number
\D anything but a #
\s space
\S anything but a space
\w any character
\W anything but a character
. any character, except newline
\b white space around words
\. period

Modiifiers
+ Match 1 or more
? Match 0 or 1
* Match 0 or more
$ Match end of string
^ Match begin of string
| either or
[] range of values
{x} expecting x amount
{x,y} expecting x to y amount

White space characters
\n new line
\s space
\t tab
\e escape
\f form feed
\r return

findall(r'<...>') r indicates regex expression
'''

import re
import json
import os


#meta data, element per ROI
PATH = []
EXTENT = []
ORIGIN = []
LABELS = []
ROI_NAME = []

#metadata text files, ONE AT A TIME
meta_files = ('Cell2_Crop3.txt',)

#read each file
for f in meta_files:
    with open('meta/'+f,'r') as text:
        lines = text.readlines()

    folder = ""
    file = ""
    path = ""
    extent = []
    origin = [None]*3
    lab_begin_line = beg = lab_end_line = end = 0
    labels = []


    #loop all lines in file
    for line in lines:

        #extract folder name
        if not folder:
            fold = re.findall(r'(?<=^#\s).+(?=\s$)', line)
            if fold:
                #TODO: if last character whitespace, remove last character
                #folder = fold[0][:-1]
                folder = fold[0]
                folder = re.sub(r'[\W]+', '_', folder)

        #extract file name
        #if not file:
        #check every line and grab the last file name
        file = re.findall(r'(?=Cell2).+(?=\.h5\s)', line)
        if file:
            file_name = file[0]+'.h5'

        #extract extent
        if len(extent) == 0:
            ext = re.findall(r'(?<=^\+\sROI\sSize\s\(pixel\)\:[ \t]).+(?=\s)',line)
            if ext:
                extent = list(map(int,re.findall(r'\d{2,3}',ext[0])))


        #extract x_origin
        if origin[0]==None:
            x_org = re.findall(r'(?<=\+\sx\:[ \t]).+(?=\s)',line)
            if x_org:
                #origin = list(map(int,[x_org[0],y_org[0],z_org[0]]))
                origin[0] = int(x_org[0])

        #extract y origin
        if origin[1] == None:
            y_org = re.findall(r'(?<=\+\sy\:[ \t]).+(?=\s)', line)
            if y_org:
                origin[1] = int(y_org[0])

        if origin[2] == None:
            z_org = re.findall(r'(?<=\+\sz\:[ \t])\S*(?=\s)', line)
            if z_org:
                origin[2] = int(z_org[0])

        #extract label section start
        if lab_begin_line == 0:
            if re.findall(r'###\sLabel\sFields\:', line):
                lab_begin_line = beg
            beg += 1

        #extract label section end
        if lab_end_line == 0:
            if re.findall(r'###\sAnnotations\:', line):
                lab_end_line = end
            end += 1

    #extract label values
    # begin reading labels at  + 1
    # index in array (from 1) = value in annotations
    x = 0
    for i in range(lab_begin_line + 1, lab_end_line - 1):
        lab = re.findall(r'(?<=\d\.\s).+(?=\s)', lines[i])
        if lab:
            labels.append(lab[0])
            labels[x] = re.sub(r'[\W]+','', labels[x]).lower()
            x += 1

    #assemble file path
    path = "/groups/funke/home/aziza/Documents/COSEM/Annotations/" + folder + "/" + file_name

    #save to global variables
    PATH.append(path)
    EXTENT.append(extent)
    ORIGIN.append(origin)
    LABELS.append(labels)
    ROI_NAME.append(f[:-4])

#combination table: combine <values> in annotated as such to compare to channel <key> in prediction
#organelles: PM, mito, golgi, vesicle, MVB, lysosome, LD, ER, NucEnv, Nuc
#keys match prediction channel names
#values match annotated label names
COMB_GT = {
    'plasma_membrane':['plasmamembrane'],
    'mito':['mitodna','mitolumen','mitomembrane'],
    'golgi':['golgilumen','golgimembrane'],
    #how is er membrane distingushed from NE
    'er':['ereslumen','erlumen','eresmembrane','ermembrane'],
    'NE':['nuclearporeout','nuclearporein','nelumen','nemembrane'],
    'nucleus':['nemembrane','nelumen','nuclearporeout','nuclearporein','hchrom','nhchrom',
               'echrom','nucleoplasm','nucleolus','chrom'],
    'vesicle':['vesiclemembrane','vesiclelumen'],
    'MVB':['mvbmembrane','mvblumen'],
    'lysosome':['lysosomemembrane','lysosomelumen'],
    'LD':['ldmembrane','ldlumen'],
    'cytosol':['cytosol']
}

#same keys as COMB_ANN but values are names of prediction channels
COMB_YHAT = {
    'plasma_membrane':['plasma_membrane'],
    'mito':['mito_DNA','mito_membrane','mito'],
    'golgi':['golgi_membrane','golgi'],
    'er':['ERES','er_membrane','er'],
    #'NE':['nuclear_pore','nuclear_pore_out','NE'],
    'nucleus':['nuclear_pore','nuclear_pore_out','NE','NHChrom','NEChrom',
               'EChrom','nucleolus','nucleus'],
    'vesicle':['vesicle_membrane','vesicle'],
    'MVB':['MVB_membrane','MVB'],
    'lysosome':['lysosome_membrane','lysosome'],
    'LD':['LD_membrane','LD'],
    'cytosol':[]
}
ORG_KEYS = list(COMB_YHAT.keys())

#Searches list <full> for values <search> and returns index
#Adds one to return value to make label
def findLabels(search,full):
    find = []
    name = []
    for s in search:
        for i in range(0,len(full)):
            if s == full[i]:
                find.append(i+1)
                name.append(s)
    return find,name


COMB_GT_FOUND = {}
COMB_GT_LABEL = {}
#loop thru all major organelles
for org in ORG_KEYS:
    #use combined dict to select required labels for given organelle
    #TODO: hard coded assuming only one ROI crop (LABELS[0])
    label,name = findLabels(COMB_GT[org],LABELS[0])
    if(label!=[]):
        COMB_GT_FOUND.update({org:name})
        COMB_GT_LABEL.update({org:label})




print('done')

print('dumping in json...')
output_path = 'meta/cosem_metadata.json'
if os.path.isfile(output_path):
    os.remove(output_path)

print('PATH:',PATH)
print('EXTENT:',EXTENT)
print('ORIGIN:',ORIGIN)
print('ROI_NAME:',ROI_NAME)
print('COMB_GT:',COMB_GT)
print('COMB_GT_FOUND:',COMB_GT_FOUND)
print('COM_GT_LABEL:',COMB_GT_LABEL)
print('COMB_YHAT:',COMB_YHAT)
print('ORG_KEYS:',ORG_KEYS)
print('LABELS:',LABELS)

with open(output_path,'w') as f:
    json.dump({'PATH':PATH,'EXTENT':EXTENT,'ORIGIN':ORIGIN,
               'LABELS':LABELS,'ROI_NAME':ROI_NAME,
               'COMB_GT':COMB_GT,'COMB_YHAT':COMB_YHAT,
               'COMB_GT_FOUND':COMB_GT_FOUND,'COMB_GT_LABEL':COMB_GT_LABEL,
               'ORG_KEYS':ORG_KEYS},f)






