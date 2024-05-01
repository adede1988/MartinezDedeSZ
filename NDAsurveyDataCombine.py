# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 09:59:00 2022

@author: pc1aod
"""

#analyzing NIMH data downloaded from the NIMH data archive

#set up the environment
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
from os import walk,getcwd
import matplotlib.pyplot as plt
import math as math
#custom helper functions for reading in data files and performing key analysis steps

pd.options.mode.chained_assignment = None  # default='warn'


path = 'G:\\Package_1204973\\txtDocs\\'


df = pd.read_csv(r'C:\Users\dtf8829\Documents\GitHub\MartinezDedeSZ\autismBiomarkersAllData2.csv') 
df = df[df['dataSet'] == 'bpSZ']


metaFiles = list(walk(path))[0][2]
metaFiles = Series([path + file for file in metaFiles])

#data frame to store the descriptions of each field in the surveys
infoKeys = pd.DataFrame()


for file in metaFiles: 
    curFile = pd.read_csv(file, sep = '\t')
    curFile.dropna(axis = 1, thresh = 50, inplace = True)
    if len(curFile) > 50: 
        numCols = len(curFile.columns)
        for col in range(numCols):
            curRow = pd.Series(data=None, index = ['itemName', 'itemDescrip', 'itemFile'], dtype = 'string')
            curRow[0] = curFile.columns[col]
            curRow[1] = curFile.iloc[0,col]
            curRow[2] = file
            infoKeys = pd.concat([infoKeys, curRow], axis = 1)
        curFile['key'] = curFile['subjectkey'] 
        df = pd.merge(df, curFile, on = ['key'], how = 'left')
        df = df.loc[:,[not s.endswith('_x') for s in df.columns]]
        df = df.loc[:,[not s.endswith('_y') for s in df.columns]]
        df.drop_duplicates(subset=['name'], keep='first', inplace=True)

    infoKeys =infoKeys.transpose()
    df.to_csv(r'C:\Users\dtf8829\Documents\GitHub\MartinezDedeSZ\combinedEEG_surveyData.csv')
    infoKeys.to_csv(r'C:\Users\dtf8829\Documents\GitHub\MartinezDedeSZ\itemKeys.csv')


