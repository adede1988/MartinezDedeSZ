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


################### EDIT HERE FOR DATA SET!!!! ########################
#Data sets: 
datSet = 'biomarkDev'  #  , socBrain, bpSZ, biomarkCon,  femaleASD  bpSZ 
#######################################################################



if datSet == 'socBrain': 
    path = r'Z:\User\pc1aod\Package_1197801' #path to the social brain dataset
elif datSet == 'biomarkCon': 
    path = r'Z:\User\pc1aod\Package_1202361' #path to the Autism biomarkers consortium for clinical trials
elif datSet == 'biomarkDev': 
    path = r'Z:\User\pc1aod\Package_1203248' #path to biomarkers of developmental trajectories and treatment in ASD
elif datSet == 'bpSZ': 
    path = r'Z:\User\pc1aod\Package_1204973' #path to Bipolar & Schizophrenia Consortium for Parsing Intermediate Phenotypes (B-SNIP 1)
elif datSet == 'femaleASD': 
    path = r'Z:\User\pc1aod\Package_1204978' #path to Multimodal Developmental Neurogenetics of Females with ASD

df = pd.read_csv(path + '\\eeg_sub_files01.txt', sep = '\t') 
df.dropna(axis = 1, thresh = 3, inplace = True)


metaFiles = list(walk(path + '\\txtFiles'))[0][2]
metaFiles = Series([path +'\\txtFiles\\'+ file for file in metaFiles])


if datSet == 'socBrain': 
    #ADOS info
    diagInfo = pd.read_csv(metaFiles[['ados4_200' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo = diagInfo.loc[:,['subjectkey', 'scoresumm_abtotal']]
    diagInfo2 = pd.read_csv(metaFiles[['ados4_201' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo2 = diagInfo2.loc[:,['subjectkey', 'scoresumm_abtotal']]
    diagout = pd.merge(diagInfo, diagInfo2, on = ['subjectkey'], how = 'outer')
    diagout = diagout.fillna(0)
    diagout['ADOS'] = 0
    diagout['ADOS'].iloc[1:] = diagout['scoresumm_abtotal_x'].iloc[1:].astype('int') + diagout['scoresumm_abtotal_y'].iloc[1:].astype('int')
    diagout = diagout.loc[:, ['subjectkey', 'ADOS']]   
    diagout = diagout.drop_duplicates(subset = ['subjectkey'])                                   
    df = pd.merge(df, diagout, on = ['subjectkey'], how = 'inner')
    df = df.loc[pd.notnull(df['data_file1']), :]
    df = df.drop_duplicates()
    df['ADOS_version'] = '4'
    #diagnosis 10 = autism 7 = ASD
    df['diag'] = 'CON'
    df['diag'].iloc[np.where(df['ADOS']>6)] = 'ASD'
    df['diag'].iloc[np.where(df['ADOS']>9)] = 'AD'
    #intelligence
    intel1 = pd.read_csv(metaFiles[['wtar' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel1 = intel1.loc[:,['subjectkey', 'wtar_ss']]
    intel1 = intel1.drop_duplicates(subset = ['subjectkey'])
    df = pd.merge(df, intel1, on = ['subjectkey'], how = 'inner')
    df = df.loc[pd.notnull(df['data_file1']), :]
    df = df.drop_duplicates()
    #make the output dataframe
    key = df['subjectkey'].iloc[1:]
    sex = df['sex'].iloc[1:]
    age= df['interview_age'].iloc[1:]
    file= df['data_file1'].iloc[1:]
    task= df['comments_misc'].iloc[1:]
    iq1 = df['wtar_ss'].iloc[1:]
    iq2 = iq1
    diag = df['diag'].iloc[1:]
    #stitching the info together
    outDF = pd.concat([sex,age,task,diag,iq1,iq2,file], axis = 1)
    outDF.columns = ['sex', 'age', 'task', 'diag', 'iq1', 'iq2', 'file' ]
    outDF['iq1_measure'] = 'WechslerReading'
    outDF['iq2_measure'] = 'WechslerReading'
    outDF['ADOS'] = df['ADOS'].iloc[1:]
    outDF['ADOS_version'] = df['ADOS_version'].iloc[1:]
    outDF['dataSet'] = datSet 
    outDF.index = key
    outDF = outDF.loc[outDF['task'].str.contains('Rest'), :]
    outDF['eyes'] = 'open'
    
    outDF.to_csv(r'Z:\User\pc1aod\CODE\GEDbounds_clusterImprove\socBrainDat.csv')
    
elif datSet == 'biomarkCon': 
    #subset to only _r.mat files for resting state data, also drop to only one per subject
    df = df.loc[['_r.mat' in file for file in df['data_file1']], :]
    df = df.drop_duplicates(subset = ['subjectkey'])
    #ADOS info
    diagInfo = pd.read_csv(metaFiles[['ados1' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo = diagInfo.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo['ADOS_version'] = 1
    diagInfo2 = pd.read_csv(metaFiles[['ados2' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo2 = diagInfo2.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo2['ADOS_version'] = 2
    diagInfo3 = pd.read_csv(metaFiles[['ados3' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo3 = diagInfo3.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo3['ADOS_version'] = 3
    diagOut = pd.merge(diagInfo, diagInfo2, on = ['subjectkey'], how = 'outer')
    diagOut = pd.merge(diagOut, diagInfo3, on = ['subjectkey'], how = 'outer')
    scoreCols = diagOut.columns[['scoresumm' in colName for colName in  diagOut.columns]]
    versionOut = diagOut.loc[:,['ADOS_version' in colName for colName in  diagOut.columns]].apply(np.nanmin, axis = 1)
    diagOut['ADOS_final'] = [diagOut[scoreCols[v-1]].iloc[ii] for ii,v in enumerate(versionOut)]
    diagOut['ADOS_version_final'] = versionOut
    diagOut = diagOut.loc[:,['subjectkey','ADOS_final','ADOS_version_final']]   
    diagOut = diagOut.drop_duplicates(subset = ['subjectkey'])  
    diagOut['ADOS_final'].iloc[0] = 0
    diagOut['ADOS_final'] = diagOut['ADOS_final'].astype('int')                             
    df = pd.merge(df, diagOut, on = ['subjectkey'], how = 'inner')
    df = df.drop_duplicates()
    df['interview_age'] = df['interview_age'].astype('int')
    df['diag'] = 'CON'
    df['diag'].iloc[np.where( (df['ADOS_final']>7) & (df['ADOS_version_final']==1) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>11) & (df['ADOS_version_final']==1) )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>6) & (df['ADOS_version_final']==2) & (df['interview_age']<60) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>9) & (df['ADOS_version_final']==2) & (df['interview_age']<60)  )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>7) & (df['ADOS_version_final']==2) & (df['interview_age']>=60) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>8) & (df['ADOS_version_final']==2) & (df['interview_age']>=60)  )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>6) & (df['ADOS_version_final']==3) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>8) & (df['ADOS_version_final']==3) )] = 'AD'
    
    #intelligence data
    intel1 = pd.read_csv(metaFiles[['das_ii_early' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel1 = intel1.loc[:,['subjectkey', 'dasii_eyr_gca_ss']]
    intel1 = intel1.drop_duplicates(subset = ['subjectkey'])
    intel2 = pd.read_csv(metaFiles[['das_ii_school' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel2 = intel2.loc[:,['subjectkey', 'dasii_sar_gca_ss']]
    intel2 = intel2.drop_duplicates(subset = ['subjectkey']) 
    intel = pd.merge(intel1, intel2, on = ['subjectkey'], how = 'outer')
    df = pd.merge(df, intel, on = ['subjectkey'], how = 'inner')
    df['dasii_sar_gca_ss'] = df['dasii_sar_gca_ss'].astype('float')
    df['dasii_eyr_gca_ss'] = df['dasii_eyr_gca_ss'].astype('float')
    
     #make the output dataframe
    
    key = df['subjectkey']
    sex = df['sex']
    age= df['interview_age']
    file= df['data_file1']
    diag = df['diag']
    iq1 = df['dasii_sar_gca_ss']
    iq2 = df['dasii_eyr_gca_ss']
    task = diag
    ADOS = df['ADOS_final']
    ADOS_version = df['ADOS_version_final']
    
    outDF = pd.concat([sex,age,task,diag,iq1,iq2,file, ADOS, ADOS_version], axis = 1)
    outDF.columns = ['sex', 'age', 'task', 'diag', 'iq1', 'iq2', 'file' , 'ADOS', 'ADOS_version']
    outDF.index = key
    
    outDF['iq1_measure'] = 'das_school'
    outDF['iq2_measure'] = 'das_early'
    outDF['dataSet'] = datSet 
    outDF['eyes'] = 'open'
    outDF['task'] = 'rest'
    
    outDF.to_csv(r'Z:\User\pc1aod\CODE\GEDbounds_clusterImprove\biomarkConDat.csv')
    
   
elif datSet == 'biomarkDev': 
     #subset to only sessions that started with rest, that way I can just take the first X minutes 
     #and be reasonably sure I've limited analysis to resting state
    #df = df.loc[['Resting' in session.split(';')[0] for session in df['experiment_notes']],:]
    df = df.drop_duplicates(subset = ['subjectkey'])
    df = df.iloc[1:,:]
    
    #ADOS: 
    diagInfo = pd.read_csv(metaFiles[['ados1' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo = diagInfo.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo['ADOS_version'] = 1
    diagInfo2 = pd.read_csv(metaFiles[['ados2' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo2 = diagInfo2.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo2['ADOS_version'] = 2
    diagInfo3 = pd.read_csv(metaFiles[['ados3' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo3 = diagInfo3.loc[:,['subjectkey', 'scoresumm_abtotal']]
    diagInfo3['ADOS_version'] = 3
    diagInfo4 = pd.read_csv(metaFiles[['adost' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo4 = diagInfo4.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo4['ADOS_version'] = 0
    diagOut = pd.merge(diagInfo, diagInfo2, on = ['subjectkey'], how = 'outer')
    diagOut = pd.merge(diagOut, diagInfo3, on = ['subjectkey'], how = 'outer')
    diagOut = pd.merge(diagOut, diagInfo4, on = ['subjectkey'], how = 'outer')
    scoreCols = diagOut.columns[['scoresumm' in colName for colName in  diagOut.columns]]
    versionOut = diagOut.loc[:,['ADOS_version' in colName for colName in  diagOut.columns]].apply(np.nanmin, axis = 1)
    diagOut['ADOS_final'] = [diagOut[scoreCols[v-1]].iloc[ii] for ii,v in enumerate(versionOut)]
    diagOut['ADOS_version_final'] = versionOut
    diagOut = diagOut.loc[:,['subjectkey','ADOS_final','ADOS_version_final']]   
    diagOut = diagOut.drop_duplicates(subset = ['subjectkey'])  
    diagOut['ADOS_final'].iloc[0] = 0
    diagOut['ADOS_final'] = diagOut['ADOS_final'].astype('int')                             
    df = pd.merge(df, diagOut, on = ['subjectkey'], how = 'inner')
    df = df.drop_duplicates()
    #diagnosis
    df['interview_age'] = df['interview_age'].astype('int')
    df['diag'] = 'CON'
    df['diag'].iloc[np.where( (df['ADOS_final']>7) & (df['ADOS_version_final']==1) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>11) & (df['ADOS_version_final']==1) )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>6) & (df['ADOS_version_final']==2) & (df['interview_age']<60) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>9) & (df['ADOS_version_final']==2) & (df['interview_age']<60)  )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>7) & (df['ADOS_version_final']==2) & (df['interview_age']>=60) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>8) & (df['ADOS_version_final']==2) & (df['interview_age']>=60)  )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>6) & (df['ADOS_version_final']==3) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>8) & (df['ADOS_version_final']==3) )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>7) & (df['ADOS_version_final']==0) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>11) & (df['ADOS_version_final']==0) )] = 'AD'    
    
    #intelligence data
    intel1 = pd.read_csv(metaFiles[['mullen' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel1 = intel1.loc[:,['subjectkey', 'scoresumm_elc_std_score']]
    intel1 = intel1.drop_duplicates(subset = ['subjectkey'])
    intel2 = pd.read_csv(metaFiles[['leiter_vr01' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel2 = intel2.loc[:,['subjectkey', 'leit_brief_iq']]
    intel2 = intel2.drop_duplicates(subset = ['subjectkey']) 
    intel = pd.merge(intel1, intel2, on = ['subjectkey'], how = 'outer')
    df = pd.merge(df, intel, on = ['subjectkey'], how = 'inner')
    df['leit_brief_iq'] = df['leit_brief_iq'].astype('float')
    df['scoresumm_elc_std_score'] = df['scoresumm_elc_std_score'].astype('float')
    
    #output dataframe:
    key = df['subjectkey']
    sex = df['sex']
    age= df['interview_age']
    file= df['data_file1']
    diag = df['diag']
    iq1 = df['leit_brief_iq']
    iq2 = df['scoresumm_elc_std_score']
    task = diag
    ADOS = df['ADOS_final']
    ADOS_version = df['ADOS_version_final']
    
    outDF = pd.concat([sex,age,task,diag,iq1,iq2,file, ADOS, ADOS_version], axis = 1)
    outDF.columns = ['sex', 'age', 'task', 'diag', 'iq1', 'iq2', 'file' , 'ADOS', 'ADOS_version']
    outDF.index = key
    
    outDF['iq1_measure'] = 'LEITER'
    outDF['iq2_measure'] = 'MullenScale'
    outDF['dataSet'] = datSet 
    outDF['eyes'] = 'open'
    outDF['task'] = 'rest'
    
    outDF.to_csv(r'Z:\User\pc1aod\CODE\GEDbounds_clusterImprove\biomarkDevDat.csv')
    
    
elif datSet == 'bpSZ': 
    #subset to only _r.mat files for resting state data, also drop to only one per subject
    df = df.loc[[('Eyes' in file) | ('eyes' in file) for file in df['comments_misc']], :]
    # add in diagnosis info
    diagInfo = pd.read_csv(metaFiles[['iec01' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo = diagInfo.loc[:,['subjectkey', 'study']]
    df = pd.merge(df, diagInfo, on = ['subjectkey'], how = 'inner')
    df = df.loc[pd.notnull(df['data_file1']), :]
    df = df.drop_duplicates()
    # add in intelligence info wechsler reading test
    intel1 = pd.read_csv(metaFiles[['wrat401' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel1 = intel1.loc[:,['subjectkey', 'wr_standardscore']]
    intel1 = intel1.drop_duplicates(subset = ['subjectkey'])
    df = pd.merge(df, intel1, on = ['subjectkey'], how = 'inner')
    df = df.loc[pd.notnull(df['data_file1']), :]
    df = df.drop_duplicates()
    # add in intelligence info wechsler spatial span scaled score
    intel2 = pd.read_csv(metaFiles[['wms_3_adult02' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel2 = intel2.loc[:,['subjectkey', 'wms3_ss_tss']]
    intel2 = intel2.drop_duplicates(subset = ['subjectkey'])
    df = pd.merge(df, intel2, on = ['subjectkey'], how = 'inner')
    df = df.loc[pd.notnull(df['data_file1']), :]
    df = df.drop_duplicates()
    #make the output dataframe
    
    key = df['subjectkey']
    sex = df['sex']
    age= df['interview_age']
    file= df['data_file1']
    task= df['comments_misc']
    iq1 = df['wr_standardscore']
    iq2 = df['wms3_ss_tss']
    diag = df['study']
    diag[['Relative of Proband with Psychotic' in name for name in diag]] = 'BPrel'
    diag[['Healthy' in name for name in diag]] = 'CON'
    diag[['Relative of Proband with Schizoaffect' in name for name in diag]] = 'SArel'
    diag[['Relative of Proband with Schizophren' in name for name in diag]] = 'SZrel'
    diag[['Psychotic' in name for name in diag]] = 'BP'
    diag[['Schizophreni' in name for name in diag]] = 'SZ'
    diag[['Schizo' in name for name in diag]] = 'SA'
    
    outDF = pd.concat([sex,age,task,diag,iq1,iq2,file], axis = 1)
    outDF.columns = ['sex', 'age', 'task', 'diag', 'iq1', 'iq2', 'file' ]
    outDF['iq1_measure'] = 'WechslerReading'
    outDF['iq2_measure'] = 'WeschlerSpatialSpan'
    outDF['ADOS'] = 999
    outDF['ADOS_version'] = 'NA'
    outDF['dataSet'] = 'bpSZ' 
   
    outDF['eyes'] = 'open'
    outDF['eyes'].loc[['Closed' in cnd for cnd in outDF['task']]] = 'closed'
    outDF.index = key
    outDF.to_csv(r'Z:\User\pc1aod\CODE\GEDbounds_clusterImprove\bpSZDat.csv')
    
elif datSet == 'femaleASD': 
    #subset to only _r.mat files for resting state data, also drop to only one per subject
    df = df.loc[[('Eyes' in file) | ('eyes' in file) for file in df['experiment_notes']], :]
    df = df.drop_duplicates(subset = ['subjectkey'])

    #ADOS info
    diagInfo4 = pd.read_csv(metaFiles[['ados4' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo4 = diagInfo4.loc[:,['subjectkey', 'scoresumm_abtotal']]
    diagInfo4['ADOS_version'] = 4
    diagInfo3 = pd.read_csv(metaFiles[['ados3' in file for file in metaFiles]].iloc[0], sep = '\t')
    diagInfo3 = diagInfo3.loc[:,['subjectkey', 'scoresumm_overalltotal']]
    diagInfo3['ADOS_version'] = 3
    diagOut = pd.merge(diagInfo3, diagInfo4, on = ['subjectkey'], how = 'outer')
    scoreCols = diagOut.columns[['scoresumm' in colName for colName in  diagOut.columns]]
    versionOut = diagOut.loc[:,['ADOS_version' in colName for colName in  diagOut.columns]].apply(np.nanmin, axis = 1)
    diagOut['ADOS_final'] = [diagOut[scoreCols[v-3]].iloc[ii] for ii,v in enumerate(versionOut)]
    diagOut['ADOS_version_final'] = versionOut
    diagOut = diagOut.loc[:,['subjectkey','ADOS_final','ADOS_version_final']]   
    diagOut = diagOut.drop_duplicates(subset = ['subjectkey'])  
    diagOut['ADOS_final'].iloc[0] = 0
    diagOut['ADOS_final'] = diagOut['ADOS_final'].astype('int')   
    #many controls didn't do the ADOS, so diagnosis for them is qualitative 
    conDiag = pd.read_csv(metaFiles[['ndar' in file for file in metaFiles]].iloc[0], sep = '\t')
    conDiag['diag2'] = 'CON'
    conDiag['diag2'].loc[[('Proband' in x) | ('999' in x) | ('ASD' in x) for x in conDiag['phenotype']]] = '999'
    conDiag = conDiag.loc[:, ['subjectkey', 'diag2']]
    conDiag = conDiag.drop_duplicates(subset = ['subjectkey'])
    
    diagOut = pd.merge(diagOut,conDiag, on = ['subjectkey'], how = 'outer')              
    df = pd.merge(df, diagOut, on = ['subjectkey'], how = 'inner')
    df = df.drop_duplicates()
    df['interview_age'] = df['interview_age'].astype('int')
    df['diag'] = df['diag2']
    df['diag'].iloc[np.where( (df['ADOS_final']>6) & (df['ADOS_version_final']==3) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>8) & (df['ADOS_version_final']==3) )] = 'AD'
    df['diag'].iloc[np.where( (df['ADOS_final']>6) & (df['ADOS_version_final']==4) )] = 'ASD'
    df['diag'].iloc[np.where( (df['ADOS_final']>9) & (df['ADOS_version_final']==4) )] = 'AD'
    df = df.loc[df['diag'] != '999',:]

    #intelligence
    intel2 = pd.read_csv(metaFiles[['das_ii_school' in file for file in metaFiles]].iloc[0], sep = '\t')
    intel2 = intel2.loc[:,['subjectkey', 'dasii_sar_gca_ss']]
    intel2 = intel2.drop_duplicates(subset = ['subjectkey']) 
    df = pd.merge(df, intel2, on = ['subjectkey'], how = 'inner')
    df['dasii_sar_gca_ss'] = df['dasii_sar_gca_ss'].astype('float')
    
    #make the output dataframe
    
    key = df['subjectkey']
    sex = df['sex']
    age= df['interview_age']
    file= df['data_file1']
    diag = df['diag']
    iq1 = df['dasii_sar_gca_ss']
    iq2 = df['dasii_sar_gca_ss']
    task = diag
    ADOS = df['ADOS_final']
    ADOS_version = df['ADOS_version_final']
    
    outDF = pd.concat([sex,age,task,diag,iq1,iq2,file, ADOS, ADOS_version], axis = 1)
    outDF.columns = ['sex', 'age', 'task', 'diag', 'iq1', 'iq2', 'file' , 'ADOS', 'ADOS_version']
    outDF.index = key
    
    outDF['iq1_measure'] = 'das_school'
    outDF['iq2_measure'] = 'das_school'
    outDF['dataSet'] = datSet 
    outDF['eyes'] = 'open'
    outDF['task'] = 'rest'
    
    outDF.to_csv(r'Z:\User\pc1aod\CODE\GEDbounds_clusterImprove\femaleASDDat.csv')



