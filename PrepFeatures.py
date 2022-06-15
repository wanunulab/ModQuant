#!/usr/bin/env python3
#Similar to hypermodDistance
import sys
import os
import csv
import re
import itertools
import argparse
import numpy as np
import pandas as pd
import time
import datetime
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from scipy import signal
from scipy import fft as spfft
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import Normalizer
from sklearn.model_selection import cross_val_score
from sklearn.decomposition import PCA
from sklearn.metrics import plot_confusion_matrix
from sklearn.ensemble import StackingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier,ExtraTreesClassifier,AdaBoostClassifier,GradientBoostingClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from matplotlib import colors
import sklearn as sk
from sklearn.model_selection import train_test_split

print('The python version is {}.'.format(sys.version))
print('The scikit-learn version is {}.'.format(sk.__version__))
print('The numpy version is {}.'.format(np.__version__))
print('The pandas version is {}.'.format(pd.__version__))

class TrainingPrep:
    ''' 
    Prepare labeled training data, one with modified reads, 
    and the other with control unmodified reads. Both contain
    the same feature space and are prepared for ML model training.
    '''
    def __init__ (self, fname_control, fname_modified, control_label, modified_label, fourierCoeffs):    
        '''Contructor: saves attribute fname '''
        self.fname_control = fname_control #can be a list of files
        self.fname_modified = fname_modified #can be a list of files
        self.modifiedLabel = modified_label
        self.controlLabel = control_label
        self.fourierCoeffs = int(fourierCoeffs)
        
    def openf (self):
        '''Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname) 

    def numericalBases (self, training_dataframe):
        '''Convert basecalls to numerical values.'''
        #le = preprocessing.LabelEncoder()
        base_to_num_map = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
        substring = "B"
        #updated_training_dataframe = self.training_dataframe
        for (columnName, columnData) in training_dataframe.iteritems():
            if substring in columnName:
                #print("Hit Base Change")
                training_dataframe = training_dataframe.replace({columnName:base_to_num_map})
        return training_dataframe
    
    def preparetraining (self):
        '''Convert training data into pandas dataframe from csv files'''
        control_file_absolute_path_first = os.path.abspath(self.fname_control[0])
        print(f'Control file path: {control_file_absolute_path_first}')
        #Control_features_df = pd.read_csv(self.fname_control[0])
        Control_features_df = pd.read_csv(control_file_absolute_path_first)
        if len(self.fname_control) > 1:
            for i in range(1, len(self.fname_control)):
                additional_control_path = os.path.abspath(self.fname_control[i])
                additional_control_df = pd.read_csv(additional_control_path)
                #additional_control_df = pd.read_csv(self.fname_control[i])
                ##Control_features_df = Control_features_df.append(additional_control_df, ignore_index = True)
                Control_features_df = pd.concat([Control_features_df, additional_control_df])
        Control_features_df = Control_features_df.dropna()

        modified_file_absolute_path_first = os.path.abspath(self.fname_modified[0])
        print(f'Modified file path: {modified_file_absolute_path_first}')
        #Modified_features_df = pd.read_csv(self.fname_modified[0])
        Modified_features_df = pd.read_csv(modified_file_absolute_path_first)
        if len(self.fname_modified) > 1:
            for i in range(1, len(self.fname_modified)):
                additional_modified_path = os.path.abspath(self.fname_modified[i])
                additional_modified_df = pd.read_csv(additional_modified_path)
                #additional_modified_df = pd.read_csv(self.fname_modified[i])
                #Modified_features_df = Modified_features_df.append(additional_control_df, ignore_index = True)
                Modified_features_df = pd.concat([Modified_features_df, additional_modified_df])
        Modified_features_df = Modified_features_df.dropna()

        control_type_label = [self.controlLabel] * len(Control_features_df)
        mod_type_label = [self.modifiedLabel] * len(Modified_features_df)
        Control_features_df["type"] = control_type_label
        Modified_features_df["type"] = mod_type_label

        #frames = [Control_features_df, Modified_features_df]
        frames = [Modified_features_df, Control_features_df]
        combined_df = pd.concat(frames, ignore_index=True)
        prior_NA = len(combined_df)
        #Combined_df = Combined_df.dropna()
        Combined_df = self.numericalBases(combined_df) #Basecalls should be numeric values instead of A,C,G,T,- 

        #substring = "raw.current.samples"
        substring = "raw"
        for (columnName, columnData) in Combined_df.iteritems():
            if substring in columnName:
                Combined_df[columnName] = Combined_df[columnName].str.split(",")
                for idx, current in enumerate(Combined_df[columnName]):
                    Combined_df[columnName][idx] = list(np.float_(current))

        post_NA = len(Combined_df)
        return Combined_df, prior_NA, post_NA

    def fourierCoefficients (self, events, cutoff_freq = None):
        n = self.fourierCoeffs
        removed = 0
        allfourierCoeffs = []
        allfourierPhasors = []
        if n == -1:
            for event in events:
                event_fourierCoeffs = np.fft.fft(event)
                allfourierCoeffs.append(event_fourierCoeffs)
        else:
            for event in events:
                if (len(event) % 2 == 0 and (len(event)/2) + 1 >= n) or (len(event) % 2 != 0 and (len(event)+1)/2 >= n):
                    cfft=np.fft.rfft(event)/len(event)
                    fft_phase = np.angle(cfft)
                    cfft=np.real(cfft) #keeping the sign of the magnitude
                    event_fourierCoeffs=np.array(cfft, copy=True)
                    event_fourierPhasors=np.array(fft_phase, copy=True)    
                    allfourierCoeffs.append(event_fourierCoeffs)
                    allfourierPhasors.append(event_fourierPhasors) #not using currently
                else:
                    removed += 1
                    event_fourierCoeffs = [float("NaN")]
                    event_fourierPhasors = [float("NaN")]
                    allfourierCoeffs.append(event_fourierCoeffs)
                    allfourierPhasors.append(event_fourierPhasors)
        return allfourierCoeffs, removed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-cf' ,'--fnameControl', nargs = '+', type = str, required = True,
                        help = 'Path to control file you want to prepare.')
    parser.add_argument('-mf' ,'--fnameModified', nargs = '+', type = str, required = True,
                        help = 'Path to control file you want to prepare.')
    parser.add_argument('-lc', '--controlLabel', nargs = '?', type = str, default = 'unmodified',
                        help = 'What do you wish to label modified training data?')
    parser.add_argument('-lm', '--modifiedLabel', nargs = '?', type = str, default = 'modified',
                        help = 'What do you wish to label modified training data?')
    parser.add_argument('-fc', '--FourierCoeffs', nargs = '?', type = str, default = '3',
                        help = 'What Fourier coefficients do you want to extract from the signal?')
    parser.add_argument('-o', '--outputDataframe', nargs = '?', type = str, default = 'training_data',
                        help = 'What directory and name do you want for your prepared output training data (pandas datframe)?')
    
    args = parser.parse_args()
    main_prep = TrainingPrep(args.fnameControl, args.fnameModified, args.controlLabel, args.modifiedLabel, args.FourierCoeffs)
    #main_prep = TrainingPrep(control_file_absolute_path, modified_file_absolute_path, args.controlLabel, args.modifiedLabel, args.FourierCoeffs)
    
    training_dataframe, priorNANs, postNANs = main_prep.preparetraining()
    print(len(training_dataframe))
    #print(f"Size of original dataframe before NaN removal: {priorNANs}")
    print(f"Size of original dataframe after NaN removal: {postNANs}")

    df_NaN = training_dataframe[training_dataframe.isna().any(axis=1)]

    print(f"Size of NaN portion of training dataframe: {len(df_NaN)}")


    ##########Adding Fourier coefficients to final dataframe#################
    substring = "raw"
    #substring = "raw.current.samples"

    for (columnName, columnData) in training_dataframe.iteritems():
        if substring in columnName:
            print(f"Preparing Fourier coefficients for {columnName} samples")
            #print(training_dataframe[columnName])
            Current_Fourier_Coefficients, removed_events = main_prep.fourierCoefficients(events = training_dataframe[columnName], cutoff_freq = None)
            #print(len(Current_Fourier_Coefficients))
            try:
                fft_col_name = columnName.split("_")
            except ValueError:
                fft_col_name = columnName.split(".")

            if "upstream" in fft_col_name:
                fft_col_name = f"kmer_{fft_col_name[0]}_{fft_col_name[-1]}"
            else:
                fft_col_name = f"kmer_{fft_col_name[-1]}"
            Coeffs_df = pd.DataFrame()
            for event_Coeffs in Current_Fourier_Coefficients:
                d = {}
                for i in range(1, int(args.FourierCoeffs)): #delete 1, to preserve first coefficient which equals current mean
                    if not np.isnan(event_Coeffs[0]):
                        d[f"fc{i+1}_{fft_col_name}"] = round(np.absolute(event_Coeffs[i]), 2)
                    else:
                        d[f"fc{i+1}_{fft_col_name}"] = float("NaN")
                #Coeffs_df = Coeffs_df.append(d, ignore_index = True) #pd.append is deprecated
                Coeffs_df = pd.concat([Coeffs_df, pd.DataFrame([d])], ignore_index = True)
            Coeffs_df.reset_index(drop = True, inplace = True)
            #print(Coeffs_df)
            training_dataframe.reset_index(drop = True, inplace = True)
            training_dataframe = pd.concat( [training_dataframe, Coeffs_df], axis = 1 )

    training_dataframe = training_dataframe.dropna()
    print(f"Number of dropped events that did not meet the minimum FC value: {removed_events}")
    print(training_dataframe)

    now = datetime.datetime.now()
    training_data_directory = "prepared_training_data"
    if not os.path.isdir(training_data_directory):
        os.system( f"mkdir {training_data_directory}" )
    os.chdir(training_data_directory)
    training_dataframe.to_pickle(f"{args.outputDataframe}_{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}_{now.second}.pkl")


if __name__ == '__main__':
    main()

