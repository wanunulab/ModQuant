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
    def __init__ (self, fname_direct, direct_label, fourierCoeffs):    
        '''Contructor: saves attribute fname '''
        self.fname_direct = fname_direct #can be a list of files
        self.directLabel = direct_label
        self.fourierCoeffs = int(fourierCoeffs)
        
    def openf (self):
        '''Handle file opens, allowing STDIN.'''
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname) 

    def numericalBases (self, direct_dataframe:pd.DataFrame):
        '''Convert basecalls to numerical values.'''
        #le = preprocessing.LabelEncoder()
        base_to_num_map = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
        substring = "B"
        #updated_training_dataframe = self.training_dataframe
        for (columnName, columnData) in direct_dataframe.iteritems():
            if substring in columnName:
                #print("Hit Base Change")
                direct_dataframe = direct_dataframe.replace({columnName:base_to_num_map})
        return direct_dataframe

    def quartileCompute (self, kmers_):
        '''Compute and add first quartile (25%), median (50%), and third quartile (75%) ionic current for every K-mer'''
        #event.first_quartiles=[]
        #event.medians=[]
        #event.third_quartiles=[]
        quartiles = []
        for kmer_current in kmers_:
            #print(f"Kmer Current: {kmer_current}")
            #print(f"Sample type: {type(kmer_current[0])}")
            quant_array = np.quantile(kmer_current, [0.25, 0.5, 0.75], interpolation="midpoint")
            quartiles.append([quant_array[0], quant_array[1], quant_array[2]])
        return quartiles
            #event.first_quartiles.append(quant_array[0])
            #event.medians.append(quant_array[1])
            #event.third_quartiles.append(quant_array[2])

    def quantileSignal (self, direct_dataframe:pd.DataFrame):
        '''Compute and add first quartile (25%), median (50%), and third quartile (75%) ionic current for every K-mer'''
        substring = "raw"
        for (columnName, columnData) in direct_dataframe.iteritems():
            if substring in columnName:
                print(f"Preparing quartile current values for {columnName} samples")
                #print(direct_dataframe[columnName])
                kmer_current_quant = self.quartileCompute(kmers_ = direct_dataframe[columnName])#, cutoff_freq = None)
                try:
                    quant_col_name = columnName.split(".")
                except ValueError:
                    quant_col_name = columnName.split("_")

                if "upstream" in quant_col_name:
                    quant_col_name = f"kmer_{quant_col_name[0]}_{quant_col_name[-1]}"
                else:
                    quant_col_name = f"kmer_{quant_col_name[-1]}"
                Quant_df = pd.DataFrame()
                for kmer_quant in kmer_current_quant:
                    d = {}
                    for i in range(0, len(kmer_quant)): 
                        d[f"Quartile_{i+1}_{quant_col_name}"] = round(kmer_quant[i], 3)
                    #else:
                    #    d[f"fc{i+1}_{fft_col_name}"] = float("NaN")
                    #Coeffs_df = Coeffs_df.append(d, ignore_index = True) #pd.append is deprecated
                    Quant_df = pd.concat([Quant_df, pd.DataFrame([d])], ignore_index = True)
                Quant_df.reset_index(drop = True, inplace = True)
                #print(Coeffs_df)
                direct_dataframe.reset_index(drop = True, inplace = True)
                direct_dataframe = pd.concat( [direct_dataframe, Quant_df], axis = 1 )

        return direct_dataframe

    
    def preparedirect (self):
        '''Convert training data into pandas dataframe from csv files'''
        direct_file_absolute_path_first = os.path.abspath(self.fname_direct[0])
        print(f'Direct file path: {direct_file_absolute_path_first}')
        #Control_features_df = pd.read_csv(self.fname_control[0])
        Direct_features_df = pd.read_csv(direct_file_absolute_path_first)
        if len(self.fname_direct) > 1:
            for i in range(1, len(self.fname_control)):
                additional_direct_path = os.path.abspath(self.fname_direct[i])
                additional_direct_df = pd.read_csv(additional_direct_path)
                Direct_features_df = pd.concat([Direct_features_df, additional_direct_df])
        print(f"Size of Direct dataframe prior to dropping NaN reads: {len(Direct_features_df)}")
        Direct_features_df = Direct_features_df.dropna()
        print(f"Size of Direct dataframe prior after dropping NaN reads: {len(Direct_features_df)}")


        direct_type_label = [self.directLabel] * len(Direct_features_df)
        
        Direct_features_df["type"] = direct_type_label
        

        
        ##frames = [Modified_features_df, Control_features_df]
        ##combined_df = pd.concat(frames, ignore_index=True)
        #Combined_df = Combined_df.dropna()
        Combined_df = self.numericalBases(Direct_features_df) #Basecalls should be numeric values instead of A,C,G,T,- 
        Combined_df = Combined_df.reset_index()

        #substring = "raw.current.samples"
        substring = "raw"
        for (columnName, columnData) in Combined_df.iteritems():
            if substring in columnName:
                Combined_df[columnName] = Combined_df[columnName].str.split(",")
                Combined_df[columnName] = Combined_df[columnName].apply(lambda x: [float(i) for i in x])
                ##for idx, current in enumerate(Combined_df[columnName]):
                    ##print(idx)
                    ##print(current)
                    
                    #Combined_df[columnName][idx] = list(np.float_(current))
                    #print("after")
                    #print(current)
                    #print(Combined_df[columnName][idx])


        return Combined_df

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
    parser.add_argument('-df' ,'--fnameDirect', nargs = '+', type = str, required = True,
                        help = 'Path to control file you want to prepare.')
    parser.add_argument('-dl', '--directLabel', nargs = '?', type = str, default = 'unspecified',
                        help = 'What do you wish to label modified training data?')
    parser.add_argument('-fc', '--FourierCoeffs', nargs = '?', type = str, default = '0',
                        help = 'What Fourier coefficients do you want to extract from the signal?')
    parser.add_argument('-o', '--outputDataframe', nargs = '?', type = str, default = 'training_data',
                        help = 'What directory and name do you want for your prepared output training data (pandas datframe)?')
    
    args = parser.parse_args()
    main_prep = TrainingPrep(args.fnameDirect, args.directLabel, args.FourierCoeffs)
    #main_prep = TrainingPrep(control_file_absolute_path, modified_file_absolute_path, args.controlLabel, args.modifiedLabel, args.FourierCoeffs)
    
    direct_dataframe = main_prep.preparedirect()
    print(direct_dataframe)
    #print(f"Size of original dataframe before NaN removal: {priorNANs}")
    print(f"Size of original dataframe after NaN removal: {len(direct_dataframe)}")

    df_NaN = direct_dataframe[direct_dataframe.isna().any(axis=1)]

    print(f"Size of NaN portion of direct dataframe: {len(df_NaN)}")


    ##########Adding Fourier coefficients to final dataframe#################
    substring = "raw"
    #substring = "raw.current.samples"
    if int(args.FourierCoeffs) > 0:
        print("Adding Fourier Coefficients to the feature space")
        for (columnName, columnData) in direct_dataframe.iteritems():
            if substring in columnName:
                print(f"Preparing Fourier coefficients for {columnName} samples")
                #print(training_dataframe[columnName])
                Current_Fourier_Coefficients, removed_events = main_prep.fourierCoefficients(events = direct_dataframe[columnName], cutoff_freq = None)
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
                direct_dataframe.reset_index(drop = True, inplace = True)
                direct_dataframe = pd.concat( [direct_dataframe, Coeffs_df], axis = 1 )

    direct_dataframe = main_prep.quantileSignal(direct_dataframe) #adding quartiles signal
    direct_dataframe = direct_dataframe.dropna()
    #print(f"Number of dropped events that did not meet the minimum FC value: {removed_events}")
    print(direct_dataframe)

    now = datetime.datetime.now()
    direct_data_directory = "prepared_direct_data"
    #direct_data_directory = "prepared_SGNEx_data"
    if not os.path.isdir(direct_data_directory):
        os.system( f"mkdir {direct_data_directory}" )
    os.chdir(direct_data_directory)
    direct_dataframe.to_pickle(f"{args.outputDataframe}_{now.month}{now.day}{now.year}_{now.hour}{now.minute}{now.second}.pkl")


if __name__ == '__main__':
    main()

