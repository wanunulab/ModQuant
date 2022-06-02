#!/usr/bin/env python3
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
from sklearn.preprocessing import OrdinalEncoder
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

class ModelTrain:
    ''' 
    Assign, train, and test different supervised machine learning models.
    '''
    def __init__ (self, file_traing_dataframe, model_type, normalize_data, train_test_ratio, excluded_features, class_ratio):    
        '''Contructor: saves attribute fname '''
        self.training_dataframe = pd.read_pickle(file_traing_dataframe) #unpickled pandas dataframe
        self.model_type = model_type #default will be gradient boosting
        self.normalize_data = normalize_data #default is set to True
        self.train_test_ratio = train_test_ratio #default will be 75-25
        self.excluded_features = excluded_features #list of features to drop prior to training model, set to None
        self.class_ratio = class_ratio #automatically sets the dataframe to have a 50-50 sample size for each class (mod and unmod)
    
    def numericalBases (self):
        '''Convert basecalls to numerical values.'''
        #le = preprocessing.LabelEncoder()
        base_to_num_map = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
        substring = "B"
        updated_training_dataframe = self.training_dataframe
        for (columnName, columnData) in updated_training_dataframe.iteritems():
            if substring in columnName:
            	#print("Hit Base Change")
            	updated_training_dataframe = updated_training_dataframe.replace({columnName:base_to_num_map})
        return updated_training_dataframe
    
    def loadModel (self):
        '''Load model.'''
        if self.model_type == "gbc":
        	model = GradientBoostingClassifier(n_estimators=200, max_depth=5)
        elif self.model_type == "rfc":
        	model = RandomForestClassifier(n_estimators=200, max_depth=5)
        elif self.model_type == "logisticReg":
        	model = LogisticRegression(solver='lbfgs')
        elif self.model_type == "knn":
        	model = KNeighborsClassifier(n_neighbors=5)
        elif self.model_type == "svm":
        	model = SVC(decision_function_shape='ovr', probability=True, kernel='rbf', C=200, verbose=True, gamma='scale', class_weight=None)
        return model

    def normalizeData (self, training_data, testing_data):
        '''Normalize input data.'''
        scaler = sk.preprocessing.StandardScaler()
        x_train = scaler.fit_transform( training_data )
        x_test = scaler.transform( testing_data )
        return x_train, x_test, scaler

    def classRatio (self, input_dataframe):
        '''Sample size of different classes in input data.'''
        #self.class_ratio[0] = mod_fraction
        #self.class_ratio[1] = unmod_fraction
        ###Only does a 50-50 sample balance for now###
        n_mod = len(input_dataframe.loc[input_dataframe["type"] == "modified"])
        n_unmod = len(input_dataframe.loc[input_dataframe["type"] == "unmodified"])
        if n_unmod >= n_mod:
        	resamp_unmod = input_dataframe.loc[input_dataframe["type"] == "unmodified"].sample(n_mod)
        	resamp_mod = input_dataframe.loc[input_dataframe["type"] == "modified"].sample(n_mod)
        else:
        	resamp_unmod = input_dataframe.loc[input_dataframe["type"] == "unmodified"].sample(n_unmod)
        	resamp_mod = input_dataframe.loc[input_dataframe["type"] == "modified"].sample(n_unmod)

        updated_training_dataframe = pd.concat([resamp_unmod, resamp_mod])

        return updated_training_dataframe

    def removeFeatures (self, input_dataframe):
        '''Feature dropout function.'''
        updated_training_dataframe = input_dataframe.drop(columns = self.excluded_features)
        return updated_training_dataframe


    def train_test_model (self):
        '''Feature dropout function.'''
        tt_data = self.numericalBases()
        tt_data = self.classRatio(tt_data)
      
        if self.excluded_features == '0':
        	print("No features will be dropped out from the dataframe prior to model training.")
        else:
        	tt_data = self.removeFeatures(tt_data)
        feature_space = []
        removal_ID = "read_ID"
        remove_samples = "raw_current"
        remove_type = "type"

        for (columnName, columnData) in tt_data.iteritems():
        	if removal_ID in columnName or remove_samples in columnName or remove_type in columnName:
        		pass
        	else:
        		feature_space.append(columnName)
        x = tt_data.loc[:, feature_space]
        y = tt_data.loc[:, ["type"]]
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=self.train_test_ratio/100)#, random_state=i+70)
        if self.normalize_data == 'T':
        	x_train, x_test, scaler = self.normalizeData(x_train, x_test)
        model = self.loadModel()
        model.fit(x_train, y_train.to_numpy().flatten())
        model_prediction = model.predict(x_test)
        return len(tt_data), len(x_train), len(x_test), feature_space, accuracy_score(y_test, model_prediction)



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f' ,'--file', nargs = '?', type = str, required = True,
                        help = 'Path to prepared dataframe.')
    parser.add_argument('-m' ,'--model', nargs = '?', type = str, default = 'gbc', const = 'gbc', choices = ['gbc', 'rfc', 'knn', 'svm', 'logisticReg'],
                        help = 'ML model to fit your data to (default is gradient boosting).')
    parser.add_argument('-n', '--normalize', nargs = '?', type = str, default = 'T', const = 'T', choices = ['T', 'F'],
                        help = 'Do you wish to normalize your data prior to model training and testing?')
    parser.add_argument('-t', '--testsize', nargs = '?', type = int, default = 25, const = 25, choices = [5, 10, 15, 20, 25, 30, 35, 40],
                        help = 'Fraction of data used for testing set?')
    parser.add_argument('-e', '--exclude', nargs = '*', type = str, default = '0',
                        help = 'What features do you wish to drop from your data?')
    parser.add_argument('-r', '--ratio', nargs = '?', type = int, default = 50,
                        help = 'Set a class ratio (set to 50-50 by default, Not currently available)?')
    
    args = parser.parse_args()
    main_prep = ModelTrain(args.file, args.model, args.normalize, args.testsize, args.exclude, args.ratio)
    dataframe_size, training_size, testing_size, features, model_accuracy = main_prep.train_test_model()
    print(f"Size of your dataframe used to fit model prior splitting into a training and testing set: {dataframe_size}")
    print(f"Size of your training set: {training_size}")
    print(f"Size of your testing set: {testing_size}")
    print(f"Features used to train model: {features}")
    print(f"Accuracy of your model: {round(model_accuracy, 3)}")

if __name__ == '__main__':
	main()




        






#if not self.excluded_features:
#        	print("No features will be dropped out from the dataframe prior to model training.")
#        	return False