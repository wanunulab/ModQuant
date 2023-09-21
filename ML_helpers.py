#!/usr/bin/env python3
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import sys
import os
import pandas as pd
from sklearn import preprocessing
import sklearn as sk
import scipy
print('The numpy version is {}.'.format(np.__version__))
print('The scikit-learn version is {}.'.format(sk.__version__))
print('The pandas version is {}.'.format(pd.__version__))
print('The scipy version is {}.'.format(scipy.__version__))
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier,ExtraTreesClassifier,AdaBoostClassifier,GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, classification_report, confusion_matrix
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import RocCurveDisplay

import seaborn as sns
from matplotlib import gridspec
import matplotlib.patches as mpatches
import matplotlib.font_manager
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from matplotlib import colors


print('The python version is {}.'.format(sys.version))
print('The scikit-learn version is {}.'.format(sk.__version__))
print('The numpy version is {}.'.format(np.__version__))
print('The pandas version is {}.'.format(pd.__version__))

class ModelLoad:
    '''Load ML model.'''
    def __init__(self,model_type,seed_model:int):
        '''
        Initialize the ModelLoad class.

        Parameters:
        - model_type (str): The type of ML model to load.
        - seed_model (int): The seed value for randomization (if applicable).
        '''
        self.model_type=model_type
        self.seed_model=seed_model
    
    def load_model(self):
        '''
        Load and return the specified ML model.

        Returns:
        - model: The loaded ML model based on the specified model_type.
        '''
        if self.model_type=="gbc":
            model=GradientBoostingClassifier(n_estimators=200, max_depth=5, random_state=self.seed_model)
        elif self.model_type=="rfc":
            model=RandomForestClassifier(n_estimators=200, max_depth=5, random_state = self.seed_model)
        elif self.model_type == "logisticReg":
            model=LogisticRegression(solver='lbfgs', random_state = self.seed_model)
        elif self.model_type=="knn":
            model=KNeighborsClassifier(n_neighbors=5)
        elif self.model_type=="svm":
            model=SVC(decision_function_shape='ovr',probability=True,kernel='rbf',C=200,verbose=True,gamma='scale',class_weight=None,random_state=self.seed_model)
        else:
            model=self.model_type
        return model
    

class ModelFitTest:
    '''Fit and train model'''
    def __init__(self,x_train:pd.DataFrame,x_test:pd.DataFrame,y_train:pd.DataFrame,y_test:pd.DataFrame,
                 labels,seed_sample:int,model_type:str,seed_model:int,test_size=0.2):
        self.x_train=x_train
        self.x_test=x_test
        self.y_train=y_train
        self.y_test=y_test
        self.labels=labels
        self.seed_sample=seed_sample
        self.model_type=model_type
        self.seed_model=seed_model
        self.test_size=test_size

    def model_fit(self):
        '''
        Train and return the specified ML model.

        Returns:
        - trained_model: Dictionary containing the trained model and related information.
        '''
        # Scale the parameters using standard scaler, fitted to the train data statistics and applied to both train and test sets
        scaler=sk.preprocessing.StandardScaler() #initialize normalization function
        x_train_norm=scaler.fit_transform(self.x_train) #normalize training data and save normalization parameters
        initializer=ModelLoad(self.model_type,self.seed_model)
        # define GBC model and train it
        model=initializer.load_model()
        model.fit(x_train_norm,self.y_train) #train GBC 
        trained_model={
            'model':model,
            'seed_model':self.seed_model,
            'scaler':scaler,
        }
        return trained_model
    
    def model_pred(self,model,scaler,x_test,y_test):
        '''
        Make predictions using the trained model and evaluate its performance.

        Parameters:
        - model: The trained ML model.
        - scaler: The scaler used for normalization.
        - x_test: Test data features.
        - y_test: Test data labels.

        Returns:
        - result: Dictionary containing model predictions and evaluation metrics.
        '''
        x_test_normalized=scaler.transform(x_test)
        # apply to test set to evaluate performance
        y_pred=model.predict(x_test_normalized) #GBC predictions
        labs,true=np.unique(y_test,return_counts=True)
        unique,counts=np.unique(y_pred, return_counts=True) #GBC prediction results
        x_test_results_df=self.output_test_df_results(x_test,y_test,y_pred)
        probability_scores=self.prediction_score(model,x_test_normalized)
        feature_imp=self.extract_feature_importance(model)

        
        result={
            'model':model,
            'seed_sample':self.seed_sample,
            'seed_model':self.seed_model,
            'scaler':scaler,
            'x_test':x_test,
            'x_test_normalized':x_test_normalized,
            'y_test':y_test,
            'y_pred':y_pred,
            'x_test_results_df':x_test_results_df,
            'counts_true':dict(zip(labs,true)),
            'counts_pred':dict(zip(unique,counts)),
            'accuracy':accuracy_score(y_test,y_pred),
            "probability_scores":probability_scores,
            'conf_mat':confusion_matrix(y_true=y_test,y_pred=y_pred,labels=self.labels,normalize="true"),
            'conf_mat_raw':confusion_matrix(y_true=y_test,y_pred=y_pred,labels=self.labels,normalize=None),
            'feature_imp':feature_imp
        }
        return result

    def model_fit_pred(self):
        '''
        Train the specified ML model and make predictions on the test set.

        Returns:
        - result: Dictionary containing trained model, predictions, and evaluation metrics.
        '''
        scaler=sk.preprocessing.StandardScaler() #initialize normalization function
        x_train_normalized=scaler.fit_transform(self.x_train)
        x_test_normalized=scaler.transform(self.x_test)
        initializer=ModelLoad(self.model_type,self.seed_model)
        model=initializer.load_model()
        model.fit(x_train_normalized,self.y_train)
        # apply to test set to evaluate performance
        y_pred=model.predict(x_test_normalized) #predict
        labs,true=np.unique(self.y_test,return_counts=True)
        unique,counts=np.unique(y_pred,return_counts=True)
        x_test_results_df=self.output_test_df_results(self.x_test,self.y_test,y_pred)
        probability_scores=self.prediction_score(model,x_test_normalized)
        feature_imp=self.extract_feature_importance(model)

        result={
        'model':model,
        'seed_sample':self.seed_sample,
        'seed_model':self.seed_model,
        'scaler':scaler,
        'test_size':self.test_size,
        'labels':self.labels,
        'x_train':self.x_train,
        'x_train_normalized':x_train_normalized,
        'x_test':self.x_test,
        'x_test_normalized':x_test_normalized,
        'y_train':self.y_train,
        'y_test':self.y_test,
        'y_pred':y_pred,
        'x_test_results_df':x_test_results_df,
        'counts_true':dict(zip(labs,true)),
        'counts_pred':dict(zip(unique,counts)),
        'accuracy':accuracy_score(self.y_test,y_pred),
        "probability_scores":probability_scores,
        'conf_mat':confusion_matrix(y_true=self.y_test,y_pred=y_pred,labels=self.labels,normalize="true"),
        'conf_mat_raw':confusion_matrix(y_true=self.y_test,y_pred=y_pred,labels=self.labels,normalize=None),
        'feature_imp':feature_imp
        }
        return result
        

    def output_test_df_results(self,x_test:pd.DataFrame,y_test:pd.DataFrame,y_pred):
        '''
        Create a DataFrame with test set results.

        Parameters:
        - x_test (pd.DataFrame): Test data features.
        - y_test (pd.DataFrame): Test data labels.
        - y_pred: Predicted labels.

        Returns:
        - x_test_results_df: DataFrame containing test set results.
        '''
        x_test_reset=x_test.reset_index(drop=True)
        y_test_reset=y_test.reset_index(drop=True)
        classification_results=[]
        for i,row in x_test_reset.iterrows():
            true_label=y_test_reset.loc[i]
            pred_label=y_pred[i]
            if true_label==pred_label:
                classification_results.append('Correctly Classified')
            else:
                classification_results.append('Misclassified')
        x_test_with_results=x_test_reset.copy()
        x_test_with_results['type']=y_test_reset
        x_test_with_results['Classification']=classification_results
        return x_test_with_results
    
    def prediction_score(self,model,x_test_normalized):
        '''
        Calculate prediction scores (if applicable).

        Parameters:
        - model: The trained ML model.
        - x_test_normalized: Normalized test data features.

        Returns:
        - probability_scores: Probability scores (if applicable, None otherwise).
        '''
        if self.model_type=='gbc':
            probability_scores=model.predict_proba(x_test_normalized)
        else:
            probability_scores=None
        return probability_scores
    
    def extract_feature_importance(self,model):
        if self.model_type=='gbc':
            ##feature_imp=pd.DataFrame(model.feature_importances_,index=df.drop('type',axis='columns').columns)
            feature_importance=pd.DataFrame(model.feature_importances_,index=self.x_train.columns)
        else:
            feature_importance=None
        return feature_importance
    

class IterateModel:
    def __init__(self,df:pd.DataFrame,feature_space,labels,model_type,seed=1234,test_size=0.2,n_splits=5,n_repeats=10,train_sample_sz=None):
        self.df=df
        self.feature_space=feature_space
        self.labels=labels
        self.model_type=model_type
        self.seed=seed
        self.test_size=test_size
        self.n_splits=n_splits
        self.n_repeats=n_repeats
        self.n=self.n_splits*self.n_repeats
        self.train_sample_sz=train_sample_sz

    def balance_dataframe(self):
        limiting_sample_size=min(self.df["type"].value_counts())
        df_balanced=self.df.groupby("type").sample(n=limiting_sample_size,random_state=self.seed)
        df_balanced.reset_index(drop=True,inplace=True)
        return df_balanced
    
    def model_seed_generator(self):
        rng=np.random.default_rng(self.seed)
        seeds_model=rng.integers(low=1000,high=10000,size=self.n)
        return seeds_model
    
    def sample_seed_generator(self):
        rng=np.random.default_rng(self.seed)
        seeds_sample=rng.integers(low=1000,high=10000,size=self.n*21)
        return seeds_sample
    
    def rsfk(self,df:pd.DataFrame,seeds_model):
        results=[]
        train_indices=[]
        test_indices=[]
        model_train_test_instances=[]
        rskf=RepeatedStratifiedKFold(n_splits=self.n_splits,n_repeats=self.n_repeats,random_state=self.seed)
        for ((train_index,test_index),seed_model) in zip(rskf.split(df.loc[:,self.feature_space],df["type"]),seeds_model):
            if self.train_sample_sz is not None:
                x_train_df=df.iloc[train_index][self.feature_space+["type"]].groupby("type").sample(n=self.train_sample_sz,random_state=self.seed)
            else:
                x_train_df=df.iloc[train_index][self.feature_space+["type"]]
            model_train_test=ModelFitTest(x_train_df.loc[:,self.feature_space],
                                          df.iloc[test_index][self.feature_space],
                                          x_train_df["type"],
                                          df.iloc[test_index]["type"],
                                          labels=self.labels,
                                          seed_sample=self.seed,
                                          model_type=self.model_type,
                                          seed_model=seed_model,
                                          test_size=0.2)
            model_train_test_instances.append(model_train_test)
            result=model_train_test.model_fit_pred()
            results.append(result)
            train_indices.append(train_index)
            test_indices.append(test_index)
            
        summary={"fit_test_instances":model_train_test_instances,
                 "iteration_results":results,
                 "train_indices":train_indices,
                 "test_indices":test_indices
                }
        return summary
    
    def feature_importance_prep(self,results):
        imp_iters=100*sum([result['feature_imp'] for result in results])/(len(results))
        imp_iters.columns=['Weight']
        imp_iters['Feature']=imp_iters.index
        imp_iters=imp_iters.reset_index(drop=True)
        return imp_iters
            
    def ratio_feature_condition_iter(self,test_df,true_positive_test_df,true_negative_test_df,feature_name,condition_lambda):
        if len(true_positive_test_df)>0:
            n_true_postive_condition=len(true_positive_test_df[true_positive_test_df[feature_name].apply(condition_lambda)])
            rate_true_postive_condition=n_true_postive_condition/len(test_df)
        else:
            n_true_postive_condition=0
            rate_true_postive_condition=0
        if len(true_negative_test_df)>0:
            n_true_negative_condition=len(true_negative_test_df[true_negative_test_df[feature_name].apply(condition_lambda)])
            rate_true_negative_condition=n_true_negative_condition/len(test_df)
        else:
            n_true_negative_condition=0
            rate_true_negative_condition=0
        return n_true_postive_condition,rate_true_postive_condition,n_true_negative_condition,rate_true_negative_condition
    
    def model_iter(self,balanced=True):
        seeds_model=self.model_seed_generator()
        if balanced==True:
            input_df=self.balance_dataframe()
        else:
            input_df=self.df
        summary=self.rsfk(input_df,seeds_model)
        
        # get the average of feature importances
        if self.model_type=='gbc': 
            imp_iters=self.feature_importance_prep(summary["iteration_results"])
        else:
            imp_iters=None
        
        summary["input_dataframe"]=input_df
        summary["seed_value"]=self.seed
        summary["avg_accuracy"]=np.mean([result['accuracy'] for result in summary["iteration_results"]])
        summary["std_accuracy"]=np.std([result['accuracy'] for result in summary["iteration_results"]])
        summary["avg_conf_mat_mean"]=np.mean([result['conf_mat'] for result in summary["iteration_results"]],axis=0)
        summary["avg_conf_mat_std"]=np.std([result['conf_mat'] for result in summary["iteration_results"]],axis=0)
        summary["avg_feature_importance"]=imp_iters
        return summary
    
    def ratio_model_iter(self,true_positive_label,true_negative_label,feature_name=None,condition_lambda=None):
        seeds_sample=self.sample_seed_generator()
        summary=self.model_iter()
        input_df=summary["input_dataframe"]
        ratio_sweep_results=[]

        for (trained_tested_model_object,train_index,test_index,trained_tested_model_result) in zip(summary["fit_test_instances"],summary["train_indices"],summary["test_indices"],summary["iteration_results"]):
            model=trained_tested_model_result["model"]
            scaler=trained_tested_model_result["scaler"]
            seed_model=trained_tested_model_result["seed_model"]
            ##test_set_sz=trained_tested_model_result["test_size"]
            test_set_sz=len(test_index)

            looper=np.arange(0.0,1.05,0.05) #TODO: make 0.05 step a variable that the user can define, and make the 1.05 endpoint 1.0000001 
            for (percentage1,percentage2,seed_sample) in zip(looper[0::1],looper[-1::-1],seeds_sample):
                percentage1=round(percentage1,4) #percent of psi (100% syn)
                percentage2=round(percentage2,4) #percent of control (0% syn)
                n_true_negative=int(percentage2*(test_set_sz/2)) #[100%, 95%, 90%,...,10%, 5%, 0%]
                n_true_positive=int(percentage1*(test_set_sz/2)) #[0%, 5%, 10%,...,90%, 95%, 100%]
                true_negative_test_df=input_df.iloc[test_index].loc[input_df.iloc[test_index]["type"]==true_negative_label].sample(n_true_negative,random_state=seed_sample)
                true_positive_test_df=input_df.iloc[test_index].loc[input_df.iloc[test_index]["type"]==true_positive_label].sample(n_true_positive,random_state=seed_sample) 
                test_df_updated=pd.concat([true_negative_test_df,true_positive_test_df])

                assert (feature_name is None and condition_lambda is None) or (feature_name is not None and condition_lambda is not None), "Both feature_name and condition_lambda must be provided together."
                
                if feature_name is not None and condition_lambda is not None:
                    n_true_postive_condition,rate_true_postive_condition,n_true_negative_condition,rate_true_negative_condition=self.ratio_feature_condition_iter(test_df_updated,
                                                                                                                                                                    true_positive_test_df,
                                                                                                                                                                    true_negative_test_df,
                                                                                                                                                                    feature_name,
                                                                                                                                                                    condition_lambda)
                else:
                    n_true_postive_condition=None
                    rate_true_postive_condition=None
                    n_true_negative_condition=None
                    rate_true_negative_condition=None

                result=trained_tested_model_object.model_pred(model,scaler,test_df_updated.loc[:,self.feature_space],test_df_updated["type"])
                result['true_positive_percentage']=round(percentage1*100, 2)
                result['true_positive_condition_sample_sz']=n_true_postive_condition
                result['true_positive_condition_rate']=rate_true_postive_condition
                result['true_negative_percentage']=round(percentage2*100,2)
                result['true_negative_condition_sample_sz']=n_true_negative_condition
                result['true_negative_condition_rate']=rate_true_negative_condition
                ratio_sweep_results.append(result)
        summary["ratio_sweep_results"]=ratio_sweep_results
        return summary



class ModelStats:
    def __init__(self,ratio_sweep_results:list,true_positive_label:str,true_negative_label:str):
        self.ratio_sweep_results=ratio_sweep_results
        self.true_positive_label=true_positive_label
        self.true_negative_label=true_negative_label

    def model_performance(self):
        #ntc,ptc,nfc,pfc,ntr,ptr,nfr,pfr
        key_names=['negative_true_calls','positive_true_calls','negative_false_calls','positive_false_calls',
                     'negative_true_rate','positive_true_rate','negative_false_rate','positive_false_rate']
        d={self.true_negative_label:{key:{'means':[],'std_devs':[]} for key in key_names},
           self.true_positive_label:{key:{'means':[],'std_devs':[]} for key in key_names}}
        unique_percentages=set(value['true_positive_percentage'] for value in self.ratio_sweep_results if 'true_positive_percentage' in value)
        unique_percentages_list=list(unique_percentages)
        negative_percentages_list=sorted(unique_percentages_list,reverse=True)
        positive_percentages_list=sorted(unique_percentages_list,reverse=False)
        
        for positive_percentage,negative_percentage in zip(positive_percentages_list,negative_percentages_list):
            positive_percentage_model_results=[result for result in self.ratio_sweep_results if 'true_positive_percentage' in result and result['true_positive_percentage']==positive_percentage]
            negative_percentage_model_results=[result for result in self.ratio_sweep_results if 'true_negative_percentage' in result and result['true_negative_percentage']==negative_percentage]
            positive_results_array=np.empty((0,8),dtype=float)
            negative_results_array=np.empty((0,8),dtype=float)
            for positive_ratio_pred,negative_ratio_pred in zip(positive_percentage_model_results,negative_percentage_model_results):
                positive_ratio_output=self.output_stats(positive_ratio_pred)
                positive_results_array=np.vstack((positive_results_array,positive_ratio_output))
                negative_ratio_output=self.output_stats(negative_ratio_pred)
                negative_results_array=np.vstack((negative_results_array,negative_ratio_output))
            positive_percentage_means=np.mean(positive_results_array,axis=0)
            positive_percentage_stds=np.std(positive_results_array,axis=0)
            for key,mean_value,std_value in zip(key_names,positive_percentage_means,positive_percentage_stds):
                d[self.true_positive_label][key]['means'].append(mean_value)
                d[self.true_positive_label][key]['std_devs'].append(std_value)
            negative_percentage_means=np.mean(negative_results_array,axis=0)
            negative_percentage_stds=np.std(negative_results_array,axis=0)
            for key,mean_value,std_value in zip(key_names,negative_percentage_means,negative_percentage_stds):
                d[self.true_negative_label][key]['means'].append(mean_value)
                d[self.true_negative_label][key]['std_devs'].append(std_value)
        return d

    def output_stats(self,ratio_pred):
        if self.true_positive_label in ratio_pred['counts_true']:
            n_true_positive=ratio_pred['counts_true'][self.true_positive_label]
        else:
            n_true_positive=0
        if self.true_negative_label in ratio_pred['counts_true']:
            n_true_negative=ratio_pred['counts_true'][self.true_negative_label]
        else:
            n_true_negative=0
        stat_results=self.stats_report(ratio_pred['conf_mat_raw'],
                                      ratio_pred['true_positive_percentage'], 
                                      ratio_pred['true_negative_percentage'], 
                                      n_true_positive, 
                                      n_true_negative)
        return stat_results

    
    def calculate_calls(self,true_negative_calls,true_positive_calls,false_negative_calls,false_positive_calls,n_true_negative,n_true_positive):
        '''Compute the True and False calls for unmodified and modified data in the test set'''
        ntc=round(true_negative_calls/(n_true_negative+n_true_positive),5)
        ptc=round(true_positive_calls/(n_true_negative+n_true_positive),5)
        nfc=round(false_negative_calls/(n_true_negative+n_true_positive),5)
        pfc=round(false_positive_calls/(n_true_negative+n_true_positive),5)
        return ntc,ptc,nfc,pfc
    
    def calculate_rates(self,true_negative_calls,true_positive_calls,false_negative_calls,false_positive_calls,n_true_negative,n_true_positive):
        '''Compute the True and False rates for unmodified and modified data in the test set'''
        if n_true_negative==0:
            ntr=0
            pfr=0
        else:
            ntr=round(true_negative_calls/n_true_negative,5)
            pfr=round(false_positive_calls/n_true_negative,5)
        if n_true_positive==0:
            ptr=0
            nfr=0
        else:
            ptr=round(true_positive_calls/n_true_positive,5)
            nfr=round(false_negative_calls/n_true_positive, 5)
        return ntr,ptr,nfr,pfr

    def stats_report (self,conf_mat,true_positive_percentage,true_negative_percentage,n_true_positive,n_true_negative):
        if len(conf_mat)==1 and n_true_positive==0:
            true_negative_calls=conf_mat[0][0] #Control reads that are correctly called
            true_positive_calls=0 #Modified reads that are correctly called
            false_negative_calls=0 #Modified reads miscalled as a control site
            false_positive_calls=0 #Control reads miscalled as a modified site
        elif len(conf_mat)==1 and n_true_negative==0:
            true_negative_calls=0 #Control reads that are correctly called
            true_positive_calls=conf_mat[0][0] #Modified reads that are correctly called
            false_negative_calls=0 #Modified reads miscalled as a control site
            false_positive_calls=0 #Control reads miscalled as a modified site
        else:
            true_negative_calls = conf_mat[0][0] #Control reads that are correctly called
            true_positive_calls = conf_mat[1][1] #Modified reads that are correctly called
            false_negative_calls = conf_mat[1][0] #Modified reads miscalled as a control site
            false_positive_calls = conf_mat[0][1] #Control reads miscalled as a modified site

        ntc,ptc,nfc,pfc=self.calculate_calls(true_negative_calls,true_positive_calls,false_negative_calls,false_positive_calls,n_true_negative,n_true_positive)
        ntr,ptr,nfr,pfr=self.calculate_rates(true_negative_calls,true_positive_calls,false_negative_calls,false_positive_calls,n_true_negative,n_true_positive)
        
        if ntc>=true_negative_percentage:
            ntc=true_negative_percentage

        if ptc>=true_positive_percentage:
            ptc=true_positive_percentage
            
        return [ntc,ptc,nfc,pfc,ntr,ptr,nfr,pfr]











