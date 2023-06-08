# -*- coding: utf-8 -*-
"""
Created on Mon May 24 15:08:21 2021

@author: Mark Barbet
"""


import os
#from typing import final
import numpy as np
import re
import cantera as ct
import copy
import pandas as pd
#import progressbar, time, sys
import time, sys
import enlighten
import yaml
import doe_object as dobj

class ranking():
    
    def __init__(self,doe_obj:dobj.doe_object):
        self.inputs=doe_obj.input_options
    
    # def __init__(self,module0=None,module1=None,settings={'batch-size':10,
    #                                                       'total_new_exp':10,
    #                                                       'output-csv':'output.csv'}):
        
        self.subloopBool=True
        #self.module0=module0
        #self.module1=module1
        self.doe_obj=doe_obj
        self.settings={}
        self.settings['batch-size']=doe_obj.input_options['batch-size']
        self.settings['total_new_exp']=doe_obj.input_options['total_new_exp']
        self.settings['output-csv']=doe_obj.input_options['output-csv']
        self.excluded_yamls=[]
        total_iters=int(np.ceil(self.settings['total_new_exp']/self.settings['batch-size']))
        #print(total_iters)
        #total_exps=0
        final_exp_dataframe=pd.DataFrame(columns=['experiment','ratio'])
        #current_yamls=copy.deepcopy(self.module1.yaml_file_list)
        S_countH=0
        S_countV=0
        new_Z=copy.deepcopy(self.doe_obj.Z_original)
        new_Y=copy.deepcopy(self.doe_obj.Y_original)
        new_X_list=list(self.doe_obj.X_original['value'])
        self.updated_S=copy.deepcopy(self.doe_obj.S_original)
        
        #self.printProgressBar(0, total_iters, prefix = 'Finding Best Experiments:', suffix = 'Complete', length = 70)
        #self.down()
        #total = progressbar.ProgressBar(maxval=total_iters)
        #total.start()
        self.manager = enlighten.get_manager()
        self.mainloop=self.manager.counter(total=total_iters,desc='Overall Progress',unit='batches',color='green')
        
        for i in np.arange(total_iters):
            #print(i,"wtf")
            if i==0:
                self.ranking=self.get_rankings(self.excluded_yamls,self.updated_S,i,countH=S_countH, countV=S_countV)
            elif i>0:
                self.ranking=self.get_rankings(self.excluded_yamls,self.updated_S,i,
                                               countH=S_countH, countV=S_countV,
                                               Z_prev=new_Z,Y_prev=new_Y,X_prev=new_X_list)
            if i+1<total_iters or np.remainder(self.settings['total_new_exp'],self.settings['batch-size'])==0:
                #print("Inside first statement")
                final_exp_dataframe=pd.concat([final_exp_dataframe,self.ranking.head(self.settings['batch-size'])],sort=False)
            elif np.remainder(self.settings['total_new_exp'],self.settings['batch-size'])==0 and i+1==total_iters:
                final_exp_dataframe=pd.concat([final_exp_dataframe,
                           self.ranking.head(np.remainder(self.settings['total_new_exp'],
                                                          self.settings['batch-size']))],
                          sort=False)
                #print("Inside second statement")
        #print(final_exp_dataframe)
        #final_exp_dataframe.to_csv(os.path.join(self.module0.startup_data['working_dir'],
        #                                        self.settings['output-csv']),index=False)  

            self.excluded_yamls=list(final_exp_dataframe['experiment'])
            
            #print('We are in iteration: '+str(i) )
            self.updated_S,new_Y,new_Z,new_X_list,S_countH,S_countV=self.get_updated_S(self.updated_S,self.excluded_yamls,new_Z,new_Y,new_X_list, countH=S_countH,countV=S_countV)
            #print('shape')
            #print(np.shape(self.updated_S),np.shape(new_Z))
            current_yamls=[]
            for j,item in enumerate(self.doe_obj.experiment_yamls):
                if item not in self.excluded_yamls:
                    current_yamls=current_yamls+[item]
            
            #total.update(i)
            self.mainloop.update()
            if i==total_iters-1:
                np.savetxt(os.path.join(self.inputs['working_dir'],'final_S.csv'),self.updated_S,delimiter=",")
                new_Z.to_csv(os.path.join(self.inputs['working_dir'],'final_Z.csv'),index=False)
                new_Y.to_csv(os.path.join(self.inputs['working_dir'],'final_Y.csv'),index=False)
            #self.printProgressBar(i + 1, total_iters, prefix = 'Finding Best Experiments:', suffix = 'Complete', length = 70)
        #print(np.shape(self.updated_S),np.shape(new_X_list))
        #total.finish()
        self.subloop.close()
        self.manager.stop()
        matrices=pd.DataFrame(data=self.updated_S,columns=new_X_list)
        matrices['rows']=new_Z['value']
        #np.savetxt(os.path.join(self.module0.startup_data['working_dir'],
        #                                       'S_updated.csv'),self.updated_S)
        new_Z.to_csv(os.path.join(self.inputs['working_dir'],
                                                'Z.csv'),index=False)
        #cov=pd.DataFrame(data=self.updated_c,columns=new_X_list)
        #cov['rows']=new_X_list
        #cov.to_csv(os.path.join(self.module0.startup_data['working_dir'],
        #                                        'cov-debugging.csv'),index=False)
        matrices.to_csv(os.path.join(self.inputs['working_dir'],
                                                'matrices-debugging.csv'),index=False)

        final_exp_dataframe.to_csv(os.path.join(self.inputs['working_dir'],
                                                self.settings['output-csv']),index=False)
            
    def up(self):
        # My terminal breaks if we don't flush after the escape-code
        sys.stdout.write('\x1b[1A')
        sys.stdout.flush()

    def down(self):
        # I could use '\x1b[1B' here, but newline is faster and easier
        sys.stdout.write('\n')
        sys.stdout.flush()


    # Print iterations progress
    def printProgressBar(self,iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total: 
            print()      
    
    def get_updated_S(self,S,yamls,new_Z,new_Y,new_X_list,countH=0,countV=0):
        
        
        S_countV=countV
        S_countH=countH
        
        S_proposed=copy.deepcopy(S)
        for i,file in enumerate([yamls[-1]]):
            #print('Shit '+str(i))
            parametersZ=self.get_Z(os.path.join(self.inputs['working_dir'],file),self.doe_obj.experiment_yamls.index(file))
            new_Z=self.construct_Z_new(parametersZ,new_Z)
            parametersY=self.get_Y(os.path.join(self.inputs['working_dir'],file),self.doe_obj.experiment_yamls.index(file))
            new_Y=self.construct_Y_new(parametersY,new_Y)
            self.experiment_length=self.get_exp_length(os.path.join(self.inputs['working_dir'],file))
            X_to_add=self.get_X_names(parametersZ)
            #print(len(X_to_add))
            S1_new,S2_new,S3_new,S4_new,S5_new=self.get_new_S_chunks(self.num_rxns,X_to_add,parametersY,parametersZ,
                                           self.doe_obj.S_original,
                                           self.doe_obj.Y_original,
                                           self.doe_obj.Z_original,
                                           self.doe_obj.experiment_matrices[self.doe_obj.experiment_yamls.index(file)]['S'])
            new_X_list=new_X_list+X_to_add
            #print(i,new_Y)
            #print('Count H: '+str(S_countH)+', Count V: '+str(S_countV)+', Exp Length: '+str(self.experiment_length)) 
            S_proposed=self.build_S(S1_new,S2_new,S3_new,S4_new,S5_new,S_proposed,countV=S_countV,countH=S_countH)
            #S_countH=S_countH+ self.get_prior_phys_param_len(parametersY) 
            #print(np.shape(S_proposed),np.shape(new_Z))
            S_countH=S_countH+len(X_to_add)
            #print(X_to_add)
            #if i>0:
            S_countV=self.experiment_length+S_countV+len(X_to_add)
            #elif i==0:
            #    S_countV=S_countV+S_countH
            #print('New Count H: '+str(S_countH)+', New Count V: '+str(S_countV)+', Exp Length: '+str(self.experiment_length))
        return (S_proposed,new_Y,new_Z,new_X_list,S_countH,S_countV)
            
        
    def get_S_current_columns(self,doe_obj=None):
        colnames=list(doe_obj.X_original['value'])
        rownames=list(doe_obj.Y_original['value'])
        return (colnames,rownames)
    
    
    
    def get_Z(self,file,indexer):
        data=self.load_to_obj(file)
        previous_exp_index=int(self.rownames_nominal[-1].split('_')[-1])
        current_exp_index=previous_exp_index+1
        exclude_list=[]
        A=[]
        N=[]
        Ea=[]
        for i in np.arange(self.num_rxns):
            A=A+['A_'+str(i)]
        for i in np.arange(self.num_rxns): 
            N=N+['n_'+str(i)]
        for i in np.arange(self.num_rxns):
            Ea=Ea+['Ea_'+str(i)]
        exclude_list=A+N+Ea
        #print(exclude_list)
        Z_to_add=self.doe_obj.experiment_matrices[indexer]['Z'][~self.doe_obj.experiment_matrices[indexer]['Z']['value'].isin(exclude_list)].copy()
        new_Z_names=[]
        for j,item in enumerate(list(Z_to_add['value'])):
            temp=item.split('_')
            if 'experiment' in temp[-1]:
                temp2=re.split('(\d+)',temp[-1])
                #print(temp2)
                temp2[-2]=str(current_exp_index)
                temp[-1]=''.join(temp2)
                
            else:
                temp[-1]=str(current_exp_index)
            new_Z_names=new_Z_names+['_'.join(temp)]
        Z_to_add['value']=new_Z_names
        #Y_to_add=self.module1.matrices[indexer]['Z'][~self.module1.matrices[indexer]['Z']['value'].isin(exclude_list)]
        
        
        return Z_to_add
    
    def get_Y(self,file,indexer):
        data=self.load_to_obj(file)
        previous_exp_index=int(self.rownames_nominal[-1].split('_')[-1])
        current_exp_index=previous_exp_index+1
        exclude_list=[]
        A=[]
        N=[]
        Ea=[]
        for i in np.arange(self.num_rxns):
            A=A+['A_'+str(i)]
        for i in np.arange(self.num_rxns): 
            N=N+['n_'+str(i)]
        for i in np.arange(self.num_rxns):
            Ea=Ea+['Ea_'+str(i)]
        exclude_list=A+N+Ea
        #print(exclude_list)
        Y_to_add=self.doe_obj.experiment_matrices[indexer]['Y'][~self.doe_obj.experiment_matrices[indexer]['Y']['value'].isin(exclude_list)].copy()
        new_Y_names=[]
        for j,item in enumerate(list(Y_to_add['value'])):
            temp=item.split('_')
            if 'experiment' in temp[-1]:
                temp2=re.split('(\d+)',temp[-1])
                #print(temp2)
                temp2[-2]=str(current_exp_index)
                temp[-1]=''.join(temp2)
                
            else:
                temp[-1]=str(current_exp_index)
            new_Y_names=new_Y_names+['_'.join(temp)]
        Y_to_add['value']=new_Y_names
        #Y_to_add=self.module1.matrices[indexer]['Z'][~self.module1.matrices[indexer]['Z']['value'].isin(exclude_list)]
        
        
        return Y_to_add
    
    def get_X_names(self, Z):
        for i, item in enumerate(list(Z['value'])):
            if 'T_experiment' in item:
                start_exp=i
                break
            
        X_values=list(Z['value'])[start_exp:]
        return X_values
    
    def construct_Z_new(self,add_Z,Z):
        new_Z=pd.concat([Z,add_Z],ignore_index=True)
        return new_Z
    
    def construct_Y_new(self,add_Y,Y):
        new_Y=pd.concat([Y,add_Y],ignore_index=True)
        return new_Y
    
    # def get_S_block(self,col_min,col_max,row_min,row_max, S):
    #     #print(np.shape(S))
    #     #return S[row_min:row_max,col_min:col_max]
    #     return S[col_min:col_max,row_min:row_max]
    def get_S_block(self,col_min,col_max,row_min,row_max, S,flag=False,X_add=None):
        #print(np.shape(S))
        #return S[row_min:row_max,col_min:col_max]
        y_len=3*self.num_rxns+len(self.inputs['observables'])+X_add
        if flag:
            #y_len=len(self.doe_obj.Y_original)
            #print("Fuckshit is : "+str(y_len))
            temp_S=np.zeros((3*self.num_rxns+len(self.inputs['observables'])+X_add,3*self.num_rxns+X_add))
            indices=np.arange(y_len-len(self.inputs['observables']))
            temp_S[indices+len(self.inputs['observables']),indices]=1.0
            temp_S[0:len(self.inputs['observables']),0:np.shape(S)[1]]=S
            #print(temp_S)
            S=temp_S
            #print(col_min,col_max,row_min,row_max)
            #print(np.shape(temp_S))
        #print(np.shape(S))
        return S[col_min:col_max,row_min:row_max]

    # def get_new_S_chunks(self,rxn_count,X_add,Y_add,add_Z,S_current,Y,Z,S_new):
    #     #print(np.shape(S_current),np.shape(S_new))
        
        
    #     S1_new=self.get_S_block(0,self.experiment_length,0,3*self.num_rxns,S_new)
    #     S2_new=self.get_S_block(0, self.experiment_length, 3*self.num_rxns, 3*self.num_rxns+len(X_add)+1, S_new)
    #     #print(np.shape(S1_new))
    #     #print(np.shape(S2_new))
    #     S3_new=self.get_S_block(self.experiment_length,self.experiment_length+3*self.num_rxns,
    #                             3*self.num_rxns,3*self.num_rxns+len(X_add)+1,S_new)
    #     #print(np.shape(S3_new))
    #     S4_new=self.get_S_block(self.experiment_length+3*self.num_rxns,
    #                             self.experiment_length+3*self.num_rxns+len(X_add)+1,
    #                             0,3*self.num_rxns,S_new)
    #     #print(np.shape(S4_new))
    #     S5_new=self.get_S_block(self.experiment_length+3*self.num_rxns,
    #                             self.experiment_length+3*self.num_rxns+len(X_add)+1,
    #                             3*self.num_rxns,3*self.num_rxns+len(X_add)+1,S_new)
    #     #print(np.shape(S5_new))
    #     #print(S5_new)
    #     return (S1_new,S2_new,S3_new,S4_new,S5_new)
    def get_new_S_chunks(self,rxn_count,X_add,Y_add,add_Z,S_current,Y,Z,S_new):
        #print(np.shape(S_current),np.shape(S_new))
        #time.sleep(5)
        #print(np.shape(S_new))
        with open(os.path.join(os.getcwd(),'log.txt'),'a') as f:
            f.write(str(np.shape(S_new)))
            f.write('\n')
        #print(S_new)
        #print("Length is: "+str(self.experiment_length))
        S1_new=self.get_S_block(0,self.experiment_length,0,3*self.num_rxns,S_new,X_add=len(X_add))
        #print('EXP Length',self.experiment_length)
        S2_new=self.get_S_block(0, self.experiment_length, 3*self.num_rxns, 3*self.num_rxns+len(X_add)+1, S_new,X_add=len(X_add))
        #print(np.shape(S1_new))
        #print(np.shape(S2_new))
        #print(self.experiment_length,self.experiment_length+3*self.num_rxns,
        #                        3*self.num_rxns,3*self.num_rxns+len(X_add)+1)
        S3_new=self.get_S_block(self.experiment_length,self.experiment_length+3*self.num_rxns,
                                3*self.num_rxns,3*self.num_rxns+len(X_add)+1,S_new,flag=True,X_add=len(X_add))
        #print(np.shape(S3_new))
        S4_new=self.get_S_block(self.experiment_length+3*self.num_rxns,
                                self.experiment_length+3*self.num_rxns+len(X_add)+1,
                                0,3*self.num_rxns,S_new,flag=True,X_add=len(X_add))
        #print(np.shape(S4_new))
        S5_new=self.get_S_block(self.experiment_length+3*self.num_rxns,
                                self.experiment_length+3*self.num_rxns+len(X_add)+1,
                                3*self.num_rxns,3*self.num_rxns+len(X_add)+1,S_new,flag=True,X_add=len(X_add))
        #print(np.shape(S5_new))
        #print(S5_new)
        return (S1_new,S2_new,S3_new,S4_new,S5_new)
    
    def get_exp_length(self,file):
        data=self.load_to_obj(file)
        mol_length=0
        conc_length=0
        for i,item in enumerate(data['datapoints']['mole-fraction']):
            if item['csvfile'] != None:
                mol_length=mol_length+1
        for i,item in enumerate(data['datapoints']['concentration']):
            if item['csvfile'] != None:
                conc_length=conc_length+1
                
        return mol_length+conc_length
    
    def get_prior_exp_len(self,priorY):
        #print(len(priorY))
        for i,item in enumerate(list(priorY['value'])):
            if 'A_' in item:
                end_exp_index=i
                break
        return len(list(priorY['value'])[0:end_exp_index])
        
    def get_prior_phys_param_len(self,priorY):
        #print(priorY)
        for i,item in enumerate(list(priorY['value'])):
            if 'Ea_' in item:
                final_Ea=i
        if self.inputs['MSI_settings']['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.inputs['working_dir'],
                                             self.inputs['MSI_settings']['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        elif self.inputs['MSI_settings']['rate_constant_targets'] == '':
            target_length=0
        return len(list(priorY['value'])[final_Ea+1:])-target_length
        
        
    def construct_zeroes(self,S3,S1,S4,countV=0,countH=0):
        prior_exp_len=self.get_prior_exp_len(self.doe_obj.Y_original)
        Z1=np.zeros((prior_exp_len,np.shape(S3)[1]))
        prior_phys_params=self.get_prior_phys_param_len(self.doe_obj.Y_original)
        #Need to get length of targets for Z2
        if self.inputs['MSI_settings']['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.inputs['working_dir'],
                                             self.inputs['MSI_settings']['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        elif self.inputs['MSI_settings']['rate_constant_targets'] == '':
            target_length=0
            
        #print(prior_phys_params,target_length)
        Z2=np.zeros((prior_phys_params+target_length+countV,np.shape(S3)[1]))
        Z3=np.zeros((np.shape(S1)[0],prior_phys_params+countH))
        Z4=np.zeros((np.shape(S4)[0],prior_phys_params+countH))
        
        return (Z1,Z2,Z3,Z4)
        
        
        
    def build_S(self,S1,S2,S3,S4,S5,S_old,countV=0,countH=0):
        
        Z1,Z2,Z3,Z4=self.construct_zeroes(S3,S1,S4,countV=countV,countH=countH)
        #print(np.shape(Z1))
        #print(np.shape(S3))
        #print(np.shape(Z2))
        block1=np.block([[Z1],[S3],[Z2]])
        block2=np.block([[S1,Z3,S2],[S4,Z4,S5]])
        #print(np.shape(Z2),np.shape(Z3),np.shape(Z4))
        
        #print('S',np.shape(S_old))
        #print('Z1',np.shape(Z1))
        # print('S3',np.shape(S3))
        # print('Z2',np.shape(Z2))
        # print('S1',np.shape(S1))
        # print('Z3',np.shape(Z3))
        # print('S2',np.shape(S2))
        # print('S4',np.shape(S4))
        # print('Z4',np.shape(Z4))
        # print('S5',np.shape(S5))
        #print(countV)
        #print(np.shape(block1))
        #print(np.shape(block2),'s')
        #print(np.shape(S_old))
        S=np.block([[S_old,block1],[block2]])
        return S
        
    def get_covariance(self,s):
        c = np.dot(np.transpose(s),s)
        c = np.linalg.inv(c)
        
    
        return (c)
    def return_posteriors(self,c,Z):
        covariance_posterior_df = pd.DataFrame(c)
        covariance_posterior_df.columns = Z
        covariance_posterior_df.reindex(labels = Z)
        posterior_diag = np.diag(c)
        posterior_sigmas = np.sqrt(posterior_diag)
        posterior_sigmas_df = pd.DataFrame({'parameter': Z,'value': posterior_sigmas.reshape((posterior_sigmas.shape[0],))})
        posterior_diag_df =  pd.DataFrame({'parameter': Z,'value': posterior_diag.reshape((posterior_diag.shape[0],))})
        sorted_posterior_diag  = posterior_diag_df.sort_values(by=['value'])
        
    def get_normalized_S(self,S,Z,Y):
        
        s=S/np.array(Z['Uncertainty'])[:,None]
        return s
    
    def calculate_sigmas_for_rate_constants(self,k_target_value_S_matrix,
                                            k_target_values_parsed_csv,
                                            unique_reactions,
                                            gas,
                                            covariance):

            
            reaction_list_from_mechanism = gas.reaction_equations()
            sigma_list_for_target_ks = [[] for reaction in range(len(unique_reactions))]
            shape = k_target_value_S_matrix.shape
            #print(unique_reactions)
            for row in range(shape[0]):
                #print(k_target_value_S_matrix,'blahasfjhbjhb')
                SC = np.dot(k_target_value_S_matrix[row,:],covariance)
                sigma_k = np.dot(SC,np.transpose(k_target_value_S_matrix[row,:]))
                sigma_k = np.sqrt(sigma_k)
                #print(row)
                #print(k_target_values_parsed_csv['Reaction'][row])
                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
                
            return sigma_list_for_target_ks
        
    def get_k_block(self,S):
        if self.inputs['MSI_settings']['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.inputs['working_dir'],
                                             self.inputs['MSI_settings']['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        return S[len(S[:,0])-target_length:,:]
    
    def get_k_block_proposed(self,S,S_old):
        if self.inputs['MSI_settings']['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.inputs['working_dir'],
                                             self.inputs['MSI_settings']['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        return S[len(S_old[:,0])-target_length:len(S_old[:,0]),:]
    
    def get_unique_elements(self,l,gas):
        
        uniques=[]
        for i,item in enumerate(l):
            index=list(gas.reaction_equations()).index(item)
            if index not in uniques:
                uniques=uniques+[index]
        return uniques




    def calculate_sigma(self,S_row,C):
        SC = np.dot(S_row,C)
        sigma = np.dot(SC,np.transpose(S_row))
        sigma = np.sqrt(sigma)
        return sigma


    def get_ignition_block(self,S):
        exp_len=self.get_prior_exp_len(self.doe_obj.Y_original)

        ig_block=S[exp_len-1,:]

        return ig_block


    def load_to_obj(self, path:str = ''):
        """
        Takes in a file path for a yaml file and returns a dictionary of 
        simulation information.

        Parameters
        ----------
        path : str, optional
            The path to where the yaml file is stored. The default is ''.

        Returns
        -------
        config : dictionary
            An unorganized dictionary that contains information reguarding
            the experiment the yaml input file was written for.

        """
        with open(path) as f:
            config = yaml.load(f,Loader=yaml.FullLoader)
        return config


    def get_yaml_conditions(self,fname):

        template=self.load_to_obj(fname)
        outputs={}
        outputs['temperature']=float(template['common-properties']['temperature']['value-list'][0])
        outputs['pressure']=float(template['common-properties']['pressure']['value'])
        outputs['residence-time']=float(template['apparatus']['residence-time']['value'])
        for i,spec in enumerate(template['common-properties']['composition']):
            outputs[spec['species']]=spec['mole-fraction']
        
        return outputs
    
    
    def load_to_obj(self, path:str = ''):
        """
        Takes in a file path for a yaml file and returns a dictionary of 
        simulation information.

        Parameters
        ----------
        path : str, optional
            The path to where the yaml file is stored. The default is ''.

        Returns
        -------
        config : dictionary
            An unorganized dictionary that contains information reguarding
            the experiment the yaml input file was written for.

        """
        with open(path) as f:
            config = yaml.load(f,Loader=yaml.FullLoader)
        return config
    
    
    def get_rankings(self,excluded_yamls,S,iteration,countH=0,countV=0,Z_prev=None,Y_prev=None,X_prev=None):
        #self.up()
        #sub = progressbar.ProgressBar(maxval=len(self.module1.yaml_file_list))
        #sub.start()
        #if self.subloopBool:
        self.subloop=self.manager.counter(total=len(self.doe_obj.experiment_matrices),desc='Possible experiments',unit='Experiments',color='red',leave=True)
        #self.subloopBool=False
        ranking_list=[]
        self.rownames_nominal,self.colnames_nominal=self.get_S_current_columns(doe_obj=self.doe_obj)
        self.mech=os.path.join(self.inputs['working_dir'],self.inputs['MSI_settings']['chemical_model'])
        gas=ct.Solution(self.mech)
        self.num_rxns=len(gas.reactions())
        if re.match('[Rr]ate[-_ ][Cc]onstant',self.inputs['quantity_of_interest']):
            k_block_og=self.get_k_block(self.doe_obj.S_original)
            targets=pd.read_csv(os.path.join(self.inputs['working_dir'],
                                             self.inputs['MSI_settings']['rate_constant_targets']))
            sigma_list_og=self.calculate_sigmas_for_rate_constants(k_block_og,
                                                                    targets,
                                                                    self.get_unique_elements(list(targets['Reaction']),gas),
                                                                    gas,
                                                                    self.doe_obj.covar_original)
            target_reaction=self.inputs['target_reaction']['equation']
            target_index=list(gas.reaction_equations()).index(target_reaction)
            sigma_list_index=self.get_unique_elements(list(targets['Reaction']),gas).index(target_index)
            #print(sigma_list[sigma_list_index][0])
            original_posterior=sigma_list_og[sigma_list_index][0]

        elif re.match('[Ii]gnition[-_ ][Dd]elay',self.inputs['quantity_of_interest']):
            '''This block collects the original uncertainty in the ignition delay quantity of interest. Gets the final experiment line
               and estimates the uncertainty in the ignition delay.'''

            ig_block_og=self.get_ignition_block(self.doe_obj.S_original)
            print(ig_block_og)
            #print(self.module0.)
            sigma_ig=self.calculate_sigma(ig_block_og,self.doe_obj.covar_original)
            original_posterior=sigma_ig

        final_yamls=[]
        for i,file in enumerate(self.doe_obj.experiment_yamls):
            
            if file not in excluded_yamls:
                final_yamls.append(file)
                print('File is:'+str(file))
                data=self.load_to_obj(os.path.join(self.inputs['working_dir'],file))
                #parametersX,parametersY,parametersZ=self.get_experiment_columns(os.path.join(self.module1.input_options['working_dir'],file),i)
                parametersZ=self.get_Z(os.path.join(self.inputs['working_dir'],file),i)
                parametersY=self.get_Y(os.path.join(self.inputs['working_dir'],file),i)
                if iteration==0:
                    
                    #print(parametersZ)
                    self.experiment_length=self.get_exp_length(os.path.join(self.inputs['working_dir'],file))
                    #print(self.experiment_length)
                    X_to_add=self.get_X_names(parametersZ)
                    
                    new_Z=self.construct_Z_new(parametersZ,self.doe_obj.Z_original)
                    
                    new_Y=self.construct_Y_new(parametersY,self.doe_obj.Y_original)
                    
                elif iteration>0:
                    new_Z=self.construct_Z_new(parametersZ,Z_prev)
                    X_to_add=self.get_X_names(parametersZ)
                    new_Y=self.construct_Y_new(parametersY,Y_prev)
                    self.experiment_length=self.get_exp_length(os.path.join(self.inputs['working_dir'],file))
                S1_new,S2_new,S3_new,S4_new,S5_new=self.get_new_S_chunks(self.num_rxns,X_to_add,parametersY,parametersZ,
                                           self.doe_obj.S_original,
                                           self.doe_obj.Y_original,
                                           self.doe_obj.Z_original,
                                           self.doe_obj.experiment_matrices[i]['S'])
                #print(self.module0.initial_optimization.Y_data_frame['value'][635:])
                #print('CountH: '+str(countH)+', countV: '+str(countV))
                S_proposed=self.build_S(S1_new,S2_new,S3_new,S4_new,S5_new,S,countH=countH,countV=countV)

                if iteration==0:
                    new_X_list=list(self.doe_obj.X_original['value'])+X_to_add
                elif iteration>0:
                    new_X_list=X_prev+X_to_add
                #print(list(self.module0.initial_optimization.z_data_frame['value']))
                #print(np.shape(S_proposed))
                #print(self.module0.initial_optimization.z_data_frame)
                #print(len(new_Z),np.shape(S_proposed))
                
                s=self.get_normalized_S(S_proposed, new_Z, new_Y)
                c=self.get_covariance(s)
                #print(c)
                self.updated_c=c
                #print(np.shape(c))
                if re.match('[Rr]ate[-_ ][Cc]onstant',self.inputs['quantity_of_interest']):
                    k_block=self.get_k_block_proposed(S_proposed,self.doe_obj.S_original)
                    targets=pd.read_csv(os.path.join(self.inputs['working_dir'],
                                                 self.inputs['MSI_settings']['rate_constant_targets']))
                    sigma_list=self.calculate_sigmas_for_rate_constants(k_block,
                                                                        targets,
                                                                        self.get_unique_elements(list(targets['Reaction']),gas),
                                                                        gas,
                                                                        c)
                    #print(sigma_list)
                    target_reaction=self.inputs['target_reaction']['equation']
                    target_index=list(gas.reaction_equations()).index(target_reaction)
                    sigma_list_index=self.get_unique_elements(list(targets['Reaction']),gas).index(target_index)
                    #print(sigma_list[sigma_list_index][0],'text')
                    ranking_list=ranking_list+[sigma_list[sigma_list_index][0]/original_posterior]
                    #print(sigma_list,'poo')
                elif re.match('[Ii]gnition[-_ ][Dd]elay',self.inputs['quantity_of_interest']):
                    print('Entered correct statement')
                    ig_block=self.get_ignition_block(S_proposed)
                    if i==len(self.doe_obj.experiment_yamls)-1:
                        print(ig_block)
                    sigma=self.calculate_sigma(ig_block,c)
                    ranking_list=ranking_list+[sigma/original_posterior]
                #elif re.match('[Ii]gnition[_ -][Dd]elay',self.module0.startup_data['quantity_of_interest']):
            #if np.remainder(i,10)==0:
                #sub.update(i)
            self.subloop.update()
        #sub.finish()
        #print(ranking_list)
        output_ranking=pd.DataFrame(columns=['experiment','ratio'])
        output_ranking['experiment']=final_yamls
        temps=[]
        pres=[]
        restime=[]
        conds=self.get_yaml_conditions(os.path.join(self.inputs['working_dir'],final_yamls[0]))
        species_index=copy.deepcopy(list(conds.keys()))
        #print(list(conds.keys()).remove('temperature'))
        species_index.remove('temperature')
        species_index.remove('pressure')
        species_index.remove('residence-time')
        specs={}
        for i,spec in enumerate(species_index):
            specs[spec]=[]
        for i,fname in enumerate(final_yamls):
            conds=self.get_yaml_conditions(os.path.join(self.inputs['working_dir'],fname))
            temps.append(conds['temperature'])
            pres.append(conds['pressure'])
            restime.append(conds['residence-time'])
            for k,spec in enumerate(species_index):
                specs[spec].append(conds[spec])
        #print(ranking_list)
        output_ranking['ratio']=ranking_list
        output_ranking['temperature']=temps
        output_ranking['pressure']=pres
        output_ranking['residence-time']=restime
        for i,spec in enumerate(species_index):
            output_ranking[spec]=specs[spec]
        output_ranking.sort_values(by='ratio',ascending=True,inplace=True)

        #output_ranking.to_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             #'output_rankings.csv'),index=False)
        #print(output_ranking)
        #print(output_ranking)
        return output_ranking
            #posteriors=self.return_posteriors(c,new_Z)
            #posteriors.to_csv('test.csv')
            #print(c)
            #print(new_S[687:])
            #print(X_to_add)
            #print(parametersY)
            #print('x',parametersX)
            #print('y',parametersY)
            #print('z',parametersZ)