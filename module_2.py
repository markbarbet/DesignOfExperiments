# -*- coding: utf-8 -*-
"""
Created on Mon May 24 15:08:21 2021

@author: Mark Barbet
"""


import os
import numpy as np
import re
import cantera as ct
import copy
import pandas as pd

class ranking():
    
    def __init__(self,module0=None,module1=None,settings={'batch-size':10,
                                                          'total_new_exp':10,
                                                          'output-csv':'output.csv'}):


        self.module0=module0
        self.module1=module1
        self.settings=settings
        self.excluded_yamls=[]
        total_iters=int(np.ceil(self.settings['total_new_exp']/self.settings['batch-size']))
        #print(total_iters)
        #total_exps=0
        final_exp_dataframe=pd.DataFrame(columns=['experiment','ratio'])
        #current_yamls=copy.deepcopy(self.module1.yaml_file_list)
        self.updated_S=copy.deepcopy(self.module0.S_original)
        for i in np.arange(total_iters):
            #print(i,"wtf")
            self.ranking=self.get_rankings(self.excluded_yamls)
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
        final_exp_dataframe.to_csv(os.path.join(self.module0.startup_data['working_dir'],
        #                                        self.settings['output-csv']),index=False)    
            #self.excluded_yamls=list(final_exp_dataframe['experiment'])
            #self.updated_S=self.get_updated_S(self.updated_S,self.excluded_yamls)
            #current_yamls=[]
            #for j,item in enumerate(self.module1.yaml_file_list):
                #if item not in self.excluded_yamls:
                    #current_yamls=current_yamls+[item]
        #final_exp_dataframe.to_csv(os.path.join(self.module0.startup_data['working_dir'],
                                                #self.settings['output-csv']),index=False)
            
        
        
    def get_updated_S(self,S,yamls):
        
        
        S_countV=0
        S_countH=0
        new_Z=copy.deepcopy(self.module0.initial_optimization.z_data_frame)
        new_Y=copy.deepcopy(self.module0.initial_optimization.Y_data_frame)
        for i,file in enumerate(yamls):
            parametersZ=self.get_Z(os.path.join(self.module1.input_options['working_dir'],file),self.module1.yaml_file_list.index(file))
            new_Z=self.construct_Z_new(parametersZ,new_Z)
            parametersY=self.get_Y(os.path.join(self.module1.input_options['working_dir'],file),self.module1.yaml_file_list.index(file))
            new_Y=self.construct_Y_new(parametersY,new_Y)
            self.experiment_length=self.get_exp_length(os.path.join(self.module1.input_options['working_dir'],file))
            X_to_add=self.get_X_names(parametersZ)
            S1_new,S2_new,S3_new,S4_new,S5_new=self.get_new_S_chunks(self.num_rxns,X_to_add,parametersY,parametersZ,
                                           self.module0.S_original,
                                           self.module0.initial_optimization.Y_data_frame,
                                           self.module0.initial_optimization.z_data_frame,
                                           self.module1.matrices[self.module1.yaml_file_list.index(file)]['S'])
            new_X_list=list(self.module0.initial_optimization.X_data_frame['value'])+X_to_add
            S_proposed=self.build_S(S1_new,S2_new,S3_new,S4_new,S5_new,S,countV=S_countV,countH=S_countH)
            S_countH=S_countH+ self.get_prior_phys_param_len(parametersY)  
            S_countV=self.experiment_length+S_countV+S_countH
            
            
        
    def get_S_current_columns(self,msi_instance=None):
        colnames=list(msi_instance.X_data_frame['value'])
        rownames=list(msi_instance.Y_data_frame['value'])
        return (colnames,rownames)
    
    
    
    def get_Z(self,file,indexer):
        data=self.module1.load_to_obj(file)
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
        Z_to_add=self.module1.matrices[indexer]['Z'][~self.module1.matrices[indexer]['Z']['value'].isin(exclude_list)].copy()
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
        data=self.module1.load_to_obj(file)
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
        Y_to_add=self.module1.matrices[indexer]['Y'][~self.module1.matrices[indexer]['Y']['value'].isin(exclude_list)].copy()
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
    
    def get_S_block(self,col_min,col_max,row_min,row_max, S):
        #print(np.shape(S))
        #return S[row_min:row_max,col_min:col_max]
        return S[col_min:col_max,row_min:row_max]
    
    def get_new_S_chunks(self,rxn_count,X_add,Y_add,add_Z,S_current,Y,Z,S_new):
        #print(np.shape(S_current),np.shape(S_new))
        
        
        S1_new=self.get_S_block(0,self.experiment_length,0,3*self.num_rxns,S_new)
        S2_new=self.get_S_block(0, self.experiment_length, 3*self.num_rxns, 3*self.num_rxns+len(X_add)+1, S_new)
        #print(np.shape(S1_new))
        #print(np.shape(S2_new))
        S3_new=self.get_S_block(self.experiment_length,self.experiment_length+3*self.num_rxns,
                                3*self.num_rxns,3*self.num_rxns+len(X_add)+1,S_new)
        #print(np.shape(S3_new))
        S4_new=self.get_S_block(self.experiment_length+3*self.num_rxns,
                                self.experiment_length+3*self.num_rxns+len(X_add)+1,
                                0,3*self.num_rxns,S_new)
        #print(np.shape(S4_new))
        S5_new=self.get_S_block(self.experiment_length+3*self.num_rxns,
                                self.experiment_length+3*self.num_rxns+len(X_add)+1,
                                3*self.num_rxns,3*self.num_rxns+len(X_add)+1,S_new)
        #print(np.shape(S5_new))
        #print(S5_new)
        return (S1_new,S2_new,S3_new,S4_new,S5_new)
    
    def get_exp_length(self,file):
        data=self.module1.load_to_obj(file)
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
        
        for i,item in enumerate(list(priorY['value'])):
            if 'Ea_' in item:
                final_Ea=i
        if self.module0.MSI_settings['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             self.module0.MSI_settings['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        elif self.module0.MSI_settings['rate_constant_targets'] == '':
            target_length=0
        return len(list(priorY['value'])[final_Ea+1:])-target_length
        
        
    def construct_zeroes(self,S3,S1,S4,countV=0,countH=0):
        prior_exp_len=self.get_prior_exp_len(self.module0.initial_optimization.Y_data_frame)
        Z1=np.zeros((prior_exp_len,np.shape(S3)[1]))
        prior_phys_params=self.get_prior_phys_param_len(self.module0.initial_optimization.Y_data_frame)
        #Need to get length of targets for Z2
        if self.module0.MSI_settings['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             self.module0.MSI_settings['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        elif self.module0.MSI_settings['rate_constant_targets'] == '':
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
        # print(np.shape(block1))
        # print(np.shape(block2))
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
                                            covarience):

            
            reaction_list_from_mechanism = gas.reaction_equations()
            sigma_list_for_target_ks = [[] for reaction in range(len(unique_reactions))]
            shape = k_target_value_S_matrix.shape
            #print(unique_reactions)
            for row in range(shape[0]):
                #print(row)
                SC = np.dot(k_target_value_S_matrix[row,:],covarience)
                sigma_k = np.dot(SC,np.transpose(k_target_value_S_matrix[row,:]))
                sigma_k = np.sqrt(sigma_k)
                #print(row)
                #print(k_target_values_parsed_csv['Reaction'][row])
                indx = reaction_list_from_mechanism.index(k_target_values_parsed_csv['Reaction'][row])
                sigma_list_for_target_ks[unique_reactions.index(indx)].append(sigma_k)
                
            return sigma_list_for_target_ks
        
    def get_k_block(self,S):
        if self.module0.MSI_settings['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             self.module0.MSI_settings['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        return S[len(S[:,0])-target_length:,:]
    
    def get_k_block_proposed(self,S,S_old):
        if self.module0.MSI_settings['rate_constant_targets'] != '':
            targets=pd.read_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             self.module0.MSI_settings['rate_constant_targets']))
            target_length=len(targets['Reaction'])
            
        return S[len(S_old[:,0])-target_length:len(S_old[:,0]),:]
    
    def get_unique_elements(self,l,gas):
        
        uniques=[]
        for i,item in enumerate(l):
            index=list(gas.reaction_equations()).index(item)
            if index not in uniques:
                uniques=uniques+[index]
        return uniques
    
    def get_rankings(self,excluded_yamls):
        ranking_list=[]
        self.rownames_nominal,self.colnames_nominal=self.get_S_current_columns(msi_instance=self.module0.initial_optimization)
        self.mech=os.path.join(self.module1.input_options['working_dir'],self.module0.MSI_settings['chemical_model'])
        gas=ct.Solution(self.mech)
        self.num_rxns=len(gas.reactions())
        if re.match('[Rr]ate[-_ ][Cc]onstant',self.module0.startup_data['quantity_of_interest']):
            k_block_og=self.get_k_block(self.module0.S_original)
            targets=pd.read_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             self.module0.MSI_settings['rate_constant_targets']))
            sigma_list_og=self.calculate_sigmas_for_rate_constants(k_block_og,
                                                                    targets,
                                                                    self.get_unique_elements(list(targets['Reaction']),gas),
                                                                    gas,
                                                                    self.module0.C_original)
            target_reaction=self.module0.startup_data['target_reaction']['equation']
            target_index=list(gas.reaction_equations()).index(target_reaction)
            sigma_list_index=self.get_unique_elements(list(targets['Reaction']),gas).index(target_index)
            #print(sigma_list[sigma_list_index][0])
            original_posterior=sigma_list_og[sigma_list_index][0]
        
        for i,file in enumerate(self.module1.yaml_file_list):
            if file not in excluded_yamls:
                data=self.module1.load_to_obj(os.path.join(self.module1.input_options['working_dir'],file))
                #parametersX,parametersY,parametersZ=self.get_experiment_columns(os.path.join(self.module1.input_options['working_dir'],file),i)
                parametersZ=self.get_Z(os.path.join(self.module1.input_options['working_dir'],file),i)
                parametersY=self.get_Y(os.path.join(self.module1.input_options['working_dir'],file),i)
                #print(parametersZ)
                self.experiment_length=self.get_exp_length(os.path.join(self.module1.input_options['working_dir'],file))
                #print(self.experiment_length)
                X_to_add=self.get_X_names(parametersZ)
                
                new_Z=self.construct_Z_new(parametersZ,self.module0.initial_optimization.z_data_frame)
                
                new_Y=self.construct_Y_new(parametersY,self.module0.initial_optimization.Y_data_frame)
                
                S1_new,S2_new,S3_new,S4_new,S5_new=self.get_new_S_chunks(self.num_rxns,X_to_add,parametersY,parametersZ,
                                           self.module0.S_original,
                                           self.module0.initial_optimization.Y_data_frame,
                                           self.module0.initial_optimization.z_data_frame,
                                           self.module1.matrices[i]['S'])
                #print(self.module0.initial_optimization.Y_data_frame['value'][635:])
                S_proposed=self.build_S(S1_new,S2_new,S3_new,S4_new,S5_new,self.module0.S_original)
                new_X_list=list(self.module0.initial_optimization.X_data_frame['value'])+X_to_add
                #print(list(self.module0.initial_optimization.z_data_frame['value']))
                #print(np.shape(S_proposed))
                #print(self.module0.initial_optimization.z_data_frame)
                #print(len(new_Z))
                s=self.get_normalized_S(S_proposed, new_Z, new_Y)
                c=self.get_covariance(s)
                #print(np.shape(c))
                if re.match('[Rr]ate[-_ ][Cc]onstant',self.module0.startup_data['quantity_of_interest']):
                    k_block=self.get_k_block_proposed(S_proposed,self.module0.S_original)
                    targets=pd.read_csv(os.path.join(self.module0.startup_data['working_dir'],
                                                 self.module0.MSI_settings['rate_constant_targets']))
                    sigma_list=self.calculate_sigmas_for_rate_constants(k_block,
                                                                        targets,
                                                                        self.get_unique_elements(list(targets['Reaction']),gas),
                                                                        gas,
                                                                        c)
                    target_reaction=self.module0.startup_data['target_reaction']['equation']
                    target_index=list(gas.reaction_equations()).index(target_reaction)
                    sigma_list_index=self.get_unique_elements(list(targets['Reaction']),gas).index(target_index)
                    #print(sigma_list[sigma_list_index][0])
                    ranking_list=ranking_list+[sigma_list[sigma_list_index][0]/original_posterior]
                #elif re.match('[Ii]gnition[_ -][Dd]elay',self.module0.startup_data['quantity_of_interest']):
                
        #print(ranking_list)
        output_ranking=pd.DataFrame(columns=['experiment','ratio'])
        output_ranking['experiment']=self.module1.yaml_file_list
        output_ranking['ratio']=ranking_list
        output_ranking.sort_values(by='ratio',ascending=True,inplace=True)
        #output_ranking.to_csv(os.path.join(self.module0.startup_data['working_dir'],
                                             #'output_rankings.csv'),index=False)
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