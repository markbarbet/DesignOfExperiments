# -*- coding: utf-8 -*-
"""
Created on Fri May  7 12:17:15 2021

@author: Mark Barbet
"""

import os
import cantera as ct
import MSI
import pandas as pd
import re
import numpy as np
import yaml
import copy
import doe_object as obj


class DoE():
    
    
    def __init__(self,doe_obj:obj.doe_object):
    
        self.experiments=doe_obj.input_options['experiments']
        self.quantity_of_interest=doe_obj.input_options['quantity_of_interest']
        self.yaml_template=doe_obj.input_options['yaml_template']
        self.MSI_settings=doe_obj.input_options['MSI_settings']
        self.startup_data=doe_obj.input_options
    
    
    # def __init__(self,startup_data:dict={'experiments':[],
    #                                         'quantity_of_interest':'',
    #                                         'qoi_exp':False,
    #                                         'yaml_template':'',
    #                                         'working_dir':'',
    #                                         'target_reaction':{'temperatures':[1000],
    #                                                         'pressure':ct.one_atm,
    #                                                         'equation':'',
    #                                                         'mixture':{}}},
    #                 MSI_settings={'chemical_model':'','reaction_uncertainty':'',
    #                             'rate_constant_targets':''}):
    #     self.experiments=startup_data['experiments']
    #     self.quantity_of_interest=startup_data['quantity_of_interest']
    #     self.yaml_template=startup_data['yaml_template']
    #     self.MSI_settings=MSI_settings
    #     self.startup_data=startup_data
        
        if self.experiments:
            if  re.match('[rR]ate[- ][Cc]onstant',self.quantity_of_interest):
                self.exp_data_rate_constant()
            elif re.match('[Ii]gnition[-_ ][Dd]elay[-_ ][Tt]ime',self.quantity_of_interest):
                self.exp_data_no_rate_constant()
            else: 
                self.exp_data_no_rate_constant()
        elif not self.experiments:
            if  re.match('[rR]ate[- ][Cc]onstant',self.quantity_of_interest):
                self.startup_rate_constant()
            elif re.match('[Ii]gnition[-_ ][Dd]elay[-_ ][Tt]ime',self.quantity_of_interest):
                self.startup_no_rate_constant()
            else: 
                self.startup_no_rate_constant()
                
        self.initial_optimization = self.run_msi(rate_constant_target_value_data=self.MSI_settings['rate_constant_targets'])   
        print(self.initial_optimization.z_data_frame)     
        # wait=input("Press enter to continue")
        self.S_original=self.initial_optimization.S_matrix
        old_X_list=list(self.initial_optimization.X_data_frame['value'])
        Z=copy.deepcopy(self.initial_optimization.z_data_frame)
        matrices=pd.DataFrame(data=self.S_original,columns=old_X_list)
        matrices['rows']=Z['value']
        matrices.to_csv(os.path.join(self.startup_data['working_dir'],
                                                'matrices-original.csv'),index=False)
        self.C_original=self.initial_optimization.covarience        
        covar=pd.DataFrame(data=self.C_original,columns=old_X_list)
        covar.to_csv(os.path.join(self.startup_data['working_dir'],
                                                'covar-original.csv'),index=False)
        Z.to_csv(os.path.join(self.startup_data['working_dir'],
                                                'Z_old.csv'),index=False)                                        
        self.uncertainties=self.calculate_uncertainty(self.S_original,self.C_original)
        doe_obj.set_priors(self.S_original,covar,Z,self.initial_optimization.X_data_frame,self.initial_optimization.Y_data_frame)
        #return doe_obj
    
    def generate_rate_targets(self):
        tempgas=ct.Solution(os.path.join(self.startup_data['working_dir'],self.MSI_settings['chemical_model']))
        reaction_equation=self.startup_data['target_reaction']['equation']
        tempgas.TPX=self.startup_data['target_reaction']['temperatures'][0],self.startup_data['target_reaction']['pressure'],{'Ar':1.0}
        #rate=tempgas.forward_rate_constants[tempgas.reaction_equations().index(reaction_equation)]
        units_reaction_types=['ElementaryReaction',
                                  'PlogReaction',
                                  'ChebyshevReaction']
        index=tempgas.reaction_equations().index(reaction_equation)
        if type(tempgas.reaction(index)).__name__ in units_reaction_types:
                coeff_sum = sum(tempgas.reaction(index).reactants.values())
                    #ask about the mixture composition
                if self.startup_data['target_reaction']['pressure'] == 0:
                    pressure = 1e-9
                else:
                    pressure = self.startup_data['target_reaction']['pressure']
                if self.startup_data['target_reaction']['mixture'] !={}:
                    #print(self.startup_data['target_reaction']['mixture'])
                    tempgas.TPX = self.startup_data['target_reaction']['temperatures'][0],pressure,self.startup_data['target_reaction']['mixture']
    
                else:
                    tempgas.TPX = self.startup_data['target_reaction']['temperatures'][0],pressure,{'Ar':1}
                reaction_number_in_cti = tempgas.reaction_equations().index(reaction_equation)
                k = tempgas.forward_rate_constants[reaction_number_in_cti]
                if coeff_sum==1:
                    k=k
                if coeff_sum==2:
                    k = k*1000
                if coeff_sum==3:
                    k=k*1000000
        elif type(tempgas.reaction(index)).__name__ == 'ThreeBodyReaction' or type(tempgas.reaction(index)).__name__ =='FalloffReaction':
                coeff_sum = sum(tempgas.reaction(index).reactants.values())
                    #ask about the mixture composition
                if self.startup_data['target_reaction']['pressure'] == 0:
                    pressure = 1e-9
                else:
                    pressure = self.startup_data['target_reaction']['pressure']
                if self.startup_data['target_reaction']['mixture'] !={}:
                    tempgas.TPX = self.startup_data['target_reaction']['temperatures'][0],pressure,self.startup_data['target_reaction']['mixture']
    
                else:
                    tempgas.TPX = self.startup_data['target_reaction']['temperatures'][0],pressure,{'Ar':1}
                reaction_number_in_cti = tempgas.reaction_equations().index(reaction_equation)
                k = tempgas.forward_rate_constants[reaction_number_in_cti]
                if coeff_sum==1:
                    k=k
                if coeff_sum==2:
                    k = k*1000
                if coeff_sum==3:
                    k=k*1000000
                    
        if self.MSI_settings['rate_constant_targets']=='':
            targets=pd.DataFrame(columns=['Reaction','temperature','pressure','M','k','ln_unc_k','W'])
            targets['Reaction']=[reaction_equation]
            targets['temperature']=self.startup_data['target_reaction']['temperatures']
            targets['pressure']=[self.startup_data['target_reaction']['pressure']]
            targets['M']=[0]
            targets['k']=[k]
            targets['ln_unc_k']=[10**6]
            targets['W']=[1]
            targets.to_csv(os.path.join(self.startup_data['working_dir'],'rate_target_QoI.csv'),index=False)
            self.MSI_settings['rate_constant_targets']='rate_target_QoI.csv'
            
        elif self.MSI_settings['rate_constant_targets']!= '':
            targets=pd.read_csv(os.path.join(self.startup_data['working_dir'],self.MSI_settings['rate_constant_targets']))
            targets['Reaction']=list(targets['Reaction'])+[reaction_equation]
            targets['temperature']=list(targets['temperature'])+self.startup_data['target_reaction']['temperatures']
            targets['pressure']=list(targets['pressure'])+[self.startup_data['target_reaction']['pressure']]
            targets['M']=list(targets['M'])+[0]
            targets['k']=list(targets['k'])+[k]
            targets['ln_unc_k']=list(targets['ln_unc_k'])+[10**6]
            targets['W']=list(targets['W'])+[1]
            targets.to_csv(os.path.join(self.startup_data['working_dir'],self.MSI_settings['rate_constant_targets'],index=False))
            
            
            
        
    
    def startup_rate_constant(self):
        self.generate_rate_targets()
        
        self.yaml_file_list=self.experiments+[self.yaml_template]
        
    def run_ignition_delay(self,conditions:dict):
        p = MSI.cti_core.cti_processor.Processor(os.path.join(self.startup_data['working_dir'],self.MSI_settings['chemical_model']))
        print(conditions['pressures'])
        print(conditions['temperatures'])
        print(conditions['observables'])
        print(conditions['conditions'])
        print(conditions['target_type'])
        print(conditions['target'])
        sim = MSI.simulations.instruments.ignition_delay.ignition_delay_wrapper(pressures=conditions['pressures'],
                                                                                temperatures=conditions['temperatures'],
                                                                                observables=conditions['observables'], 
                                                                                kineticSens=0, 
                                                                                physicalSens=0, 
                                                                                conditions=conditions['conditions_to_run'], 
                                                                                thermalBoundary=conditions['thermalBoundary'], 
                                                                                mechanicalBoundary=conditions['mechanicalBoundary'],
                                                                                processor = p,
                                                                                finalTime=conditions['finalTime'],
                                                                                target=conditions['target'],
                                                                                target_type=conditions['target_type'],
                                                                                save_physSensHistories=0,
                                                                                save_timeHistories=1,fullParsedYamlFile=conditions)
        solution,trash=sim.run()
        return solution
        
    def write_fake_csv(self,filename,data):
        
        data.to_csv(filename,index=False)
        
        return filename
        
    def startup_no_rate_constant(self):
        yaml_class_inst = MSI.simulations.yaml_parser.Parser()
        yaml_object=yaml_class_inst.load_to_obj(path=os.path.join(self.startup_data['working_dir'],self.yaml_template[0]))
        yaml_dict = yaml_class_inst.parse_ignition_delay_obj(loaded_exp=yaml_object)

        if re.match('[Ii]gnition[- ][Dd]elay', self.quantity_of_interest):
            #if not yaml_dict['csvFiles']:
                #print('Template yaml csv file is empty: Assuming data for quantity of interest must be generated.')
                
                '''Include code here to run the ignition delay at the conditions of the QoI, then write csv files
                   and edit the yaml file to include it'''

                solution=self.run_ignition_delay(yaml_dict)
                solution.rename(columns={'delay':'tau_s'},inplace=True)
                solution=solution[['temperature','tau_s','pressure']]
                outfile=self.write_fake_csv(os.path.join(self.startup_data['working_dir'],'temp_data.csv'),solution)
                yaml_dict['csvFiles'].append(outfile)
                yaml_dict['ignitionDelayCsvFiles'].append(outfile)
                yaml_dict['ignitionDelayRelativeUncertainty']=[10000]
                yaml_dict['ignitionDelayAbsoluteUncertainty']=[1.0]
                with open(os.path.join(self.startup_data['working_dir'],self.yaml_template[0])) as f:
                    yaml_file = yaml.load(f,Loader=yaml.FullLoader)
                
                #print(yaml_file)
                yaml_file['datapoints']['ignition-delay'][0]['csvfile']=outfile
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['name']='tau'
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['observable']='tau'
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['relative-uncertainty']= 10000
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['absolute-uncertainty']=1.0
                with open(os.path.join(self.startup_data['working_dir'],self.yaml_template[0]),'w') as f:
                    yaml.safe_dump(yaml_file, f,default_flow_style=False)
                
                
            #else:


                #pass
            
            
        self.yaml_file_list=self.experiments+[self.yaml_template]
        
            
    def exp_data_rate_constant(self):
        self.generate_rate_targets()
        
        self.yaml_file_list=self.experiments
    def exp_data_no_rate_constant(self):
        '''Simply provide startup_data['experiments']' with the appropriate
           yaml files to run msi'''
        if re.match('[Ii]gnition[- ][Dd]elay', self.quantity_of_interest):
            if not self.startup_data['qoi_exp']:
                '''Code enters this block if there is experimental data but the qoi is not among them'''
                yaml_class_inst = MSI.simulations.yaml_parser.Parser()
                #print(self.yaml_template)
                yaml_object=yaml_class_inst.load_to_obj(path=os.path.join(self.startup_data['working_dir'],self.yaml_template[0]))
                yaml_dict = yaml_class_inst.parse_ignition_delay_obj(loaded_exp=yaml_object)
                solution=self.run_ignition_delay(yaml_dict)
                solution.rename(columns={'delay':'tau_s'},inplace=True)
                solution=solution[['temperature','tau_s','pressure']]
                outfile=self.write_fake_csv(os.path.join(self.startup_data['working_dir'],'temp_data.csv'),solution)
                yaml_dict['csvFiles'].append(outfile)
                yaml_dict['ignitionDelayCsvFiles'].append(outfile)
                yaml_dict['ignitionDelayRelativeUncertainty']=[10000]
                yaml_dict['ignitionDelayAbsoluteUncertainty']=[1.0]
                with open(os.path.join(self.startup_data['working_dir'],self.yaml_template[0])) as f:
                    yaml_file = yaml.load(f,Loader=yaml.FullLoader)
                
                #print(yaml_file)
                yaml_file['datapoints']['ignition-delay'][0]['csvfile']=outfile
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['name']='tau'
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['observable']='tau'
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['relative-uncertainty']= 10000
                yaml_file['datapoints']['ignition-delay'][0]['targets'][0]['absolute-uncertainty']=1.0
                with open(os.path.join(self.startup_data['working_dir'],self.yaml_template[0]),'w') as f:
                    yaml.safe_dump(yaml_file, f,default_flow_style=False)
                self.yaml_file_list=self.experiments+[self.yaml_template]
                #print(self.yaml_file_list)
            elif self.startup_data['qoi_exp']:
                '''Code enters here if the quantity of interest is among the experimental data'''
                yaml_class_inst = MSI.simulations.yaml_parser.Parser()
                #print(self.experiments)
                yaml_object=yaml_class_inst.load_to_obj(path=os.path.join(self.startup_data['working_dir'],self.experiments[0][0]))
                yaml_dict = yaml_class_inst.parse_ignition_delay_obj(loaded_exp=yaml_object)
                solution=self.run_ignition_delay(yaml_dict)
                outfile=self.write_fake_csv(os.path.join(self.startup_data['working_dir'],'temp_data.csv'),solution)
                yaml_dict['csvFiles'].append(outfile)
                yaml_dict['ignitionDelayCsvFiles'].append(outfile)
                yaml_dict['ignitionDelayRelativeUncertainty']=[10000]
                yaml_dict['ignitionDelayAbsoluteUncertainty']=[1.0]
                print(self.yaml_template)
                #with open(self.yaml_template,'w') as f:
                #   yaml.safe_dump(yaml_dict, f,default_flow_style=False)
                self.yaml_file_list=self.experiments

                #self.yaml_file_list=self.experiments
    
        
    def calculate_uncertainty(self,S,C):
        
        uncertainties=np.matmul(np.matmul(S,C),np.transpose(S))
        
        return uncertainties
        
        
    def run_msi(self,rate_constant_target_value_data=''):
        
        files_to_include = self.yaml_file_list
        #print(self.yaml_file_list)
        
        working_directory = self.startup_data['working_dir']
        
        cti_file = self.MSI_settings['chemical_model']
        
        reaction_uncertainty_csv = self.MSI_settings['reaction_uncertainty']
        
        rate_constant_target_value_data = rate_constant_target_value_data
        #print(rate_constant_target_value_data)
        
        print(files_to_include)
        MSI_instance = MSI.optimization.optimization_shell.MSI_optimization(cti_file,
                                                                            0.01,
                                                                            1,
                                                                            1,working_directory,
                                                                            files_to_include,
                                                                            reaction_uncertainty_csv,
                                                                            rate_constant_target_value_data)
        MSI_instance.one_run_optimization()
        #print(MSI_instance.X_data_frame['value'])
        return MSI_instance
    
        
        
        
        
