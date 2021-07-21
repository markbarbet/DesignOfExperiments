# -*- coding: utf-8 -*-
"""
Created on Sat May 15 15:19:44 2021

@author: Mark Barbet
"""

import cantera as ct
import MSI
import os
import pandas as pd
import yaml
import re
import numpy as np
import module_0 as m0
import multiprocessing
from shutil import copyfile

def get_matrices_parallel(arg):
        file=arg[0]
        working_dir=arg[1]
        cti_file=arg[2]
        reaction_uncertainty_csv=arg[3]
        #print(arg)
        return run_simulation_parallel([[file,]],working_dir,cti_file,reaction_uncertainty_csv)

def run_simulation_parallel(yaml_list,working_dir,cti_file,reaction_uncertainty_csv):
        #print(yaml_list,working_dir,cti_file,reaction_uncertainty_csv)
        files_to_include = yaml_list
        
        working_directory = working_dir
        
        cti_file = cti_file
        
        reaction_uncertainty_csv = reaction_uncertainty_csv
        
        rate_constant_target_value_data = ''
        
        original_cti=os.path.join(working_directory,cti_file)
        temp_cti=os.path.splitext(yaml_list[0][0])[0]
        temp_cti=cti_file.split('.')[0]+'_'+temp_cti+'_temp.cti'
        temp_cti=os.path.join(working_directory,temp_cti)
        copyfile(original_cti,temp_cti)
        newfile=os.path.split(temp_cti)[-1]
        
        MSI_instance = MSI.optimization.optimization_shell.MSI_optimization(newfile,
                                                                            0.01,
                                                                            1,
                                                                            1,
                                                                            working_directory,
                                                                            files_to_include,
                                                                            reaction_uncertainty_csv,
                                                                            rate_constant_target_value_data)
        MSI_instance.run_simulations()
        
        
        #self.add_yaml_data()
        
        MSI_instance.get_matrices()
        
        S=MSI_instance.S_matrix
        #X=MSI_instance.X
        Z=MSI_instance.z_data_frame
        Y=MSI_instance.Y_data_frame
        #print(Z['value'][630:])
        #print(Y)
        os.remove(temp_cti)
        os.remove(os.path.splitext(temp_cti)[0]+'_updated.cti')
        return {'S':S,'Y':Y,'Z':Z}

class potential_experiments():
    
    
    def __init__(self,input_options:dict={'initialization':None,
                                          'experiment_type':'JSR',
                                          'observables':['H2','O2'],
                                          'temperature_range':[600,1150], 
                                          'pressure_range':[1.0],
                                          'residence_time':[1.2],
                                          'mixture_species':{'H2':[0.001,0.02],'O2':[0.001,0.05]},
                                          'diluent':'Ar',
                                          'constructor_settings':{'intervals':{'species':{'all'},
                                                                               'pressure':[1.0],
                                                                               'temperature':50,
                                                                               'residence_time':0.5},
                                                                  'method':'random_sample',
                                                                  'random_sample_settings':{'n':[10,5]}},
                                          'experiment_yaml_template':'',
                                          'working_dir':'',
                                          'thermal-boundary':'isothermal',
                                          'mechanical-boundary':'constant pressure',
                                          'volume':0.000085,
                                          'yaml_output_name':'DoE_yaml_',
                                          'parallel-computing':True}):
        self.input_options=input_options
        self.constructor_settings=self.input_options['constructor_settings']
        if re.match('[Jj][Ss][Rr]',self.input_options['experiment_type']) or re.match('[Jj]et[- ][Ss]tirred[- ][Rr]eactor',self.input_options['experiment_type']):
            self.exp_type='JSR'
            self.conditions_list=self.JSR_conditions_constructor(self.input_options['constructor_settings'])

            self.yaml_file_list=self.JSR_yaml_constructor(self.conditions_list)
            
            if not self.input_options['parallel-computing']:
                self.matrices=self.get_matrices()
                
            elif self.input_options['parallel-computing']:
                self.cores=input_options['cores']
                args=self.get_args()
                with multiprocessing.Pool(processes=self.cores) as pool:
                    self.matrices=pool.map(get_matrices_parallel,args)
                    
                    
            
    def get_args(self):
        
        args=[]
        for i,file in enumerate(self.yaml_file_list):
            args=args+[[file,
                       self.input_options['initialization'].startup_data['working_dir'],
                       self.input_options['initialization'].MSI_settings['chemical_model'],
                       self.input_options['initialization'].MSI_settings['reaction_uncertainty']]]
        #print(args)
        return args
            
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
    
    def get_matrices(self):
        matrices=[]
        #print(self.yaml_file_list)
        for i,file in enumerate(self.yaml_file_list):
            matrices.append(self.run_simulation([[file,]]))
            
        return(matrices)
        
    
    
    
    
    def run_simulation(self,yaml_list):
        files_to_include = yaml_list
        
        working_directory = self.input_options['initialization'].startup_data['working_dir']
        
        cti_file = self.input_options['initialization'].MSI_settings['chemical_model']
        
        reaction_uncertainty_csv = self.input_options['initialization'].MSI_settings['reaction_uncertainty']
        
        rate_constant_target_value_data = ''
        
        MSI_instance = MSI.optimization.optimization_shell.MSI_optimization(cti_file,
                                                                            0.01,
                                                                            1,
                                                                            1,
                                                                            working_directory,
                                                                            files_to_include,
                                                                            reaction_uncertainty_csv,
                                                                            rate_constant_target_value_data)
        MSI_instance.run_simulations()
        
        
        #self.add_yaml_data()
        
        MSI_instance.get_matrices()
        
        S=MSI_instance.S_matrix
        #X=MSI_instance.X
        Z=MSI_instance.z_data_frame
        Y=MSI_instance.Y_data_frame
        #print(Z['value'][630:])
        #print(Y)
        
        return {'S':S,'Y':Y,'Z':Z}
    
    
        
        
    
    def JSR_yaml_constructor(self,conditions):
        
        template=self.load_to_obj(path=os.path.join(self.input_options['working_dir'],
                                                    self.input_options['experiment_yaml_template']))
        outfilenames=[]
        #print(conditions)
        #print(conditions['index'])
        for i,index in enumerate(conditions['index']):
            
            template['common-properties']['temperature']['value-list']=[float(round(conditions['temperatures'][i],9))]
            template['common-properties']['pressure']['value']=float(round(conditions['pressures'][i],9))
            #print(template['apparatus'])
            template['apparatus']['residence-time']['value']=float(round(conditions['restimes'][i],9))
            template['common-properties']['composition']=[]
            mole_sum=0
            for j,species in enumerate(list(self.input_options['mixture_species'].keys())): 
                template['common-properties']['composition'].append({'species':species,
                                                                     'mole-fraction':float(round(conditions[species][i],9)),
                                                                     'relative-uncertainty':0.05})
                mole_sum=mole_sum+conditions[species][j]
                #Above value of relative uncertainty is not necessary for this procedure, leave unchanged
            diluent_value=1.0-mole_sum
            template['common-properties']['composition'].append({'species':self.input_options['diluent'],
                                                                 'mole-fraction':float(round(diluent_value,9)),
                                                                 'relative-uncertainty':0.05})
            template['common-properties']['assumptions']['thermal-boundary']=self.input_options['thermal-boundary']
            template['common-properties']['assumptions']['mechanical-boundary']=self.input_options['mechanical-boundary']
            template['apparatus']['volume']=float(round(self.input_options['volume'],9))
            template['file-author']={'name':'Design of Experiments Code'}
            template['apparatus']['kind']='JSR'
            template['apparatus']['facility']='Computer Simulated'
            template['apparatus'].pop('institution',None)
            template['common-properties']['assumptions']['equation-of-state']='ideal gas'
            
            outfilenames.append(self.input_options['yaml_output_name']+str(i+1)+'.yaml')
            template['datapoints']={'mole-fraction':[]}
            for j,species in enumerate(self.input_options['observables']):
                template['datapoints']['mole-fraction'].append({'csvfile':os.path.join(self.input_options['working_dir'],
                                                                                       'jsr_simulation_data_'+species+'.csv'),
                                                                'targets':[{'name':species,
                                                                            'species':species,
                                                                            'absolute-uncertainty':5e-5,
                                                                            'relative-uncertainty':0.05}]})
            template['datapoints']['concentration']=[]
            template['datapoints']['concentration'].append({'csvfile':None,
                                                     'targets':[{'name':None,
                                                                 'species':None,
                                                                 'absolute-uncertainty':None,
                                                                 'relative-uncertainty':None}]})
            
            with open(os.path.join(self.input_options['working_dir'],self.input_options['yaml_output_name']+str(i+1)+'.yaml'),'w') as f:
                    yaml.safe_dump(template, f,default_flow_style=False)
        return outfilenames
    def JSR_conditions_constructor(self,settings):
        conditions=pd.DataFrame()
        
        if re.match('[Ii]ntervals',settings['method']):
            if 'all' in list(settings['intervals']['species'].keys()):
                #Assumes that the same interval is applied to all non-diluent species
                fractions={}
                for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                    fractions[species]=np.arange(self.input_options['mixture_species'][species][0],
                                                 self.input_options['mixture_species'][species][1]+settings['intervals']['species']['all'],
                                                 settings['intervals']['species']['all'])
                    
            if len(self.input_options['pressure_range'])==1:
                    
                pressures=self.input_options['pressure_range']
                
            if len(self.input_options['temperature_range'])==2:
                temperatures=np.arange(self.input_options['temperature_range'][0],
                                       self.input_options['temperature_range'][1]+settings['intervals']['temperature'],
                                       settings['intervals']['temperature'])
            elif len(self.input_options['temperature_range'])==1 or len(self.input_options['temperature_range'])>2:
                temperatures=self.input_options['temperature_range']
                
                
                
            if len(self.input_options['residence_time'])==1:
                restimes=self.input_options['residence_time']
                
                
                
        elif re.match('[Rr]andom[- _][Ss]ample',settings['method']):
            
            fractions={}
            
            
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                
                fractions[species]=np.random.uniform(low=self.input_options['mixture_species'][species][0],
                                                     high=self.input_options['mixture_species'][species][1],
                                                     size=settings['random_sample_settings']['n'])
                
            temperatures=np.random.uniform(low=self.input_options['temperature_range'][0],
                                           high=self.input_options['temperature_range'][1],
                                           size=settings['random_sample_settings']['n'])
            
            pressures=np.random.uniform(low=self.input_options['pressure_range'][0],
                                           high=self.input_options['pressure_range'][1],
                                           size=settings['random_sample_settings']['n'])
            
            
            restimes=np.random.uniform(low=self.input_options['residence_time'][0],
                                       high=self.input_options['residence_time'][1],
                                       size=settings['random_sample_settings']['n'])
            
            conditions['temperatures']=temperatures
            conditions['pressures']=pressures
            conditions['restimes']=restimes
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                conditions[species]=fractions[species]
            conditions['index']=np.arange(1,settings['random_sample_settings']['n']+1)
            #print(conditions)
            conditions.to_csv(os.path.join(os.getcwd(),'test_conditions.csv'),index=False)
        return conditions
                
            
            
                
                
                
            
                    
                    
                    
                    
                
                