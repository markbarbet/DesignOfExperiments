# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:46:37 2021

@author: Mark Barbet
"""


import module_0
import module_1
import module_2
import sys
import os
import fire
import yaml
import re


class design_of_experiments():
    
    
    '''
    Design-of-experiments container class.  Contains all the attributes of an experimental design procedure
    and determines which sub-functions to run based on information provided in an input file.
    
    
    '''
    
    
    def __init__(self,input_options):
        '''
        Initializes the 'multiscale_informatics' class, and calls 'self.get_parameters'
        to seperate and sort parameters from the parsed input file.

        Parameters
        ----------
        input_options : List
            List of strings of the parsed lines of the MSI input file.

        Returns
        -------
        None.

        '''
        self.input_options=input_options
        self.get_parameters()
        
    def get_prior_experiments(self):
        self.module_0_input['experiments']=[]
        if 'prior_experiments' not in list(self.input_options.keys()):
            self.module_0_input['experiments']=[]
        elif self.input_options['prior_experiments']==None:
            self.module_0_input['experiments']=[]
        elif self.input_options['prior_experiments']['list']==None:
            self.module_0_input['experiments']=[]
        elif 'prior_experiments' in list(self.input_options.keys()) or self.input_options['prior_experiments']['list']!=[]:
            if self.input_options['prior_experiments']['list'][0]['file1']==None:
                self.module_0_input['experiments']=[]
            else:
                for i,dic in enumerate(self.input_options['prior_experiments']['list']):
                    self.module_0_input['experiments'].append([dic['file1']])
                    if dic['file2']!=None:
                        self.module_0_input['experiments'][-1].append(dic['file2'])
            
            
    def get_QoI(self):
        
        self.module_0_input['quantity_of_interest']=self.input_options['quantity_of_interest']
        
        
    def get_yaml_template(self):
        self.module_0_input['yaml_template']=[self.input_options['experiment_yaml_template']]
        
    def get_working_dir(self):
        self.working_dir=self.input_options['working_directory']
        self.module_0_input['working_dir']=self.working_dir
        
    def get_target_reaction(self):
        self.module_0_input['target_reaction']={}
        self.module_0_input['target_reaction']['equation']=self.input_options['target_reaction']['equation']
        self.module_0_input['target_reaction']['pressure']=self.input_options['target_reaction']['pressure']['value']*101325
        tempdic={}
        for i,el in enumerate(self.input_options['target_reaction']['mixture']):
            tempdic[el['species']]=el['fraction']
        if None in list(tempdic.keys()):
            tempdic={}    
        self.module_0_input['target_reaction']['mixture']=tempdic
        
        self.module_0_input['target_reaction']['temperatures']=self.input_options['target_reaction']['temperature']['values']
    
    def get_MSI_settings(self):
        self.MSI_settings={}
        self.MSI_settings['chemical_model']=self.input_options['MSI_settings']['chemical_model']
        self.MSI_settings['reaction_uncertainty']=self.input_options['MSI_settings']['reaction_uncertainty']
        if self.input_options['MSI_settings']['rate_constant_targets']==None:
            self.MSI_settings['rate_constant_targets']=''
        else:
            self.MSI_settings['rate_constant_targets']=self.input_options['MSI_settings']['rate_constant_targets']
     
    def get_experiment_attributes(self):
        self.module_1_input['experiment_type']=self.input_options['new_experiments']['experiment_type']
        self.module_1_input['diluent']=self.input_options['new_experiments']['diluent']
        self.module_1_input['mechanical-boundary']=self.input_options['new_experiments']['mechanical-boundary']
        self.module_1_input['thermal-boundary']=self.input_options['new_experiments']['thermal-boundary']
        self.module_1_input['volume']=self.input_options['new_experiments']['volume']
        self.module_1_input['yaml_output_name']=self.input_options['new_experiments']['yaml_output_name']
        #self.module_1_input['yaml_template']=self.module_0_input['yaml_template']
        self.module_1_input['working_dir']=self.module_0_input['working_dir']
        self.module_1_input['observables']=self.input_options['new_experiments']['observables']
        self.module_1_input['observables_abs_uncertainties']=self.input_options['new_experiments']['observables_abs_uncertainties']
        self.module_1_input['observables_rel_uncertainties']=self.input_options['new_experiments']['observables_rel_uncertainties']
        for i,item in enumerate(self.module_1_input['observables']):
            if not item:
                self.module_1_input['observables'][i]='NO'
        
    def conditions_constructor(self):
        self.module_1_input['constructor_settings']={'intervals':{'species':{'all'},
                                                                               'pressure':[1.0],
                                                                               'temperature':50,
                                                                               'residence_time':0.5},
                                                                  'method':'random_sample',
                                                                  'random_sample_settings':{'n':10}}
        if re.match('[Rr]andom[- _][Ss]ample',self.input_options['new_experiments']['constructor_settings']['method']):
            self.module_1_input['constructor_settings']['method']='random_sample'
            
            self.module_1_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.module_1_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.module_1_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.module_1_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
            self.module_1_input['constructor_settings']['random_sample_settings']['n']=self.input_options['new_experiments']['constructor_settings']['random_sample_size']
            self.module_1_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            
    def get_experiment_yaml_template(self):
        self.module_1_input['experiment_yaml_template']=self.input_options['new_experiments']['yaml_template']
    
    def get_mod2_settings(self):
        self.mod2_settings={}
        self.mod2_settings['batch-size']=self.input_options['new_experiments']['batch-size']
        self.mod2_settings['total_new_exp']=self.input_options['new_experiments']['total_new_experiments']
        self.mod2_settings['output-csv']=self.input_options['new_experiments']['output-csv']
    
    
    def get_parallel_setting(self):
        self.module_1_input['parallel-computing']=self.input_options['new_experiments']['parallel-computation']
        #print(type(self.module_1_input['parallel-computing']))
        self.module_1_input['cores']=int(self.input_options['new_experiments']['cores'])
    
    def get_parameters(self):
        self.module_0_input={}
        
        self.get_prior_experiments()
        self.get_QoI()
        self.get_yaml_template()
        
        self.get_working_dir()
        self.get_target_reaction()
        self.get_MSI_settings()
        
        self.module_1_input={}
        self.module_1_input['initialization']=None
        self.get_experiment_attributes()
        self.get_experiment_yaml_template()
        self.conditions_constructor()
        self.get_parallel_setting()
        self.get_mod2_settings()
        
        
    def run_DoE(self):
        #print(self.module_0_input)
        self.mod0=module_0.DoE(startup_data=self.module_0_input,MSI_settings=self.MSI_settings)
        self.module_1_input['initialization']=self.mod0
        self.mod1=module_1.potential_experiments(input_options=self.module_1_input)
        self.mod2=module_2.ranking(module0=self.mod0,module1=self.mod1,settings=self.mod2_settings)
        
def parser(input_file):
    
    with open(input_file) as f:    
            config = yaml.load(f,Loader=yaml.FullLoader)
            for i,dic in enumerate(config['new_experiments']['mixture_species']):
                if not dic['species']:
                    config['new_experiments']['mixture_species'][i]['species']='NO'
                    
            
    return config

def main(input_file=''):
    '''
    A design-of-experiments code for use primarily on the Burke lab jet-stirred-reactor.  
    In principle, could be extended for use with other systems.

    Parameters
    ----------
    input_file : String
        The default is ''.  Enter the path to the directory where the design of experiments input file is located.

    Returns
    -------
    simulation : TYPE
        DESCRIPTION.

    '''
    if input_file=='':
        print('Please run program with defined input file using --input_file=FILEPATH')

    elif input_file !='':
        
        input_options=parser(input_file)
        simulation=design_of_experiments(input_options)
        simulation.run_DoE()
        #simulation.write_convergence()
        
        return simulation

if __name__ == '__main__':
    a=fire.Fire(main)

