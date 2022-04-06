# -*- coding: utf-8 -*-
"""
Created on Thu May 20 14:46:37 2021

@author: Mark Barbet
"""


from numpy import exp
from pandas.core.algorithms import rank
import module_0
import module_1
import module_2
import sys
import os
import fire
import yaml
import re
import doe_object as dobj
import importlib as imp



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
        self.parsed_input['experiments']=[]
        if 'prior_experiments' not in list(self.input_options.keys()):
            self.parsed_input['experiments']=[]
        elif self.input_options['prior_experiments']==None:
            self.parsed_input['experiments']=[]
        elif self.input_options['prior_experiments']['list']==None:
            self.parsed_input['experiments']=[]
        elif 'prior_experiments' in list(self.input_options.keys()) or self.input_options['prior_experiments']['list']!=[]:
            if self.input_options['prior_experiments']['list'][0]['file1']==None:
                self.parsed_input['experiments']=[]
            else:
                for i,dic in enumerate(self.input_options['prior_experiments']['list']):
                    self.parsed_input['experiments'].append([dic['file1']])
                    if dic['file2']!=None:
                        self.parsed_input['experiments'][-1].append(dic['file2'])
            
            
    def get_QoI(self):
        
        self.parsed_input['quantity_of_interest']=self.input_options['quantity_of_interest']
        if ',' in self.parsed_input['quantity_of_interest']:
            self.parsed_input['quantity_of_interest']=self.parsed_input['quantity_of_interest'].split(',')
            
        
    def get_yaml_template(self):
        self.parsed_input['yaml_template']=[self.input_options['experiment_yaml_template']]
        
    def get_working_dir(self):
        self.working_dir=self.input_options['working_directory']
        self.parsed_input['working_dir']=self.working_dir
        
    def get_target_reaction(self):
        self.parsed_input['target_reaction']={}
        self.parsed_input['target_reaction']['equation']=self.input_options['target_reaction']['equation']
        self.parsed_input['target_reaction']['pressure']=self.input_options['target_reaction']['pressure']['value']*101325
        tempdic={}
        for i,el in enumerate(self.input_options['target_reaction']['mixture']):
            tempdic[el['species']]=el['fraction']
        if None in list(tempdic.keys()):
            tempdic={}    
        self.parsed_input['target_reaction']['mixture']=tempdic
        
        self.parsed_input['target_reaction']['temperatures']=self.input_options['target_reaction']['temperature']['values']
    
    def get_MSI_settings(self):
        self.parsed_input['MSI_settings']={}
        self.parsed_input['MSI_settings']['msi_iterations']=self.input_options['MSI_settings']['iterations']
        self.parsed_input['MSI_settings']['chemical_model']=self.input_options['MSI_settings']['chemical_model']
        self.parsed_input['MSI_settings']['reaction_uncertainty']=self.input_options['MSI_settings']['reaction_uncertainty']
        if self.input_options['MSI_settings']['rate_constant_targets']==None:
            self.parsed_input['MSI_settings']['rate_constant_targets']=''
        else:
            self.parsed_input['MSI_settings']['rate_constant_targets']=self.input_options['MSI_settings']['rate_constant_targets']
     
    def get_experiment_attributes(self):
        self.parsed_input['experiment_type']=self.input_options['new_experiments']['experiment_type']
        self.parsed_input['diluent']=self.input_options['new_experiments']['diluent']
        self.parsed_input['mechanical-boundary']=self.input_options['new_experiments']['mechanical-boundary']
        self.parsed_input['thermal-boundary']=self.input_options['new_experiments']['thermal-boundary']
        self.parsed_input['volume']=self.input_options['new_experiments']['volume']
        self.parsed_input['yaml_output_name']=self.input_options['new_experiments']['yaml_output_name']
        #self.module_1_input['yaml_template']=self.module_0_input['yaml_template']
        #self.parsed_input['working_dir']=self.parsed_input['working_dir']
        self.parsed_input['observables']=self.input_options['new_experiments']['observables']
        self.parsed_input['observables_abs_uncertainties']=self.input_options['new_experiments']['observables_abs_uncertainties']
        self.parsed_input['observables_rel_uncertainties']=self.input_options['new_experiments']['observables_rel_uncertainties']
        for i,item in enumerate(self.parsed_input['observables']):
            if not item:
                self.parsed_input['observables'][i]='NO'
        
    def conditions_constructor(self):
        self.parsed_input['constructor_settings']={'intervals':{'species':{'all'},
                                                                               'pressure':[1.0],
                                                                               'temperature':50,
                                                                               'residence_time':0.5},
                                                                  'method':'random_sample',
                                                                  'random_sample_settings':{'n':10}}
        if re.match('[Rr]andom[- _][Ss]ample',self.input_options['new_experiments']['constructor_settings']['method']):
            self.parsed_input['constructor_settings']['method']='random_sample'
            
            self.parsed_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.parsed_input['temperature-uncertainty']=self.input_options['new_experiments']['temperature_range']['uncertainty']
            self.parsed_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.parsed_input['pressure-uncertainty']=self.input_options['new_experiments']['pressure_range']['uncertainty']
            self.parsed_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.parsed_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
            self.parsed_input['constructor_settings']['random_sample_settings']['n']=self.input_options['new_experiments']['constructor_settings']['random_sample_size']
            self.parsed_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            self.parsed_input['restime-uncertainty']=self.input_options['new_experiments']['residence_time']['uncertainty']
        
        elif re.match('[Gg]rid[-_ ][Ll]og',self.input_options['new_experiments']['constructor_settings']['method']):
            self.parsed_input['gridpoints']={}
            self.parsed_input['constructor_settings']['method']='grid_log'
            self.parsed_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.parsed_input['gridpoints']['temperature']=self.input_options['new_experiments']['temperature_range']['gridpoints']
            self.parsed_input['temperature-uncertainty']=self.input_options['new_experiments']['temperature_range']['uncertainty']
            self.parsed_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.parsed_input['gridpoints']['pressure']=self.input_options['new_experiments']['pressure_range']['gridpoints']
            self.parsed_input['pressure-uncertainty']=self.input_options['new_experiments']['pressure_range']['uncertainty']
            
            self.parsed_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.parsed_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
                self.parsed_input['gridpoints'][dic['species']]=dic['gridpoints']
                
            
            self.parsed_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            self.parsed_input['gridpoints']['restime']=self.input_options['new_experiments']['residence_time']['gridpoints']
            self.parsed_input['restime-uncertainty']=self.input_options['new_experiments']['residence_time']['uncertainty']
            
            

        elif re.match('[Hh]alton[-_ ][Ss]ampling',self.input_options['new_experiments']['constructor_settings']['method']):
            self.parsed_input['constructor_settings']['method']='halton_sampling'
            
            self.parsed_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.parsed_input['temperature-uncertainty']=self.input_options['new_experiments']['temperature_range']['uncertainty']
            self.parsed_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.parsed_input['pressure-uncertainty']=self.input_options['new_experiments']['pressure_range']['uncertainty']
            self.parsed_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.parsed_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
            self.parsed_input['constructor_settings']['random_sample_settings']['n']=self.input_options['new_experiments']['constructor_settings']['random_sample_size']
            self.parsed_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            self.parsed_input['restime-uncertainty']=self.input_options['new_experiments']['residence_time']['uncertainty']


        elif re.match('[Ll]og[-_ ][Hh]alton[-_ ][Ss]ampling',self.input_options['new_experiments']['constructor_settings']['method']):
            self.parsed_input['constructor_settings']['method']='log_halton_sampling'
            
            self.parsed_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.parsed_input['temperature-uncertainty']=self.input_options['new_experiments']['temperature_range']['uncertainty']
            self.parsed_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.parsed_input['pressure-uncertainty']=self.input_options['new_experiments']['pressure_range']['uncertainty']
            self.parsed_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.parsed_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
            self.parsed_input['constructor_settings']['random_sample_settings']['n']=self.input_options['new_experiments']['constructor_settings']['random_sample_size']
            self.parsed_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            self.parsed_input['restime-uncertainty']=self.input_options['new_experiments']['residence_time']['uncertainty']

        elif re.match('[Ss]obol[-_ ][Ss]ampling',self.input_options['new_experiments']['constructor_settings']['method']):
            self.parsed_input['constructor_settings']['method']='sobol_sampling'
            
            self.parsed_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.parsed_input['temperature-uncertainty']=self.input_options['new_experiments']['temperature_range']['uncertainty']
            self.parsed_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.parsed_input['pressure-uncertainty']=self.input_options['new_experiments']['pressure_range']['uncertainty']
            self.parsed_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.parsed_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
            self.parsed_input['constructor_settings']['random_sample_settings']['n']=self.input_options['new_experiments']['constructor_settings']['random_sample_size']
            self.parsed_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            self.parsed_input['restime-uncertainty']=self.input_options['new_experiments']['residence_time']['uncertainty']

        elif re.match('[Ll]og[_ -][Ss]obol[-_ ][Ss]ampling',self.input_options['new_experiments']['constructor_settings']['method']):
            self.parsed_input['constructor_settings']['method']='log_sobol_sampling'
            
            self.parsed_input['temperature_range']=[self.input_options['new_experiments']['temperature_range']['low'],
                                                      self.input_options['new_experiments']['temperature_range']['high']]
            self.parsed_input['temperature-uncertainty']=self.input_options['new_experiments']['temperature_range']['uncertainty']
            self.parsed_input['pressure_range']=[self.input_options['new_experiments']['pressure_range']['low'],
                                                      self.input_options['new_experiments']['pressure_range']['high']]
            self.parsed_input['pressure-uncertainty']=self.input_options['new_experiments']['pressure_range']['uncertainty']
            self.parsed_input['mixture_species']={}
            for i,dic in enumerate(self.input_options['new_experiments']['mixture_species']):
                self.parsed_input['mixture_species'][dic['species']]=[dic['low'],dic['high']]
            self.parsed_input['constructor_settings']['random_sample_settings']['n']=self.input_options['new_experiments']['constructor_settings']['random_sample_size']
            self.parsed_input['residence_time']=[self.input_options['new_experiments']['residence_time']['low'],
                                                      self.input_options['new_experiments']['residence_time']['high']]
            self.parsed_input['restime-uncertainty']=self.input_options['new_experiments']['residence_time']['uncertainty']

    def get_experiment_yaml_template(self):
        self.parsed_input['experiment_yaml_template']=self.input_options['new_experiments']['yaml_template']
    
    def get_mod2_settings(self):
        #self.mod2_settings={}
        self.parsed_input['ranking-scheme']=self.input_options['new_experiments']['ranking-scheme']
        self.parsed_input['batch-size']=self.input_options['new_experiments']['batch-size']
        self.parsed_input['total_new_exp']=self.input_options['new_experiments']['total_new_experiments']
        self.parsed_input['output-csv']=self.input_options['new_experiments']['output-csv']
        
    
    
    def get_parallel_setting(self):
        self.parsed_input['parallel-computing']=self.input_options['new_experiments']['parallel-computation']
        #print(type(self.module_1_input['parallel-computing']))
        self.parsed_input['cores']=int(self.input_options['new_experiments']['cores'])


    def get_qoi_exp(self):
        if self.parsed_input['yaml_template']==None and self.parsed_input['experiments']==[]:
            print('Error: No experiment template for predicting priors, and no prior experiments included')
        elif self.parsed_input['yaml_template']==None and not self.parsed_input['experiments']==[]:
            #print('Assume first file in the prior experiments contains quantity of interest to be calculated')
            self.parsed_input['qoi_exp']=True
        elif self.parsed_input['yaml_template']!=None:
            self.parsed_input['qoi_exp']=False
            
    def get_module_bools(self):
        self.parsed_input['priors_bool']=self.input_options['new_experiments']['priors']
        self.parsed_input['exp_gen_bool']=self.input_options['new_experiments']['experiment_generation']
        self.parsed_input['ranking_bool']=self.input_options['new_experiments']['ranking']
        
    def get_previous_doe_object(self,obj_path):
        doe_obj=dobj.doe_object.read_object(obj_path)
        return doe_obj
    
    def get_prior_doe_obj_path(self):
        self.parsed_input['prior_doe_obj']=self.input_options['new_experiments']['prior_doe_obj']
        
    def get_output_doe(self):
        self.parsed_input['final_doe_obj']=self.input_options['new_experiments']['final_doe_obj']
        
    
    def get_parameters(self):
        self.parsed_input={}
        
        self.get_prior_experiments()
        self.get_QoI()
        
        self.get_yaml_template()
        self.get_qoi_exp()
        self.get_working_dir()
        self.get_target_reaction()
        self.get_MSI_settings()
        self.get_module_bools()
        self.get_prior_doe_obj_path()
        #self.module_1_input={}
        self.parsed_input['initialization']=None
        self.get_experiment_attributes()
        self.get_experiment_yaml_template()
        self.conditions_constructor()
        self.get_parallel_setting()
        self.get_mod2_settings()
        self.get_output_doe()
        
        
    def initialize_doe_obj(self):
        
        return(dobj.doe_object(self.parsed_input)) 
        
    def run_DoE(self):
        #print(self.module_0_input)
        
        priorbool=self.parsed_input['priors_bool']
        expbool=self.parsed_input['exp_gen_bool']
        rankbool=self.parsed_input['ranking_bool']
        #print(priorbool,expbool,rankbool)
        if priorbool and not expbool and not rankbool:
            doe_obj=self.initialize_doe_obj()
            self.mod0=module_0.DoE(doe_obj)
            self.mod0=None
            
        elif not priorbool and expbool and not rankbool:
            doe_obj=dobj.doe_object.read_object(self.parsed_input['prior_doe_obj'])
            self.mod1=module_1.potential_experiments(doe_obj=doe_obj)
            self.mod1=None
            
        elif not priorbool and not expbool and rankbool:
            doe_obj=dobj.doe_object.read_object(self.parsed_input['prior_doe_obj'])
            imp.import_module(self.parsed_input['ranking-scheme'])
            self.mod2=eval(self.parsed_input['ranking-scheme']).ranking(doe_obj=doe_obj)
            self.mod2=None
            
        elif priorbool and expbool and not rankbool:
            doe_obj=self.initialize_doe_obj()
            self.mod0=module_0.DoE(doe_obj)
            self.mod0=None
            self.mod1=module_1.potential_experiments(doe_obj=doe_obj)
            self.mod1=None
            #print(doe_obj.experiment_matrices)
            
        elif not priorbool and expbool and rankbool:
            doe_obj=dobj.doe_object.read_object(self.parsed_input['prior_doe_obj'])
            self.mod1=module_1.potential_experiments(doe_obj=doe_obj)
            self.mod1=None
            imp.import_module(self.parsed_input['ranking-scheme'])
            self.mod2=eval(self.parsed_input['ranking-scheme']).ranking(doe_obj=doe_obj)
            self.mod2=None
            
        elif priorbool and expbool and rankbool:
            doe_obj=self.initialize_doe_obj()
            self.mod0=module_0.DoE(doe_obj)
            self.mod0=None
            self.mod1=module_1.potential_experiments(doe_obj=doe_obj)
            self.mod1=None
            mod=imp.import_module(self.parsed_input['ranking-scheme'])
            print(dir())
            self.mod2=mod.ranking(doe_obj=doe_obj)
            #self.mod2=eval(self.parsed_input['ranking-scheme']).ranking(doe_obj=doe_obj)
            self.mod2=None
            
        elif not priorbool and not expbool and not rankbool:
            print("Error:  cannot run code without at least one module set to run.")
        
        #print('doe',self.parsed_input['final_doe_obj'])
        if self.parsed_input['final_doe_obj'] != None:
            doe_obj_dir=os.path.join(self.working_dir,self.parsed_input['final_doe_obj'])
            doe_obj.write_object(doe_obj_dir)   
        
        #print(doe_obj.covar_original)
        #self.module_1_input['initialization']=self.mod0
        #self.mod1=module_1.potential_experiments(doe_obj=doe_obj)
        #self.mod2=module_2.ranking(module0=self.mod0,module1=self.mod1,settings=self.mod2_settings)
        
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

