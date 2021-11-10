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
import traceback

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
        try:
            MSI_instance.run_simulations()
            exclude=''
            MSI_instance.get_matrices()
        
            S=MSI_instance.S_matrix
            #X=MSI_instance.X
            Z=MSI_instance.z_data_frame
            Y=MSI_instance.Y_data_frame
            #print(Z['value'][630:])
            #print(Y)
            os.remove(temp_cti)
            os.remove(os.path.splitext(temp_cti)[0]+'_updated.cti')
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            print('Simulation failed to converge: removing this experiment from YAMLs to try.')
            exclude=files_to_include[0][0]
            S=None
            Y=None
            Z=None
            os.remove(temp_cti)
            #os.remove(os.path.splitext(temp_cti)[0]+'_updated.cti')
        #self.add_yaml_data()
        newfile=None
        temp_cti=None
        original_cti=None
        working_directory=None
        files_to_include=None
        cti_file=None
        yaml_list=None
        working_dir=None
        reaction_uncertainty_csv=None
        MSI_instance=None
        return {'S':S,'Y':Y,'Z':Z,'excluded_yaml':exclude}

class potential_experiments():
    
    
    def __init__(self,input_options:dict={'initialization':None,
                                          'experiment_type':'JSR',
                                          'observables':['H2','O2'],
                                          'observables_abs_uncertainties':[0.00005,0.00005],
                                          'observables_rel_uncertainties':[0.05,0.05],
                                          'temperature_range':[600,1150], 
                                          'temperature-uncertainty': 0.01,
                                          'pressure_range':[1.0],
                                          'pressure-uncertainty': 0.01,
                                          'residence_time':[1.2],
                                          'restime-uncertainty': 0.02,
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
                div_args=list(self.divide_list(args,10000))
                self.matrices=[]
                for count, item in enumerate(div_args):
                    with multiprocessing.Pool(processes=self.cores) as pool:
                        temp_mat=pool.map(get_matrices_parallel,item)
                    self.matrices=self.matrices+temp_mat
                #with multiprocessing.Pool(processes=self.cores) as pool:
                #    self.matrices=pool.map(get_matrices_parallel,args,chunksize=10)

                included_files=[]
                included_matrices=[]
                for i,mat in enumerate(self.matrices):
                    if mat['excluded_yaml']=='':
                        included_files.append(self.yaml_file_list[i])
                        included_matrices.append(mat)
                    
                self.yaml_file_list=included_files
                self.matrices=included_matrices
                    #print(self.yaml_file_list)
                    
                    
    def divide_list(self,arglist,n):
        for i in range(0,len(arglist),n):
            yield arglist[i:i+n]

    
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
            template['common-properties']['temperature']['relative-uncertainty']=float(self.input_options['temperature-uncertainty'])
            template['common-properties']['pressure']['value']=float(round(conditions['pressures'][i],9))
            template['common-properties']['pressure']['relative-uncertainty']=float(self.input_options['pressure-uncertainty'])
            #print(template['apparatus'])
            template['apparatus']['residence-time']['value']=float(round(conditions['restimes'][i],9))
            template['apparatus']['residence-time']['relative-uncertainty']=float(self.input_options['restime-uncertainty'])
            template['common-properties']['composition']=[]
            mole_sum=0
            for j,species in enumerate(list(self.input_options['mixture_species'].keys())): 
                template['common-properties']['composition'].append({'species':species,
                                                                     'mole-fraction':float(round(conditions[species][i],9)),
                                                                     'relative-uncertainty':0.05})
                mole_sum=mole_sum+conditions[species][i]
                #Above value of relative uncertainty is not necessary for this procedure, leave unchanged
            diluent_value=1.0-mole_sum
            #print(diluent_value)
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
                                                                                       self.input_options['yaml_output_name']+str(i+1)+'_'+species+'.csv'),
                                                                'targets':[{'name':species,
                                                                            'species':species,
                                                                            'absolute-uncertainty':self.input_options['observables_abs_uncertainties'][j],
                                                                            'relative-uncertainty':self.input_options['observables_rel_uncertainties'][j]}]})
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

        elif re.match('[Ll]og[-_ ][Hh]alton[-_ ][Ss]ampling',settings['method']):
            dimension=3+len(list(self.input_options['mixture_species'].keys()))
            data_gen=self.halton(dim=dimension, nbpts=settings['random_sample_settings']['n'])
            #Normalizing across temperature range
            data_gen[:,0]=data_gen[:,0]*(self.input_options['temperature_range'][1]-self.input_options['temperature_range'][0])+self.input_options['temperature_range'][0]

            #Normalizing across pressure range
            data_gen[:,1]=data_gen[:,1]*(self.input_options['pressure_range'][1]-self.input_options['pressure_range'][0])+self.input_options['pressure_range'][0]

            #Normalize across residence time range
            data_gen[:,2]=data_gen[:,2]*(self.input_options['residence_time'][1]-self.input_options['residence_time'][0])+self.input_options['residence_time'][0]

            #Normalize across species ranges
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                low_limit=np.log10(self.input_options['mixture_species'][species][0])
                if low_limit==0.0:
                    low_limit=-6.0
                high_limit=np.log10(self.input_options['mixture_species'][species][1])
                
                data_gen[:,i+3]=(data_gen[:,i+3]*(high_limit-low_limit)+low_limit)
                data_gen[:,i+3]=np.power(10,data_gen[:,i+3])

            conditions['temperatures']=data_gen[:,0]
            conditions['pressures']=data_gen[:,1]
            conditions['restimes']=data_gen[:,2]
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                conditions[species]=data_gen[:,i+3]
            conditions['index']=np.arange(1,settings['random_sample_settings']['n']+1)
            conditions.to_csv(os.path.join(os.getcwd(),'test_conditions.csv'),index=False)

        elif re.match('[Ll]og[-_ ][Ss]obol[-_ ][Ss]ampling',settings['method']):
            import sobol_seq
            dimension=3+len(list(self.input_options['mixture_species'].keys()))
            data_gen=sobol_seq.i4_sobol_generate(dimension,settings['random_sample_settings']['n'])
            #Normalizing across temperature range
            data_gen[:,0]=data_gen[:,0]*(self.input_options['temperature_range'][1]-self.input_options['temperature_range'][0])+self.input_options['temperature_range'][0]

            #Normalizing across pressure range
            data_gen[:,1]=data_gen[:,1]*(self.input_options['pressure_range'][1]-self.input_options['pressure_range'][0])+self.input_options['pressure_range'][0]

            #Normalize across residence time range
            data_gen[:,2]=data_gen[:,2]*(self.input_options['residence_time'][1]-self.input_options['residence_time'][0])+self.input_options['residence_time'][0]

            #Normalize across species ranges
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                low_limit=np.log10(self.input_options['mixture_species'][species][0])
                high_limit=np.log10(self.input_options['mixture_species'][species][1])
                if low_limit==0.0:
                    low_limit=-6.0
                data_gen[:,i+3]=(data_gen[:,i+3]*(high_limit-low_limit)+low_limit)
                data_gen[:,i+3]=np.power(10,data_gen[:,i+3])

            conditions['temperatures']=data_gen[:,0]
            conditions['pressures']=data_gen[:,1]
            conditions['restimes']=data_gen[:,2]
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                conditions[species]=data_gen[:,i+3]
            conditions['index']=np.arange(1,settings['random_sample_settings']['n']+1)
            conditions.to_csv(os.path.join(os.getcwd(),'test_conditions.csv'),index=False)


        elif re.match('[Hh]alton[-_ ][Ss]ampling',settings['method']):
            #first need to get dimension count
            dimension=3+len(list(self.input_options['mixture_species'].keys()))
            data_gen=self.halton(dim=dimension, nbpts=settings['random_sample_settings']['n'])

            #Now need to project the data along the species, temperature, pressure, residence time ranges
            
            #Normalizing across temperature range
            data_gen[:,0]=data_gen[:,0]*(self.input_options['temperature_range'][1]-self.input_options['temperature_range'][0])+self.input_options['temperature_range'][0]

            #Normalizing across pressure range
            data_gen[:,1]=data_gen[:,1]*(self.input_options['pressure_range'][1]-self.input_options['pressure_range'][0])+self.input_options['pressure_range'][0]

            #Normalize across residence time range
            data_gen[:,2]=data_gen[:,2]*(self.input_options['residence_time'][1]-self.input_options['residence_time'][0])+self.input_options['residence_time'][0]

            #Normalize across species ranges
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                data_gen[:,i+3]=data_gen[:,i+3]*(self.input_options['mixture_species'][species][1]-self.input_options['mixture_species'][species][0])+self.input_options['mixture_species'][species][0]

            conditions['temperatures']=data_gen[:,0]
            conditions['pressures']=data_gen[:,1]
            conditions['restimes']=data_gen[:,2]
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                conditions[species]=data_gen[:,i+3]
            conditions['index']=np.arange(1,settings['random_sample_settings']['n']+1)
            conditions.to_csv(os.path.join(os.getcwd(),'test_conditions.csv'),index=False)

        elif re.match('[Ss]obol[-_ ][Ss]ampling',settings['method']):
            import sobol_seq
            dimension=3+len(list(self.input_options['mixture_species'].keys()))
            data_gen=sobol_seq.i4_sobol_generate(dimension,settings['random_sample_settings']['n'])
            #Now need to project the data along the species, temperature, pressure, residence time ranges
            
            #Normalizing across temperature range
            data_gen[:,0]=data_gen[:,0]*(self.input_options['temperature_range'][1]-self.input_options['temperature_range'][0])+self.input_options['temperature_range'][0]

            #Normalizing across pressure range
            data_gen[:,1]=data_gen[:,1]*(self.input_options['pressure_range'][1]-self.input_options['pressure_range'][0])+self.input_options['pressure_range'][0]

            #Normalize across residence time range
            data_gen[:,2]=data_gen[:,2]*(self.input_options['residence_time'][1]-self.input_options['residence_time'][0])+self.input_options['residence_time'][0]

            #Normalize across species ranges
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                data_gen[:,i+3]=data_gen[:,i+3]*(self.input_options['mixture_species'][species][1]-self.input_options['mixture_species'][species][0])+self.input_options['mixture_species'][species][0]

            conditions['temperatures']=data_gen[:,0]
            conditions['pressures']=data_gen[:,1]
            conditions['restimes']=data_gen[:,2]
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                conditions[species]=data_gen[:,i+3]
            conditions['index']=np.arange(1,settings['random_sample_settings']['n']+1)
            conditions.to_csv(os.path.join(os.getcwd(),'test_conditions.csv'),index=False)

        return conditions

    def halton(self, dim: int, nbpts: int):
        #print('Hi Mark the TA, can you explain to me how the Second Law of Thermodynamics can be applied to my everyday life?')
        import math
        h = np.full(nbpts * dim, np.nan)
        p = np.full(nbpts, np.nan)
        P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
        lognbpts = math.log(nbpts + 1)
        for i in range(dim):
            b = P[i]
            n = int(math.ceil(lognbpts / math.log(b)))
            for t in range(n):
                p[t] = pow(b, -(t + 1))

            for j in range(nbpts):
                d = j + 1
                sum_ = math.fmod(d, b) * p[0]
                for t in range(1, n):
                    d = math.floor(d / b)
                    sum_ += math.fmod(d, b) * p[t]

                h[j*dim + i] = sum_
        return h.reshape(nbpts, dim)

    
                
            
            
                
                
                
            
                    
                    
                    
                    
                
                