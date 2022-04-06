import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
import shutil
import zipfile



class doe_object():
    
    
    def __init__(self, parsed_inputs):
        
        #Initialization for class.  Takes parsed input file.
        self.input_options=parsed_inputs
        
        
        self.total_new_exps=0
        #self.conditions=None
        self.excluded_files=None
        self.new_uncertainties={}
        self.selected_conditions=None
        self.selected_yamls=None
        self.sorted_uncertainties=None
        self.predictions=None
        self.selected_S=[]
        self.selected_Z=[]
        self.selected_Y=[]
        
        
        #Information derived from prior models
        self.S_original=None
        self.covar_original=None
        self.Z_original=None
        self.X_original=None
        self.Y_original=None
        self.prior_experiments=[]
    
    #def set_selected_matrices_all(self,slist,ylist,zlist):
    #    self.selected_S=slist
    #    self.selected_Y=ylist
    #    self.selected_Z=zlist
        
        
    def set_selected_conditions(self,selected_conds):
        self.selected_conditions=selected_conds
        
    def set_selected_yamls(self,selected_yamls):
        self.selected_yamls=selected_yamls
        
    def set_priors(self,S,C,Z,X,Y):
        self.S_original=S
        self.covar_original=C
        self.Z_original=Z
        self.X_original=X
        self.Y_original=Y   
    
    def set_selected_matrices(self,S,Z,Y):
        self.selected_S.append(S)
        self.selected_Y.append(Y)
        self.selected_Z.append(Z)
    
    def plot_convergence(self,log=True):
        if self.sorted_uncertainties == None:
            print('Cannot plot convergence until algorithm to find new uncertainties is run')
        elif self.sorted_uncertainties != None and log:
            #Placeholder code
            pass
        elif self.sorted_uncertainties != None and not log:
            #Placeholder code
            pass
    
    def set_matrices(self,matrices):
        self.experiment_matrices=matrices
    
    def set_yaml_list(self,yamls):
        self.experiment_yamls=yamls
    
    def write_exps_S_to_file(self,directoryName,yamls):
        S_mats=[]
        for i,filename in enumerate(yamls):
            S_mats.append(self.experiment_matrices[i]['S'])
        np.savez_compressed(os.path.join(directoryName,"S_matrices"),*S_mats)
    
    
    def write_exps_Z_to_file(self, directoryName,yamls):
        tempdata=pd.DataFrame(columns=['value'])
        tempdata['value']=self.experiment_matrices[0]['Z']['value']
        #tempdata.join()
        #print(self.experiment_matrices[0]['Z'])
        #self.experiment_matrices[0]['Z'].rename(columns={'Uncertainty':'1'},inplace=True)
        for i, mat in enumerate(yamls):
            tempdata=tempdata.join(self.experiment_matrices[i]['Z']['Uncertainty'])
            #print(b)
            tempdata.rename(columns={'Uncertainty':str(i+1)},inplace=True)
        tempdata.to_csv(os.path.join(directoryName,'Z_matrices.csv'),index=False)
    
    def write_exps_Y_to_file(self, directoryName,yamls):
        tempdata=pd.DataFrame(columns=['value'])
        tempdata['value']=self.experiment_matrices[0]['Y']['value']
        #tempdata.join()
        print(self.experiment_matrices[0]['Y'])
        #self.experiment_matrices[0]['Z'].rename(columns={'Uncertainty':'1'},inplace=True)
        for i, mat in enumerate(yamls):
            tempdata=tempdata.join(self.experiment_matrices[i]['Y']['ln_difference'])
            tempdata.rename(columns={'ln_difference':str(i+1)},inplace=True)
        tempdata.to_csv(os.path.join(directoryName,'Y_matrices.csv'),index=False) 
        
    def plot_sensitivities(self,reactions=[],n_rxns=0):
        #Code to plot sensitivities.  Placeholder code.
        pass
    
    def write_priors(self,directory):
        np.save(os.path.join(directory,'S_original'),self.S_original)
        np.save(os.path.join(directory,'c_original'),self.covar_original)
        self.Z_original.to_csv(os.path.join(directory,'Z_original.csv'),index=False)
        self.Y_original.to_csv(os.path.join(directory,'Y_original.csv'),index=False)
        self.X_original.to_csv(os.path.join(directory,'X_original.csv'),index=False)
    
    def write_input_options(self,filename,object_out):
        print(object_out)
        with open(filename,'w') as outfile:
            output=yaml.safe_dump(object_out,outfile,default_flow_style=False)
        
    def write_selected_exps(self,directoryName):
        for i,name in enumerate(self.selected_yamls):
            shutil.copyfile(os.path.join(self.input_options['working_dir'],name),
                            os.path.join(directoryName,name.split('.')[0]+'_selected.yaml'))
        for i, item in enumerate(self.selected_S):
            
            np.save(os.path.join(directoryName,'S_selected_'+str(i+1)),self.selected_S[i])
            self.selected_Z[i].to_csv(os.path.join(directoryName,'selected_Z_'+str(i+1)+'.csv'),index=False)
            self.selected_Y[i].to_csv(os.path.join(directoryName,'selected_Y_'+str(i+1)+'.csv'),index=False)
            
    def write_experiment_yamls(self,directory): 
        #if not os.path.isdir(os.path.join(directory,'potential_yamls')):
        #    os.mkdir(os.path.join(directory,'potential_yamls'))
        with zipfile.ZipFile(os.path.join(directory,'psuedo_exps.zip'),'w') as f:
            for i,name in enumerate(self.experiment_yamls):
               f.write(os.path.join(self.input_options['working_dir'],name),arcname=name)
            
                       
            
        
    #def write_doe_results(self,directory):
    def write_selected_conds(self,directory):
        self.selected_conditions.to_csv(os.path.join(directory,'selected_conds.csv'),index=False)
        
    def write_object(self,output_dir):
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        
            
        self.write_exps_Y_to_file(output_dir,self.experiment_yamls)
        self.write_exps_Z_to_file(output_dir,self.experiment_yamls)
        self.write_exps_S_to_file(output_dir,self.experiment_yamls)
        
        self.write_input_options(os.path.join(output_dir,'doe_inputs.yaml'),self.input_options)
        self.write_priors(output_dir)
        if self.prior_experiments:
            p_experiments=pd.DataFrame(columns=['prior_experiment_yaml'])
            p_experiments['prior_experiment_yaml']=self.prior_experiments
            p_experiments.to_csv(os.path.join(output_dir,'prior_experiments.csv'),index=False)
        
        
        #Check if ranking was done by existence of selected conds
        if hasattr(self,'selected_conditions'):
            self.write_selected_exps(output_dir)
            self.write_selected_conds(output_dir)
        self.write_experiment_yamls(output_dir)
        shutil.make_archive(os.path.join(self.input_options['working_dir'],
                                         os.path.basename(os.path.normpath(output_dir))),
                            'zip',
                            output_dir)
        #shutil.make_archive(os.path.basename(os.path.normpath(output_dir)),
        #                    'zip', self.input_options['working_dir'])
        
        shutil.rmtree(output_dir)
    
    def read_object(self, input_file):
        unzipped_file=zipfile.ZipFile(input_file,'r')
        
        if not os.path.isdir(os.path.join(os.getcwd(),'temp_doe_obj_dir')):
            os.mkdir(os.path.join(os.getcwd(),'temp_doe_obj_dir'))
        unzipped_file.extractall(path=os.path.join(self.input_options['working_dir'],'temp_doe_obj_dir'))
        file_list=os.listdir(os.path.join(self.input_options['working_dir'],'temp_doe_obj_dir'))
        file_dir=os.path.join(self.input_options['working_dir'],'temp_doe_obj_dir')
        
        #Get input options
        if 'doe_inputs.yaml' in file_list:
            with open(os.path.join(file_dir,'doe_inputs.yaml')) as f:    
                self.input_options = yaml.load(f,Loader=yaml.FullLoader)
                
        #Get Priors
        self.covar_original=np.load(os.path.join(file_dir,'c_original.npy'))
        self.S_original=np.load(os.path.join(file_dir,'S_original.npy'))
        self.X_original=pd.read_csv(os.path.join(file_dir,'X_original.csv'))
        self.Z_original=pd.read_csv(os.path.join(file_dir,'Z_original.csv'))
        self.Y_original=pd.read_csv(os.path.join(file_dir,'Y_original.csv'))
        
        #Check if experiments exist and if so, import them.
        if os.path.exists(os.path.join(file_dir,'Z_matrices.csv')):
            print('Importing pseudo-experiment matrices')
            
        #Import yaml files
            with zipfile.Zipfile(os.path.join(file_dir,'pseudo_exps.zip')) as f:
                if not os.path.isdir(self.input_options['working_dir']):
                    os.mkdir(self.input_options['working_dir'])
                self.set_yaml_list(f.namelist())
                f.extractall(path=self.input_options['working_dir'])
            
            
            #Import experiment matrices
            Z_mats=pd.read_csv(os.path.join(file_dir,'Z_matrices.csv'))
            Y_mats=pd.read_csv(os.path.join(file_dir,'Y_matrices.csv'))
            S_mats=np.load(os.path.join(file_dir,'S_matrices.npz'))
            self.experiment_matrices=[]
            for i,file in enumerate(self.experiment_yamls):
                self.experiment_matrices.append({})
                self.experiment_matrices[i]['Z']=Z_mats['value',str(i+1)]
                self.experiment_matrices[i]['Y']=Y_mats['value',str(i+1)]
                self.experiment_matrices[i]['S']=S_mats[i]
        
        
        #Import selected files if they exist
        if os.path.exists(os.path.join(file_dir,'selected_Z_1.csv')):
            print('Importing selected matrices and conditions')
            self.selected_S=np.load(os.path.join(file_dir,'S_selected_1.npy'))
            self.selected_Y=pd.read_csv(os.path.join(file_dir,'selected_Y_1.csv'))
            self.selected_Z=pd.read_csv(os.path.join(file_dir,'selected_Z_1.csv'))
            
            #Read in the selected experiment yaml files
            self.selected_yamls=[]
            for i,name in enumerate(file_list):
                if self.input_options['yaml_output_name'] in name:
                    self.selected_yamls.append(name)
            self.selected_conditions=pd.read_csv(os.path.join(file_dir,'selected_conds.csv'))
        
        
    def add_experiments_to_priors(self, yaml_list=[]):
        #Placeholder code.  Add experiments to already existing object to go another iteration.
        pass
        
        

            
        
        
        
        
        
        
        
        
        