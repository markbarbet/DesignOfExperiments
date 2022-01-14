import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt





class doe_object():
    
    
    def __init__(self, parsed_inputs):
        
        #Initialization for class.  Takes parsed input file.
        self.input_options=parsed_inputs
        
        
        self.total_new_exps=0
        self.conditions=None
        self.excluded_files=None
        self.new_uncertainties={}
        self.selected_conditions=None
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
        self.Z_original=None
        self.prior_experiments=[]
        
        
    def plot_convergence(self,log=True):
        if self.sorted_uncertainties == None:
            print('Cannot plot convergence until algorithm to find new uncertainties is run')
        elif self.sorted_uncertainties != None and log:
            #Placeholder code
            pass
        elif self.sorted_uncertainties != None and not log:
            #Placeholder code
            pass
        
    def plot_sensitivities(self,reactions=[],n_rxns=0):
        #Code to plot sensitivities.  Placeholder code.
        pass
        
    def write_object(self, output_file):
        #Placeholder code.  Function to write the object to files.
        pass
    
    def read_object(self, input_file):
        #Placeholder code.  Function to import data about DoE object from file.
        pass
        
    def add_experiments_to_priors(self, yaml_list=[]):
        #Placeholder code.  Add experiments to already existing object to go another iteration.
        pass
        
        

            
        
        
        
        
        
        
        
        
        