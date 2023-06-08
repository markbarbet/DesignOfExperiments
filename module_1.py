# -*- coding: utf-8 -*-
"""
Created on Sat May 15 15:19:44 2021

@author: Mark Barbet
"""

from asyncio.subprocess import DEVNULL, STDOUT
from cProfile import run
#from msilib.schema import Directory
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
import copy
import doe_object as dobj
from joblib import Parallel, delayed
import time
import soln2ck as w
from os.path import exists
import subprocess
import platform

def get_matrices_opensmoke(arg):
        file=arg[0]
        working_dir=arg[1]
        cti_file=arg[2]
        reaction_uncertainty_csv=arg[3]
        mix_species=arg[-2]
        #print(os.path.splitext(cti_file))
        ckfile=os.path.join(os.path.splitext(cti_file)[0],'.inp')
        num_rxns=arg[-1]
        #run_opensmoke_preprocessor(working_dir,ckfile)
        return run_simulation_opensmoke(file,working_dir,ckfile,reaction_uncertainty_csv,num_rxns)
    
def run_opensmoke_preprocessor(wdir,ckfile):
    lines=[]
    lines.append('Dictionary CHEMKIN_PreProcessor\n')
    lines.append('{\n')
    lines.append('\t\t@Thermodynamics\t'+os.path.join(wdir,'therm.dat')+';\n')
    lines.append('\t\t@Kinetics\t\t'+os.path.join(wdir,ckfile)+';\n')
    lines.append('\t\t@Output\t\t'+os.path.join(wdir,'kineticsFolder')+';\n')
    lines.append('}\n')
    with open(os.path.join(wdir,'preprocess.dic'),'w') as f:
        f.writelines(lines)
    
    lines=[]
    if platform.system()=='Windows':
        lines.append('\"%OPENSMOKEPP_EXE_FOLDER%\\OpenSMOKEpp_CHEMKIN_PreProcessor.exe\" --input "'+os.path.join(wdir,'preprocess.dic')+'"')
        with open(os.path.join(wdir,'preprocess.bat'),'w') as f:
            f.writelines(lines)
    elif platform.system()=='Linux':
        lines.append('#!/bin/sh\n')
        lines.append('OpenSMOKEpp_CHEMKIN_PreProcessor.sh --input "'+os.path.join(wdir,'preprocess.dic')+'"')
        with open(os.path.join(wdir,'preprocess.sh'),'w') as f:
            f.writelines(lines)
    #print(os.path.join(wdir,'preprocess.bat'))
    if platform.system()=='Windows':
        subprocess.call(os.path.join(wdir,'preprocess.bat'),stdout=DEVNULL,stderr=STDOUT)
    elif platform.system()=='Linux':
        subprocess.call(['chmod','-R','+x',wdir],stdout=DEVNULL,stderr=STDOUT)
        subprocess.call(os.path.join(wdir,'preprocess.sh'),stdout=DEVNULL,stderr=STDOUT)
    
    
def run_simulation_opensmoke(yaml,working_dir,ckfile,reaction_uncertainty_csv,num_rxns):
    conditions,observables,spec_list,target_unc,physical_uncertainties=get_conds_yaml(yaml,working_dir)    
    run_id=build_opensmoke_dirs(working_dir,ckfile,yaml)
    bats=[]
    ofolder=os.path.join(os.path.join(working_dir,'simulation_'+str(run_id)),'Output')
    instfolder=os.path.join(working_dir,'simulation_'+str(run_id))
    inputfile,batfile,subdir=write_opensmoke_inp(conditions,working_dir,ckfile,run_id,observables,'input.dic',ofolder,kineticsens=True)
    bats.append(batfile)
    inps=[]
    inps.append(os.path.join(subdir,inputfile))
    #print(conditions.keys(),'looooooool')
    for i,item in enumerate(conditions.keys()):
        if item=='temperature':
            temperatureconds=copy.deepcopy(conditions)
            temperatureconds['temperature']=1.01*temperatureconds['temperature']
            tempinp,tempbatfile,tempsubdir=write_opensmoke_inp(temperatureconds,working_dir,ckfile,run_id,observables,'temperature_input.dic',os.path.join(instfolder,'TempOutputs'))
            bats.append(tempbatfile)
            inps.append(os.path.join(subdir,tempinp))
            
        elif item=='pressure':
            pressureconds=copy.deepcopy(conditions)
            pressureconds['pressure']=1.01*pressureconds['pressure']
            tempinp,tempbatfile,tempsubdir=write_opensmoke_inp(pressureconds,working_dir,ckfile,run_id,observables,'pressure_input.dic',os.path.join(instfolder,'PressOutputs'))
            bats.append(tempbatfile)
            inps.append(os.path.join(subdir,tempinp))
            
        elif item not in ['volume','residence-time','temperature','pressure','AR','Ar','ar','HE','He','he']:
            xj=conditions[item]
            delxj=0.01*xj
            newxj=np.divide(np.multiply(xj+delxj,1.0-xj),1.0-xj-delxj)
            specs=[]
            vals=[]
            for j,key in enumerate(conditions.keys()):
                if key not in ['residence-time','temperature','pressure','volume'] and key!=item:
                    specs.append(key)
                    vals.append(conditions[key])
            specs.append(item)
            vals.append(newxj)
            #print(vals)
            valsum=np.sum(np.array(vals))
            vals=np.divide(np.array(vals),valsum)
            specconds=copy.deepcopy(conditions)
            for j,item in enumerate(specs):
                specconds[item]=vals[j]
                    
            tempinp,tempbatfile,tempsubdir=write_opensmoke_inp(specconds,working_dir,ckfile,run_id,observables,'X_'+item+'_input.dic',os.path.join(instfolder,'X'+str(i-3)+'Outputs'))
            bats.append(tempbatfile)
            inps.append(os.path.join(subdir,tempinp))
            
        elif item=='residence-time':
            resconds=copy.deepcopy(conditions)
            resconds['residence-time']=1.01*resconds['residence-time']
            tempinp,tempbatfile,tempsubdir=write_opensmoke_inp(resconds,working_dir,ckfile,run_id,observables,'res_input.dic',os.path.join(instfolder,'ResOutputs'))
            bats.append(tempbatfile)
            inps.append(os.path.join(subdir,tempinp))
            
    for i,file in enumerate(bats):
        run_opensmoke_jsr(file,subdir,inps[i])
            
    
    
    return build_matrices_opensmoke(working_dir,run_id,yaml,observables,num_rxns,reaction_uncertainty_csv,spec_list,target_unc,physical_uncertainties)

def build_matrices_opensmoke(wdir,run_id,yaml,observables,num_rxns,uncertainty_csv,mix_species,target_unc,physical_uncertainties):
    exclude=''
    
    var_names=[]
    folder=os.path.join(wdir,'simulation_'+str(run_id))
    folder=os.path.join(folder,'Output')
    tempfolder=os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),'TempOutputs')
    pressfolder=os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),'PressOutputs')
    resfolder=os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),'ResOutputs')
    specfolders=[]
    for i,sp in enumerate(mix_species):
        if sp not in ['Ar','AR','ar','he','He','HE']:
            specfolders.append(os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),'X'+str(i+1)+'Outputs'))
    
          
    if os.path.exists(folder):
        #
        sens_file_list={}
        for i,sp in enumerate(observables):
            sens_file_list['Sensitivities.'+sp+'.xml']=os.path.exists(os.path.join(folder,'Sensitivities.'+sp+'.xml'))
        phys_sens_file_list={}
        phys_sens_file_list['TempOutput']=os.path.exists(os.path.join(tempfolder,'Output.out'))
        phys_sens_file_list['PressOutput']=os.path.exists(os.path.join(pressfolder,'Output.out'))
        phys_sens_file_list['ResOutput']=os.path.exists(os.path.join(resfolder,'Output.out'))
          
        for i,sp in enumerate(mix_species):
            if sp not in ['Ar','AR','ar','he','He','HE']:
                phys_sens_file_list['X'+str(i+1)+'Outputs']=os.path.exists(os.path.join(specfolders[i],'Output.out'))
        #print(sens_file_list)
        #print(phys_sens_file_list)
        if os.path.exists(os.path.join(folder,'Output.out')) and all(sens_file_list.values()) and all(phys_sens_file_list.values()):
            with open(os.path.join(folder,'FinalSummary.out'),'r') as tempfile:
                 summary=tempfile.readlines()
                 if 'inf' in summary[7].split()[-1] or 'nan' in summary[7].split()[-1]:
                    exclude=yaml
                 if 'inf' in summary[8].split()[-1] or 'nan' in summary[8].split()[-1]:
                    exclude=yaml
                 if 'nan' in summary[9].split()[-1] or 'inf' in summary[9].split()[-1]:
                    exclude=yaml
                 if 'nan' in summary[10].split()[-1] or 'inf' in summary[10].split()[-1]:
                    exclude=yaml
                 if 'nan' in summary[11].split()[-1] or 'inf' in summary[11].split()[-1]:
                    exclude=yaml
                 if 'nan' in summary[18].split()[-1] or 'inf' in summary[18].split()[-1]:
                    exclude=yaml
                 #if 'nan' in summary[19].split()[-1] or 'inf' in summary[19].split()[-1]:
                  #  exclude=yaml
                 #if 'nan' in summary[20].split()[-1] or 'inf' in summary[20].split()[-1]:
                  #  exclude=yaml
                 #if 'nan' in summary[21].split()[-1] or 'inf' in summary[21].split()[-1]:
                  #  exclude=yaml
                 #if 'nan' in summary[22].split()[-1] or 'inf' in summary[22].split()[-1]:
                  #  exclude=yaml
                 #print(summary)
            #Get solutions for observables
            data=pd.read_csv(os.path.join(folder,'Output.out'),delimiter=r"\s+",header=0)
            obs_vals=np.zeros(len(observables))
            new_cols=[]
            for i,string in enumerate(data.columns):
                new_cols.append(string.split('(')[0])
            data.columns=new_cols
            
            
            ###Build Y vector
            Zval=[]
            for i,spe in enumerate(observables):
                #obs_vals[i]=data[spe+'_x']
                obs_vals[i]=0.0
                var_names.append(spe+'_experiment0')
                Zval.append(spe+'_experiment0')
            pad_length=3*num_rxns+3+len(mix_species)
            obs_vals=np.append(obs_vals,np.zeros(pad_length))
            for i in range(num_rxns):
                var_names.append('A_'+str(i))
                Zval.append('A_'+str(i))
            for i in range(num_rxns):
                var_names.append('n_'+str(i))
                Zval.append('n_'+str(i))
            for i in range(num_rxns):
                var_names.append('Ea_'+str(i))
                Zval.append('Ea_'+str(i))
            var_names.append('T_experiment_0')
            Zval.append('T_experiment_0')
            var_names.append('P_experiment_0')
            Zval.append('P_experiment_0')
            for i,sp in enumerate(mix_species):
                var_names.append('X_'+str(i)+'_experiment_0')
                Zval.append('X_'+sp+'_experiment_0')
            var_names.append('R_experiment_0')
            Zval.append('R_experiment_0')
            Y=pd.DataFrame(columns=['value','ln_difference'])
            Y['value']=var_names
            Y['ln_difference']=obs_vals  
            
            ###Build Z vector
            #z_array=np.zeros(len(obs_vals))
            z_array=[]
            Z=pd.DataFrame(columns=['value','Uncertainty'])
            z_array=z_array+get_data_uncertainties(data,target_unc,observables)
            k_unc=pd.read_csv(os.path.join(wdir,uncertainty_csv))
            z_array=z_array+list(k_unc['Uncertainty A (unit)'])
            z_array=z_array+list(k_unc['Uncertainty N (unit)'])
            z_array=z_array+list(k_unc['Uncertainty Ea (unit)'])
            #print(physical_uncertainties)
            z_array.append(physical_uncertainties[0]['temperature']['rel'])
            z_array.append(physical_uncertainties[1]['pressure']['rel'])
            for i,sp in enumerate(mix_species):
                z_array.append(physical_uncertainties[i+3][sp]['rel'])
            z_array.append(physical_uncertainties[2]['residence-time']['rel']) 
            Z['value']=Zval
            Z['Uncertainty']=z_array
            Y.to_csv(os.path.join(wdir,'Y_'+str(run_id)+'.csv'),index=False)
            Z.to_csv(os.path.join(wdir,'Z_'+str(run_id)+'.csv'),index=False)
            ###Build S matrix
            #S=np.zeros((len(Z['value']),len(Z['value'])-len(observables)))
            #indices=np.arange(len(Z['value'])-len(observables))
            #S[indices+len(observables),indices]=1.0
            return {'excluded_yaml':exclude,'y_len':len(Y['value'])}
        else:
            exclude=yaml 
    
            return {'excluded_yaml':exclude}
    else:
            exclude=yaml 
    
            return {'excluded_yaml':exclude}

def get_data_uncertainties(data,target_unc,observables):
    uncertainties=[]
    for i,spe in enumerate(observables):
        tempabs=target_unc[spe]['abs']
        temprel=target_unc[spe]['rel']
        tempunc=np.sqrt((np.divide(tempabs,data[spe+'_x'][0]))**2+(temprel)**2)
        uncertainties.append(tempunc)
    return uncertainties

  

def get_conds_yaml(file,wdir):
    #temp=yaml.load(os.path.join(wdir,yaml),Loader=yaml.FullLoader)
    with open(os.path.join(wdir,file)) as f:
            temp = yaml.load(f,Loader=yaml.FullLoader)
            
    outputs={}
    species_list=[]
    phys_uncertainties=[]
    outputs['temperature']=float(temp['common-properties']['temperature']['value-list'][0])
    phys_uncertainties.append({'temperature':{'rel':temp['common-properties']['temperature']['relative-uncertainty']}})
    outputs['pressure']=float(temp['common-properties']['pressure']['value'])
    phys_uncertainties.append({'pressure':{'rel':temp['common-properties']['pressure']['relative-uncertainty']}})
    phys_uncertainties.append({'residence-time':{'rel':temp['apparatus']['residence-time']['relative-uncertainty']}})
    outputs['residence-time']=float(temp['apparatus']['residence-time']['value'])
    outputs['volume']=np.multiply(float(temp['apparatus']['reactor-volume']['value']),1000000.0)
    for i,spec in enumerate(temp['common-properties']['composition']):
        outputs[spec['species']]=spec['mole-fraction']
        if spec['species'] not in ['Ar','AR','ar','He','HE','he']:
            species_list.append(spec['species'])
            phys_uncertainties.append({spec['species']:{'rel':spec['relative-uncertainty']}})
    observables=[]
    obs_unc_dict={}
    for i,entry in enumerate(temp['datapoints']['mole-fraction']):
        #print(entry)
        observables.append(entry['targets'][0]['species'])
        obs_unc_dict[entry['targets'][0]['species']]={'abs':entry['targets'][0]['absolute-uncertainty'],
                                                   'rel':entry['targets'][0]['relative-uncertainty']}
        
        
        
        
    return (outputs,observables,species_list,obs_unc_dict,phys_uncertainties)

def write_opensmoke_inp(condsdict,wdir,ckfile,run_id,observables,outfile,outputfolder,kineticsens=False):
    lines=[]
    lines.append('Dictionary PerfectlyStirredReactor\n')
    lines.append('{\n')
    lines.append('\t\t@KineticsFolder\t'+os.path.join(wdir,'kineticsFolder')+';\n')
    lines.append('\t\t@Type\t\t\t\t\tIsothermal-ConstantPressure;\n')
    lines.append('\t\t@InletStatus\t\t\tinlet-mixture;\n')
    lines.append('\t\t@ResidenceTime\t\t\t'+str(condsdict['residence-time'])+' s;\n')
    lines.append('\t\t@Volume\t\t\t\t\t'+str(condsdict['volume'])+' cm3;\n')
    if kineticsens:
        lines.append('\t\t@SensitivityAnalysis\tsensitivity-options;\n')
    lines.append('\t\t@OdeParameters\t\t\tode-parameters;\n')
    lines.append('\t\t@Options\t\t\toutput-options;\n')
    lines.append('}\n')
    
    #lines.append('Dictionary DoE_PreProcessor\n')
    #lines.append('{\n')
    #lines.append('\t\t@Kinetics\t\t\t')
    
    
    lines.append('Dictionary inlet-mixture\n')
    lines.append('{\n')
    lines.append('\t\t@Temperature\t\t'+str(condsdict['temperature'])+' K;\n')
    lines.append('\t\t@Pressure\t\t'+str(condsdict['pressure']*ct.one_atm)+' Pa;\n')
    tempstr='\t\t@MoleFractions\t\t'
    for i,item in enumerate(condsdict.keys()):
        if item not in ['temperature','pressure','residence-time','volume']:
            tempstr=tempstr+' '+item+' '+str(condsdict[item])
    tempstr=tempstr+';\n'
    lines.append(tempstr)
    lines.append('}\n')
    
    lines.append('Dictionary sensitivity-options\n')
    lines.append('{\n')
    lines.append('\t@Type\t\t\t\t\tkinetic-constants;\n')
    lines.append('\t@DenseSolver\t\tEigen;\n')
    lines.append('\t@DenseFullPivoting\ttrue;\n')
    lines.append('\t@SubSteps\t\t\t5;\n')
    lines.append('\t@Species\t\t\t')
    for i,item in enumerate(observables):
        lines[-1]=lines[-1]+item+' '
    lines[-1]=lines[-1]+';\n'
    lines.append('}\n')
    
    lines.append('Dictionary ode-parameters\n')
    lines.append('{\n')
    lines.append('\t@OdeSolver\t\t\tOpenSMOKE;\n')
    lines.append('\t@RelativeTolerance\t\t1e-16;\n')
    lines.append('\t@AbsoluteTolerance\t\t1e-17;\n')
    lines.append('\t@MaximumNumberOfSteps\t\t50000;\n')
    lines.append('}\n')
    lines.append('Dictionary output-options\n')
    lines.append('{\n')
    lines.append('\t@OutputFolder\t\t'+outputfolder+';\n')
    lines.append('}\n')
    
    with open(os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),outfile),'w') as f:
        f.writelines(lines)
        
    #lines=[]
    #lines.append('\"%OPENSMOKEPP_EXE_FOLDER%\\OpenSMOKEpp_PerfectlyStirredReactor.exe\" --input input.dic')
    #with open(os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),'run.bat'),'w') as f:
    #    f.writelines(lines)
    #subprocess.call(os.path.join(os.path.join(wdir,'simulation_'+str(run_id)),'run.bat'))
    
    
    if platform.system()=='Windows':
        return (outfile,'run.bat',os.path.join(wdir,'simulation_'+str(run_id)))
    elif platform.system()=='Linux':
        return (outfile,'run.sh',os.path.join(wdir,'simulation_'+str(run_id)))

def run_opensmoke_jsr(batfile,wdir,inputfile):
    lines=[]
    if platform.system()=='Windows':
        
        lines.append('\"%OPENSMOKEPP_EXE_FOLDER%\\OpenSMOKEpp_PerfectlyStirredReactor.exe\" --input "'+inputfile+'"')
    elif platform.system()=='Linux':
        lines.append('#!/bin/sh\n')
        lines.append('OpenSMOKEpp_PerfectlyStirredReactor.sh --input "'+inputfile+'"')
    with open(os.path.join(wdir,batfile),'w') as f:
        f.writelines(lines)
    subprocess.call(['chmod','-R','+x',wdir],stdout=DEVNULL,stderr=STDOUT)
    #subprocess.call(os.path.join(wdir,batfile),stdout=DEVNULL,stderr=STDOUT)
    subprocess.call(os.path.join(wdir,batfile))


    
    with open(inputfile,'r') as f:
        r=f.readlines()
        for i,line in enumerate(r):
            if 'Dictionary output-options' in line:
                line_ind=i+2
        direc=r[line_ind].split()[-1].rstrip().rstrip(';')
    os.remove(os.path.join(direc,'Output.history'))

    return

def build_opensmoke_dirs(working_dir,ckfile,yaml):
    run_id=yaml.split('.')[0].split('_')[-1]
    if not os.path.exists(os.path.join(working_dir,'simulation_'+run_id)):
        os.makedirs(os.path.join(working_dir,'simulation_'+run_id))
    return run_id

def get_matrices_parallel(arg):
        file=arg[0]
        working_dir=arg[1]
        cti_file=arg[2]
        reaction_uncertainty_csv=arg[3]
        #print(arg)
        ckfile=os.path.join(os.path.splitext(cti_file),'.inp')
        if not exists(ckfile):
            tempgas=ct.Solution(cti_file)
            w.write(tempgas,ckfile,os.path.join(os.path.dirname(ckfile),'therm.dat'))
            del(tempgas)
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
            #print(MSI_instance.S_matrix[0:4,9:13])
            #print(working_directory)
            #print(files_to_include)
            #print(reaction_uncertainty_csv)
            #print(rate_constant_target_value_data)
            #print(newfile)
            S=MSI_instance.S_matrix
            #X=MSI_instance.X
            Z=MSI_instance.z_data_frame
            Y=MSI_instance.Y_data_frame
            #print(Z['value'][630:])
            #print(Y)
            os.remove(temp_cti)
            os.remove(os.path.splitext(temp_cti)[0]+'_updated.cti')
            newfile=None
            temp_cti=None
            original_cti=None
            files_to_include=None
            cti_file=None
            reaction_uncertainty_csv=None
            MSI_instance=None
            #print(yaml_list[0][0])
            #print(yaml_list[0][0].split('_')[-1].split('.')[0])
            np.save(os.path.join(working_directory,'S_'+yaml_list[0][0].split('_')[-1].split('.')[0]),S)
            
            np.savetxt(os.path.join(working_directory,'S_'+yaml_list[0][0].split('_')[-1].split('.')[0]+'.csv'),S,delimiter=',')
            Y.to_csv(os.path.join(working_directory,'Y_'+yaml_list[0][0].split('_')[-1].split('.')[0]+'.csv'),index=False)
            Z.to_csv(os.path.join(working_directory,'Z_'+yaml_list[0][0].split('_')[-1].split('.')[0]+'.csv'),index=False)
            
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
        
        return {'excluded_yaml':exclude}


class potential_experiments():
    
    def __init__(self,doe_obj:dobj.doe_object):
        self.input_options=doe_obj.input_options
        self.constructor_settings=self.input_options['constructor_settings']
        
        
    
    # def __init__(self,input_options:dict={'initialization':None,
    #                                       'experiment_type':'JSR',
    #                                       'observables':['H2','O2'],
    #                                       'observables_abs_uncertainties':[0.00005,0.00005],
    #                                       'observables_rel_uncertainties':[0.05,0.05],
    #                                       'temperature_range':[600,1150], 
    #                                       'temperature-uncertainty': 0.01,
    #                                       'pressure_range':[1.0],
    #                                       'pressure-uncertainty': 0.01,
    #                                       'residence_time':[1.2],
    #                                       'restime-uncertainty': 0.02,
    #                                       'mixture_species':{'H2':[0.001,0.02],'O2':[0.001,0.05]},
    #                                       'diluent':'Ar',
    #                                       'constructor_settings':{'intervals':{'species':{'all'},
    #                                                                            'pressure':[1.0],
    #                                                                            'temperature':50,
    #                                                                            'residence_time':0.5},
    #                                                               'method':'random_sample',
    #                                                               'random_sample_settings':{'n':[10,5]}},
    #                                       'experiment_yaml_template':'',
    #                                       'working_dir':'',
    #                                       'thermal-boundary':'isothermal',
    #                                       'mechanical-boundary':'constant pressure',
    #                                       'volume':0.000085,
    #                                       'yaml_output_name':'DoE_yaml_',
    #                                       'parallel-computing':True,
    #                                       'gridpoints':{}}):
    #     self.input_options=input_options
    #     self.constructor_settings=self.input_options['constructor_settings']
        if re.match('[Jj][Ss][Rr]',self.input_options['experiment_type']) or re.match('[Jj]et[- ][Ss]tirred[- ][Rr]eactor',self.input_options['experiment_type']):
            self.exp_type='JSR'
            self.conditions_list=self.JSR_conditions_constructor(self.input_options['constructor_settings'])

            self.yaml_file_list=self.JSR_yaml_constructor(self.conditions_list)
            
            if not self.input_options['parallel-computing']:
                self.matrices=self.get_matrices()  
                
            elif self.input_options['parallel-computing']:
                self.cores=self.input_options['cores']
                args=self.get_args()
                div_args=list(self.divide_list(args,10000))
                self.matrices=[]
                #for count, item in enumerate(div_args):
                
                
                ckfile=os.path.splitext(args[0][2])[0]+'.inp'
                
                if not exists(os.path.join(self.input_options['working_dir'],ckfile)):
                        tempgas=ct.Solution(os.path.join(self.input_options['working_dir'],args[0][2]))
                        w.write(tempgas,os.path.join(self.input_options['working_dir'],ckfile),os.path.join(os.path.dirname(os.path.join(self.input_options['working_dir'],ckfile)),'therm.dat'))
                        num_rxns=len(tempgas.reactions())
                        del(tempgas)
                start=time.time()
                #temp_mat=Parallel(n_jobs=self.cores)(delayed(get_matrices_parallel)(arg) for i,arg in enumerate(args))
                run_opensmoke_preprocessor(self.input_options['working_dir'],ckfile)
                with multiprocessing.Pool(processes=self.cores,maxtasksperchild=1) as pool:
                        temp_mat=pool.map(get_matrices_opensmoke,args)
                        pool.close()
                        pool.join()
                        self.matrices=temp_mat
                self.write_S_mats()
                #with multiprocessing.Pool(processes=self.cores) as pool:
                #    self.matrices=pool.map(get_matrices_parallel,args,chunksize=10)
                stop=time.time()
                print('{:.4f} s'.format(stop-start))
                #self.matrices=temp_mat
                included_files=[]
                included_matrices=[]
                #print(len(self.matrices))
                for i,mat in enumerate(self.matrices):
                    if mat['excluded_yaml']=='':
                        included_files.append(self.yaml_file_list[i])
                        #included_matrices.append(mat)
                        mat['S']=np.load(os.path.join(self.input_options['working_dir'],'S_'+str(i+1)+'.npy'))
                        mat['Y']=pd.read_csv(os.path.join(self.input_options['working_dir'],'Y_'+str(i+1)+'.csv'))
                        mat['Z']=pd.read_csv(os.path.join(self.input_options['working_dir'],'Z_'+str(i+1)+'.csv'))
                        included_matrices.append(mat)
                    loop_len=i
                self.yaml_file_list=included_files
                self.matrices=included_matrices
                doe_obj.set_matrices(self.matrices)
                doe_obj.set_yaml_list(self.yaml_file_list)
                #print(len(self.matrices))
                self.garbage_collection(self.input_options['working_dir'],loop_len)
                #return self.input_options
                    #print(self.yaml_file_list)
    def write_S_mats(self):
        #print(self.matrices)
        tempgas=ct.Solution(os.path.join(self.input_options['working_dir'],self.input_options['MSI_settings']['chemical_model']))
        for i,mat in enumerate(self.matrices):
            if mat['excluded_yaml']=='':
                y_len=mat['y_len']
                folder=os.path.join(self.input_options['working_dir'],'simulation_'+str(i+1))
                S=np.zeros((len(self.input_options['observables']),y_len-len(self.input_options['observables'])))
                #S=np.zeros((y_len,y_len-len(self.input_options['observables'])))
                indices=np.arange(y_len-len(self.input_options['observables']))
                #S[indices+len(self.input_options['observables']),indices]=1.0
                with open(os.path.join(self.input_options['working_dir'],self.input_options['yaml_output_name']+str(i+1)+'.yaml'),'r') as f:
                        config = yaml.load(f,Loader=yaml.FullLoader)
                        temp=float(config['common-properties']['temperature']['value-list'][0])
                        pres=float(config['common-properties']['pressure']['value'])
                        Xfrac={}
                        for k,specie in enumerate(config['common-properties']['composition']):
                            Xfrac[specie['species']]=specie['mole-fraction']
                tempgas.TPX=temp,pres,Xfrac
                forward_rates=tempgas.forward_rate_constants
                constants=[]
                with open(os.path.join(os.path.join(folder,'Output'),'Sensitivities.xml'),'r') as f:
                    ratedata=f.readlines()
                    constparam=[False,False]
                    for j,line in enumerate(ratedata):
                        if 'constant-parameters' in line and constparam[0] and not constparam[1]:
                            constparam[1]=True
                        
                        if constparam[0] and not constparam[1]:
                            constants.append(float(line.rstrip('\n')))
                            
                        if 'constant-parameters' in line and not constparam[0] and not constparam[1]:
                            constparam[0]=True
                            
                constants=np.array(constants)         
                #Get solutions for observables
                specdata=pd.read_csv(os.path.join(os.path.join(folder,'Output'),'Output.out'),delimiter=r"\s+",header=0)
                #obs_vals=np.zeros(len(observables))
                new_cols=[]
                for j,string in enumerate(specdata.columns):
                    new_cols.append(string.split('(')[0])
                specdata.columns=new_cols
                for j,sp in enumerate(self.input_options['observables']):
                    with open(os.path.join(os.path.join(folder,'Output'),'Sensitivities.'+sp+'.xml'),'r') as f:
                        data=f.readlines()
                    ksens=data[3].rstrip('\n').rstrip(' ').split(' ')
                    #print(ksens)
                    spec_X=float(specdata[sp+'_x'])
                    ksens=np.array(list(map(float,ksens)))
                    ksens=np.multiply(ksens,constants)
                    ksens=np.divide(ksens,spec_X)
                    #if i==258:
                    #    print(ksens)
                    Asens=ksens
                    nsens=ksens*np.log(temp)
                    Easens=np.divide(ksens,-temp)
                    #finalrow=np.empty((Asens.size+nsens.size+Easens.size),dtype=Asens.dtype)
                    #finalrow[0::3]=Asens
                    #finalrow[1::3]=nsens
                    #finalrow[2::3]=Easens
                    S[j,:len(ksens)]=Asens
                    S[j,len(ksens):2*len(ksens)]=nsens
                    S[j,2*len(ksens):3*len(ksens)]=Easens
                    #S[j,:3*len(ksens)]=finalrow
                physSensBlock=np.zeros((len(self.input_options['observables']),np.shape(S)[1]-3*len(ksens)))
                physLen=3+len(self.input_options['mixture_species'].keys())
                nomdata=pd.read_csv(os.path.join(os.path.join(os.path.join(self.input_options['working_dir'],'simulation_'+str(i+1)),'Output'),'Output.out'),delimiter=r"\s+",header=0)
                #obs_vals=np.zeros(len(self.input_options['observables']))
                directory=os.path.join(self.input_options['working_dir'],'simulation_'+str(i+1))
                new_cols=[]
                for j,string in enumerate(nomdata.columns):
                    new_cols.append(string.split('(')[0])
                nomdata.columns=new_cols
                for j in range(physLen):
                    ####Temp sen
                    if j==0:
                        tempdata=pd.read_csv(os.path.join(os.path.join(directory,'TempOutputs'),'Output.out'),delimiter=r"\s+",header=0)
                        tempcols=[]
                        
                        
                    
                    ####Pres sen
                    elif j==1:
                        tempdata=pd.read_csv(os.path.join(os.path.join(directory,'PressOutputs'),'Output.out'),delimiter=r"\s+",header=0)
                        tempcols=[]
                    ####Spec sens
                    elif j>1 and j<physLen-1:
                        tempdata=pd.read_csv(os.path.join(os.path.join(directory,'X'+str(j-1)+'Outputs'),'Output.out'),delimiter=r"\s+",header=0)
                        tempcols=[]
                    ####Res sen
                    elif j==physLen-1:
                        tempdata=pd.read_csv(os.path.join(os.path.join(directory,'ResOutputs'),'Output.out'),delimiter=r"\s+",header=0)
                        tempcols=[]
                        
                    for k,string in enumerate(tempdata.columns):
                        #print(string.split('(')[0])
                        tempcols.append(string.split('(')[0])
                    tempdata.columns=tempcols 
                    for k,sp in enumerate(self.input_options['observables']):
                            xk=np.log(nomdata[sp+'_x'])
                            xkprime=np.log(tempdata[sp+'_x'])
                            dif=np.subtract(xkprime,xk)
                            dif=dif/0.01
                            physSensBlock[k,j]=dif  
                S[0:len(self.input_options['observables']),3*len(ksens):]=physSensBlock
                np.save(os.path.join(self.input_options['working_dir'],'S_'+str(i+1)),S)
                np.savetxt(os.path.join(self.input_options['working_dir'],'S_'+str(i+1)+'.csv'),S,delimiter=',')

    def garbage_collection(self,working_dir,num):
        #Function to delete extraneous S,Y,Z files along with csv files and yaml files 
        #leftover from the "fake" experiment generation
        for i in range(num):

            #Delete S
            if os.path.exists(os.path.join(working_dir,'S_'+str(i+1)+'.npy')):
                os.remove(os.path.join(working_dir,'S_'+str(i+1)+'.npy'))
            #if os.path.exists(os.path.join(working_dir,'Y_'+str(i+1)+'.csv')):
            #    os.remove(os.path.join(working_dir,'Y_'+str(i+1)+'.csv'))
            #if os.path.exists(os.path.join(working_dir,'Z_'+str(i+1)+'.csv')):
            #    os.remove(os.path.join(working_dir,'Z_'+str(i+1)+'.csv'))


            #if os.path.exists(os.path.join(working_dir,'DoE_yaml_'+str(i+1)+'.yaml')):
            #    os.remove(os.path.join(working_dir,'DoE_yaml_'+str(i+1)+'.yaml'))
        
        for item in os.listdir(working_dir):
            if 'DoE_yaml_' in item and '.csv' in item:
                os.remove(os.path.join(working_dir,item))

                    
                    
    def divide_list(self,arglist,n):
        for i in range(0,len(arglist),n):
            yield arglist[i:i+n]

    
    def get_args(self):
        
        args=[]
        tempgas=ct.Solution(os.path.join(self.input_options['working_dir'],self.input_options['MSI_settings']['chemical_model']))
        numrxns=len(tempgas.reactions())
        for i,file in enumerate(self.yaml_file_list):
            args=args+[[file,
                       self.input_options['working_dir'],
                       self.input_options['MSI_settings']['chemical_model'],
                       self.input_options['MSI_settings']['reaction_uncertainty'],
                       list(self.input_options['mixture_species'].keys()),
                       numrxns]]
        
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
        
        working_directory = self.input_options['working_dir']
        
        cti_file = self.input_options['chemical_model']
        
        reaction_uncertainty_csv = self.input_options['reaction_uncertainty']
        
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
        
        if re.match('[Gg]rid[-_ ][Ll]og',settings['method']):
            dimension=0
            conds_dict={}
            var_list=[]
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                dimension=dimension+1
                points=self.input_options['gridpoints'][species]
                temparray=np.zeros(points)
                if self.input_options['mixture_species'][species][0]==0:
                    temparray[1:]=np.logspace(np.log10(1e-5),np.log10(self.input_options['mixture_species'][species][1]),num=points-1)
                    #temparray[1:]=np.power(10,temparray[1:])
                elif self.input_options['mixture_species'][species][0]!=0:
                    temparray=np.logspace(np.log10(self.input_options['mixture_species'][species][0]),np.log10(self.input_options['mixture_species'][species][1]),num=points)
                    #temparray=np.power(10,temparray)
                conds_dict[species]=copy.deepcopy(temparray)
            
            conds_dict['temperatures']=np.linspace(self.input_options['temperature_range'][0],self.input_options['temperature_range'][1],
                                                   num=self.input_options['gridpoints']['temperature'])
            conds_dict['pressures']=np.linspace(self.input_options['pressure_range'][0],self.input_options['pressure_range'][1],
                                                   num=self.input_options['gridpoints']['pressure'])
            conds_dict['restimes']=np.linspace(self.input_options['residence_time'][0],self.input_options['residence_time'][1],
                                                   num=self.input_options['gridpoints']['restime'])
            
            keys=[]
            for i,species in enumerate(list(self.input_options['mixture_species'].keys())):
                var_list.append(conds_dict[species])
                keys.append(species)
            keys.append('temperatures')
            keys.append('pressures')
            keys.append('restimes')
            var_list.append(conds_dict['temperatures'])
            var_list.append(conds_dict['pressures'])
            var_list.append(conds_dict['restimes'])
            temparray=np.array(np.meshgrid(*var_list)).T.reshape(-1,len(var_list)).T
            for i,value in enumerate(keys):
                conditions[value]=temparray[i]
            conditions['index']=np.arange(1,len(temparray[0])+1)
            conditions.to_csv(os.path.join(os.getcwd(),'test_conditions.csv'),index=False)  
            dimension=dimension+3
            
            
            
                
        
        elif re.match('[Ii]ntervals',settings['method']):
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
                if np.isneginf(low_limit):
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
                if np.isneginf(low_limit):
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

    
                
            
            
                
                
                
            
                    
                    
                    
                    
                
                
