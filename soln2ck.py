import cantera as ct
import os
import numpy as np
from decimal import Decimal

def write(gas,output_file_chem,outputfile_therm):
    def write_elements(output_file,gas):
        with open(output_file,'w') as f:
            f.write('ELEMENTS\n')
            for i,el in enumerate(gas.element_names):
                if i!=(len(gas.element_names)-1):
                    f.write(el+'  ')
                elif i==(len(gas.element_names)-1):
                    f.write(el+'\n')
            f.write('END\n')
        
        return

    def write_species(output_file,gas):
        with open(output_file,'a') as f:
            f.write('SPECIES\n')
            for i,sp in enumerate(gas.species_names):
                if (i+1)%5!=0 and i!=(len(gas.species_names)-1):
                    f.write(sp+'  ')
                elif (i+1)%5!=0 and i==(len(gas.species_names)-1):
                    f.write(sp+'\n')
                elif (i+1)%5==0:
                    f.write(sp+'\n')
            f.write('END\n')
        
        return

    def write_thermo(output_file,gas):
        with open(output_file,'w') as f:
            f.write('THERMO\n')
            f.write('0200.00  1000.00  6000.00\n')
            for i,sp in enumerate(gas.species()):
                line1=''
                line1=line1+sp.name
                line1=line1+(24-len(sp.name))*' '
                compstr=''
                for j,el in enumerate(sp.composition.keys()):
                    #print(str(int(sp.composition[el])))
                    compstr=compstr+el+(5-len(el)-len(str(int(sp.composition[el]))))*' '+str(int(sp.composition[el]))
                line1=line1+compstr
                print(sp)
                line1=line1+(44-len(line1))*' '+'G'
                line1=line1+(48-len(line1))*' '+str(round(sp.thermo.min_temp,2))
                line1=line1+(57-len(line1))*' '+str(round(sp.thermo.max_temp,2))
                #if len(sp.thermo.input_data['temperature-ranges'])==3:
                line1=line1+(66-len(line1))*' '+str(round(sp.thermo.coeffs[0],2))
                line1=line1+(79-len(line1))*' '+'1\n'   
                    
                lines=[line1]
                #if len(sp.thermo.coeffs)==14:
                #    tempcoeffs=list(sp.thermo.coeffs)+[0.00]
                #elif len(sp.thermo.coeffs)==15:
                #    tempcoeffs=list(sp.thermo.coeffs)
                for k in [0,1,2]:
                    templine=''
                    for l in [0,1,2,3,4]:
                        if l+5*k<14:
                            tempstr=(15-len("{:2.8E}".format(sp.thermo.coeffs[l+k*5+1])))*' '+"{:2.8E}".format(sp.thermo.coeffs[l+k*5+1])
                            templine=templine+tempstr
                            if l==4:
                                templine=templine+(79-len(templine))*' '+str(int(k+2))+'\n'
                        elif l+5*k==14:
                            templine=templine+(79-len(templine))*' '+str(int(k+2))+'\n'
                    if k==2:
                        templine=templine+'\n\n'
                    lines=lines+[templine]
                f.writelines(lines) 
        
        return

    def build_arrhenius(equation_object, equation_type):
                """
                Builds Arrhenius coefficient string
    
                :param equation_objects
                    cantera equation object1
                :param equation_type:
                    string of equation type
                """
                calories_constant = 4184.0
                coeff_sum = sum(equation_object.reactants.values())
                pre_exponential_factor = equation_object.rate.pre_exponential_factor
                temperature_exponent = equation_object.rate.temperature_exponent
                activation_energy = (equation_object.rate.activation_energy /
                                    calories_constant)
    
                if equation_type == 'ElementaryReaction':
                    if coeff_sum == 1:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor))
                    if coeff_sum == 2:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor*10**3))
                    if coeff_sum == 3:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor*10**6))
                if equation_type == 'ThreeBodyReaction':
                    if coeff_sum == 1:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor*10**3))
                    if coeff_sum == 2:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor*10**6))
    
                if (equation_type != 'ElementaryReaction'
                    and equation_type != 'ThreeBodyReaction'):
                    pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor))
    
                arrhenius = [pre_exponential_factor,
                            temperature_exponent,
                            activation_energy
                            ]
                return arrhenius

    def build_modified_arrhenius(equation_object, t_range):
                """
                Builds Arrhenius coefficient strings for high and low temperature ranges
    
                :param equation_objects
                    cantera equation object
                :param t_range:
                    simple string ('high' or 'low') to designate temperature range
                """
                calories_constant = 4184.0
                if t_range == 'high':
                    pre_exponential_factor = equation_object.high_rate.pre_exponential_factor
                    temperature_exponent = equation_object.high_rate.temperature_exponent
                    activation_energy = (equation_object.high_rate.activation_energy /
                                        calories_constant)
                    dictlength=0
                    try:
                        dictlength=list(equation_object.products.values())[0]
                        #reaclength=equation_object.reactants.values()[0]
                    except:
                        pass          
                    
                    prodcount=0
                    reaccount=0
                    for num in np.arange(len(list(equation_object.products.values()))):
                        prodcount=prodcount+list(equation_object.products.values())[num]
                    for num in np.arange(len(list(equation_object.reactants.values()))):
                        reaccount=reaccount+list(equation_object.reactants.values())[num]
                        
                        
                    if prodcount==1.0 and reaccount==1.0:
                        not_single_and_equal=False
                    else:
                        not_single_and_equal=True
                    if len(equation_object.products) == 1 and dictlength==1.0 and not_single_and_equal:
                        
                        
                            pre_exponential_factor = str(
                                            '{:.5E}'.format(pre_exponential_factor*10**3))
                        
                    else:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor))
                    arrhenius_high = [pre_exponential_factor,
                                        temperature_exponent,
                                        activation_energy
                                        ]
                    return arrhenius_high
    
                if t_range == 'low':
                    pre_exponential_factor = equation_object.low_rate.pre_exponential_factor
                    temperature_exponent = equation_object.low_rate.temperature_exponent
                    activation_energy = (equation_object.low_rate.activation_energy /
                                        calories_constant)
                    dictlength=0
                    try:
                        dictlength=list(equation_object.products.values())[0]
                    except:
                        pass
                    
                    prodcount=0
                    reaccount=0
                    for num in np.arange(len(list(equation_object.products.values()))):
                        prodcount=prodcount+list(equation_object.products.values())[num]
                    for num in np.arange(len(list(equation_object.reactants.values()))):
                        reaccount=reaccount+list(equation_object.reactants.values())[num]
                        
                        
                    if prodcount==1.0 and reaccount==1.0:
                        not_single_and_equal=False
                    else:
                        not_single_and_equal=True
                    
                    
                    if len(equation_object.products) == 1 and dictlength==1.0 and not_single_and_equal:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor*10**6))
                    else:
                        pre_exponential_factor = str(
                                        '{:.5E}'.format(pre_exponential_factor*10**3))
                    arrhenius_low = [pre_exponential_factor,
                                    temperature_exponent,
                                    activation_energy
                                    ]
                    return arrhenius_low
            
            
    def write_reactions(output_file,gas):
        
        rxn_lines=['REACTIONS\n']
        for i,rxn in enumerate(gas.reactions()):
            if type(rxn).__name__=='ElementaryReaction':
                rxn_lines=rxn_lines+build_elementary_rxn(rxn,i)
            elif type(rxn).__name__=='ThreeBodyReaction':
                rxn_lines=rxn_lines+build_threebody_rxn(rxn,i)           
            elif type(rxn).__name__=='FalloffReaction':
                rxn_lines=rxn_lines+build_falloff_rxn(rxn,i)        
            elif type(rxn).__name__=='PlogReaction':
                rxn_lines=rxn_lines+build_plog_rxn(rxn,i)
        rxn_lines=rxn_lines+['END\n']
        with open(output_file,'a') as f:
            f.writelines(rxn_lines)
        
        return

    def build_elementary_rxn(rxn,i):
        rxn_lines=[]
        data=build_arrhenius(rxn,type(rxn).__name__)
        #print(data)
        for i,item in enumerate(data):
            data[i]=float(item)
            if i==0:
                data[i]="{:2.6E}".format(float(item))
            elif i==1:
                data[i]=str(round(data[i],6))
            elif i==2:
                data[i]=str(round(data[i],6))
        line=rxn.equation.replace(' ','')
        line=line+(60-len(line))*' '
        line=line+data[0]+'   '+data[1]+'   '+data[2]+'\n'
        rxn_lines.append(line)
        if rxn.duplicate:
            rxn_lines.append('DUPLICATE\n')
        return rxn_lines

    def build_threebody_rxn(rxn,i):
        rxn_lines=[]
        data=build_arrhenius(rxn,type(rxn).__name__)
        for i,item in enumerate(data):
            data[i]=float(item)
            if i==0:
                data[i]="{:2.6E}".format(float(item))
            elif i==1:
                data[i]=str(round(data[i],6))
            elif i==2:
                data[i]=str(round(data[i],6))
        line=rxn.equation.replace(' ','')
        line=line+(60-len(line))*' '
        line=line+data[0]+'   '+data[1]+'   '+data[2]+'\n'
        rxn_lines.append(line)
        if hasattr(rxn,'efficiencies'):
            if rxn.efficiencies.keys():
                eff=' '
                for i,el in enumerate(rxn.efficiencies.keys()):
                    eff=eff+el+'/'+str(rxn.efficiencies[el])+'/'
                    if i==(len(rxn.efficiencies.keys())-1):
                        eff=eff+'\n'
                rxn_lines=rxn_lines+[eff]
        if rxn.duplicate:
            rxn_lines.append('DUPLICATE\n')
        return rxn_lines

    def build_falloff_rxn(rxn,i):
        rxn_lines=[]
        data_high=build_modified_arrhenius(rxn,'high')
        for i,item in enumerate(data_high):
            data_high[i]=float(item)
            if i==0:
                data_high[i]="{:2.6E}".format(float(item))
            elif i==1:
                data_high[i]=str(round(data_high[i],6))
            elif i==2:
                data_high[i]=str(round(data_high[i],6))
        if hasattr(rxn,'low_rate'):
            data_low=build_modified_arrhenius(rxn,'low')
            for i,item in enumerate(data_low):
                data_low[i]=float(item)
                if i==0:
                    data_low[i]="{:2.6E}".format(float(item))
                elif i==1:
                    data_low[i]=str(round(data_low[i],6))
                elif i==2:
                    data_low[i]=str(round(data_low[i],6))
        line=rxn.equation.replace(' ','')
        line=line+(60-len(line))*' '
        line=line+data_high[0]+'   '+data_high[1]+'   '+data_high[2]+'\n'
        rxn_lines=rxn_lines+[line]
        if hasattr(rxn,'low_rate'):
            low=' LOW/'
            low=low+(60-len(low))*' '
            low=low+data_low[0]+'   '+data_low[1]+'   '+data_low[2]+'\n'
            rxn_lines=rxn_lines+[low]
        if hasattr(rxn,'falloff'):
            if rxn.falloff.falloff_type=='Troe':
                troe_line=' TROE/  '+"{:.4E}".format(rxn.falloff.parameters[0])
                troe_line=troe_line+'  '+"{:.4E}".format(rxn.falloff.parameters[1])
                troe_line=troe_line+'  '+"{:.4E}".format(rxn.falloff.parameters[2])
                if rxn.falloff.parameters[2]!=0:
                    troe_line=troe_line+'  '+"{:.4E}".format(rxn.falloff.parameters[3])
                troe_line=troe_line+'/\n'
                rxn_lines=rxn_lines+[troe_line]
        if hasattr(rxn,'efficiencies'):
            if rxn.efficiencies.keys():
                eff=' '
                for i,el in enumerate(rxn.efficiencies.keys()):
                    eff=eff+el+'/'+str(rxn.efficiencies[el])+'/'
                    if i==(len(rxn.efficiencies.keys())-1):
                        eff=eff+'\n'
                rxn_lines=rxn_lines+[eff]
        if rxn.duplicate:
            rxn_lines.append('DUPLICATE\n')
        return rxn_lines
        

    def build_plog_rxn(rxn,i):
        rxn_lines=[]
        #f.write('#  Reaction '+str(m)+'\n')
        line=rxn.equation.replace(' ','')
        rxn_lines=rxn_lines+[line]
        #f.write('pdep_arrhenius(\''+equation_object.equation+'\',\n')
        n=len(rxn.rates)
        calories_constant = 4184.0
        count=1
        
        for i,plog in enumerate(rxn.rates):
            plog_line='PLOG / '
            #print(plog)
            #pre_exponential_factor=plog[1].pre_exponential_factor
            #temperature_exponent=plog[1].temperature_exponent
            activation_energy=plog[1].activation_energy/calories_constant
            convertedPressure=float(plog[0])*9.86923e-6
            if sum(rxn.reactants.values())==1.0:
                #f.write('               [('+'%s' % float('%.5g' % convertedPressure)+', \'atm\'), ')
                #f.write('%.7e' % Decimal(str(plog[1].pre_exponential_factor)))
                plog_line=plog_line+"{:.6E}".format(convertedPressure)+'   '
                plog_line=plog_line+"{:.6E}".format(plog[1].pre_exponential_factor)+'   '
            elif sum(rxn.reactants.values())==2.0:
                #f.write('               [('+'%s' % float('%.5g' % convertedPressure)+', \'atm\'), ')
                #f.write('%.7e' % Decimal(str(plog[1].pre_exponential_factor*1000.0)))
                plog_line=plog_line+"{:.6E}".format(convertedPressure)+'   '
                plog_line=plog_line+"{:.6E}".format(plog[1].pre_exponential_factor*1000.0)+'   '
            #f.write(', '+'%s' % float('%.3g' % plog[1].temperature_exponent))
            #f.write(', '+'%s' % float('%.5g' % activation_energy))
            #f.write(']')
            plog_line=plog_line+'%s' % float('%.5g' % plog[1].temperature_exponent)+'   '
            plog_line=plog_line+'%s' % float('%.6g' % activation_energy)+' /\n'
            if count==n:
                rxn_lines[0]=rxn_lines[0]+(60-len(rxn_lines[0]))*' '
                if sum(rxn.reactants.values())==1.0:
                    rxn_lines[0]=rxn_lines[0]+"{:.3E}".format(plog[1].pre_exponential_factor)+'   '
                elif sum(rxn.reactants.values())==2.0:
                    rxn_lines[0]=rxn_lines[0]+"{:.3E}".format(plog[1].pre_exponential_factor*1000.0)+'   '
                rxn_lines[0]=rxn_lines[0]+'%s' % float('%.5g' % plog[1].temperature_exponent)+'   '
                rxn_lines[0]=rxn_lines[0]+'%s' % float('%.6g' % activation_energy)+'\n'
            #    if rxn.duplicate is True:
            #        f.write(', options = \'duplicate\')\n\n')
            #    else:
            #        f.write(')\n\n')
            #elif count<n:
            #    f.write(',\n')
            count=count+1
            rxn_lines=rxn_lines+[plog_line]
        
        if rxn.duplicate:
            rxn_lines.append('DUPLICATE\n')
            
        return rxn_lines


    
    write_elements(output_file_chem,gas)
    write_species(output_file_chem,gas)
    write_thermo(outputfile_therm,gas)
    write_reactions(output_file_chem,gas)
    

    return