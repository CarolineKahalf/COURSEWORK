# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 04:24:31 2019

@author: Kahalf
"""

#This program reports information about DNA nucleotide sequences that may encode proteins.
from collections import Counter

nucleotide_mass={
        "A":135.128,
        "C":111.103,
        "G":151.128,
        "T":125.107,
        "-":100.00,
        }
#Above is a dictonary that assigns the mass value to each nucleotide
def read_input_file(file_name):
    z=open(file_name,"r")
    file_lines = z.readlines()
    z.close()
    data_file= [(file_lines[x].strip('\n'),file_lines[x+1].strip('\n')) for x in range(0,len(file_lines),2)]
    return data_file
#the function below counts the occurence of each nucleotide
def nucleotide_occurence(nucleotides):
    each_nucleotide= ["A","C","G","T"]
    count =[]
    
    for nuc in each_nucleotide:
        count.append(nucleotides.count(nuc))
        
    return count
    
#below is a function that comes up with a list of codons
def codon_list(nucleotide,nucleotides):
    list_of_codons=[]
    list_of_codons = [nucleotide.replace("-","")[i:i+3] for i in range(0, len(nucleotide.replace("-","")), 3)]
    return list_of_codons
#below is a function that calculates total mass of a DNA structure
def nucleotide_mass_function(nucleotides):
    total=0
    each = ["A","C","G","T","-"]
    for nuc in each:
        count = nucleotides.count(nuc)
        total = total +(count* nucleotide_mass[nuc])
    return total 

#below is a function that computes the mass percentage of each nucleotide in a codon
def computation_of_mass_percentages(nucleotide,nucleotides,total):
    count=nucleotides.Counter(nucleotide)#counts the number of times a nucleotide occurs in a structure
    total_mass_per_nucleotide=nucleotide_mass[nucleotide]*count
    mass_percentage=(total_mass_per_nucleotide/total)*100

    return round(mass_percentage,1)
    
#below function helps determine whether a given structure of nucleotides satisfies the necessary conditions to be a protein    
def protein(codon_list,total_mass_per_nucleotide):
    if len(codon_list)>4:
        if codon_list[0]=="ATG":
            if codon_list[-1]=="TAA" or"TAG" or"TGA":
                if ((total_mass_per_nucleotide[1]+total_mass_per_nucleotide[2]))>30:
                    return"YES"
                else:
                    return"NO"
            else:
                return"NO"
        else:
            return"NO"
    else:
        return"NO"
#this is the main function of python
if __name__=='__main__':
    print("This program reports information about DNA and nucleotide sequences that may encode proteins")
    In_file=input("Please enter Input filename")
    Out_file=input("Please enter Output filename")
    f=open(Out_file,'w')
    input_data=read_input_file(In_file)
    for x, y in input_data:
        f.write('Region name: {}\n'.format(x))
        f.write('Nucleotides: {}\n'.format(y.upper()))
        count=nucleotide_occurence(y)
        f.write('Nuc.Counts: {}\n'.format(count[:4]))
        total=nucleotide_mass_function(y)
        f.write('Total Mass %: {} of {}\n'.format(round(total,1)))
        list_of_codons = codon_list()
        f.write('Codons List: {}\n'.format(codon_list))
        f.write('Is Protein?: {} \n\n'.format(protein))
        
    f.close()