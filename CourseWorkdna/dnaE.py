from collections import Counter
mass_percentages = {'a': 135.128,'c': 111.103,'g': 151.128,'t': 125.107,'-': 100.00}

def read_input_file(filename):
    f = open(filename,'r')
    data = f.readlines()
    f.close() 
    list_of_tuples = [(data[r].strip('\n'),data[r+1].strip('\n')) for r in range(len(data)) if r%2 == 0]
    return list_of_tuples
    
def is_protein(codon_list,masses):
    if len(codon_list) > 4:
        if codon_list[0] == "ATG":
            if (masses[1]+masses[2]) >= 30 :
                if codon_list[-1] in ["TAA","TGA","TAG"]:
                    return "YES"
                else:
                    return "NO"
            else:
                return "NO"
        else:
            return "NO"
    else:
        return "NO"


def build_codon_list(nucleotide):
    list_of_codons = [nucleotide.replace("-","").upper()[i:i+3] for i in range(0, len(nucleotide.replace("-","")), 3)]
    
    return list_of_codons

def compute_mass_percentage(nucleotide, nucleotides, total):
    count = nucleotides.count(nucleotide)
    num = mass_percentages[nucleotide]*count
    percent = (num/total) *100
    return round(percent,1)

def generate_nuc_count_and_parecentage(y):
    total = 0.0  
    compute_mass_list = []
    nuc_occur = []
    for i in ['a', 'c', 'g', 't', '-']:
        new_y = y.lower()
        count =  new_y.count(i)
        nuc_occur.append(count)
        total = total + (mass_percentages[i]*count)
    
    for x in ['a', 'c', 'g', 't']:
            compute_mass_list.append(compute_mass_percentage(x, new_y, total))
    
    return total,compute_mass_list,nuc_occur

if __name__ == '__main__':
    print('This Program reports information about DNA \n'
          + 'nucleotides sequences that may encode proteins. \n')
    inputfile = input('Input filename?')
    outputfile = input('Output file name?')
    input_data = read_input_file(inputfile)
    f = open(outputfile, 'w')
    for x, y in input_data:
        total,compute_mass_list,nuc_occur = generate_nuc_count_and_parecentage(y)
        f.write('Region name: {}\n'.format(x))
        f.write('Nucleotides: {}\n'.format(y.upper()))
        f.write('Nuc.Counts: {}\n'.format(nuc_occur[:4]))
        f.write('Total Mass %: {} of {}\n'.format(compute_mass_list,round(total,1)))
        codon_list = build_codon_list(y.upper())
        f.write('Codons List: {}\n'.format(codon_list))
        f.write('Is_Protein?: {} \n\n'.format(is_protein(codon_list,compute_mass_list)))
    f.close()
        
        
    