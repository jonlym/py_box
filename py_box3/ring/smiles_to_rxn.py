# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 11:54:14 2016

@author: Jonathan Lym
This script converts SMILES code to a readable format. Note that the entries must be present in the dictionary for it to work.
"""

# Dictionary that converts between SMILES code and a custom format. In this
# dictionary, Cu is the composite atom representing a metal surface
chem_dict = {
    'O=C=O': 'CO2',
    '[HH]': 'H2',
    'C=O': 'CO',
    '[{Cu}]': 'CU',
    '[{Cu}H]': 'H',
    'C(=O)O[{Cu}]': 'HCOO(S)',
    'C(=O)([{Cu}])O': 'COOH(S)',
    'O(CO[{Cu}])[{Cu}]': 'H2CO2(S)',
    'O=CO': 'HCOOH(g)',
    'C(O[{Cu}])([{Cu}])O': 'HCOOH(S) (?)',
    'C([{Cu}])([{Cu}])(O)O': 'C(OH)2(S)',
    'C(=O)([{Cu}])[{Cu}]': 'CO(S) (?)',
    '[{Cu}]O': 'OH(S)',
    'C(O[{Cu}])O':'H2COOH(S)',
    'C([{Cu}])(O)O':'HC(OH)2(S)',
    'O=C[{Cu}]':'HCO(S)',
    'C(O[{Cu}])([{Cu}])[{Cu}]':'HCO(S) (?)',
    'C([{Cu}])([{Cu}])([{Cu}])O':'COH(S)',
    'O':'H2O(g)',
    'OCO':'H2C(OH)2(g)',
    'O(C[{Cu}])[{Cu}]':'H2CO(S) (?)',
    'C([{Cu}])([{Cu}])O':'HCOH(S)',
    'O=C':'H2CO(g)',
    'C([{Cu}])([{Cu}])([{Cu}])[{Cu}]':'C(S)',
    'C([{Cu}])O':'H2COH(S)',
    'O([{Cu}])C':'CH3O(S)',
    'C([{Cu}])([{Cu}])[{Cu}]':'HC(S)',
    'OC':'CH3OH',
    '[{Cu}]C[{Cu}]':'CH2(S)',
    '[{Cu}]C':'CH3(S)',
    'C':'CH4',
    '[{Cu}]O[{Cu}]': 'O(S)'
}

file_path = 'C:\\Users\Jonathan Lym\Google Drive\\UDel Documents\\UDel Research\RING\RING_install\RING_program\Projects'
#Do not put on the .txt file extension
input_files = ['MethanolPathway', 'MethanolPathwayShort', 'CH4Pathway', 'HCOOHPathway', 'H2COHPathway', 'H2COPathway', 'CO2_Hydrogenation_All_Reactions']

for input_file in input_files:
    print(("-"*15))
    print(("Processing file: %s" % input_file)) 
    output_file = input_file+'_trans'
    f_in = open("%s\%s.txt" % (file_path, input_file), 'r')
    lines = f_in.readlines()
    f_in.close()
    
    f_out = open("%s\%s.txt" % (file_path, output_file), 'w')
    for line in lines:
        #The line contains a reaction
        line = line.replace('\n','')
        out_line = ""
        if ">>" in line:
            # Split the smiles string into reactants and products
            rxn_line = line.split('>>')
            react_line = rxn_line[0].split('.')
            prod_line = rxn_line[1].split('.')
            
            #Translate the smiles code for reactants and products
            nReact = len(react_line)
            nProd = len(prod_line)
            for i in range(nReact):
                out_line += chem_dict[react_line[i]]
                if i != (nReact-1):
                    out_line += ' + '
                else:
                    out_line += ' --> '
                    
            for i in range(nProd):
                out_line += chem_dict[prod_line[i]]
                if i != (nProd-1):
                    out_line += ' + '
        else:
            out_line = line
        print(out_line)
        f_out.write("%s\n" % out_line)
    f_out.close()
    print(('Completed translating SMILES reactions from %s.txt.' % input_file))
    print('Result written to:')
    print(('%s\%s.txt' % (file_path, output_file)))
print("All done translations!")