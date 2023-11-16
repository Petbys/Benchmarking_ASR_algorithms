import os
import sys
import argparse
import pandas as pd

cwd = os.getcwd()
# load data
print(cwd)
snp_info = pd.read_excel("/Users/petterbystrom/Documents/Applied Bioinformatics/snp_spreadsheet.xlsx", "Supp. table 14", header=1, index_col=0)
taxa_info = pd.read_excel("/Users/petterbystrom/Documents/Applied Bioinformatics/snp_spreadsheet.xlsx", "Supp. table 13", header=1)
headers_path="/Users/petterbystrom/Documents/Applied Bioinformatics/Headers.txt"
with open(headers_path, 'r') as file:
    file_contents = file.readlines()
Headers_from_fasta = sorted([line.replace('>', '').replace(':', '').replace('\n','') for line in file_contents])
print("header",Headers_from_fasta)
#print(Headers_from_fasta)
enteries = set(snp_info.columns.to_list()) & set(taxa_info["Genome ID"].to_list()) & set(Headers_from_fasta)
print(len(enteries))
unique_to_snpinfo = [element for element in snp_info.columns.to_list() if element not in enteries]
unique_to_taxainfo = [element for element in set(taxa_info["Genome ID"].to_list()) if element not in enteries]
unique_to_Headersfromfasta = [element for element in Headers_from_fasta if element not in enteries]

print(unique_to_snpinfo)
print(unique_to_taxainfo)
print(unique_to_Headersfromfasta)
print("\n")

snp_unique =['CO92 Ref.', 'BSK001', 'BSK003', 'BSK001-003', 'London-ES', 'Barcelona3031', 'London-6330', 'SC2-SPN7', 'CHE001', 'Y.pseudotuberculosis IP32953 (outgroup)']
taxa_unique = ['CHE1', 'London  6330', 'London 8124/8291/11972', 'Barcelona 3031', 'SC2_SPN7']
fasta_unique = ['AGU007.B0102', 'AGU010.B0102', 'AGU025.B0102', 'BED024.A0102', 'BED028.A0102', 'BED030.A0102', 'BED034.A0102', 
                'BRA001.A0101', 'BSK001-003.A0101.A0102.A0103-malt', 'Barcelona3031', 
                'CHE.C1-C2-3-C4', 'COLC1-COLC2a_COLC2b', 'ELW098_JK1548', 'LAI009.A0101', 'LBG002.A0101', 
                'London_EastSmithfield_8124_8291_11972', 'London_StMaryGraces_6330', 'MAN008.B0101', 'NAB003.B0101', 'NMS002.A0101', 
                'Reference_gi|16120353|ref|NC_003143.1| Yersinia pestis CO92 chromosome, complete genome', 'SC2_SPN7', 'STA001.A0101', 'STN002.A0101', 'STN007.A0101', 'STN008.A0101', 'STN013.A0101', 'STN014.A0101', 'STN019.A0101',
                  'STN020.A0101', 'STN021.A0101', 'TRP002.A0101-02', 'outgroup_Y.pseudo']
modern_ancient_dict={}
print(len(fasta_unique))
lineage_to_modern_ancient = dict(zip(taxa_info['Genome ID'], taxa_info['Modern/Ancient']))
#print(lineage_to_modern_ancient)
#print(snp_info.describe())
rows_without_N = snp_info[~snp_info.apply(lambda row: row.astype(str).str.contains('N', case=False)).any(axis=1)]

print(rows_without_N)
#print(snp_info)
#print(snp_info.columns)

count_of_n_df = pd.DataFrame(snp_info.apply(lambda x: x.eq('N').sum()), columns=['count_of_n'])
#print(count_of_n_df)

#print(sorted(list(count_of_n_df['count_of_n'])))
Ancient_df = count_of_n_df[count_of_n_df.apply(lambda row: row>108).any(axis=1)]
#print(sorted(list(Ancient_df['count_of_n'])))
#print(len(list(Ancient_df['count_of_n'])))
#print(Ancient_df)
Headers_df = sorted(count_of_n_df.index.to_list())
#print(Headers_df)
listofmismatch=[]
"""
for i in range(len(Headers_from_fasta)):
    print(Headers_df[i],Headers_from_fasta[i])
    if Headers_df[i][0:5]== Headers_from_fasta[i][0:5]:
        print("Yes")
    else:
        print("No")
        
        listofmismatch.append([Headers_df[i],Headers_from_fasta[i]])
"""
#print(listofmismatch)
#print(len(listofmismatch))
#common_elements = set(Headers_df) & set(Headers_from_fasta)
#print(len(common_elements))
#list1_first_6 = [element[0:7] for element in Headers_from_fasta]
#list2_first_6 = [element[0:7] for element in Headers_df]

# Find common elements based on the first 6 letters
#common_elements = set(list1_first_6) & set(list2_first_6)
#print(common_elements)
#print("\n")
#unique_to_list1 = [element for element in Headers_df if element not in Headers_from_fasta]
'''
# Find elements unique to list2
unique_to_list2 = [element for element in Headers_from_fasta if element not in Headers_df]
print(unique_to_list1)
list1_first_6 = [element[0:6] for element in unique_to_list1]
list2_first_6 = [element[0:6] for element in unique_to_list2]
common_elements2 = (set(list1_first_6) & set(list2_first_6))
all_shared_elements = common_elements2 & common_elements
print(len(all_shared_elements))
#not_in_set1 = [element[0:6] for element in Headers_from_fasta if element not in all_shared_elements]
#not_in_set2 = [element[0:6] for element in Headers_df if element not in all_shared_elements]
#print(not_in_set1)
#print(not_in_set2)
#print(len(common_elements))
#print(Headers_df)
#print(Headers_from_fasta)



"""
n=0
for i in count_of_n_df.iterrows():
    print(i[0],Headers_from_fasta[n])
    n+=1
"""
'''