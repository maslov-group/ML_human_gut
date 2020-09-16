########### Self-customized setting
import pandas as pd
import numpy as np

########### Load Chia network as net (containing information of metabolite consumption and production)
net = pd.read_csv('./data/pruned_chia_network.csv')
mean_net = net.groupby('microbes_ID').mean()
selfish_net = mean_net[mean_net.iloc[:,1] == 2]
i_selfish = selfish_net.index.values   #### i_selfish returns IDs of microbes don't generate byproducts
print(net.head())
print('###################################################################################################')

########### Load names of all nodes in the Chia network
names = pd.read_csv('./data/names_ID.txt',sep=': ')
names.set_index('IDs', inplace=True)
print(names.head())
print('###################################################################################################')

########### Load names of all nodes in the Chia network
i_intake = pd.read_csv('./data/nutrient_intake_ID.txt',sep=': ')
i_intake = i_intake['IDs'].values
print(i_intake)
print('###################################################################################################')

########### Load all gut metagenomic data of all 41 individuals
metagenome_all = pd.read_csv('./data/abundance_matched_thai.txt', sep='\t')
metagenome_all.head()
metagenome_all = metagenome_all.groupby('Chia_id').sum().iloc[1:,].reset_index()
metagenome_ID = metagenome_all['Chia_id']
#print((metagenome_ID!=0).sum())
metagenome = metagenome_all[metagenome_ID!=0].iloc[:,3:]
metagenome_ID = metagenome_ID[metagenome_ID!=0]
print(metagenome.head())
print('###################################################################################################')

########### Load all gut metabolome data of all 41 individuals
metabolome_all = pd.read_excel('./data/metabolome_matched_thai_modified_by_Tong.xlsx')
metabolome_all = metabolome_all.groupby('Chia_id').sum().iloc[1:,].reset_index()
metabolome_ID = metabolome_all['Chia_id']
#print((metabolome_ID!=0).sum())
metabolome = metabolome_all[metabolome_ID!=0].iloc[:,2:]
metabolome_ID = metabolome_ID[metabolome_ID!=0]
print(metabolome.head())
print('###################################################################################################')

intersected_names = np.intersect1d(metagenome.columns.values, metabolome.columns.values)
metagenome = metagenome[intersected_names]
metabolome = metabolome[intersected_names]
print('Intersection between metagenome and metabolome:')
print(metagenome.head())
print(metabolome.head())

########### pickle all processed data which are useful for simulations
import pickle

pickle_out = open("Chia_network.pickle","wb")
#pickle.dump([net, i_selfish, i_intake, names], pickle_out)
pickle.dump([net, i_selfish, i_intake, names], pickle_out, protocol=2)
pickle_out.close()

pickle_out = open("data.pickle","wb")
#pickle.dump([metagenome_ID, metagenome, metabolome_ID, metabolome], pickle_out)
pickle.dump([metagenome_ID, metagenome, metabolome_ID, metabolome], pickle_out, protocol=2)
pickle_out.close()