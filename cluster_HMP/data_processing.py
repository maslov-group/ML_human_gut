########### Self-customized setting
import pandas as pd
import numpy as np

########### Load the metadata of HMP
metadata = pd.read_csv('./data/HMP_metadata.txt', delimiter='\t')
i_healthy_ids = metadata.loc[metadata['diagnosis']=='nonIBD', 'External ID']
metadata.head()

########### Load the mapped metabolites and microbes from HMP to Chia network
metagenome = pd.read_csv('./data/HMP_metagenome_matched.txt', delimiter='\t')
metabolome = pd.read_csv('./data/HMP_metabolome_matched_wb.txt', delimiter='\t')
metagenome = metagenome[['Microbe_name', 'Chia_name', 'Chia_id', 'Match'] + list(i_healthy_ids.values)]
metabolome = metabolome[['Compound_name', 'Kegg_id', 'Chia_name', 'Chia_id', 'Method'] + list(i_healthy_ids.values)]

########### Select the microbial strains have a maped existence in the Chia network
#metagenome_mapped_in_Chia = metagenome[metagenome['Chia_id'] > 0]
metagenome_mapped_in_Chia = metagenome[metagenome['Match'] == 'T']
metagenome = metagenome_mapped_in_Chia.groupby('Chia_id').apply(np.sum).iloc[:,4:]
metagenome_ID = metagenome.reset_index()['Chia_id']
metagenome.head()

print('The microbial coverage in Chia is',len(metagenome_mapped_in_Chia) / len(metagenome))
print(len(metagenome),'OTUs are shrinked to ', len(metagenome))

########### Select all metabolites in the metabolome have a maped existence in the Chia network
#metabolome_mapped_in_Chia = metabolome[metabolome['Chia_id'] > 0]
metabolome_mapped_in_Chia = metabolome[metabolome['Chia_id'] > 0]
metabolome = metabolome_mapped_in_Chia.groupby('Chia_id').apply(np.sum).iloc[:,5:]
metabolome = metabolome[metagenome.columns]
metabolome_ID = metabolome.reset_index()['Chia_id']
metabolome.head()

print('The metabolite coverage in Chia is',len(metabolome_mapped_in_Chia) / len(metabolome))
print(len(metabolome),'metabolites are shrinked to', len(metabolome))

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
#i_intake = pd.read_csv('nutrient_intake_ID_infact.txt',sep=': ')
i_intake = pd.read_csv('./data/nutrient_intake_ID.txt',sep=': ')
i_intake = i_intake['IDs'].values
print(i_intake)
print('###################################################################################################')

########### pickle all processed data which are useful for simulations
import pickle

pickle_out = open("Chia_network.pickle","wb")
pickle.dump([net, i_selfish, i_intake, names], pickle_out)
pickle_out.close()

pickle_out = open("data.pickle","wb")
pickle.dump([metagenome_ID, metagenome, metabolome_ID, metabolome], pickle_out)
pickle_out.close()