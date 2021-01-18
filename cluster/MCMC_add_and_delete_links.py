########### Self-customized setting
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import sklearn
import math
from numpy import array
from scipy.sparse import csr_matrix
import numpy.matlib
from scipy.optimize import minimize, rosen, rosen_der
from scipy.stats import pearsonr
import sys

########### import the pickled file containing all processed data which are useful for simulations (the processing is
########### done in "Trophic_model_for_gut_data_processing.ipynb")
import pickle
pickle_in = open("Chia_network.pickle","rb")
net, i_selfish, i_intake, names = pickle.load(pickle_in)

pickle_in = open("data.pickle","rb")
metagenome_ID, metagenome, metabolome_ID, metabolome = pickle.load(pickle_in)

#i_nonzero_microbes = net['microbes_ID'].unique()
#i_nonzero_microbes = np.sort(i_nonzero_microbes)
i_nonzero_microbes = metagenome_ID.values.copy()
i_nonzero_metabolites = net['metabolites_ID'].unique()
i_nonzero_metabolites = np.sort(i_nonzero_metabolites)


df_metabolites = pd.DataFrame.from_dict({'oldID': i_nonzero_metabolites, 'newID':list(range(len(i_nonzero_metabolites)))})
df_metabolites.set_index('oldID', inplace=True)
df_microbes = pd.DataFrame.from_dict({'oldID': i_nonzero_microbes, 'newID':list(range(len(i_nonzero_microbes)))})
df_microbes.set_index('oldID', inplace=True)


outgoingNodes = df_metabolites.reindex(net['metabolites_ID'].values).values.flatten()
ingoingNodesTemp = df_microbes.reindex(net['microbes_ID'].values).values.flatten()
edge_types = net.iloc[~np.isnan(ingoingNodesTemp),2].values
outgoingNodes = outgoingNodes[~np.isnan(ingoingNodesTemp)]
ingoingNodes = ingoingNodesTemp[~np.isnan(ingoingNodesTemp)].astype(int)

net_reduced = pd.DataFrame.from_dict({'metabolites': outgoingNodes, 'microbes':ingoingNodes, 'edgeType':edge_types})
net = net_reduced.copy()
net_temp = net.copy()
net['edgeType'][net['edgeType']==5] = 2
net_temp['edgeType'][net_temp['edgeType']==5] = 3
net = net.append(net_temp).drop_duplicates()
net_ori = net.copy()

metagenome_ID_reduced = df_microbes.reindex(metagenome_ID).values.flatten()
metagenome_ID = metagenome_ID_reduced[~np.isnan(metagenome_ID_reduced)].astype(int)

metabolome_ID_reduced = df_metabolites.reindex(metabolome_ID).values.flatten()
metabolome_ID = metabolome_ID_reduced[~np.isnan(metabolome_ID_reduced)].astype(int)


i_selfish_reduced = df_microbes.reindex(i_selfish).values.flatten()
i_selfish = i_selfish_reduced[~np.isnan(i_selfish_reduced)].astype(int)

i_intake_reduced = df_metabolites.loc[i_intake].values.flatten()
i_intake = i_intake_reduced[~np.isnan(i_intake_reduced)].astype(int)


MAX_ID_microbes = len(i_nonzero_microbes)  # MAX_ID_microbes is the maximum of ID labels for microbes.
MAX_ID_metabolites = len(i_nonzero_metabolites)  # MAX_ID_metabolites is the maximum of ID labels for metabolites.

def Ain_out(b_real, net):
    '''
    Ain_out is a function used to create sparse matrices made of metabolites and microbes 
    where metabolite consumption and production is considered. The matrices created are "m2b" and "b2m":
    (1) m2b is a matrix determines the nutrient splitting among microbes, and
    (2) b2m is a matrix determines the byproducts generation.
    Both matrices have rows representing bacterial species and columns representing metablolites.
    Two matrices are created based on (1) the metabolite consumption and production network which is 
    encode in "net" as a dataframe, and (2) the experimentally measured bacterial abundance "b_real".
    '''
    valid_index = np.where((net['edgeType']==2) | (net['edgeType']==5))[0]
    row = net['metabolites'].iloc[valid_index]
    col = net['microbes'].iloc[valid_index]
    data = np.ones((len(valid_index),))
    m2b = csr_matrix( (data,(row,col)), shape=(MAX_ID_metabolites, MAX_ID_microbes)).toarray()#.todense()

    valid_index = np.where((net['edgeType']==3) | (net['edgeType']==5))[0]
    row = net['metabolites'].iloc[valid_index]
    col = net['microbes'].iloc[valid_index]
    data = np.ones((len(valid_index),))
    b2m = csr_matrix( (data,(row,col)), shape=(MAX_ID_metabolites, MAX_ID_microbes)).toarray()#.todense()

    ########## Normalize the b2m by out_degree
    out_degree = b2m.sum(0).copy()
    out_degree[out_degree==0]=100
    b2m = (b2m / out_degree)

    ########## Normalize the m2b by proportion of microbial abundance in each individual
    b_real_repmat = numpy.matlib.repmat(b_real[np.newaxis,:], MAX_ID_metabolites, 1)
    m2b = m2b * b_real_repmat
    in_degree = m2b.sum(1)
    in_degree[in_degree==0]=100
    m2b = m2b / numpy.matlib.repmat(in_degree[:,np.newaxis], 1, MAX_ID_microbes)
    m2b = np.float32(m2b)
    b2m = np.float32(b2m)
    return [m2b, b2m]

def m2b_multiple_levels(f, m2b, b2m, numLevels_max):
    '''
    m2b_multiple_levels is a function used to generate matrices involving the calculation of metabolite 
    byproducts and microbial biomass after several trophic levels/layers. Those matrices are:
    (1) m2m_layer is a conversion matrix from the nutrient intake to the metabolite byproducts at a trophic
    level or layer.
    (2) m2m_total is a conversion matrix from the nutrient intake to a summation of metabolite byproducts at
    all trophic levels or layers.
    (3) m2b_total is a conversion matrix from the nutrient intake to a summation of all microbial/bacterial 
    biomass gain at all trophic levels or layers.
    Those matrices are computed based on (1) metabolite consumption matrix "m2b", (2) metabolite byproduct
    generation matrix "b2m", (3) byproduct/leakage fraction "f", and (4) number of trophic levels/layers in the 
    simulation "numLevels_max".
    '''
    m2m_layer = np.zeros((MAX_ID_metabolites, MAX_ID_metabolites, numLevels_max));  
    #m2b_total = np.zeros((MAX_ID_metabolites, MAX_ID_microbes));  
    m2b_total = np.zeros((MAX_ID_metabolites, MAX_ID_metabolites));  
    
    f_mul = numpy.matlib.repmat(f, 1, MAX_ID_metabolites)
    #s_step =  np.dot(b2m, m2b.T) # s_step is the conversion matrix of each trophic level/layer
    s_step =  np.dot(b2m, f_mul*m2b.T) # s_step is the conversion matrix of each trophic level/layer
    s_step_ii = np.eye(MAX_ID_metabolites, MAX_ID_metabolites);
    #f_mul = numpy.matlib.repmat(f[np.newaxis,:], MAX_ID, 1)#numpy.matlib.repmat(f, 1, MAX_ID)
    #f_mul = numpy.matlib.repmat(f, 1, MAX_ID_metabolites)
    
    for ii in range(numLevels_max):
        # m2b_total is a series made of s_step of each layer
        #m2b_total = m2b_total + f_mul**(ii) * s_step_ii
        #m2m_layer[:,:,ii] = f_mul**(ii) * s_step_ii
        m2b_total = m2b_total + s_step_ii
        m2m_layer[:,:,ii] = s_step_ii
        s_step_ii = np.dot(s_step_ii, s_step);  
    m2m_total = m2b_total;
    m2b_total = np.dot((1 - f_mul) * m2b.T, m2b_total) # m2b_total has an extra multiplication of m2b and (1-f).
    return [m2b_total, m2m_total, m2m_layer]

def pred_error(x, b_real, m2b_total, i_intake):
    '''
    pred_error is a function used to compute the logarithmic error between experimentally measured
    metagenome and predicted metagenome computed from the model for a certain nutrient intake. It relies on 
    (1) x: the nutrient intake, (2) i_intake: IDs of the nutrient intake, (3) m2b_total: a conversion matrix 
    from thenutrient intake to the total biomass, and (4) b_real: experimentally measured metagenome. The 
    first three is used to compute the metagenome predicted by the model "ba_pred" and compare it with the 
    experimentally measured metagenome "b_real".
    '''
    x_full = np.zeros((MAX_ID_metabolites,));
    x_full[i_intake] = x;
    
    ba_pred = np.dot(m2b_total, x_full)
    pred_error = (np.log10(ba_pred + 1e-6) - np.log10(b_real +1e-6)) / np.log10(b_real +1e-6);
    pred_error = np.sqrt(np.dot(pred_error, pred_error.T)) #np.sqrt(np.sum(pred_error**2))
    return pred_error

def calc_metabolome(x, m2m_layer, numLevels_max, i_intake):
    '''
    calc_metabolome is a function used to calculate the metabolome from the fitted nutrient intake from the
    model. It relies on (1) x: the nutrient intake, (2) i_intake: IDs of the nutrient intake, (3) m2m_layer: 
    a conversion matrix from the nutrient intake to the metabolite byproducts at a trophic level or layer, 
    and (4) numLevels_max: the number of trophic levels/layers in the model. The metabolome in the model is 
    assumed to be composed of two parts: (1) met_levels: all metabolites in the final trophic level/layer 
    (which is considered to be reaching the end of the gut because of the finite gut length and gut motility.),
    and (2) met_leftover_levels: all unusable metabolites from all previous trophic levels/layers. 
    '''
    i_unused = np.where(np.sum(m2b.T,0) == 0)[0]
    met_levels = np.zeros((MAX_ID_metabolites, numLevels_max));  
    met_leftover_levels = np.zeros((MAX_ID_metabolites, numLevels_max));
    
    x_full = np.zeros((MAX_ID_metabolites,));
    x_full[i_intake] = x;
    
    for ii in range(numLevels_max):
        met_levels[:,ii] = np.dot(m2m_layer[:,:,ii], x_full)
        if ii==0:
            met_leftover_levels[i_unused,ii] = x_full[i_unused]
        else:
            met_leftover_levels[i_unused,ii] = met_levels[i_unused,ii-1]
            
    return [met_levels, met_leftover_levels]

import random
######## Select all microbes in the Thai data:
b_real = np.zeros((MAX_ID_microbes,))
b_real[metagenome_ID] = 0.05

######## Compute matrices involving the metabolite consumption and generation:
m2b, b2m = Ain_out(b_real, net)

m2b_ori = (m2b!=0).astype(int).copy()
b2m_ori = (b2m!=0).astype(int).copy()
x_ori = np.concatenate([m2b_ori.flatten(), b2m_ori.flatten()])
net_ori = net.copy()

def pred_error_addingLinks(i_train, x, m2b_ori, b2m_ori, x_ori):
    '''
    pred_error_addingLinks is a function used to compute the logarithmic error between experimentally measured
    metagenome and predicted metagenome computed from the model for a certain nutrient intake. It relies on 
    (1) x: the nutrient intake, (2) i_intake: IDs of the nutrient intake, (3) m2b_total: a conversion matrix 
    from thenutrient intake to the total biomass, and (4) b_real: experimentally measured metagenome. The 
    first three is used to compute the metagenome predicted by the model "ba_pred" and compare it with the 
    experimentally measured metagenome "b_real".
    '''
    max_links = m2b_ori.shape[0] * m2b_ori.shape[1] # maximal number of links = numberr of specis * number of metabolites
    
    ######## Convert x to net structure (convert the adjacency matrix into the edge list):
    thres = 0.1
    x = x - x_ori
    ### delete links:
    net_ori_deleted = net_ori.copy()
    m2b_deleted = x[:max_links].reshape((m2b_ori.shape[0], m2b_ori.shape[1]))
    a = df_metabolites['newID'].values[np.where(m2b_deleted <= -thres)[0]]
    b = metagenome_ID[np.where(m2b_deleted <= -thres)[1]]
    c = [2] * np.where(m2b_deleted <= -thres)[1].shape[0]
    net_deleted = pd.DataFrame({net_ori.columns[0]:list(a), net_ori.columns[1]:list(b), net_ori.columns[2]:c})
    net_ori_deleted = pd.concat([net_ori_deleted, net_deleted]).drop_duplicates(keep=False)
    ### consumption links:
    m2b_added = x[:max_links].reshape((m2b_ori.shape[0], m2b_ori.shape[1]))
    a = df_metabolites['newID'].values[np.where(m2b_added >= thres)[0]]
    b = metagenome_ID[np.where(m2b_added >= thres)[1]]
    c = [2] * np.where(m2b_added >= thres)[1].shape[0]
    net_added_consumption = pd.DataFrame({net_ori.columns[0]:list(a), net_ori.columns[1]:list(b), net_ori.columns[2]:c})
    ### production links:
    b2m_added = x[max_links:].reshape((m2b_ori.shape[0], m2b_ori.shape[1]))
    a = df_metabolites['newID'].values[np.where(b2m_added >= thres)[0]]
    b = metagenome_ID[np.where(b2m_added >= thres)[1]]
    c = [3] * np.where(b2m_added >= thres)[1].shape[0]
    net_added_production = pd.DataFrame({net_ori.columns[0]:list(a), net_ori.columns[1]:list(b), net_ori.columns[2]:c})
    ### overall edge list of the new network
    net = pd.concat([net_ori_deleted, net_added_consumption, net_added_production])
    
    f_byproduct = 0.9
    f = f_byproduct * np.ones((MAX_ID_microbes,1));
    f[i_selfish] = 0.0;  # The byproduct/leakage fraction f for microbes that don't generate byproducts is set as 0.

    numLevels_max = 4

    corr_list_aveDiet = np.zeros((i_train.shape[0], 2))
    log_list_aveDiet = np.zeros((i_train.shape[0], 2))
    order_dev_list = np.zeros((i_train.shape[0], MAX_ID_metabolites)) ## order difference in metabolome
    order_dev_metagenome_list = np.zeros((i_train.shape[0], MAX_ID_microbes)) ## order difference in metagenome
    numMetabolites_list = np.zeros((i_train.shape[0]))
    for pa in range(i_train.shape[0]):
        ######## Select an ID of a individual and pull out its experimentally measured metagenome:
        b_real = np.zeros((MAX_ID_microbes,))
        b_real[metagenome_ID] = metagenome.iloc[:, i_train[pa]] / np.sum(metagenome.iloc[:, i_train[pa]])

        ######## Compute matrices involving the metabolite consumption and generation:
        m2b, b2m = Ain_out(b_real, net)
        #m2b_total, m2m_total, m2m_layer = m2b_multiple_levels(f, m2b, b2m, numLevels_max)
        m2b_total, m2m_total, m2m_layer = m2b_multiple_levels(f, m2b, b2m, numLevels_max)

        ######## As long as the optimized nutrient intake is found by using the above optimization solver, the
        # optimized nutrient intake "res.x" is used to generate the predict metagenome and metabolome. They are
        # visually and statistically compared to the experimentally measured metagenome and metabolome.
        
        ######## Assign diet
        x0 = np.ones((len(i_intake),)) / len(i_intake)

        ######## The model is converted into an optimization problem where the nutrient intake is constantly changed to
        # minimize the logarithmic error between experimentally measured metagenome and predicted metagenome computed 
        # from the model for a certain nutrient intake.
        fun2 = lambda y: pred_error(y, b_real, m2b_total, i_intake)
        bnds = ((0, 100), ) * len(x0)
        res = minimize(fun2, x0, method='SLSQP', bounds=bnds, options={'disp': False, 'maxiter': 1000}, tol=1e-3, callback=None)

        ######## As long as the optimized nutrient intake is found by using the above optimization solver, the
        # optimized nutrient intake "res.x" is used to generate the predict metagenome and metabolome. They are
        # visually and statistically compared to the experimentally measured metagenome and metabolome.
        x_full = np.zeros((MAX_ID_metabolites,));
        x_full[i_intake] = res.x;
        ba_pred = np.dot(m2b_total, x_full)

        met_levels, met_leftover_levels = calc_metabolome(res.x, m2m_layer, numLevels_max, i_intake)
        metabolome_measured = np.zeros((MAX_ID_metabolites,))
        metabolome_measured[metabolome_ID] = metabolome.iloc[:, i_train[pa]]
        metabolome_pred = met_levels[:, -1] + met_leftover_levels.sum(1) #np.dot(m2m_total, x_full)
        i_common = np.where(metabolome_measured * metabolome_pred > 1e-5)[0]
        metabolome_pred_common = metabolome_pred[i_common] / np.sum(metabolome_pred[i_common])
        metabolome_measured_common = metabolome_measured[i_common] / np.sum(metabolome_measured[i_common])
        ######## compute the bias of the order of magnitude
        order_dev_metagenome_list[pa, :] = np.log10(b_real + 1e-5) - np.log10(ba_pred + 1e-5)
        order_dev_list[pa, i_common] = np.log10(metabolome_measured_common + 1e-5) - np.log10(metabolome_pred_common + 1e-5)
        numMetabolites_list[pa] = i_common.shape[0]

        corr_list_aveDiet[pa, 0] = pearsonr(ba_pred[b_real>0], b_real[b_real>0])[0]
        corr_list_aveDiet[pa, 1] = pearsonr(metabolome_pred_common, metabolome_measured_common)[0]
            
        if i_common.shape[0] >= 2:
            corr_list_aveDiet[pa, 0] = pearsonr(ba_pred[b_real>0], b_real[b_real>0])[0]
            corr_list_aveDiet[pa, 1] = pearsonr(metabolome_pred_common, metabolome_measured_common)[0]
        else:
            corr_list_aveDiet[pa, 0] = -1
            corr_list_aveDiet[pa, 1] = -1
        
        log_list_aveDiet[pa, 0] = np.mean(np.abs(np.log10(ba_pred[b_real>0]+1e-6) - np.log10(b_real[b_real>0]+1e-6)))
        log_list_aveDiet[pa, 1] = np.mean(np.abs(np.log10(metabolome_pred_common+1e-6) - np.log10(metabolome_measured_common+1e-6)))
        
    x_full = np.zeros((MAX_ID_metabolites,));
    x_full[i_intake] = res.x;
    
    ba_pred = np.dot(m2b_total, x_full)
    
    pred_error1 = np.mean(log_list_aveDiet[:,0])
    pred_error2 = np.mean(log_list_aveDiet[:,1])
    pred_error3 = np.sum(np.abs(x))
    pred_error4 = np.mean(numMetabolites_list)
    hyper_reg = 0.001
    pred_errorTotal = pred_error2 + hyper_reg * pred_error3 - (pred_error4 - 20) * 0.003 # with reward
    
    return [pred_errorTotal, np.abs(np.sum(order_dev_list, 0)), np.sum(order_dev_metagenome_list, 0)]

i_train = np.array(range(metabolome.shape[1])) #np.sort(np.random.choice(range(metabolome.shape[1]), 30, replace=False))
i_test = np.setdiff1d(range(metabolome.shape[1]), i_train)

error_list = []
current_step_list = []
pos_x_list = []
metID_list = []
microbeID_list = []
prior_list = []
fun = lambda x: pred_error_addingLinks(i_train, x, m2b_ori, b2m_ori, x_ori)
max_links = m2b.shape[0]*m2b.shape[1]
x = x_ori.copy()
error_before, bias_metabolome, bias_metagenome = fun(x)
bias_metabolome_ori = bias_metabolome.copy()
prior_prob1 = np.matlib.repmat(np.exp(3 * (bias_metagenome) / metabolome.shape[1])[np.newaxis, :], m2b.shape[0], 1) * np.matlib.repmat(np.exp(3 * bias_metabolome_ori / metabolome.shape[1])[:, np.newaxis], 1, m2b.shape[1])
prior_prob2 = np.matlib.repmat(np.exp(3 * bias_metabolome_ori / metabolome.shape[1])[:, np.newaxis], 1, m2b.shape[1])
prior_prob = np.concatenate([prior_prob1.flatten(), prior_prob2.flatten()]) + 0.1
prior_prob = prior_prob / np.sum(prior_prob)
print('The original error is', error_before)
error_list.append(error_before)

inverseKT = 5000
Twindow = 500
numStepsNotAdded = 0
error_window = []
# MCMC
for i in range(10000):
    if i%50==0:
        print(i)
    if np.random.uniform(0,1,1)[0] < 0.5:
        i_x = np.random.choice(np.where(x==0)[0], 1, p=prior_prob[x==0]/np.sum(prior_prob[x==0]))[0]
    else:
        p_delete = 1.0 / prior_prob[x==1]
        p_delete = p_delete / np.sum(p_delete)
        i_x = np.random.choice(np.where(x==1)[0], 1, p=p_delete)[0]
        #i_x = np.random.choice(np.where(x==1)[0], 1)[0]
    x[i_x] = 1 - x[i_x]
    error_after, bias_metabolome, bias_metagenome = fun(x)
    ## Consumption links:
    prior_prob1 = np.matlib.repmat(np.exp(3 * (bias_metagenome) / metabolome.shape[1])[np.newaxis, :], m2b.shape[0], 1) * np.matlib.repmat(np.exp(3 * bias_metabolome_ori / metabolome.shape[1])[:, np.newaxis], 1, m2b.shape[1])
    ## Production links:
    prior_prob2 = np.matlib.repmat(np.exp(3 * bias_metabolome_ori / metabolome.shape[1])[:, np.newaxis], 1, m2b.shape[1])
    prior_prob = np.concatenate([prior_prob1.flatten(), prior_prob2.flatten()]) + 0.1
    prior_prob = prior_prob / np.sum(prior_prob)
    if np.random.uniform(0,1,1)[0] <= np.min([1, np.exp((inverseKT)*(error_before - error_after))]): 
        error_before = error_after
        if x[i_x] == 1:
            print('Addition accepted, error is', error_before)
        else:
            print('Deletion accepted, error is', error_before)
        error_list.append(error_before)
        current_step_list.append(i)
        pos_x_list.append(i_x)
        prior_list.append(prior_prob[i_x])
        if i_x < max_links:
            row_num = i_x // m2b_ori.shape[1]
            col_num = i_x - row_num * m2b_ori.shape[1]
        elif i_x >= max_links:
            i_x = i_x - m2b_ori.shape[0] * m2b_ori.shape[1]
            row_num = i_x // b2m_ori.shape[1]
            col_num = i_x - row_num * b2m_ori.shape[1]
        metID_list.append(row_num)
        microbeID_list.append(col_num)
        numStepsNotAdded = 0
    else:  ## not accepted   
        x[i_x] = 1 - x[i_x]
        numStepsNotAdded += 1
    error_window.append(error_before)
    if (i > Twindow) and ((error_window[-1] - error_window[-Twindow]) > -(np.sqrt(Twindow) / inverseKT)):
        break

######## Convert x to net structure:
a = np.array(metID_list)
b = np.array(microbeID_list)
c = np.ones([len(metID_list)], dtype = int) * 3
c[np.where(np.array(pos_x_list) < max_links)[0]] = 2
net_modified = pd.DataFrame({net_ori.columns[0]:list(a), net_ori.columns[1]:list(b), net_ori.columns[2]:c})
i_added = x[pos_x_list].astype(bool)
i_deleted = ~x[pos_x_list].astype(bool)
net_added = net_modified[i_added]
net_deleted = net_modified[i_deleted]
net_new = pd.concat([net_ori, net_added, net_deleted]).drop_duplicates(keep=False)

NUMADDED = len(net_added)

df_added_metabolites = names.loc[df_metabolites.reset_index().set_index('newID').loc[net_added.iloc[-NUMADDED:,0]].values.flatten()]
df_added_metabolites.reset_index(inplace=True)
df_added_metabolites.columns = ['metabolite ID', 'metabolite names'] 

df_added_microbes = names.loc[df_microbes.reset_index().set_index('newID').loc[net_added.iloc[-NUMADDED:,1]].values.flatten()]
df_added_microbes.reset_index(inplace=True)
df_added_microbes.columns = ['microbes ID', 'microbes names'] 

df_added_edgeTypes = net_added.iloc[-NUMADDED:,2]
df_added_edgeTypes = pd.DataFrame(df_added_edgeTypes)
df_added_edgeTypes.index = range(NUMADDED)

df_added_errorReduced = -np.diff(error_list)[i_added]
df_added_errorReduced = pd.DataFrame(df_added_errorReduced)
df_added_errorReduced.index = range(NUMADDED)
df_added_errorReduced.columns = ['error reduced'] 

df_added_step = np.array(current_step_list)[i_added]
df_added_step = pd.DataFrame(df_added_step)
df_added_step.index = range(NUMADDED)
df_added_step.columns = ['step added'] 

df_added_pos = np.array(pos_x_list)[i_added]
df_added_pos = pd.DataFrame(df_added_pos)
df_added_pos.index = range(NUMADDED)
df_added_pos.columns = ['position in x'] 

df_added_tables = pd.concat([df_added_metabolites, df_added_microbes, df_added_edgeTypes,  df_added_errorReduced, df_added_step, df_added_pos], axis=1, sort=False)
df_added_tables.sort_values(by = 'error reduced', ascending=False)
df_added_tables.to_csv('added_links_prior_distribution_run' + sys.argv[1]+'.csv')

NUMdeleted = len(net_deleted)

df_deleted_metabolites = names.loc[df_metabolites.reset_index().set_index('newID').loc[net_deleted.iloc[-NUMdeleted:,0]].values.flatten()]
df_deleted_metabolites.reset_index(inplace=True)
df_deleted_metabolites.columns = ['metabolite ID', 'metabolite names'] 

df_deleted_microbes = names.loc[df_microbes.reset_index().set_index('newID').loc[net_deleted.iloc[-NUMdeleted:,1]].values.flatten()]
df_deleted_microbes.reset_index(inplace=True)
df_deleted_microbes.columns = ['microbes ID', 'microbes names'] 

df_deleted_edgeTypes = net_deleted.iloc[-NUMdeleted:,2]
df_deleted_edgeTypes = pd.DataFrame(df_deleted_edgeTypes)
df_deleted_edgeTypes.index = range(NUMdeleted)

df_deleted_errorReduced = -np.diff(error_list)[i_deleted]
df_deleted_errorReduced = pd.DataFrame(df_deleted_errorReduced)
df_deleted_errorReduced.index = range(NUMdeleted)
df_deleted_errorReduced.columns = ['error reduced'] 

df_deleted_step = np.array(current_step_list)[i_deleted]
df_deleted_step = pd.DataFrame(df_deleted_step)
df_deleted_step.index = range(NUMdeleted)
df_deleted_step.columns = ['step deleted'] 

df_deleted_pos = np.array(pos_x_list)[i_deleted]
df_deleted_pos = pd.DataFrame(df_deleted_pos)
df_deleted_pos.index = range(NUMdeleted)
df_deleted_pos.columns = ['position in x'] 

df_deleted_tables = pd.concat([df_deleted_metabolites, df_deleted_microbes, df_deleted_edgeTypes,  df_deleted_errorReduced, df_deleted_step, df_deleted_pos], axis=1, sort=False)
df_deleted_tables.sort_values(by = 'error reduced', ascending=False)
df_deleted_tables.to_csv('deleted_links_prior_distribution_run' + sys.argv[1]+'.csv')
