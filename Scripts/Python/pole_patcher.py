# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 15:08:52 2016

@author: Nicolò Lardelli @D-ERDW ETH Zürich
"""
import sys
import pandas as pd
from matplotlib import pyplot as pp
import h5py
import numpy as np

if __name__=="__main__":

    
    try:
        argv = sys.argv
        filename = argv[1]
    except RuntimeError as e:
        print(e)
        print('Supposed usage: python pole_patcher.py filename_tobepatched')
        sys.exit()

    h5_file = h5py.File(filename,'r+')

    try:
        theta_points = np.array(h5_file['mesh']['grid_theta'])
        theta_points = np.insert(theta_points,0,np.pi)
        theta_points = np.insert(theta_points,len(theta_points),0.)

        del h5_file['mesh']['grid_theta']
        h5_file.create_dataset('mesh/grid_theta',data = theta_points)

        phi_points = np.array(h5_file['mesh']['grid_phi'])
        phi_points = np.insert(phi_points,len(phi_points),2*np.pi)

        del  h5_file['mesh']['grid_phi']
        h5_file.create_dataset('mesh/grid_phi',data = phi_points)
        

        
        


    except RuntimeError as e:
        print(e)
        sys.exit()
            
    for group in list(h5_file):
        for subg in list(h5_file[group]):

            ref = h5_file[group][subg]
            
            dataset = np.array(ref)
            if len(dataset.shape)!=3:
                continue

            name = ref.name

            # patch boundaries in phi direction
            row = dataset[:,:,0]
            dataset = np.insert(dataset,dataset.shape[2],row,axis=2)

            # patch boundaries in theta direction bottom

            row = dataset[:,0,:]
            row = np.ones_like(row)*np.mean(row)
            dataset = np.insert(dataset,0,row,axis = 1)

            row = dataset[:,-1,:]
            row = np.ones_like(row)*np.mean(row)
            dataset = np.insert(dataset,dataset.shape[1],row,axis = 1)

            del h5_file[name]

            h5_file.create_dataset(name,data = dataset)
            
            
           
