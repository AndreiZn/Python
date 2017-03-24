#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 21:41:52 2017

@author: andreyzn
"""

import pandas as pd
import numpy as np

a = pd.read_table('./32Mpc/out_2.list.id', nrows = 17)

#read
halo = pd.read_table('./32Mpc/out_2.list.id', sep = ' ', skiprows=range(1,17))

#select
halo2 = halo[['X', 'Y', 'Z', 'M200c']]

#save
halo2.to_csv ('./32Mpc/halo_selected.csv')

#open
#halo3 = pd.read_csv('./64Mpc/halo_selected.csv', index_col = 0)