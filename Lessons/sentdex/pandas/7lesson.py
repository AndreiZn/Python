# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 01:35:54 2017

@author: Andrey
"""

#lesson 6
import quandl
import pandas as pd
import pickle #save any python object (saves bite code)
# pandas has its own pickling functionality (its faster for large dataframes)

api_key = 'yjxUqzGvctxnxz16ErHx'

def state_list():    
    fifty_states = pd.read_html('https://simple.wikipedia.org/wiki/List_of_U.S._states')
    return fifty_states[0][0][1:]
 
def grab_initial_state_data():
    states = state_list()    
    main_df = pd.DataFrame()

    #this is a column:
    #print (fifty_states[0][0])
    
    for abbv in states:
        query = "FMAC/HPI_"+str(abbv)
        df = quandl.get(query, authtoken = api_key)
        df.rename (columns={'Value':abbv}, inplace = True)
        
        if main_df.empty:
            main_df = df
        else:
            main_df = main_df.join(df)
    
    
    print (main_df.head())
    
    pickle_out = open('fifty_states.pickle', 'wb') # write bites
    pickle.dump(main_df, pickle_out)
    pickle_out.close()
      
grab_initial_state_data()

pickle_in = open('fifty_states.pickle', 'rb')
HPI_data = pickle.load(pickle_in)
print(HPI_data)

# pandas version of pickling (only 2 lines):

HPI_data.to_pickle('Data_lesson6.pickle')
HPI_data2 = pd.read_pickle('Data_lesson6.pickle')
print(HPI_data2)