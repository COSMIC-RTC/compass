# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:15:03 2015

@author: fvidal
"""
import pandas as pd
import os
import glob
os.environ["SHESHA_ROOT"]+"/data/simuDB/"
from subprocess import check_output
import time
#
#resAll = pd.DataFrame( columns=colnames.keys()) # res is the local dataframe for THIS data set
#resAll = resAll.append(colnames, ignore_index=True)  #Fill dataframe
#resAll = db.fillDf(resAll,h5u.params_dictionary(config))
def dfInfo(df):
    df.memory_usage(memory_usage='deep')


def createDf(colnamesList):
    if(type(colnamesList) is list):
        df = pd.DataFrame( columns=colnamesList) # res is the local dataframe for THIS data set
    else:
        df = pd.DataFrame( columns=colnamesList.keys()) # res is the local dataframe for THIS data set
    return df

def addcolumn(df, colnameList):
    nbcol = len(df.columns)
    for i in range(len(colnameList)):
        print "adding %s at col num:%d" % (colnameList[i], nbcol+i)
        df.insert(nbcol+i ,colnameList[i], None)
    return df


def fillDf(df, colnames):
    df = df.append(colnames, ignore_index=True)  #Fill dataframe
    return df


def mergeDB(df1, df2):
    print "TBDone"
    #
    dfmerged = df1.append(df2, ignore_index= True)
    return dfmerged
    #ind = resAll.duplicated("filename") find duplicates

def readDataBase(name='compassDB', dbFormat=".h5", fullpath=None):
    """
    database must be stored in $SHESHA_ROOT/data/simuDB/

    resAll = readDataBase() # reads compassDB.h5
    resAll = readDataBase("myDB") # .h5 format by default
    resAll = readDataBase(name="myDB", dbFormat=".csv")

    """

    if(name.find(".")<0):
        name = name+dbFormat
    if(fullpath is None):
        fullpath = os.environ["SHESHA_ROOT"]+"/data/simuDB/"+name

    if(not glob.glob(fullpath)):
        print "Cannot find database %s" % fullpath
        return 0
    else:
        print "reloading database from: %s" %  fullpath


    try:
        if(dbFormat == ".h5"):
            tmp = pd.read_hdf(fullpath,'resAll')
        elif(dbFormat == ".csv"):
            tmp = pd.read_csv(fullpath)
        else:
            print "Format ", dbFormat, " not recognized!"
            return
        return tmp
    except:
        print "Error Could not read database"
        return 0




def saveDataBase(df, name='compassDB', dbFormat=".h5"):
    #http://glowingpython.blogspot.fr/2014/08/quick-hdf5-with-pandas.html

    datapath = os.environ["SHESHA_ROOT"]+"/data/simuDB/"

    if(len(glob.glob(datapath+'/'+name+dbFormat))): # if the file "really exists" anyway
        currTime = time.strftime("%Y-%m-%d_%H:%M:%S", time.gmtime())
        CurrDbFile= datapath +'/'+name+dbFormat
        #check_output("mv "+ CurrDbFile  + " " +datapath +'/backup/'+name+dbFormat +".save_"+currTime+"_GMT")
        #print "Saved backup file: %s" % (CurrDbFile +".save_"+currTime+"_GMT")


    fullname = datapath +'/'+name+dbFormat
    if(dbFormat == ".h5"):
        hdf = pd.HDFStore(fullname)
        hdf.put('resAll', df, data_columns=True)
        hdf.close()
    elif(dbFormat == ".csv"):
        df.to_csv(fullname)
    else:
        print "ERROR format %s NOT recognized" % dbFormat




#resAll = readDataBase()
