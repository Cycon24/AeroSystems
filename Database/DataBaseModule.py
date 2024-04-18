# -*- coding: utf-8 -*-
"""
Created on Tue May 30 08:37:50 2023

@author: martinellibr
"""

import sqlite3
import numpy as np
from contextlib import closing

# _______________________________________________________________________________________________________
# ______________________________________ DATABASE INTERACTION FUNCTIONS _________________________________
def Connect(dbFileName):
    '''
    Parameters
    ----------
    dbFileName : String
        The file name for the database you wish to open/create. Ex: 'myDataBase.db'

    Returns
    -------
    Object
        A cursor object that is needed for any operations with the database.

    '''
    #____________ ESTABLISH CONNECTION ______________________
    dataBaseFN = dbFileName
    
    # File created automatically if not already existing
    # If want to create a RAM stored database that terminiates after
    #   python finishes, use fileName = ":memory:"
    con = sqlite3.connect(dataBaseFN)
    
    # Verify connection: -> Gives total number of rows that have been changes by connection
    if con.total_changes == 0:
        print('Connection Successful')
    else:
        print('Connection Unsuccessful')
    return con, con.cursor()

def CreateTable(cursorObj, tableName, colHeaders, colTypes):
    '''
    Parameters
    ----------
    cursorObj : Object
        Cursor object created from Connect() function.
    tableName : String
        The name of the table to create. Ex: 'MyTable2'
    colHeaders : Array of Strings
        Contains the names as strings for the table headers. Ex: ['Col1','name','numbers']
    colTypes : Array of Strings
        Contains the data types as strings for each column. Ex: ['Integer','Text','float']. Not case-sensitve.

    Returns
    -------
    Successful : Boolean
        True if creation was successful.

    '''
    
    Successful = False
    try:
        createCMD = 'CREATE TABLE {} ('.format(tableName)
        
        if len(colHeaders) == len(colTypes):
            for i, header in enumerate(colHeaders):
                if i == 0:
                    createCMD += '{} {}'.format(header,colTypes[i].upper())
                else:
                    createCMD += ', {} {}'.format(header,colTypes[i].upper())
            createCMD += ')'
          
        cursorObj.execute(createCMD)
        print(createCMD)
        Successful = True
    except:
        print('Table named "{}" already exists'.format(tableName))
    return Successful

def AddData(cursorObj, tableName, data2D, colTypes):
    '''
    Parameters
    ----------
    cursorObj : Object
        Cursor object created from Connect() function.
    tableName : String
        The name of the table to create. Ex: 'MyTable2'
    data2D : Matrix of Strings
        A 2d array of strings with the same number of columns as the table and any number of rows.
    colTypes : Array of Strings
        Contains the data types as strings for each column. Ex: ['Integer','Text','float']. Not case-sensitve.

    Returns
    -------
    bool
        True if addition of data was successful.

    '''
    # Data to add must be same num cols as table
    # Get header names
    headers = list(map(lambda x: x[0], cursorObj.description))
    hCols = len(headers)
    
    # Manage data to ensure it is the right dimensions
    data = np.array(data2D)
    if len(np.shape(data)) > 2:
        print('ERROR: Data to add has too many dimensions.')
        return False
    else:
        try:
            dRows, dCols = np.shape(data)
        except ValueError:
            dRows = 1
            dCols = np.shape(data)[0]
            
            # tempData = np.empty((dRows,dCols))
            # tempData[0,:] = data[:]
            # data = [data]
            print(data)
    
    if not hCols == dCols:
        print('ERROR: Data does not have the same number of columns as the table:\n\tTable: {}\n\tData: {}'.format(hCols,dCols))
        return False
    
    # Now we know the data matches the table size, add data
    # Add check for column types?
    for i in range(0,dRows): 
        addCMD = "INSERT INTO {} VALUES (".format(tableName)
        for j, val in enumerate(data[i,:]):
            
            # if neither of those then its a string
            if j == 0:
                # Change types to what they need to be for the table
                if colTypes[j].lower() == 'float' or colTypes[j].lower() == 'integer':
                    addCMD += "{}".format(val)
                else:
                    addCMD += "'{}'".format(val)
                
            else:
                # Change types to what they need to be for the table
                if colTypes[j].lower() == 'float' or colTypes[j].lower() == 'integer':
                    addCMD += ", {}".format(val)
                else:
                    addCMD += ", '{}'".format(val)
                
        addCMD += ")"
        # print(addCMD)
        cursorObj.execute(addCMD)
        
    
    print("All data added successfully")
    return True


def ReadData(cursorObj, tableName, dataToRetrieve=[], Filters=[]):
    '''
    Parameters
    ----------
    cursorObj : Object
        Cursor object created from Connect() function.
    tableName : String
        The name of the table to create. Ex: 'MyTable2'
    dataToRetrieve : String Array, optional
        The column headers of the data you would like to receive. The default is [].
    Filters : String Array, optional
        An array of strings or a string consisting of all the conditional arguments in SQL format to filter the code. Ex: "col1 > 5". The default is [].

    Returns
    -------
    List
        A list consisting of all the rows of the requested data. Each row is grouped by index, row 1 = List[0]. To access specific data: data1 = List[0][0].

    '''
    # dataToRetrieve will be list of the columns wanted
    # Likely set this up as a constant?
    readCMD = "SELECT "
    
    if dataToRetrieve == []:
        readCMD += "*" # Retrieve all
    else:
        for i, colH in enumerate(dataToRetrieve):
            if i==0:
                readCMD += "{}".format(colH)
            else:
                readCMD += ", {}".format(colH)
    
    readCMD += " FROM {}".format(tableName)
    
    if not Filters == []:
        readCMD += " WHERE "
        for filt in Filters:
            readCMD += filt
    
    return cursorObj.execute(readCMD).fetchall()

def DeleteData_PK(cursorObj, tableName, PK_col_hdr, PK_ids):
    # Will delete rows of data of particular tables in case a mistake was made.
    # To limit mistakes on data entry, ask the user to confirm all entries before committing the changes to database with connection.commit()
    # Deleting Data
    
    if type(PK_ids) == list or  isinstance(PK_ids, np.ndarray) :
        for pk in PK_ids:
            com = "DELETE FROM {} WHERE {} = {}".format(tableName, PK_col_hdr, pk)
            cursorObj.execute(com)
    else:
        com = "DELETE FROM {} WHERE {} = {}".format(tableName, PK_col_hdr, PK_ids)
        cursorObj.execute(com)
    # show update
    rows = cursorObj.execute("SELECT * FROM {}".format(tableName)).fetchall()
    # print(rows)
    
    return None

def DeleteTableData(cursorObj, tableName):
        com = "DELETE FROM {}".format(tableName)
        cursorObj.execute(com)
        return None
        
def Disconnect(dbFileName):
    with closing(sqlite3.connect(dbFileName)) as con:
        with closing(con.cursor()) as cursor:
            rows = cursor.execute("SELECT 1").fetchall()
            if rows == [(1,)]:
                print('Connection Closed Successfully')
            else:
                print('Connection Closed Unuccessfully')



def GetUniqueProperties(cursorObj, tableName, targetHeader):
    '''
    Parameters
    ----------
    cursorObj : Object
        DESCRIPTION.
    tableName : String
        DESCRIPTION.
    targetHeader : String
        DESCRIPTION.

    Returns
    -------
    UniqueProps : List
        All of the unique properties within the column targetHeader in tableName of database.

    '''
    UniqueProps = []
    allProps = ReadData(cursorObj, tableName, targetHeader)
    
    for i, prop in enumerate(allProps):
        if not(prop in UniqueProps) and not (prop==''):
            UniqueProps.append(prop)
            
    return UniqueProps




