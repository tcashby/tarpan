#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Database Processor

    All methods needed to create and update the sqlite database.

"""

import os

import sqlite3 as sql
from sqlite3 import Error

from modules.file_processor import file_processor

class db_processor(object):

    def __init__(self, args):
#        print(args.database)
#        print(args.refGenome)
#        print(args.pipeline)        
#        print(args.targetBed)
#        print(args.groupBed)
#        print(args.bedFile)
#        print(args.bedType)        
#        print(args.current_dir)
#        print(args.edit_mode)        

        self.__args = args
        self.__db_id = None

        self.__vars = {}
        self.__files = {}
        self.__invalid_files = {}

        try:
            self.__createConnection()
            
            self.__organizeVars()
            self.__organizeFiles()            
            
            if args.edit_mode == 'INSERT':
                self.__processInsert()

            elif args.edit_mode == 'UPDATE':
                self.__processUpdate()                          

            self.__conn.commit() 
            self.__closeConnection()

        except Exception as e:

            if args.edit_mode == 'INSERT':
                self.__closeConnection()
                self.__deleteDB()
                
            elif args.edit_mode == 'UPDATE':
                self.__conn.rollback()
                self.__closeConnection()
                
            print('EXCEPTION: {}'.format(e))

    #--------------------------------------
    # Database connections
    #--------------------------------------

    def __createConnection(self):        
        # creates database if it doesn't exist
        
        try:        
            self.__conn = sql.connect(self.__args.database)
        except Error as e:
            print('Error Creating Connection: {}'.format(e)) 
            raise Exception('Database Error: Rollback Proceeding') 
            
    def __closeConnection(self):        
        # close database connecton
        
        try:        
            self.__conn.close()
        except Error as e:
            print('Error Closing Connection: {}'.format(e))   
            raise Exception('Database Error: Rollback Proceeding')  

    #--------------------------------------
    # Organize Data
    #--------------------------------------

    def __organizeVars(self):
        # organize vars in dict
        
        if self.__args.refGenome:
            var = 'ref_genome'
            self.__vars[var] = self.__args.refGenome
            
        if self.__args.pipeline:
            var = 'pipeline_name'
            self.__vars[var] = self.__args.pipeline
            
    
    def __organizeFiles(self):
        # confirm and add files to dict
        
        processor = file_processor(self.__args.current_dir)
        
        if self.__args.targetBed:

            bedType = 'targeted_regions'

            target_data = processor.processBed(bedType, self.__args.targetBed)
            
            if target_data is None:
                self.__invalid_files[bedType] = self.__args.targetBed
            else:
                self.__files[bedType] = target_data

        if self.__args.groupBed:

            bedType = 'groups'

            group_data = processor.processBed(bedType, self.__args.groupBed)
            
            if group_data is None:
                self.__invalid_files[bedType] = self.__args.groupBed
            else:
                self.__files[bedType] = group_data
        
        if self.__args.bedFile:
            
            for bedFile, bedType in zip(self.__args.bedFile,self.__args.bedType):
                
                bed_data = processor.processBed(bedType, bedFile)
                
                if bed_data is None:
                        self.__invalid_files[bedType] = bedFile
                else:
                    self.__files[bedType] = bed_data       
        
    #--------------------------------------
    # Process Methods
    #--------------------------------------        

    def __processInsert(self):
        # process initial creation of database.
        
        if not self.__invalid_files:
                            
            self.__createTables()        
            self.__insertDB()
            self.__insertFiles()
        
        else:
            for bedType in self.__invalid_files:
                print('Invalid File Format: {}'.format(bedType))
                
            raise Exception('Database Error: Rollback Proceeding')        
        
    def __processUpdate(self):
        # process update to database information
        
        if not self.__invalid_files:
        
            db_sql = """SELECT db_id FROM GENOM_DB LIMIT 1; """        
    
            try:  
                cursor = self.__conn.cursor()
                cursor.execute(db_sql)
                self.__db_id = cursor.fetchone()[0]
                cursor.close()
        
                print('GENOM_DB: Variables retrieved. (db_id = {})'.format(self.__db_id))
                
            except Error as e:
                print('Error Retrieving DB Info: {}'.format(e)) 
                raise Exception('Database Error: Rollback Proceeding')             
            
            if self.__vars:
                self.__updateDB()
                
            if self.__files:
                self.__deleteFiles()
                self.__insertFiles()

        else:
            for bedType in self.__invalid_files:
                print('Invalid File Format: {}'.format(bedType))
                
            raise Exception('Database Error: Rollback Proceeding')                

    #--------------------------------------
    # Database Methods
    #--------------------------------------

    def __createTables(self):       
        # create the initial tables
        
        sql = ''
        
        sql_file = self.__args.current_dir + '/db/generate_script.sql'
        
        with open(sql_file, 'r') as file:
            sql = file.read()
            #print(sql)
        
        try:           
            cursor = self.__conn.cursor()
            cursor.executescript(sql)
            cursor.close()
            
            print('DB: Tables Generated')
        except Error as e:
            print('Error Generating Tables: {}'.format(e)) 
            raise Exception('Database Error: Rollback Proceeding')         

            
    def __insertDB(self):
        # insert initial database information 
        
        db_sql = """INSERT INTO GENOM_DB (ref_genome,pipeline_name) VALUES (?,?); """
        db_vars = [self.__vars.get('ref_genome'), self.__vars.get('pipeline_name')]

        try:  
            cursor = self.__conn.cursor()
            cursor.execute(db_sql,db_vars)
            self.__db_id = cursor.lastrowid 
            cursor.close()
    
            print('GENOM_DB: Variables inserted. (db_id = {})'.format(self.__db_id))
        except Error as e:
            print('Error Inserting DB Info: {}'.format(e)) 
            raise Exception('Database Error: Rollback Proceeding')  
            
        
    def __updateDB(self):
        # update database information
        
        db_sql = """UPDATE GENOM_DB SET {} = ? WHERE db_id = ?; """    
        
        for var in self.__vars:
            sql = db_sql.format(var)
            sql_vars = [self.__vars[var], self.__db_id]
                
            try:  
                cursor = self.__conn.cursor()
                cursor.execute(sql,sql_vars)

                cursor.close()
        
                print('GENOM_DB: {} updated.'.format(var))
            except Error as e:
                print('Error Updating DB Info: {}'.format(e)) 
                raise Exception('Database Error: Rollback Proceeding')         


    def __deleteDB(self):
        # delete database
        if os.path.exists(self.__args.database):
            os.remove(self.__args.database)
        
    #--------------------------------------
    # File Methods
    #--------------------------------------        
        
    def __insertFiles(self):
        # insert bed files
        
        db_sql = """INSERT INTO GENOM_BEDFILES ({}) VALUES ({}); """

        for file_type in self.__files:

            file_data = self.__files[file_type]

            file_data.insert(0,'db_id',self.__db_id,True)
            file_data.insert(len(file_data.columns),'type',file_type,True)
            
            headers = str(','.join(file_data.columns))
            values = str(('?,' * len(file_data.columns)).strip(','))
            
            sql = db_sql.format(headers,values)
            
            try:              
                insert_count = 0                
                cursor = self.__conn.cursor()
                insert_count = cursor.executemany(sql,file_data.values.tolist()).rowcount
                cursor.close()
                print('GENOM_BEDFILES: ({}) File Processed. {} rows inserted'.format(file_type,insert_count))
            except Error as e:
                print('Error Inserting Bed Files: {}'.format(e)) 
                raise Exception('Database Error: Rollback Proceeding')
                   
    def __deleteFiles(self):
        # delete bed files
        
        db_sql = """DELETE FROM GENOM_BEDFILES WHERE db_id = ? AND type = ?; """
        
        for file_type in self.__files:
            
            db_vars = [self.__db_id, file_type]
            
            try:
                cursor = self.__conn.cursor()
                delete_count = cursor.execute(db_sql,db_vars)
                cursor.close()
                print('GENOM_BEDFILES: Deleted ({}) File: {} rows deleted'.format(file_type, str(delete_count.rowcount)))        
            except Error as e:
                print('Error Deleting Bed File: ({})'.format(e))
                raise Exception('Database Error: Rollback Proceeding')
