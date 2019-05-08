#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Entry Processor

    All methods needed to create, update and delete entries.

"""

import sqlite3 as sql
from sqlite3 import Error
#from sqlite3 import Binary
#
#import pandas as pd
#from pandas import DataFrame

import datetime
from os.path import basename

from modules.file_processor import file_processor

class entry_processor(object):

    def __init__(self, args):
#        print(args.database)
#        print(args.refGenome)
#        print(args.bedFile)
#        print(args.bedType)
#        print(args.sampleID)
#        print(args.normalID)
#        print(args.mutationFile)
#        print(args.mutationTool)
#        print(args.structVarFile)
#        print(args.structVarTool)
#        print(args.snp)
#        print(args.metrics)
#        print(args.depth)
#        print(args.log)
#        print(args.patientID)
#        print(args.sex)
        
        try:
            self.__args = args
            
            self.__tables = ['METRICS','DEPTH','SNPDIFF','STRUCTVAR','MUTATIONS','FILES','SAMPLE']
            
            self.__samp_id = args.sampleID

            self.__vars = {}
            self.__files = {}
            self.__invalid_files = {}
            
            self.__createConnection()

            self.__organizeFiles()     
            self.__organizeVars()
            
            if args.edit_mode == 'INSERT':
                self.__processInsert()

            elif args.edit_mode == 'UPDATE':
                self.__processUpdate()    
                
            elif args.edit_mode == 'DELETE':
                self.__processDelete()
            
            self.__conn.commit() 
            self.__closeConnection()
        except Exception as e:
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
        # organize non-file variables.
        
        if self.__args.normalID:
            var = 'normal_sample_id'
            self.__vars[var] = self.__args.normalID
            
        if self.__args.patientID:
            var = 'pat_id'
            self.__vars[var] = self.__args.patientID
            
        if self.__args.sex:
            var = 'pat_sex'            
            self.__vars[var] = self.__args.sex
            
    def __organizeFiles(self):
        # organize import files
        
        processor = file_processor(self.__args.current_dir)
        
        if self.__args.log:
            suffix = '_LOG'
            
            log_data = processor.processLog(self.__args.log)        
        
            if log_data is None:
                self.__invalid_files[suffix] = self.__args.log
            else:
                self.__vars['pipeline_log'] = log_data
                self.__files[suffix] = [suffix, self.__args.log, log_data, None, None, None]
    
        if self.__args.mutationFile:
            suffix = '_MUTATIONS'
            
            self.__files[suffix] = {}
            
            for mutFile, mutTool in zip(self.__args.mutationFile,self.__args.mutationTool):

                mut_data = processor.processVcf(suffix, mutFile)
                
                if mut_data is None:
                    if not self.__invalid_files.get('_MUTATIONS'):
                        self.__invalid_files[suffix] = {}
                        self.__invalid_files[suffix][mutTool] = mutFile
                else:
                    frame = mut_data[0]
                    header = mut_data[1]
                    frame.insert(0, 'mut_tool', mutTool, True)
                    self.__files[suffix][mutTool] = [suffix, mutFile, frame, header, mutTool, 'mut_tool']
                    #print(frame)
            
        if self.__args.structVarFile:
            suffix = '_STRUCTVAR'

            self.__files[suffix] = {}

            for svFile, svTool in zip(self.__args.structVarFile,self.__args.structVarTool):

                sv_data = processor.processVcf(suffix, svFile)
            
                if sv_data is None:
                    if not self.__invalid_files.get('_STRUCTVAR'):
                        self.__invalid_files[suffix] = {}                    
                        self.__invalid_files[suffix][svTool] = svFile
                else:
                    frame = sv_data[0]
                    header = sv_data[1]                    
                    frame.insert(0, 'sv_tool', svTool, True)
                    self.__files[suffix][svTool] = [suffix, svFile, frame, header, svTool, 'sv_tool']
            
        if self.__args.snp:
            suffix = '_SNPDIFF'

            snpdiff_data = processor.processTsv(suffix, self.__args.snp)
            
            if snpdiff_data is None:
                self.__invalid_files[suffix] = self.__args.snp
            else:
                self.__files[suffix] = [suffix, self.__args.snp, snpdiff_data, None, None, None]
                        
        if self.__args.metrics:
            suffix = '_METRICS'

            metrics_data = processor.processTsv(suffix, self.__args.metrics)
            
            if metrics_data is None:
                self.__invalid_files[suffix] = self.__args.metrics
            else:
                self.__files[suffix] = [suffix, self.__args.metrics, metrics_data, None, None, None]
            
        if self.__args.depth:
            suffix = '_DEPTH'

            depth_data = processor.processTsv(suffix, self.__args.depth)
            
            if depth_data is None:
                self.__invalid_files[suffix] = self.__args.depth
            else:
                self.__files[suffix] = [suffix, self.__args.depth, depth_data, None, None, None]
                
                
    #--------------------------------------
    # Process Methods
    #--------------------------------------        

    def __processInsert(self):
        # process entry insert
        
        if not self.__invalid_files:
                            
             self.__insertSample()
             self.__insertFiles()
        
        else:
            for bedType in self.__invalid_files:
                print('Invalid File Format: {}'.format(bedType))
                
            raise Exception('Database Error: Rollback Proceeding')        
        
    def __processUpdate(self):                
        # process entry update
        
        if not self.__invalid_files:
                       
            if self.__vars:
                self.__updateSample()
            
            if self.__files:
                self.__deleteFiles()
                self.__insertFiles()        
        else:
            for bedType in self.__invalid_files:
                print('Invalid File Format: {}'.format(bedType))
                
            raise Exception('Database Error: Rollback Proceeding')          
        
        
    def __processDelete(self):
        # process entry delete
        
        self.__deleteSample()

    #--------------------------------------
    # Sample Methods
    #--------------------------------------

    def __insertSample(self):
        # insert sample
        
        samp_sql = """INSERT INTO GENOM_SAMPLE 
                        (sample_id, normal_sample_id, pipeline_log, 
                         pat_id, pat_sex, date_processed) 
                      VALUES (?,?,?,?,?,?); """  
               
        samp_vars = [self.__samp_id, self.__vars.get('normal_sample_id'),
                     self.__vars.get('pipeline_log'), self.__vars.get('pat_id'),
                     self.__vars.get('pat_sex'), datetime.datetime.now()]

        try:
            cursor = self.__conn.cursor()
            cursor.execute(samp_sql,samp_vars)
            self.__samp_id = self.__args.sampleID 
            #print(self.__samp_id)
            cursor.close()
    
            print('GENOM_SAMPLE: Sample inserted. (samp_id = {})'.format(self.__samp_id))
        except Error as e:
            print('Error Inserting Sample: {}'.format(e))
            raise Exception('Database Error: Rollback Proceeding')     
            

    def __updateSample(self):
        # update sample
        
        samp_sql = """UPDATE GENOM_SAMPLE 
                         SET {} = ?
                       WHERE sample_id = ?; """  
        try:
            
            for var in self.__vars:
                sql = samp_sql.format(var)

                sql_vars = [self.__vars[var], self.__samp_id]

                cursor = self.__conn.cursor()
                cursor.execute(sql,sql_vars)

                cursor.close()
                print('GENOM_SAMPLE: {} updated.'.format(var))
        except Error as e:
            print('Error Updating Sample: {}'.format(e))
            raise Exception('Database Error: Rollback Proceeding') 
            

    def __deleteSample(self):
        # delete sample
        
        delete_sql = "DELETE FROM GENOM_{} WHERE sample_id = ?;"

        for table in self.__tables:
            delete_vars = [self.__args.sampleID]
            
            try:
                cursor = self.__conn.cursor()
                delete_count = cursor.execute(delete_sql.format(table),delete_vars)
                cursor.close()
                print('GENOM_{}: Deleted {} rows'.format(table, str(delete_count.rowcount)))        
            except Error as e:
                print('Error Deleting Entry: {}'.format(e))
                raise Exception('Database Error: Rollback Proceeding')  

    #--------------------------------------
    # File Methods
    #--------------------------------------  

    def __insertFiles(self):
        # insert files
        
        file_data = []
        
        for suffix in self.__files:
            
            if suffix == '_MUTATIONS' or suffix == '_STRUCTVAR':
                for tool in self.__files[suffix]:
                    file_data.append([suffix, self.__files[suffix][tool]])
            else:
                if suffix != '_LOG':
                    file_data.append([suffix, self.__files[suffix]])
                    
        for suffix, file_info in file_data:
            
            insert_sql = "INSERT INTO GENOM{} ({}) VALUES ({});"
            file_sql = "INSERT INTO GENOM_FILES (sample_id, file_type, file_tool, file_header, file_name, file_path) VALUES (?,?,?,?,?,?);"
            
            file_suffix = str(file_info[0])
            file_path = str(file_info[1])
            file_frame = file_info[2]
            file_header = str(file_info[3])
            file_tool = file_info[4]
    
            file_frame.insert(0, 'sample_id', self.__samp_id, True)
    
            headers = str(','.join(file_frame.columns))
            values = str(('?,' * len(file_frame.columns)).strip(','))
        
            try:    
                sql = insert_sql.format(file_suffix,headers,values)
                cursor = self.__conn.cursor()
                insert_count = cursor.executemany(sql,file_frame.values.tolist())
                cursor.close()
        
                file_vars = [self.__samp_id, file_suffix, file_tool, file_header, basename(file_path), file_path]
        
                cursor = self.__conn.cursor()
                cursor.execute(file_sql, file_vars)
                cursor.close()
                
                if file_tool:
                    print('GENOM{}: ({}) Inserted {} rows'.format(file_suffix, file_tool, str(insert_count.rowcount)))    
                else:
                    print('GENOM{}: Inserted {} rows'.format(file_suffix, str(insert_count.rowcount)))        
                
            except Error as e:
                print('Error Inserting Into {}: {}'.format('GENOM_{}'.format(file_suffix), e))
                raise Exception('Database Error: Rollback Proceeding')              

                
    def __deleteFiles(self):
        # delete files        
        
        file_data = []
        
        for suffix in self.__files:
            
            if suffix == '_MUTATIONS' or suffix == '_STRUCTVAR':
                for tool in self.__files[suffix]:
                    file_data.append([suffix, self.__files[suffix][tool], tool])
            else:
                if suffix != '_LOG':
                    file_data.append([suffix, self.__files[suffix], None])
                    
        for suffix, file_info, file_tool in file_data:
            del_sql = "DELETE FROM GENOM{} WHERE sample_id = ?"
            del_file_sql = "DELETE FROM GENOM_FILES WHERE sample_id = ? AND file_type = ?"

            tool_col = file_info[5]   
            
            del_vars = [self.__samp_id]
            del_file_vars = [self.__samp_id, suffix]
            
            if file_tool:
                del_sql += ' AND {} = ?;'
                del_sql = del_sql.format(suffix, tool_col)
                
                del_file_sql += ' AND file_tool = ?;'
                
                del_vars.append(file_tool)
                del_file_vars.append(file_tool)
            else:
                del_sql += ';'
                del_sql = del_sql.format(suffix)
                del_file_sql += ';'                
                 
            try:
                cursor = self.__conn.cursor()
                delete_count = cursor.execute(del_sql,del_vars)
                cursor.close()
                
                if file_tool:                
                    print('GENOM_{}: Deleted ({}) File: {} rows deleted'.format(suffix, file_tool, str(delete_count.rowcount)))
                else:
                    print('GENOM_{}: {} rows deleted'.format(suffix, str(delete_count.rowcount)))
            except Error as e:
                print('Error Deleting File: ({})'.format(e))
                raise Exception('Database Error: Rollback Proceeding')        
        
        
        
        
        
        
        
        
        
        
        