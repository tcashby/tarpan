#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Argument parser
    
    Parse arguments
    Contains methods to validate for intended purpose (create, update, delete)

"""

# DB Create
# create_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<BED_FILE>" -groupbed "<BED_FILE>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"

# DB Update (only one change arg required)
# update_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<BED_FILE>" -groupbed "<BED_FILE2>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"

# Entry Create
# create_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>" -mutfile "<MUTATION file>" -muttool "<MUTATION tool>" -mutfile "<MUTATION2 file>" -muttool "<MUTATION2 tool>" -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>" -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>" -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>" 

# Entry Update (only one arg change required)
# update_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>" -mutfile "<MUTATION file>" -muttool "<MUTATION tool>" -mutfile "<MUTATION2 file>" -muttool "<MUTATION2 tool>" -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>" -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>" -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>" 

# Entry Delete
# delete_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" 

import argparse
from os.path import isfile, getsize

import sqlite3 as sql
from sqlite3 import Error

class arg_parser(object):

    def __init__(self, argv):
        self.__parser = argparse.ArgumentParser(description=(""))

        self.parseArgs()

        self.__curr_dir = ''
        self.__edit_mode = ''        


    def parseArgs(self):        

        #------------------------------------

        self.bicond_args = [['bedfile','bedtype'],['mutfile','muttool'],['svfile','svtool']]
        self.file_args = ['bedfile','targetbed','groupbed','mutfile','svfile','snp','met','depth','log']
        
        #------------------------------------        
        self.__parser.add_argument("-db", "--db")
        
        self.__parser.add_argument("-refgen", "--refgen")
        self.__parser.add_argument("-pipeline", "--pipeline")
        self.__parser.add_argument("-targetbed", "--targetbed")
        self.__parser.add_argument("-groupbed", "--groupbed")
        self.__parser.add_argument("-bedfile", "--bedfile", action='append')
        self.__parser.add_argument("-bedtype", "--bedtype", action='append')
        
        self.__parser.add_argument("-sampid", "--sampid")
        self.__parser.add_argument("-normid", "--normid")
        
        self.__parser.add_argument("-mutfile", "--mutfile", action='append')
        self.__parser.add_argument("-muttool", "--muttool", action='append')
        self.__parser.add_argument("-svfile", "--svfile", action='append')
        self.__parser.add_argument("-svtool", "--svtool", action='append')
        
        self.__parser.add_argument("-snp", "--snp")
        self.__parser.add_argument("-met", "--met")
        self.__parser.add_argument("-depth", "--depth")
        self.__parser.add_argument("-log", "--log")
        
        self.__parser.add_argument("-patid", "--patid")
        self.__parser.add_argument("-sex", "--sex")
        
        self._args = vars(self.__parser.parse_args())        
        

    #---------------------------------------
    # Arg Check Functions
    #---------------------------------------    
    
    def checkDbCreateArgs(self):
        # check arguments for database creation
        
        required_args = ['db','refgen','pipeline','targetbed','groupbed']

        # check required arguments (false = fail)
        if not self.__checkArgsExist(required_args):
            raise Exception('Required Args Not Provided. DB Creation Aborted')            
        
        # check database exists (true = fail)  
        if self.__checkDatabaseExists():
            raise Exception('Database Already Exists. DB Creation Aborted')
        
        # check files exist (false = fail)  
        if not self.__checkFilesExist():
            raise Exception('File Errors. DB Creation Aborted')
                
        return True
        
    
    def checkDbUpdateArgs(self):
        # check arguments for database update
        
        required_args = ['db']
        possible_args = ['refgen','pipeline','targetbed','groupbed','bedfile']
        
        # check required arguments (false = fail)
        if not self.__checkArgsExist(required_args):
            raise Exception('Required Args Not Provided. DB Update Aborted')  
            
        # check at least one argument added
        if not self.__checkPossibleArgExists(possible_args):
            raise Exception('No Args Added to Update. DB Update Aborted')
        
        # check database exists (false = fail)
        if not self.__checkDatabaseExists():
            raise Exception('Database Does Not Exist. DB Update Aborted')
        
        # check files exist (false = fail)  
        if not self.__checkFilesExist():
            raise Exception('File Errors. DB Update Aborted')        
        
        return True

    def checkEntryCreateArgs(self):
        # check args for entry insert
        
        #required_args = ['db','sampid','mutfile','muttool','svfile','svtool','snp','met','depth']        
        #required_args = ['db','sampid','mutfile','muttool']  
        required_args = ['db','sampid','mutfile','muttool','svfile','svtool','snp','depth']  
        
        # check required arguments (false = fail)
        if not self.__checkArgsExist(required_args):
            raise Exception('Required Args Not Provided. Entry Creation Aborted')          
        
        # check database exists (false = fail)
        if not self.__checkDatabaseExists():
            raise Exception('Database Does Not Exist. Entry Creation Aborted')
        
        # check sample id exists (true = fail)
        if self.__checkEntryExists():
            raise Exception('Sample ID Already Exists. Entry Creation Aborted')
        
        # check files exist (false = fail)  
        if not self.__checkFilesExist():
            raise Exception('File Errors. Entry Creation Aborted')
        
        return True

    def checkEntryUpdateArgs(self):
        # check args for entry update
        
        required_args = ['db','sampid'] 
        possible_args = ['normid','mutfile','svfile','snp','met','depth','log','patid','sex']
        
        # check required arguments (false = fail)
        if not self.__checkArgsExist(required_args):
            raise Exception('Required Args Not Provided. Entry Update Aborted')         
        
        # check at least one argument added
        if not self.__checkPossibleArgExists(possible_args):
            raise Exception('No Args Added to Update. Entry Update Aborted')        
        
        # check database exists (false = fail)
        if not self.__checkDatabaseExists():
            raise Exception('Database Does Not Exist. Entry Update Aborted')        

        # check sample id exists (false = fail)
        if not self.__checkEntryExists():
            raise Exception('Sample ID Does Not Exist. Entry Update Aborted')
        
        # check files exist (false = fail)  
        if not self.__checkFilesExist():
            raise Exception('File Errors. Entry Update Aborted')
        
        return True

    def checkEntryDeleteArgs(self):
        # check args for entry delete
        
        required_args = ['db','sampid']        
        
        # check required arguments (false = fail)
        if not self.__checkArgsExist(required_args):
            raise Exception('Required Args Not Provided. Entry Deletion Aborted')         
        
        # check database exists (false = fail)
        if not self.__checkDatabaseExists():
            raise Exception('Database Does Not Exist. Entry Deletion Aborted')
        
        # check sample id exists (false = fail)
        if not self.__checkEntryExists():
            raise Exception('Sample ID Does Not Exist. Entry Deletion Aborted')        
        
        return True        
        
    #---------------------------------------

    def __checkDatabaseExists(self):
        # check arg database exists
        
        db = self._args['db']
        
        if not isfile(db):
            return False
        if getsize(db) < 100: # SQLite database file header is 100 bytes
            return False
    
        with open(db, 'rb') as fd:
            header = fd.read(100)

        return header[:16].decode("utf-8") == 'SQLite format 3\x00'

    def __checkEntryExists(self):

        lu_sql = """SELECT * FROM GENOM_SAMPLE WHERE sample_id = ?; """
        lu_var = [self._args['sampid']]
        
        try:  
            conn = sql.connect(self._args['db'])
            
            cursor = conn.cursor()
            cursor.execute(lu_sql,lu_var)
            rows = cursor.fetchall()

            cursor.close()
            conn.close()
    
            if rows:
                return True
            else:
                return False
            
        except Error as e:
            conn.close()
            print('Error Checking Sample ID: {}'.format(e)) 
            raise Exception('Aborting Execution') 


    def __checkArgsExist(self, req_list):
        # check args exist against required list
        
        args_exist = True
    
        for arg in req_list:
            
            if not self._args[arg]:
                
                print('Argument ({}) Is Required'.format('-{}'.format(arg)))
                args_exist = False
        
        for pair in self.bicond_args:
            
            if self._args[pair[0]]:
            
                if not self._args[pair[1]] or len(self._args[pair[0]]) != len(self._args[pair[1]]):
                
                    print('Need ({}) For Every ({})'.format('-{}'.format(pair[1]),'-{}'.format(pair[0])))
                    args_exist = False
        
        return args_exist


    def __checkPossibleArgExists(self, possible_list):
        # check at least one of the possible args exists
        
        arg_exists = False
        
        for arg in possible_list:
            
            if self._args[arg]:
                
                arg_exists = True        
        
        return arg_exists 

    
    def __checkFilesExist(self):
        # check files exist
        
        files_exist = True
        
        for arg in self.file_args:

            if self._args[arg]:

                files = []

                if type(self._args[arg]) is list:
                    for file in self._args[arg]:
                        files.append(file)
                else:
                    files.append(self._args[arg])
                    
                for file in files:
                    if not isfile(file):
                        print('Invalid File: arg: ({}) path: ({})'.format('-{}'.format(arg), file))
                        files_exist = False                    
                    
        return files_exist        


    
    #---------------------------------------    
    # Properties
    #---------------------------------------    
        
    @property    
    def database(self):
        return self._args['db']

    @property    
    def refGenome(self):
        return self._args['refgen']
    
    @property    
    def pipeline(self):
        return self._args['pipeline']    

    @property    
    def bedFile(self):
        return self._args['bedfile']
    
    @property    
    def bedType(self):
        return self._args['bedtype']

    @property    
    def targetBed(self):
        return self._args['targetbed']
    
    @property    
    def groupBed(self):
        return self._args['groupbed']    

    @property    
    def sampleID(self):
        return self._args['sampid']

    @property    
    def normalID(self):
        return self._args['normid']

    @property    
    def mutationFile(self):
        return self._args['mutfile']

    @property    
    def mutationTool(self):
        return self._args['muttool']

    @property    
    def structVarFile(self):
        return self._args['svfile']

    @property    
    def structVarTool(self):
        return self._args['svtool']

    @property    
    def snp(self):
        return self._args['snp']

    @property    
    def metrics(self):
        return self._args['met']

    @property    
    def depth(self):
        return self._args['depth']

    @property    
    def log(self):
        return self._args['log']   
    
    @property    
    def patientID(self):
        return self._args['patid']
    
    @property    
    def sex(self):
        return self._args['sex']      
        
    
    @property    
    def current_dir(self):
        return self.__curr_dir
    @current_dir.setter    
    def current_dir(self,curr_dir):
        self.__curr_dir = curr_dir  
        
    @property    
    def edit_mode(self):
        return self.__edit_mode
    @edit_mode.setter    
    def edit_mode(self,edit_mode):
        self.__edit_mode = edit_mode         
        
#prog = arg_parser(sys.argv[1:])