#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    File Processor
    
    Logic needed to confirm file formats and prepare for import

"""

from os.path import isfile
from pandas import read_csv
from sqlite3 import Binary

class file_processor(object):

    def __init__(self, current_dir):
        self.__current_dir = current_dir
    
    def processLog(self, file_path):
        # process sequence log
        
        if isfile(file_path):
            with open(file_path, 'rb') as log_file:
                ablob = log_file.read()
                return Binary(ablob)
        else:
            return None    
    
    def processTsv(self, suffix, file_path):
        # process tsv
        
        template_path = '{}/templates/TEMPLATE{}.txt'.format(self.__current_dir, suffix)
        
        frame = read_csv(file_path,delimiter='\t',encoding='utf-8',header=(0))
        template = read_csv(template_path,delimiter='\t',encoding='utf-8',header=(0))       
        
        if [x.lower() for x in list(frame.columns.values)] == [x.lower() for x in list(template.columns.values)]:
            return frame
        else:
            return None 
   
        
    def processVcf(self, suffix, file_path):
        # process vcf
        
        template_path = '{}/templates/TEMPLATE{}.txt'.format(self.__current_dir, suffix)
        
        header = None        
        lines = 0
        
        with open(file_path, 'r') as file:
            for num, line in enumerate(file, 1):
                if line.startswith("#CHROM"):
                    lines = num - 1

        with open(file_path, 'r') as file:                    
            header = [l for l in file if l.startswith('#')]

        header =  "".join(header)

        frame = read_csv(file_path,delimiter='\t',encoding='utf-8',header=(0),skiprows=lines)
        template = read_csv(template_path,delimiter='\t',encoding='utf-8',header=(0)) 

        frame.columns = frame.columns.str.replace('#CHROM','CHROM')
        
                                                  
        if len(frame.columns) == 11:
            frame.columns.values[9] = 'NORMAL'
            frame.columns.values[10] = 'TUMOR'
        elif len(frame.columns) == 10:
            frame.columns.values[9] = 'TUMOR'
            frame.columns.insert(9,'NORMAL',None,True)
            
        if [x.lower() for x in list(frame.columns.values)] == [x.lower() for x in list(template.columns.values)]:
            return [frame, header]
        else:
            return None
        
    def processBed(self, bedtype, file_path):       
        # process bed file
        
        header =['chrom','start','end','id']
        
        template_path = '{}/templates/TEMPLATE_BED.txt'.format(self.__current_dir)

        frame = read_csv(file_path,delimiter='\t',encoding='utf-8',usecols=[0,1,2,3],names=header)
        
        template = read_csv(template_path,sep='\t',encoding='utf-8',header=(0))       
        
        if [x.lower() for x in list(frame.columns.values)] == [x.lower() for x in list(template.columns.values)]:
            return frame
        else:
            return None 
        
       
        
        