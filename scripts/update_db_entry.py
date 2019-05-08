#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Update Entry

"""

# Entry Update (only one arg change required)
# update_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" -normid "<SAMPLE_ID_NORMAL>" -mutfile "<MUTATION file>" -muttool "<MUTATION tool>" -mutfile "<MUTATION2 file>" -muttool "<MUTATION2 tool>" -svfile "<STRUCTVAR file>" -svtool "<STRUCTVAR tool>" -snp "<SNPDIFF file>" -met "<METRICS file>" -depth "<DEPTH file>" -log "<LOG file>" -patid "<PAT ID>" -sex "<SEX>" 

import sys
from inspect import getsourcefile
from os.path import abspath,dirname
from modules.arg_parser import arg_parser
from modules.entry_processor import entry_processor

args = arg_parser(sys.argv[1:])

try:
    
    if args.checkEntryUpdateArgs():
        
        current_dir = dirname(abspath(getsourcefile(lambda:0)))    
        
        args.current_dir = current_dir
        args.edit_mode = 'UPDATE'        
        
        entry = entry_processor(args)

        
except Exception as e:
    print('Exception: {}'.format(e))