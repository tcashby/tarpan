#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Delete Entry

"""

# Entry Delete
# delete_db_entry.py -db "<WHATEVER.db>" -sampid "<SAMPLE_ID>" 

import sys
from inspect import getsourcefile
from os.path import abspath,dirname
from modules.arg_parser import arg_parser
from modules.entry_processor import entry_processor

args = arg_parser(sys.argv[1:])

try:
    
    if args.checkEntryDeleteArgs():
        
        current_dir = dirname(abspath(getsourcefile(lambda:0)))    
        
        args.current_dir = current_dir
        args.edit_mode = 'DELETE'        
        
        entry = entry_processor(args)

        
except Exception as e:
    print('Exception: {}'.format(e))