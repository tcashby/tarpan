#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    Update database

"""

# DB Update (only one change arg required)
# update_db.py -db "<WHATEVER.db>" -refgen "<REF_GENOME>" -pipeline "<PIPELINE>" -targetbed "<TARGET_BED>" -groupbed "<GROUP_BED>" -blacklist "<BLACKLIST_BED>" -bedfile "<BED_FILE>" -bedtype "<BED_TYPE>"


import sys
from inspect import getsourcefile
from os.path import abspath,dirname
from modules.arg_parser import arg_parser
from modules.db_processor import db_processor

args = arg_parser(sys.argv[1:])

try:
    
    if args.checkDbUpdateArgs():
        
        current_dir = dirname(abspath(getsourcefile(lambda:0)))    
        
        args.current_dir = current_dir
        args.edit_mode = 'UPDATE'
        
        db = db_processor(args)
        
except Exception as e:
    print('Exception: {}'.format(e))
