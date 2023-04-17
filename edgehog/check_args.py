#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import sys
import shutil


def check_out_dir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def check_path_exists(file):
    if not os.path.exists(file):
        sys.exit('error: Unable to locate \'%s\'' % file)
        
        
def check_file_extension(file, file_type):
    if file_type == 'orthoxml':
        if not re.match(r'.*?\.ortho(xml|XML)$', file):
            sys.exit('error: File \'%s\' has an incorrect extension: \'.orthoxml\', is expected for a HOGs file' % file)
    elif file_type == 'gff':
        if not re.match(r'.*?\.gff(3)?$', file):
            sys.exit('error: File \'%s\' has an incorrect extension: \'.gff\', \'.gff3\' are expected for a gff annotation file' % file)


def print_args(args, out_dir):
    print('output directory: \'%s\'' % out_dir)
    print('species tree:     \'%s\'' % args.species_tree)
    print('HOGs.orthoxml:    \'%s\'' % args.hogs)
    print('gff directory:    \'%s\'' % args.gff_directory)
    

def check_args(args):
    print('###################################')
    print('Checking arguments ...')
    # check cpu is integer
    # try:
    #     int(args.cpu)
    # except:
    #     sys.exit('error: Please provide an integer for the number of cpu to be used')  
    # check output directory exists
    out_dir = check_out_dir(os.path.abspath(args.output_directory))
    # check mandatory input files exist and have the correct extension
    check_path_exists(args.species_tree)
    check_path_exists(args.hogs)
    check_file_extension(args.hogs, 'orthoxml')
    if args.gff_directory:
        check_path_exists(args.gff_directory)
    elif args.hdf5:
        check_path_exists(args.hdf5)
    else:
        sys.exit('error: no annotations provided (e.g. gff directory)')
    # print config
    print_args(args, out_dir)
    return out_dir

