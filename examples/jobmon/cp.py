#!/usr/bin/env python
'''
This script copies job file to each folder and change $name to appropriate one
'''

from dlePy.vasp.job_monitor import *
import sys
import os

if __name__ == "__main__":
    # 
    # It is best to use the following file/folder structures
    # Top dir: POSCAR_name, where `name` are names of the jobs
    # Sub dirs: all subdirs named `name`. Inside each subdir, INCAR, POTCAR, KPOINTS, POSCAR and job are expected.

    jobs = {}
    list_dir=os.listdir('.')
    for rootdir in list_dir:
        if os.path.isdir(rootdir):
            for root, subFolders, files in os.walk(rootdir):
                for file_name in files:
                    if file_name == 'POSCAR' and 'RUN' not in root.strip():
                        loc = root.strip() 
                        name = loc.replace('/','.' )
                        jobs[ name ] = { 'dir'    : loc,
                                         'status' : "NOT FINISHED" }
    print ( jobs )

    loc = os.getcwd()

    for name in jobs.keys():
            print ( name, jobs[ name ][ 'dir' ] )
            os.system( 'cat job | sed \'s/$name/' + name + '/\' >  '+ jobs[ name ][ 'dir' ] + '/job' )
