#!/usr/bin/env python

from dlePy.vasp.job_monitor import *
import sys
import os

if __name__ == "__main__":
    # Submitting command
    submit_cmd           = 'sbatch job'
    #Time (s) from last write to OUTCAR for considering job is not running  
    time_from_last_write = 600
    #Time (s) between queue checks
    waiting_time         = 1500

    # Provide a list of jobs, in dictionary
    # jobs = { job_name1:{ info }, job_name2: { info }, ...}
    # where { info } are dictionaries with  
    # {  'dir' : directory,
    #    'status': "NOT FINISHED", "FINISHED", "WAITING" }
    # 
    # It is best to use the following file/folder structures
    # Top dir: POSCAR_name, where `name` are names of the jobs
    # Sub dirs: all subdirs named `name`. Inside each subdir, INCAR, POTCAR, KPOINTS, POSCAR and job are expected.
    confs = sorted( [ x for x in os.listdir( '.' ) if os.path.isfile( x ) and 'POSCAR' in x ] )
    jobs = {}
    for conf in confs:
        jobs[ conf.replace( 'POSCAR_', '' ) ] = { 'dir'    : conf.replace( 'POSCAR_', '' ),
                                                  'status'  : "NOT FINISHED" } 

    # For email notification (set send_email = True )
    send_email = False 
    # For SMTP setting, you can use free setting from mailjet.com
    # and set environment variable for using the following code, 
    # or manually edit the setting 
    email_setting = {
    'from' : os.environ[ 'from_address' ],
    'to'   : os.environ[ 'to_address' ],
    'SMTP' : {
             'server'   : os.environ[ "mail_server" ],
             'port'     : int( os.environ[ "mail_port" ] ),
             'username' : os.environ[ "mail_login" ] ,
             'password' : os.environ[ "mail_password" ],
              }
    }

    run_job_monitor( jobs, time_from_last_write, waiting_time, submit_cmd, send_email = send_email, email_setting = email_setting )

