#!/usr/bin/env python

from dlePy.jobmon.monitor import *
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

    jobs = {}
    list_dir=os.listdir('.')
    for rootdir in list_dir:
        if os.path.isdir(rootdir):
            for root, subFolders, files in os.walk(rootdir):
                for file_name in files:
                    if file_name == 'POSCAR' and 'RUN' not in root.strip() :
                        loc = root.strip() 
                        name = 'Pt_' + loc.replace('/','.' ).replace( 'CENTER','C').replace ( 'SUB','S' ).replace( 'VAC', 'V' )
                        jobs[ name ] = { 'dir'    : loc,
                                         'status' : "NOT FINISHED" }
    monitorID = 'NAME OF YOUR WORK'

    # For email notification (set send_email = True )
    send_email = True
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

    run_job_monitor( jobs, time_from_last_write = 600, time_to_check= 1500, submit_cmd = "sbatch job", send_email = True, email_setting = email_setting,
                     monitorID = monitorID,
                     if_job_done = if_vasp_done, output = 'OUTCAR', backup_files = [ 'POSCAR', 'OUTCAR', 'INCAR', 'KPOINTS', 'vasprun.xml', 'CONTCAR' ],
                     create_new_run = create_new_run_vasp  )

