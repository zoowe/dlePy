import os
import sys
import time
import getpass
import subprocess
from ..send_email import send_email
def get_mtime( file ):
    ( mode, ino, dev, nlink, uid, gid, size, atime, mtime, ctime ) = os.stat( file )
    return mtime

def if_vasp_done( outcar ):
    job_done = False
    niter    = 0
    NWS      = 10000
    mtime    = 0
    if os.path.isfile( outcar ):
        with open( outcar, 'r' ) as f:
            lines = f.readlines( )
        if 'Voluntary context switches' in lines[ -1 ]:
            job_done = True
        i = len( lines )
        while i > 0:
            i -= 1
            if "Iteration " in lines[ i ]:
                niter = lines[ i ].replace( '-','' ).replace( '(',' ' ).replace( ')','' ).split( )[ -2 ]
                niter = int( niter ) - 1
                for j in range( i, len( lines ) ):
                    if "free  energy   TOTEN" in lines[ j ]:
                        niter += 1
                        break
                i = 0
        while i < len( lines ):
            if "number of steps for IOM" in lines[ i ]:
                NWS = lines[ i ].replace( '-','' ).replace( '(',' ' ).replace( ')','' ).split()[ 2 ]
                NWS = int( NWS )
                i = len( lines )
                break
            i += 1

        mtime = get_mtime( outcar ) - time.time( )
    return job_done, niter, NWS, mtime

def get_jobs_in_queue( max_name_length = 20 ):
    os.environ["SQUEUE_FORMAT"]="%.18i %.9P %." + str( max_name_length )+ "j %.8u %.2t %.10M %.6D %.20R %q"
    username     = getpass.getuser( )
    proc         = subprocess.Popen( [ 'squeue', '-u', username ],stdout=subprocess.PIPE )
    lines        = proc.stdout.readlines( )
    job_in_queue = {}
    for i in range( 1, len( lines ) ):
        job_name                 = lines[ i ].split()[ 2 ]
        job_status               = lines[ i ].split()[ 4 ]
        job_in_queue[ job_name ] = job_status 
    return job_in_queue

def get_run_num( ):
    runs = [ int( x.replace( 'RUN', '' ) ) for x in os.listdir( '.' ) if 'RUN' in x and os.path.isdir( x ) ]
    return len( runs )

def backup( run, backup_files = [ 'POSCAR', 'OUTCAR', 'INCAR', 'KPOINTS', 'vasprun.xml', 'CONTCAR' ]  ):
    folder = 'RUN' + run
    print ( "%s" %( 'BACKUP TO FOLDER: ' + folder ) )
    os.system( 'mkdir -p ' + folder )
    for item in backup_files:
        os.system( 'cp -a ' + item + ' ' + folder + '/' )

def create_new_run( run ):
    folder = 'RUN' + run
    os.system( 'cp -a ' + folder + '/CONTCAR POSCAR' )

def submitjob( cmd = 'sbatch -q scavenger --time-min=4:0:00 job' ):
    os.system( cmd )

def resubmit( jobs, jobs_in_queue, time_from_last_write, submit_cmd, send_email = False, email_setting = {} ):
    top_dir = os.getcwd( )
    for ID, name in enumerate( sorted( jobs.keys( ) ) ):
        print ( "----------------------------------" )
        print ( "### %s: " %( name ) )
        if jobs[ name ][ 'status' ] == "FINISHED":
            print ( 'FINISHED' )
        else:
            os.chdir( jobs[ name ][ 'dir' ] )
            job_done, niter, nmaxiters, mtime = if_vasp_done( 'OUTCAR' )
            if not job_done:
                print ( 'OUTCAR is not done' )
            else:
                print ( 'OUTCAR is done' )
            print ( 'Number of iterations %i' %( niter ) )
            print ( 'OUTCAR was last updated %i s ago' %( -mtime ) )
            if name not in jobs_in_queue:
                print ( 'Job is not in queue or running' )
            else:
                print ( 'Job was submitted and is %s' %( jobs_in_queue[ name ] ) )
            if ( not job_done or niter == nmaxiters ) and ( mtime < -time_from_last_write ) and ( name not in jobs_in_queue.keys( ) ) and niter > 0:
                run =  get_run_num( ) + 1
                run = str( run )
                backup( run )
                create_new_run( run )
                print ( "%s" %( 'RESUBMITING... ') )
                sys.stdout.flush( )
                submitjob( submit_cmd )
                jobs[ name ][ 'status' ] = "SUBMITTED"
            elif job_done and niter < nmaxiters:
                print ( "%s" %( 'JOB DONE!!!' ) )
                jobs[ name ][ 'status' ] = "FINISHED"
            elif ( name not in jobs_in_queue ) and ( niter == 0 ):
                print ( "%s" %( 'JOB ERROR. Please check and resubmit this job!!!' ) )
                if jobs[ name ][ 'status' ] != "ERROR" and send_email:
                    try:
                        send_error_email( top_dir, name, email_setting )
                        print ( "Email notification was sent" )
                    except:
                        print ( "Email can not be sent" )
                jobs[ name ][ 'status' ] = "ERROR"
            else:
                print ( "%s" %( 'WAITING...' ) )
                jobs[ name ][ 'status' ] = "WAITING"

            os.chdir( top_dir )

    # Update status of WAITING job
    jobs_in_queue = get_jobs_in_queue( max_name_length = get_name_length( jobs ) )
    for name in jobs.keys( ):
        if jobs[ name ][ 'status' ] in [ "WAITING", "SUBMITTED" ]:
            if name in jobs_in_queue.keys( ):
                 jobs[ name ][ 'status' ] = jobs_in_queue[ name ] 
                
    sys.stdout.flush()

    return jobs

def count_jobs( jobs ):
    nfinished = 0
    ntotal    = 0
    nerror    = 0
    for name in sorted( jobs.keys( ) ):
        ntotal += 1
        if jobs[ name ][ 'status' ] == "FINISHED":
            nfinished += 1
        if jobs[ name ][ 'status' ] == "ERROR":
            nerror    += 1
    return nfinished, ntotal, nerror

def check_jobs( jobs, top_dir ):
    for ID, name in enumerate( sorted( jobs.keys( ) ) ):
        if len( name ) > 50:
            print ( "%s: ERROR!!! Job name must have no more than 20 characters"  %( name ) )
            sys.exit( )

        os.chdir( top_dir )
        if os.path.isdir( jobs[ name ][ 'dir' ] ):
            if os.path.isfile( jobs[ name ][ 'dir' ] + '/job' ):
                with open( jobs[ name ][ 'dir' ] + '/job' ) as f:
                    lines = f.readlines( )
                name_match = False
                for line in lines:
                    if "#SBATCH -J" in line:
                        if line.split( )[ 2 ] == name:
                             name_match = True
                             break
                if not name_match:
                    print ( "%s: ERROR!!! Expecting a line \"#SBATCH -J %s\" in %s" %( name, name, jobs[ name ][ 'dir' ] + '/job' ) )
                    sys.exit( )
            else:
                print ( "%s: ERROR!!! Expecting a job script file %s" %( name, jobs[ name ][ 'dir' ] + '/job' ) )
                sys.exit( )
        else:
            print ( "%s: ERROR!!! Expecting a directory %s" %( name, jobs[ name ][ 'dir' ] ) )
            sys.exit( )

def print_info( jobs, top_dir, jobs_in_queue ):
    print ( "==========================================" )
    print ( "Automatic monitoring and resubmitting jobs" )
    print ( "==========================================" )
    print ( " " )
    print ( "List of jobs" )
    print ( "%-10s %-20s %-15s %-s" %( 'ID', 'JOB_NAME', 'STATUS', 'DIRECTORY' ) )
    print ( "---------------------------------------------------------------" )
    n_jobs_in_queue = 0
    for ID, name in enumerate( sorted( jobs.keys( ) ) ):
        print ( "%-10s %-20s %-15s %-s" %( ID, name,
                                             jobs[ name ][ 'status' ],
                                             jobs[ name ][ 'dir' ]      ) )
        if name in jobs_in_queue.keys( ):
             n_jobs_in_queue += 1 
    print ( " " )
    print ( "Top dir: %s" %( top_dir ) )
    print ( " " )
    print ( "Number of jobs currently in queue/running: %s" %( n_jobs_in_queue ) )
    print ( " " )

def send_error_email( top_dir, name, email_setting ):
    message  = 'Subject: ERROR ' + name
    message += '\n'
    message += '\nThis job has an error. No iteration'
    message += '\nJob location: %s' %( top_dir )
    send_email ( email_setting[ 'from' ], email_setting[ 'to' ], message, email_setting[ 'SMTP' ])

def send_finish_email( top_dir, monitorID, jobs, email_setting ):
    message  = 'Subject: FINSIHED: ' + monitorID 
    message += '\n'
    message +="\n=========================================="
    message +="\nAutomatic monitoring and resubmitting jobs"
    message +="\n==========================================\n\n"
    message +="\nList of jobs"
    message +="\n%-10s %-20s %-15s %-s" %( 'ID', 'JOB_NAME', 'STATUS', 'DIRECTORY' )
    message +="\n---------------------------------------------------------------"
    for ID, name in enumerate( sorted( jobs.keys( ) ) ):
        message +="\n%-10s %-20s %-15s %-s" %( ID, name,
                                             jobs[ name ][ 'status' ],
                                             jobs[ name ][ 'dir' ]      )
    message +="\n\nTop dir: %s" %( top_dir )
    nfinished, ntotal, nerror = count_jobs( jobs )

    message += '\nThis job is finished. %i/%i are succesfully done.' % ( nfinished, ntotal )
    send_email ( email_setting[ 'from' ], email_setting[ 'to' ], message, email_setting[ 'SMTP' ] )

def get_name_length( jobs ):
    l = [ len( name ) for name in jobs.keys() ]
    return sorted( l )[ -1 ]

def run_job_monitor( jobs, time_from_last_write = 600, time_to_check= 1500, submit_cmd = "sbatch job", send_email = False, email_setting = { }, monitorID = 'Automatic Monitoring' ):
    nfinished, ntotal, nerror = count_jobs( jobs )
    jobs_in_queue             = get_jobs_in_queue( max_name_length = get_name_length( jobs ) )
    top_dir                   = os.getcwd( )
    print_info( jobs, top_dir, jobs_in_queue )
    check_jobs( jobs, top_dir )

    
    for ID, name in enumerate( sorted( jobs.keys( ) ) ):
        os.chdir( jobs[ name ][ 'dir' ] )
        if not os.path.isfile( 'OUTCAR' ) and name not in jobs_in_queue.keys( ) and jobs[ name ][ 'status' ] != 'FINISHED':
            print ( "%s was not submitted, Submitting..." %( name ) )
            sys.stdout.flush( )
            submitjob( submit_cmd )
        os.chdir( top_dir )
    while nfinished < ntotal:
        print ( time.ctime( ) )
        jobs_in_queue             = get_jobs_in_queue( max_name_length = get_name_length( jobs ) )
        jobs                      = resubmit( jobs, jobs_in_queue, time_from_last_write, submit_cmd, send_email, email_setting )
        sys.stdout.flush()
        nfinished, ntotal, nerror = count_jobs( jobs )
        print ( " " )
        print ( "Number of finished jobs: %s out of %s" %( nfinished, ntotal) )
        if nfinished + nerror == ntotal:
            print ( "Please check the jobs with errors, resubmit it and rerun this monitoring program" )
            break 
        if nfinished < ntotal:
            jobs_in_queue = get_jobs_in_queue( max_name_length = get_name_length( jobs ) )
            print_info( jobs, top_dir, jobs_in_queue )
            print ( " " )
            print ( "Waiting for %s s" %( time_to_check ) )
            sys.stdout.flush( )
            os.system( 'sleep ' + str( time_to_check ))

    # FINISHING MONITORING
    if send_email:
        send_finish_email( top_dir, monitorID, jobs, email_setting )
    print_info( jobs, top_dir, jobs_in_queue )
    print ( "AUTOMATIC MONITORING FINISHED!!!" )
