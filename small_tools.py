def file_exist( FILENAME ):
    """
    Check if FILENAME exists, is readable
    """
    from os import path, access, R_OK

    if path.exists( FILENAME ) and path.isfile( FILENAME ) and access( FILENAME, R_OK ):
        file_exists = True
        txt = 'File ' + FILENAME + ' exists and readable. You are good to go'
    else:
        file_exists = False
        txt= 'File ' + FILENAME + ' does not exist or is not readable.'
    
    return file_exists, txt  

def str_decode( s, decode = True ):
    if decode:
        return s.decode('utf8')
    else:
        return s

