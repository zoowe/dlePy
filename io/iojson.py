import json
import gzip as gz
import os

def loadjson( fjson ):
    try:
        with gz.open( fjson, 'r' ) as f:
            print ( "Load  data from %s" %( fjson ) )
            data = json.load( f )
    except IOError:
        try:
            with open( fjson, 'r' ) as f:
                print ( "Load  data from %s" %( fjson ) )
                data = json.load( f )
        except IOError:
            raise IOError( "IOError: Can not open %s" %( fjson ) )

    return data

def writejson( fjson, data ):
    try:
        with open( fjson, 'w' ) as f:
            print ( "Write data to %s" %( fjson ) )
            f.write( json.dumps( data ) )
            print ( "Compress %s" %( fjson ) )
        os.system( 'gzip ' + fjson )
    except IOError:
        raise IOError( "IOError: Can not open %s" %( fjson ) )


