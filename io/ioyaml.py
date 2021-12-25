from yaml import load
import gzip as gz

def loadyaml( fyaml ):
    try:
        with gz.open( fyaml, 'r' ) as f:
            data = load( f )
    except IOError:
        try:
            with open( fyaml, 'r' ) as f:
                data = load( f )
        except IOError:
            raise IOError( "IOError: Can not open %s" %( fyaml ) )

    return data
