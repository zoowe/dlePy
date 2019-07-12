import os
import time

import sys
import smtplib

def get_SMTP( ):
    SMTP = { }
    SMTP[ 'server' ] = os.environ[ "mail_server" ] 
    SMTP[ 'port' ]   = int( os.environ[ "mail_port" ] )
    SMTP[ 'username' ] = os.environ[ "mail_login" ] 
    SMTP[ 'password' ] = os.environ[ "mail_password" ]
    return SMTP
 
def send_email ( fromaddr, toaddrs, message, SMTP = get_SMTP( )):
    server = smtplib.SMTP( SMTP[ 'server' ], SMTP[ 'port' ] )
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login( SMTP[ 'username' ], SMTP[ 'password' ] )
    server.sendmail(fromaddr, toaddrs,  message)
    server.quit()

