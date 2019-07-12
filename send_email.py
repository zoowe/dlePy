import os
import time

import sys
import smtplib

def get_SMTP( ):
    SMTP = { }
    try:
        SMTP[ 'server' ] = os.environ[ "mail_server" ] 
    except:
        SMTP[ 'server' ] = "in-v3.mailjet.com"
    try: 
        SMTP[ 'port' ]   = int( os.environ[ "mail_port" ] )
    except:
        SMTP[ 'port' ]   = 587
    try:
        SMTP[ 'username' ] = os.environ[ "mail_login" ] 
    except:
        SMTP[ 'username' ] = "" 
    try:
        SMTP[ 'password' ] = os.environ[ "mail_password" ]
    except:
        SMTP[ 'password' ] = ""
    return SMTP
 
def send_email ( fromaddr, toaddrs, message, SMTP = get_SMTP( )):
    server = smtplib.SMTP( SMTP[ 'server' ], SMTP[ 'port' ] )
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login( SMTP[ 'username' ], SMTP[ 'password' ] )
    server.sendmail(fromaddr, toaddrs,  message)
    server.quit()

