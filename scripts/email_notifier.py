#!/usr/bin/env python3
__author__ = 'tomarovsky'
import smtplib
from email.mime.text import MIMEText
from argparse import ArgumentParser

LABEL = 'Message from server'

def main():
    me = 'atomarovsky@mcb.nsc.ru'
    you = 'st079639@student.spbu.ru'
    smtp_server = 'mcb.nsc.ru'
    msg = MIMEText(args.text)
    msg['Subject'] = LABEL
    msg['From'] = me
    msg['To'] = you
    s = smtplib.SMTP(smtp_server)
    s.sendmail(me, [you], msg.as_string())
    s.quit()

if __name__ == "__main__":
    parser = ArgumentParser(description="for sending messages")
    group_required = parser.add_argument_group('Options')
    group_required.add_argument('-t', '--text', type=str,
                                help="message text")
    args = parser.parse_args()
    main()