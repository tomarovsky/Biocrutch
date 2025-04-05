#!/usr/bin/env python3
__author__ = 'tomarovsky'
import smtplib
from email.mime.text import MIMEText
from argparse import ArgumentParser


LABEL = 'Message from server'


def main():
    me = args.sender
    you = args.recipient
    smtp_server = 'mcb.nsc.ru'
    msg = MIMEText(args.text)
    msg['Subject'] = LABEL
    msg['From'] = me
    msg['To'] = you
    s = smtplib.SMTP(smtp_server)
    s.sendmail(me, [you], msg.as_string())
    s.quit()


if __name__ == '__main__':
    parser = ArgumentParser(description='for sending messages')
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-t', '--text', type=str, help='message text')
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-s', '--sender', type=str, help='server mail')
    group_additional.add_argument('-r', '--recipient', type=str, help='mail address')
    args = parser.parse_args()
    main()
