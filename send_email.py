#!/usr/bin/env python
# -*- coding: utf-8 -*-

# import modules
import os
import sys
import argparse
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

def main():
  parser = argparse.ArgumentParser(description = "Script to send an email from a Snakefile")
  parser.add_argument("-t", "--to-address", dest = "address", help = "Email address to send to")
  parser.add_argument("-s", "--subject", dest = "subject", help = "Subject line", default = "Automated email from Robo-Phil!")
  parser.add_argument("-b", "--body", dest = "body", help = "File with body text in it", default = "Dat body doe")
  if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

  args = parser.parse_args()

  fromaddr = "robophilross@gmail.com"
  toaddr = args.address
  msg = MIMEMultipart()
  msg['From'] = fromaddr
  msg['To'] = toaddr
  msg['Subject'] = args.subject

  if os.path.isfile(args.body):
    with open(args.body) as f:
      body = f.read()
  else:
    body = args.body

  msg.attach(MIMEText(body, 'plain'))

  server = smtplib.SMTP('smtp.gmail.com', 587)
  server.starttls()
  server.login(fromaddr, "wutANicePassword1")
  text = msg.as_string()
  server.sendmail(fromaddr, toaddr, text)
  server.quit()
  
if __name__ == "__main__":
  try:
    main()
  except KeyboardInterrupt:
    sys.stderr.write("User interrupt! Ciao!\n")
    sys.exit(0)
