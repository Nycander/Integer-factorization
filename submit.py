#!/usr/bin/env python
import smtplib
import ConfigParser
import getopt
import logging
import os
import sys

from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText

_DEFAULT_CONFIG='/afs/nada.kth.se/home/3/u1rs4bi3/Public/etc/kattisrc'
_VERSION = 'Version: $Id: submit.py,v 1.3 2006/11/02 11:13:27 gkreitz Exp $'
_HELP_MSG = \
'''Usage: %s [-f] [-p <pid>] [-m <mainclsss>] [-l <lang>] [-L <libpath>] [-t <tag>] [-h?] files ...

  -p specifies problem-id. Will override default guess (first part of first filename)
  -m specifies mainclass. Will override default guess (first part of first filename)
  -l specifies language. Will override default guess (bases on suffix of first filename)
  -L specifies library directory to upload files to
  -f force, no confirmation prompt
  -h/-? displays this help'''

_RC_HELP = \
'''I failed to read in a config file from your home directory. You need to put a file 
named .kattisrc in your home directory. It should contain the following (with values
replaced in obvious ways):
[user]
username: gkreitz
password: *********
email: gkreitz@nada.kth.se'''

_LANGUAGE_GUESS = { 'java' : 'Java', 'c' : 'C', 'cpp' : 'C++', 'h' : 'C++', 'cc' : 'C++', 'cxx' : 'C++', 'c++' : 'C++' }

def makemail(cfg, files, bodyfields):
	ret = MIMEMultipart()
	bodystr = ''.join(map(lambda (k,v) : '%s : %s\n' % (k,v), bodyfields.iteritems()))
	ret.attach(MIMEText(bodystr))
	seen = { }

	for filename in files:
		if seen.has_key(filename):
			continue
		seen[filename] = True
		f = open(filename, 'r')
		msg = MIMEText(f.read(), 'plain', 'iso8859-1')
		f.close()
		msg.add_header('Content-Disposition', 'attachment', filename=os.path.split(filename)[1])
		ret.attach(msg)
	ret['Subject'] = 'Automatically generated Marvin submission'
	ret['From'] = cfg.get('user', 'email')
	ret['To'] = cfg.get('kattis', 'email')
	# Guarantees ending newline
	ret.epilogue=''
	return ret

def sendmail(cfg, msg):
	s = smtplib.SMTP(cfg.get('kattis', 'hostname'), cfg.getint('kattis', 'smtpport'))
	s.sendmail(cfg.get('user', 'email'), cfg.get('kattis', 'email'), msg.as_string())
	s.close()

def showhelp():
	print _VERSION
	print _HELP_MSG % ('submit.py')
	sys.exit(0)

if __name__ == '__main__':
	cfg = ConfigParser.ConfigParser()
	try:
		cfg.readfp(open(os.path.join(os.getenv('HOME'), '.kattisrc')));
	except IOError, (errno, strerror):
		print _RC_HELP
		sys.exit(1)

	problem=None
	mainclass=None
	librarypath=None
	force=False
	tag=None
	(opts, args) = getopt.getopt(sys.argv[1:], 'p:m:l:L:t:h?f')

	if(len(args) == 0):
		showhelp()

	filename = os.path.split(args[0])[1]
	if filename == '':
		print "Expected filename, got directory"
	spl = filename.split('.')
	problem = spl[0]
	if len(spl) > 1:
		if _LANGUAGE_GUESS.has_key(spl[-1]):
			language = _LANGUAGE_GUESS[spl[-1]]

	for (opt, arg) in opts:
		if opt == '-p':
			problem = arg
		elif opt == '-m':
			mainclass = arg
		elif opt == '-l':
			language = arg
		elif opt == '-L':
			librarypath = arg
		elif opt == '-t':
			tag = arg
		elif opt == '-h' or opt == '-?':
			showhelp()
		elif opt == '-f':
			force = True

	if mainclass == None and language.upper() == 'JAVA':
		mainclass = spl[0]
	
	bodyfields = { 'author' : cfg.get('user', 'username'), 'password' : cfg.get('user', 'password') }
	if librarypath != None:
		bodyfields['librarypath'] = librarypath
	else:
		bodyfields['problem'] = problem
		bodyfields['language'] = language
		if mainclass != None:
			bodyfields['mainclass'] = mainclass
		if tag != None:
			bodyfields['tag'] = tag
		if bodyfields['language'] == None:
			print 'No language specified and I failed to guess! Please specify a language.'
			sys.exit(1)
	mail = makemail(cfg, args, bodyfields)
	if not force:
		print 'About to send mail:'
		bodyfields['password'] = '********'
		keys = bodyfields.keys()
		keys.sort()
		for k in keys:
			print '%s: %s' % (k, bodyfields[k])
		print "Files:"
		for f in args:
			print '\t%s' % (f)

		print 'Ok? (y/N)'
		if sys.stdin.readline().upper()[:-1] != 'Y':
			print 'Cancelling'
			sys.exit(1)
	sendmail(cfg, mail)
	sys.exit(0)
