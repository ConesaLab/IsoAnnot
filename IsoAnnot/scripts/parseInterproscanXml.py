#!/usr/bin/env python2
# -*- coding: utf-8 -*-
""" 
Written by H. del Risco
hdelrisco@ufl.edu
Modified by Carlos Martínez and Alessandra Martínez
"""

from __future__ import division
from __future__ import print_function

import os, argparse
import xml.etree.ElementTree as ET


def genip(files, outFilepath):

	if os.path.isfile(outFilepath):
		os.remove(outFilepath)

	# loop through all files
	totalfiles = 0
	fout = open(outFilepath, 'w')

	for fpath in files:
		totalfiles += 1
		fnddata = False
		tree = ET.parse(fpath)
		root = tree.getroot()
		namespaces = {'ipro' : "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5"}
		n = 0
		for protein in root.findall('ipro:protein', namespaces):
			xrefs = protein.findall('ipro:xref', namespaces)
			if xrefs != None:
				for xref in xrefs: 
					idvalue = xref.get('id')
					for matches in protein.findall('ipro:matches', namespaces):
						for match in matches:
							for signature in match.findall('ipro:signature', namespaces):
								lib = signature.find('ipro:signature-library-release', namespaces)
								if lib != None:
									library = lib.get('library')
									if library == "PFAM" or library == "TMHMM" or library == "SIGNALP_EUK" or library == "MOBIDB_LITE" or library == "COILS":  # COILS??
										for locations in match.findall('ipro:locations', namespaces):
											for location in locations:
												n +=1
												if library == "PFAM":
													strpos = location.get('start')
													score = location.get('score')
													evalue = location.get('evalue')
													endpos = location.get('end')
													name = signature.get('name')
													desc = signature.get('desc')
													ac = signature.get('ac')
												elif library == "TMHMM":
													strpos = location.get('start')
													score = 'NA'
													evalue = 'NA'
													endpos = location.get('end')
													name = 'Transmembrane Helix'
													desc = signature.get('desc')
													ac = signature.get('ac')
												elif library == "SIGNALP_EUK":
													strpos = location.get('start')
													score = location.get('score')
													evalue = 'NA'
													endpos = location.get('end')
													name = signature.get('name')
													desc = 'Signal Peptide'
													ac = signature.get('ac')
												elif library == "MOBIDB_LITE":
													strpos = location.get('start')
													score = 'NA'
													evalue = 'NA'
													endpos = location.get('end')
													name = signature.get('name')
													desc = signature.get('desc')
													ac = signature.get('ac')
												elif library == "COILS":
													strpos = location.get('start')
													score = 'NA'
													evalue = 'NA'
													endpos = location.get('end')
													name = signature.get('name')
													desc = 'Coiled coil'
													ac = 'Coiled coil'
												
												fout.write(idvalue+"\t"+library+"\t"+"DOMAIN"+"\t"+ac+"\t"+name+"\t"+desc+"\t"+strpos+"\t"+endpos+"\t"+evalue+"\t"+score+"\n")

def main():
	parser = argparse.ArgumentParser(description="parseInterproScan")
	parser.add_argument('--output', required=True)
	parser.add_argument('--interproscan_files', nargs='+', required=True)

	args = parser.parse_args()

	errflg = genip(args.interproscan_files, args.output)
	if errflg:
		raise SystemExit(1)

if __name__ == "__main__":
   main()

