#!/usr/bin/env python3
# -*- coding: utf8 -*-

'''
Author: Lorena de la Fuente Lorente
Modified by: Carlos Martínez and Alessandra Martinez

Script: Creating SQANTI-like files for Reference annotations.

Important: check that your chromosomes are called equal in gtf and fasta file!
For refseq we use the gff3 but after removing partial alignments!!! Why? becuase they don't fit when transfering uniprot motifs
'''

import argparse, gzip
from IsoAnnot import openfile, argparse_nullable, read_chr_ref_acc

# TODO: adapt SQANTI to python3 and use those definitions from IsoAnnot

class myTranscripts:

	def __init__(self, transcript , chrom, strand, length, exons,  refGeneName, refGeneId, refTranscriptName, db, primary_class = "NA", coding = "non_coding", refProtein = "NA", exons_coord = [], firstProtPos ="NA", lastProtPos = "NA"):

		self.transcript = transcript
		self.length = length
		self.exons = exons
		self.refTranscriptName = refTranscriptName
		self.refGeneName = refGeneName
		self.refGeneId = refGeneId
		self.db = db
		self.refProtein = refProtein
		self.strand = strand
		self.chrom = chrom
		self.primary_class = primary_class
		self.coding = coding
		self.exons_coord = []
		self.cds_coord = [] # Added by Ale
		self.firstProtPos = firstProtPos
		self.lastProtPos = lastProtPos


	def cds_transcript_coord(self):
		if self.firstProtPos == "NA":
			return ["NA", "NA", "NA", "NA"]
		else:
			cdsStart = 0
			cdsEnd = 0
			cds_len = 0
			v = []

			for exon in self.exons_coord:
				v.extend(range(exon[0], exon[1] + 1))

			if self.strand == "+":
				v = sorted(v)
				cdsStart = v.index(self.firstProtPos) +1 
				cdsEnd = v.index(self.lastProtPos) + 1
			else:
				v = sorted(v, reverse=True)
				cdsEnd = v.index(self.firstProtPos) +1
				cdsStart = v.index(self.lastProtPos) + 1

			# To calculate correctly the CDS length
			for CDS in self.cds_coord:
				cds_len += CDS[1] - CDS[0] + 1
			# cds_len = cdsEnd - cdsStart + 1

			orf_length = cds_len/3

			return [cdsStart, cdsEnd, cds_len, orf_length]


	def NMD_calculation(self):
		if self.firstProtPos == "NA":
			return "NA"
		else:
			exonsStart = [x[0] for x in self.exons_coord]
			exonsEnd = [x[1] for x in self.exons_coord]
			exonsStart = sorted(exonsStart)
			exonsEnd = sorted(exonsEnd)
			nmd = "nonNMD"
			if self.strand == "+":

				if self.lastProtPos >= exonsStart[-1]:
					nmd = "nonNMD"

				else:
					if (exonsEnd[-2] - self.lastProtPos + 1) < 55:
						nmd = "nonNMD"
					else:
						nmd = "NMD"

			if self.strand == "-":

				if self.firstProtPos <= exonsEnd[0]:
					nmd = "nonNMD"

				else:
					if (self.firstProtPos - exonsStart[1] + 1) < 55:
						nmd = "nonNMD"
					else:
						nmd = "NMD"
			return nmd


	def __str__(self):
		na_string = "\tNA" * 16
		return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.strand, str(self.length), str(self.exons), str(self.primary_class), str(self.refGeneName), self.refTranscriptName, str(self.length), str(self.exons), na_string, str(self.coding), str(self.cds_transcript_coord()[3]), str(self.cds_transcript_coord()[2]), str(self.cds_transcript_coord()[0]), str(self.cds_transcript_coord()[1]), str(self.refProtein))



class myQueryJunctions:

	def __init__(self, transcript, junctionNumber, chrom, strand, junc_class, donor, acceptor, diff_start, diff_end, genomCoordStart, genomCoordEnd, transcriptCoord, spliceSite, canonical, indel="NA", coverage=None):

		self.transcript      = transcript
		self.junctionNumber  = junctionNumber
		self.chrom           = chrom
		self.strand          = strand
		self.junc_class      = junc_class
		self.donor           = donor
		self.acceptor        = acceptor
		self.diff_start      = diff_start
		self.diff_end        = diff_end
		self.genomCoordStart = genomCoordStart
		self.genomCoordEnd   = genomCoordEnd
		self.transcriptCoord = transcriptCoord
		self.spliceSite      = spliceSite
		self.canonical       = canonical
		self.coverage        = []
		self.indel           = indel


	def totalCov(self):
		if "NA" in self.coverage:
			totalCov = "NA"
		else:
			totalCov = sum([int(i) for i in self.coverage])
		return(totalCov)


	def coverage_files(self):
		if len(self.coverage)==0:
			return("")
		else:
			return(("\t").join(self.coverage))


	def biteDef(self):
		#if self.diff_donor != 0 and self.diff_acceptor != 0 and self.diff_donor%3 != 0 and self.diff_donor%3 != 0:
		if self.diff_start not in  [0, "NA"] and self.diff_end  not in  [0, "NA"]:
			return ("TRUE")
		else:
			return ("FALSE")


	def __str__(self):
		return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.transcript, self.junctionNumber, self.chrom, self.strand, self.genomCoordStart, self.genomCoordEnd, str(self.transcriptCoord), self.junc_class, self.donor, self.acceptor, self.diff_start, self.diff_end, self.biteDef(), self.spliceSite, self.canonical, self.indel, self.totalCov(), self.coverage_files())


def ReverseComplement(seq):
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	seq = seq.upper()
	RCseq = ""
	for base in seq:
		if base not in basecomplement:
			print("Error: NOT a DNA sequence")
			return None
			break
		RCseq = basecomplement[base]+RCseq
	return (RCseq)


def gtf_parsing(gtf_file, ref, refseqChrom={}):
	transcripts_dicc = {}

	if ref == "ensembl":
		gtf = openfile(gtf_file, "rt")

		for line in gtf:
			if line[0] != "#":
				strand = line.split("\t")[6]
				chrom = line.split("\t")[0]
				# Filter transcripts that are not present in the main chromosomes (eg: mitochondrial)
				if chrom in refseqChrom.values():
					if line.split("\t")[2]=="transcript":
						transcript_id = line.split("transcript_id")[1].split("\"")[1].strip()

						if transcript_id not in transcripts_dicc:
							if len(line.split("gene_name")) > 1:
								gene_name = line.split("gene_name")[1].split("\"")[1].strip()
							else:
								gene_name = line.split("gene_id")[1].split("\"")[1].strip()
							gene_id = line.split("gene_id")[1].split("\"")[1].strip()					
							transcripts_dicc[transcript_id] = myTranscripts(transcript_id, chrom, strand, 0, 0, gene_name, gene_id, transcript_id, ref)
						
						if ref == "ensembl":
							biotype = ""
							try:
								biotype = line.split("transcript_biotype")[1].split("\"")[1].strip()
							except: 
								pass
							transcripts_dicc[transcript_id].primary_class = biotype


					if line.split("\t")[2]=="exon":
						transcript_id = line.split("transcript_id")[1].split("\"")[1].strip()
						if transcript_id not in transcripts_dicc:
							if len(line.split("gene_name")) > 1:
								gene_name = line.split("gene_name")[1].split("\"")[1].strip()
							else:
								gene_name = line.split("gene_id")[1].split("\"")[1].strip()

							gene_id = line.split("gene_id")[1].split("\"")[1].strip()
							transcripts_dicc[transcript_id] = myTranscripts(transcript_id, chrom, strand, 0, 0, gene_name, gene_id, transcript_id, ref)
						
						exon_e = int(line.split("\t")[4])
						exon_s = int(line.split("\t")[3])
						exon_length = exon_e - exon_s + 1
						transcripts_dicc[transcript_id].length += exon_length
						transcripts_dicc[transcript_id].exons += 1
						transcripts_dicc[transcript_id].exons_coord.append([exon_s, exon_e])

					if line.split("\t")[2]=="CDS":
						transcript_id = line.split("transcript_id")[1].split("\"")[1].strip()
						if transcript_id not in transcripts_dicc:
							gene_id = line.split("gene_id")[1].split("\"")[1].strip()

							if len(line.split("gene_name")) > 1:
								gene_name = line.split("gene_name")[1].split("\"")[1].strip()
							else:
								gene_name = line.split("gene_id")[1].split("\"")[1].strip()

							transcripts_dicc[transcript_id] = myTranscripts(transcript_id, chrom, strand, 0, 0, gene_name, gene_id, transcript_id, ref)
						protein_id = line.split("protein_id")[1].split("\"")[1].strip()

						transcripts_dicc[transcript_id].refProtein = protein_id
						transcripts_dicc[transcript_id].coding = "coding"

						CDS_e = int(line.split("\t")[4])
						CDS_s = int(line.split("\t")[3])
						transcripts_dicc[transcript_id].cds_coord.append([CDS_s, CDS_e])

						exon_e = int(line.split("\t")[4])
						exon_s = int(line.split("\t")[3])
						if transcripts_dicc[transcript_id].firstProtPos == "NA":
							transcripts_dicc[transcript_id].firstProtPos = exon_s
							transcripts_dicc[transcript_id].lastProtPos =  exon_e
						else:
							if exon_s < transcripts_dicc[transcript_id].firstProtPos:
								transcripts_dicc[transcript_id].firstProtPos = exon_s
							if exon_e > transcripts_dicc[transcript_id].lastProtPos:
								transcripts_dicc[transcript_id].lastProtPos = exon_e
				else:
					continue				

	if ref == "refseq":
		gtf = openfile(gtf_file, "rt")

		for line in gtf:
			if line[0] != "#":
				strand = line.split("\t")[6]
				chr = line.split("\t")[0]
				#  Filter transcripts that are not present in the main chromosomes (eg: mitochondria)
				if chr in refseqChrom: 
					chrom = refseqChrom[line.split("\t")[0]]

					if line.split("\t")[2]=="transcript":
						transcript_id = line.split("transcript_id")[1].split(";")[0].replace('"', '').strip()
					
						if transcript_id not in transcripts_dicc:
							gene_name = line.split("gene_id")[1].split(";")[0].replace('"', '').strip()
							gene_id = line.split('db_xref "GeneID:')[1].split(";")[0].replace('"', '').strip()
							transcripts_dicc[transcript_id] = myTranscripts(transcript_id, chrom, strand, 0, 0, gene_name, gene_id, transcript_id, ref)

					if line.split("\t")[2]=="exon":
						transcript_id = line.split("transcript_id")[1].split(";")[0].replace('"', '').strip()
						if transcript_id not in transcripts_dicc:
							gene_name = line.split("gene_id")[1].split(";")[0].replace('"', '').strip()
							gene_id = line.split('db_xref "GeneID:')[1].split(";")[0].replace('"', '').strip()
							transcripts_dicc[transcript_id] = myTranscripts(transcript_id, chrom, strand, 0, 0, gene_name, gene_id, transcript_id, ref)

						exon_e = int(line.split("\t")[4])
						exon_s = int(line.split("\t")[3])
						exon_length = exon_e - exon_s + 1
						transcripts_dicc[transcript_id].length += exon_length
						transcripts_dicc[transcript_id].exons += 1
						transcripts_dicc[transcript_id].exons_coord.append([exon_s, exon_e])

					if line.split("\t")[2]=="CDS" and "protein_id" in line.split("\t")[8]:
						transcript_id = line.split("transcript_id")[1].split(";")[0].replace('"', '').strip()
						protein_id = line.split("protein_id")[1].split(";")[0].replace('"', '').strip()
						
						if transcript_id not in transcripts_dicc:
							gene_id = line.split('db_xref "GeneID:')[1].split(";")[0].replace('"', '').strip()
							gene_name = line.split("gene_id")[1].split(";")[0].strip()
							transcripts_dicc[transcript_id] = myTranscripts(transcript_id, chrom, strand, 0, 0, gene_name, gene_id, transcript_id, ref)
						
						transcripts_dicc[transcript_id].refProtein = protein_id
						transcripts_dicc[transcript_id].coding = "coding"

						CDS_e = int(line.split("\t")[4])
						CDS_s = int(line.split("\t")[3])
						transcripts_dicc[transcript_id].cds_coord.append([CDS_s, CDS_e])

						exon_e = int(line.split("\t")[4])
						exon_s = int(line.split("\t")[3])
						if transcripts_dicc[transcript_id].firstProtPos == "NA":
							transcripts_dicc[transcript_id].firstProtPos = exon_s
							transcripts_dicc[transcript_id].lastProtPos =  exon_e
						else:
							if exon_s < transcripts_dicc[transcript_id].firstProtPos:
								transcripts_dicc[transcript_id].firstProtPos = exon_s
							if exon_e > transcripts_dicc[transcript_id].lastProtPos:
								transcripts_dicc[transcript_id].lastProtPos = exon_e
				else:
					continue

	return(transcripts_dicc)



def fasta_parser(fastaFile):
	try:
		fasta = openfile(fastaFile, "rt")
	except IOError:
		print ('ERROR: Unable to read %s file\n' % fastaFile)
		raise SystemExit(1)
	try:
		seqDicc = {}
		index = 0
		for line in fasta:
			if line.startswith(">"):
				if index > 0:
					seqDicc[name] = seq
				index+=1
				name = line[1:].split()[0].rstrip()
				seq = ''
			else:
				seq += line.rstrip()
		seqDicc[name] = seq

	except IOError:
		print ('File %s without fasta format' % fastaFile)

	return(seqDicc)


def main():
	parser = argparse.ArgumentParser(description="ReferenceSQANTI")
	parser.add_argument('--output_classification', required=True)
	parser.add_argument('--output_junctions', required=True)
	parser.add_argument('--output_nmd', required=True)

	parser.add_argument('--gtf_file', required=True)
	parser.add_argument('--reference_file', required=True)
	parser.add_argument('--database', required=True)
	parser.add_argument('--chr_ref', required=True, nargs='?', const=None)

	args = parser.parse_args()

	gtf_file = args.gtf_file
	genome = args.reference_file
	ref_type = args.database # refseq or ensembl

	output_file = args.output_classification
	output_junc_file = args.output_junctions
	output_NMD_file = args.output_nmd

	genome_dicc = fasta_parser(genome)
	refseqChrom = read_chr_ref_acc(args.chr_ref)

	transcript_dicc = gtf_parsing(gtf_file, ref_type, refseqChrom)


	output_NMD = open(output_NMD_file, "w")

	with open(output_file, "w") as output:

		output.write("isoform\tchrom\tstrand\tlength\texons\tstructural_category\tassociated_gene\tassociated_transcript\tref_length\tref_exons\tdiff_to_TSS\tdiff_to_TTS\tsubcategory\tRTS_stage\tall_canonical\tmin_cov\tmin_cov_pos\tsd_cov\tFL\tn_indels\tn_indels_junc\tbite\tiso_exp\tgene_exp\tratio_exp\tFSM_class\tcoding\tORF_length\tCDS_length\tCDS_start\tCDS_end\trefProt\n")

		for i in transcript_dicc:
			if transcript_dicc[i].chrom in genome_dicc and transcript_dicc[i].exons_coord != []:
				output.write(i+"\t")
				output_NMD.write(i+"\t")
				output_NMD.write(transcript_dicc[i].NMD_calculation())
				output.write(str(transcript_dicc[i]))
				output.write("\n")
				output_NMD.write("\n")

	sites =  ["ATAC", "GCAG", "GTAG"]

	output_junc = open(output_junc_file, "w")
	output_junc.write("isoform\tjunction_number\tchrom\tstrand\tgenomic_start_coord\tgenomic_end_coord\ttranscript_coord\tjunction_category\tstart_site_category\tend_site_category\tdiff_to_Ref_start_site\tdiff_to_Ref_end_site\tbite\tsplice_site\tcanonical\tindel_near_junct\ttotal_coverage\n")

	for transcript in transcript_dicc:
		chrom = transcript_dicc[transcript].chrom
		strand = transcript_dicc[transcript].strand
		if len(transcript_dicc[transcript].exons_coord) > 1 and chrom in genome_dicc and transcript_dicc[transcript].exons_coord != []:
			exons_s = []
			exons_e = []
			for exon in transcript_dicc[transcript].exons_coord:
				exons_s.append(exon[0]-1)
				exons_e.append(exon[1]+1)
			junctions_e	= sorted(exons_s)[1:]
			junctions_s	= sorted(exons_e)[:-1]

			juncN = 0

			for i in range(len(junctions_s)):
				site = genome_dicc[chrom][(junctions_s[i]-1):(junctions_s[i]+1)] + genome_dicc[chrom][(junctions_e[i]-2):junctions_e[i]]
				if strand =="-":
					site = ReverseComplement(site)
				else:
					site = site.upper()
					if any([base not in ["A","T","C","G","N"] for base in site]):
						print("Error: NOT a DNA sequence")
						return None
						break

				if site in sites:
					canonical = "canonical"
				else:
					canonical= "non_canonical"

				juncN += 1

				output_junc.write("%s\tjunction_%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (transcript, str(juncN), chrom, strand, str(junctions_s[i]), str(junctions_e[i]), ".", "known", "known", "known", "0", "0", "NA", site, canonical, "NA", "NA"))

	output_junc.close()

if __name__ == "__main__":
	main()
