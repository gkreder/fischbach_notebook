##Modules
from string import ascii_letters
from string import ascii_letters
import pickle as pkl

def get_sequence(fasta):

		"""get the description and trimmed dna sequence"""
		in_file = open(fasta, 'r')
		content = in_file.readlines()
		in_file.close()
		content2 = []

		for i in content:
			if i != "":
				content2.append(i)
		
		content = content2
		while content[0] == "" or content[0] == "\n":
			content = content[1:]

		header = content[0]
		content = content[1:]
		content = [x.rstrip() for x in content]
		seq = "".join(content)

		if ">" not in header or ">" in seq:
			print(sys.stderr, "FASTA file not properly formatted; should be single sequence starting with '>' and sequence name.")
			# print >> sys.stderr, "FASTA file not properly formatted; should be single sequence starting with '>' and sequence name."
			# logfile.write("FASTA file not properly formatted; should started with '>' and sequence name on first line.\n")
			# logfile.close()
			sys.exit(1)

		return seq


def complement(seq):  
	complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}  
	complseq = []
	for base in seq:
		if base in complement.keys():
			complbase = complement[str(base)]
			complseq.append(complbase)
		else:
			complbase = 'n'
			complseq.append(complbase)

	return complseq 

def reverse_complement(seq):  
	seq = list(seq)  
	seq.reverse()  
	revcompl = complement(seq)
	revcomplstr = str()

	for i in revcompl:
		revcomplstr = revcomplstr + str(i)

	return  revcomplstr

def fastaseqlengths(proteins):
	names = proteins[0]
	seqs = proteins[1]
	seqlengths = {}
	a = 0

	for i in names:
		seq = seqs[a]
		seqlength = len(seq)
		seqlengths[i] = seqlength
		a += 1

	return seqlengths



# Function that reads the fasta file into a dictionary
def fastadict(fasta):
	file = open(fasta,"r")
	filetext = file.read()
	filetext = filetext.replace("\r","\n")
	filetext = filetext.strip()
	#Replaces all spaces with "_" to avoid problems
	filetext = filetext.replace(' ','_')
	filetext = filetext.split()
	dictseq = {}

	for a in filetext:
		if ">" in a[0]:
			f = str()
			d = a[1:68]

		else:
			e = a
			f += e
			dictseq[d] = f

	return dictseq



# Function that extracts all sequence names from the fasta dictionary
def lnames(fastadict):
	items = fastadict.items()
	items.sort()
	return [names for names, seqs in items]

# Function that extracts all sequences from the fasta dictionary
def lseqs(fastadict):
	items = fastadict.items()
	items.sort()
	return [seqs for names, seqs in items]



def parsegenes(genes):
	genedict = {}
	genelist = []
	joinlist = []
	joindict = {}
	accessiondict = {}
	genenr = 0

	for i in genes:
		i = i.split("     gene            ")[0]
		join = "no"
		genenr += 1

		#Find gene location info for each gene
		if "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
			location = i.split("\n")[0]
		elif "complement" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
			location = i.split("   /")[0]
			while ")" not in location.replace(" ","")[-3:]:
				locationlist = location.split("\n")
				locationlist = locationlist[:-1]
				location = ""
				for i in locationlist:
					location = location + "i"
			location = location.replace("\n","")
			location = location.replace(" ","")
		elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] == ")":
			location = i.split("\n")[0]
		elif "join" in i.split("\n")[0].lower() and i.split("\n")[0][-1] != ")":
			location = i.split("/")[0]
			while ")" not in location.replace(" ","")[-3:]:
				locationlist = location.split("\n")
				locationlist = locationlist[:-1]
				location = ""
				for i in locationlist:
					location = location + "i"
			location = location.replace("\n","")
			location = location.replace(" ","")
		else:
			location = i.split("\n")[0]

		#location info found in embl file, now extract start and end positions
		if "complement" in location.lower():
			location = location.lower()
			location = location.split("complement(")[1][:-1]
			if "join(" in location.lower():
				join = "yes"
				location = location.lower()
				location2 = location.split("join(")[1][:-1]
				start = location2.split(",")[0]
				start = start.split("..")[0]
				start = start.replace("<","")
				end = location2.split(",")[-1]

				if ".." in end:
					end = end.split("..")[1]
				end = end.replace(">","")
				joinedparts = location2.split(",")
				joinedparts2 = []

				for j in joinedparts:
					newjoinedpart = j.replace("<","")
					newjoinedpart = newjoinedpart.replace(">","")
					joinedparts2.append(newjoinedpart)

			else:
				start = location.split("..")[0]
				start = start.replace("<","")
				end = location.split("..")[1]
				end = end.replace(">","")

			strand = "-"

		else:
			if "join(" in location.lower():
				join = "yes"
				location = location.lower()
				location2 = location.split("join(")[1][:-1]
				start = location2.split(",")[0]
				start = start.split("..")[0]
				start = start.replace("<","")
				end = location2.split(",")[-1]
				if ".." in end:
					end = end.split("..")[1]
				end = end.replace(">","")
				joinedparts = location2.split(",")
				joinedparts2 = []

				for j in joinedparts:
					newjoinedpart = j.replace("<","")
					newjoinedpart = newjoinedpart.replace(">","")
					joinedparts2.append(newjoinedpart)

			else:
				start = location.split("..")[0]
				start = start.replace("<","")
				end = location.split("..")[1]
				end = end.replace(">","")

			strand = "+"

		if int(start) > int(end):
			start2 = end
			end2 = start
			start = start2
			end = end2

		#Correct for alternative codon start positions
		if "codon_start=" in i.lower():
			codonstart = i.lower().split("codon_start=")[1][0]
			if strand == "+":
				start = str(int(start) +  (int(codonstart) - 1))
			elif strand == "-":
				end = str(int(end) - (int(codonstart) - 1))

		#Find gene name for each gene, preferably locus_tag, than gene, than protein_ID
		a = 0
		b = 0
		genename = ""
		nrlines = len(i.split("\n"))
		while b == 0:
			line = i.split("\n")[a]
			if "protein_id=" in line:
				genename = (line.split("protein_id=")[1][1:-1]).replace(" ","_")
				genename = genename.replace("\\","_")
				genename = genename.replace("/","_")
				b += 1

			elif "protein_id=" in line.lower():
				genename = (line.lower().split("protein_id=")[1][1:-1]).replace(" ","_")
				genename = genename.replace("\\","_")
				genename = genename.replace("/","_")
				b += 1

			elif a == (nrlines - 1):
				genename = ""
				b += 1

			else:
				a += 1

		if len(genename) > 1:
			accnr = genename

		else:
			accnr = "no_accession_number_found"

		a = 0
		b = 0
		nrlines = len(i.split("\n"))
		while b == 0:
			line = i.split("\n")[a]
			if "gene=" in line:
				genename = (line.split("gene=")[1][1:-1]).replace(" ","_")
				genename = genename.replace("\\","_")
				genename = genename.replace("/","_")
				b += 1

			elif "gene=" in line.lower():
				genename = (line.lower().split("gene=")[1][1:-1]).replace(" ","_")
				genename = genename.replace("\\","_")
				genename = genename.replace("/","_")
				b += 1

			elif a == (nrlines - 1):
				b += 1

			else:
				a += 1

		a = 0
		b = 0
		nrlines = len(i.split("\n"))
		while b == 0:
			line = i.split("\n")[a]
			if "locus_tag=" in line:
				genename = (line.split("locus_tag=")[1][1:-1]).replace(" ","_")
				genename = genename.replace("\\","_")
				genename = genename.replace("/","_")
				b += 1
			elif "locus_tag=" in line.lower():
				genename = (line.lower().split("locus_tag=")[1][1:-1]).replace(" ","_")
				genename = genename.replace("\\","_")
				genename = genename.replace("/","_")
				b += 1
			elif a == (nrlines - 1):
				if genename == "":
					genename = "prot_ID_" + str(genenr)
				b += 1
			else:
				a += 1

		#Find sequence for each gene

		a = 0                                             ###Not all gbks contain protein sequences as translations, therefore sequences from gene clusters are now extracted from the database at a later stage if sequence is not in gbk
		b = 0
		sequence = ""
		while b < 2:
			line = i.split("\n")[a]
			if "translation=" in line:
				sequence = line.split("translation=")[1][1:]
				b += 1
				a += 1
				if line.count('"') > 1:
					sequence = line.split("translation=")[1][1:-1]
					b = 2
			elif "translation=" in line.lower():
				sequence = line.lower().split("translation=")[1][1:]
				b += 1
				a += 1
				if line.count('"') > 1:
					sequence = line.lower().split("translation=")[1][1:-1]
					b = 2
			elif a == (nrlines - 2) or a == (nrlines - 1):
				sequence = ""
				b = 2
			elif b == 1:
				if '"' in line:
					seqline = line.replace(" ","")
					seqline = seqline.split('"')[0]
					sequence = sequence + seqline
					b += 1
				else:
					seqline = line.replace(" ","")
					sequence = sequence + seqline
				a += 1
			else:
				a += 1

		sequence = sequence.upper()

		#Quality-check sequence
		forbiddencharacters = ["'",'"','=',';',':','[',']','>','<','|','\\',"/",'*','-','_','.',',','?',')','(','^','#','!','`','~','+','{','}','@','$','%','&']
		for z in forbiddencharacters:
			if z in sequence:
				sequence = ""

		#Find annotation for each gene
		a = 0
		b = 0
		while b == 0:
			line = i.split("\n")[a]
			if "product=" in line:
				annotation = line.split("product=")[1][1:]
				annotation = annotation.replace(" ","_")
				if annotation[-1] == '"':
					annotation = annotation[:-1]
				b += 1
			elif "product=" in line.lower():
				annotation = line.lower().split("product=")[1][1:]
				annotation = annotation.replace(" ","_")
				if annotation[-1] == '"':
					annotation = annotation[:-1]
				b += 1
			elif a == (nrlines - 1):
				annotation = "not_annotated"
				b += 1
			else:
				a += 1
		accessiondict[genename] = accnr
		if join == "yes":
			joinlist.append(genename)
			joindict[genename] = joinedparts2
		
		#Save data to dictionary
		if len(genename) > 1:
			genedict[genename] = [start,end,strand,annotation,sequence]
		genelist.append(genename)
	return [genelist, genedict, joinlist, joindict, accessiondict]

def cleandnaseq(dnaseq):
	dnaseq = dnaseq.replace(" ","")
	dnaseq = dnaseq.replace("\t","")
	dnaseq = dnaseq.replace("\n","")
	dnaseq = dnaseq.replace("0","")
	dnaseq = dnaseq.replace("1","")
	dnaseq = dnaseq.replace("2","")
	dnaseq = dnaseq.replace("3","")
	dnaseq = dnaseq.replace("4","")
	dnaseq = dnaseq.replace("5","")
	dnaseq = dnaseq.replace("6","")
	dnaseq = dnaseq.replace("7","")
	dnaseq = dnaseq.replace("8","")
	dnaseq = dnaseq.replace("9","")
	dnaseq = dnaseq.replace("/","")
	return dnaseq

def extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict):
	names = []
	seqs = []
	for i in genelist:
		genename = i
		#If suitable translation found in gbk, use that
		if len(genedict[i][4]) > 5:
			protseq = genedict[i][4]
			i = genedict[i]
		#If no suitable translation found in gbk, extract from DNA sequence
		else:
			i = genedict[i]
			y = int(i[0])
			z = int(i[1])
			if i[2] == "+":
				if genename in joinlist:
					geneseq = ""
					for j in joindict[genename]:
						partstart = int(j.split("..")[0])
						if ".." in j:
							partend = int(j.split("..")[1])
						else:
							partend = int(j)
						geneseqpart = dnaseq[(partstart - 1):partend]
						geneseq = geneseq + geneseqpart
				else:
					geneseq = dnaseq[(y - 1):z]
				protseq = translate(geneseq)
			elif i[2] == "-":
				if genename in joinlist:
					geneseq = ""
					joinlistrev = joindict[genename]
					joinlistrev.reverse()
					for j in joinlistrev:
						partstart = int(j.split("..")[0])
						if ".." in j:
							partend = int(j.split("..")[1])
						else:
							partend = int(j)
						geneseqpart = rc_dnaseq[(len(rc_dnaseq) - partend):(len(rc_dnaseq) - partstart + 1)]
						geneseq = geneseq + geneseqpart
				else:
					geneseq = rc_dnaseq[(len(rc_dnaseq) - z):(len(rc_dnaseq) - y + 1)]
				protseq = translate(geneseq)
		name = "input" + "|" + "c1" + "|" + i[0] + "-" + i[1] + "|" + i[2] + "|" + genename + "|" + i[3]
		seqs.append(protseq)
		names.append(name)
	proteins = [names,seqs,genelist,genedict,accessiondict]
	return proteins



def gbk2proteins(gbkfile):
	file = open(gbkfile,"r")
	filetext_all = file.read()
	filetext_all = filetext_all.replace("\r","\n")

	all_text_blocks = filetext_all.split("//")
	all_results = []


	if "     CDS             " not in filetext_all or "\nORIGIN" not in filetext_all:
			print(sys.stderr, "Exit: GBK file not properly formatted, no sequence found")
			# print >> sys.stderr, "Exit: GBK file not properly formatted, no sequence found"
			# logfile.write("Exit: GBK file not properly formatted, no sequence found or no CDS annotation found.\n")
			# logfile.close()
			sys.exit(1)


	for filetext in all_text_blocks:		
		if "     CDS             " not in filetext or "\nORIGIN" not in filetext:
			continue


		cdspart = filetext.split("\nORIGIN")[0]
		#Extract DNA sequence and calculate reverse complement of it
		dnaseq = filetext.split("\nORIGIN")[1]
		dnaseq = cleandnaseq(dnaseq)
		dnaseqlength = len(dnaseq)
		rc_dnaseq = reverse_complement(dnaseq)

		#Extract genes
		genes = cdspart.split("     CDS             ")
		genes = genes[1:]

		# print('::' in genes[1])
		for gene_index, gene in enumerate(genes):
			if '::' in gene:
				genes.pop(gene_index)

		genesdetails = parsegenes(genes)
		genelist = genesdetails[0]
		genedict = genesdetails[1]
		joinlist = genesdetails[2]
		joindict = genesdetails[3]
		accessiondict = genesdetails[4]

		#Locate all genes on DNA sequence and translate to protein sequence
		proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict)
		textlines = filetext.split("\n//")[0]
		textlines = textlines.split("\n")
		accession = ""

		for i in textlines:
			if accession == "":
				if "LOCUS       " in i:
					j = i.split("LOCUS       ")[1]
					accession = j.split(" ")[0]
					if len(accession) < 4:
						accession = ""

		#Test if accession number is probably real GenBank/RefSeq acc nr
		numbers = range(0,10)
		# from string import ascii_letters
		letters = []
		for i in ascii_letters:
			letters.append(i)

		nrnumbers = 0
		nrletters = 0
		for i in accession:
			if i in letters:
				nrletters += 1
			try:
				j = int(i)
				if j in numbers:
					nrnumbers += 1
			except:
				pass
		if nrnumbers < 3 or nrletters < 1:
			accession = ""

		all_results.append([proteins,accession,dnaseqlength])

	return(all_results)
	# return [proteins,accession,dnaseqlength]



def embl2proteins(emblfile,sequence):
	file = open(emblfile,"r")
	filetext = file.read()
	filetext = filetext.replace("\r","\n")
	file.close()

	if "FT   CDS " not in filetext or ("\nSQ" not in filetext and len(sequence) < 1):
		# logfile.write("Exit: EMBL file not properly formatted, no sequence found or no CDS annotation found.\n")
		print(sys.stderr, "Exit: EMBL file not properly formatted, no sequence found or no CDS annotation found.\n")
		# print >> sys.stderr, "Exit: EMBL file not properly formatted, no sequence found or no CDS annotation found.\n"
		# logfile.close()
		sys.exit(1)

	if len(sequence) < 1:
		cdspart = filetext.split("\nSQ  ")[0]
		#Extract DNA sequence and calculate reverse complement of it
		seqpart = filetext.split("\nSQ  ")[1]
		seqlines = seqpart.split("\n")[1:]
		dnaseq = ""

		for i in seqlines:
			dnaseq = dnaseq + i
		dnaseq = cleandnaseq(dnaseq)

	else:
		dnaseq = sequence
		cdspart = filetext
	dnaseqlength = len(dnaseq)
	rc_dnaseq = reverse_complement(dnaseq)

	#Extract genes
	genes = cdspart.split("FT   CDS             ")
	genes = genes[1:]
	genesdetails = parsegenes(genes)
	genelist = genesdetails[0]
	genedict = genesdetails[1]
	joinlist = genesdetails[2]
	joindict = genesdetails[3]
	accessiondict = genesdetails[4]

	#Locate all genes on DNA sequence and translate to protein sequence
	proteins = extractprotfasta(genelist,genedict,dnaseq,rc_dnaseq,joinlist,joindict,accessiondict)
	textlines = filetext.split("SQ   ")[0]
	textlines = textlines.split("\n")
	accession = ""

	for i in textlines:
		if accession == "":
			if "AC   " in i:
				j = i.split("AC   ")[1]
				j = j.replace(" ","")
				accession = j.split(";")[0]
				if len(accession) < 4:
					accession = ""

	#Test if accession number is probably real GenBank/RefSeq acc nr
	numbers = range(0,10)
	# from string import ascii_letters
	letters = []
	for i in ascii_letters:
		letters.append(i)

	nrnumbers = 0
	nrletters = 0

	for i in accession:
		if i in letters:
			nrletters += 1
		try:
			j = int(i)
			if j in numbers:
				nrnumbers += 1
		except:
			pass

	if nrnumbers < 3 or nrletters < 1:
		accession = ""

	return [proteins,accession,dnaseqlength]



def translate(sequence):

	#Translation table standard genetic code; according to http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
	transldict = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
								 'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C', 
								 'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*', 
								 'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W', 
								 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R', 
								 'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 
								 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 
								 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 
								 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 
								 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 
								 'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 
								 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 
								 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 
								 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 
								 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 
								 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G',
								 'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
								 'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
								 'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
								 'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
								 'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
								 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
								 'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
								 'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
								 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
								 'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
								 'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
								 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
								 'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
								 'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
								 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
								 'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'}

	triplets = []
	triplet = ""
	a = 0

	for i in sequence:
		if a < 2:
			a += 1
			triplet = triplet + i
		elif a == 2:
			triplet = triplet + i
			triplets.append(triplet)
			triplet = ""
			a = 0
	protseq = ""
	aanr = 0

	for i in triplets:
		aanr += 1
		if aanr == 1:
			protseq = protseq + "M"
		else:
			if "n" in i or "N" in i or i not in transldict.keys():
				protseq = protseq + "X"
			else:
				protseq = protseq + transldict[i]

	if  len(protseq) > 0 and protseq[-1] == "*":
		protseq = protseq[:-1]
	return protseq



def writefasta(names,seqs,file):
	e = 0
	f = len(names) - 1

	try:
		out_file = open(file,"a")
		while e <= f:
			out_file.write(">")
			out_file.write(names[e])
			out_file.write("\n")
			out_file.write(seqs[e])
			out_file.write("\n")
			e += 1
		out_file.close()

	except(IOError,OSError,NotImplementedError):
		print(sys.stderr, "FASTA file not created.")
		# print >> sys.stderr, "FASTA file not created."
		# logfile.write("FASTA file not created.\n")

	
def hmmlengths(hmmfile):
	# print('-------------------------')
	# print(hmmfile)
	# print('-------------------------')
	hmmlengthsdict = {}
	file = open(hmmfile,"r")
	filetext = file.read()
	filetext = filetext.replace("\r","\n")
	hmms = filetext.split("//")[:-1]

	for i in hmms:
		namepart = i.split("NAME  ")[1]
		name = namepart.split("\n")[0]
		lengthpart = i.split("LENG  ")[1]
		#print(lengthline)
		#tabs = lengthline.split(" ")
		#tabs2 = []
		#for j in tabs:
		#  if j != "":
		#    tabs2.append(j)
		#print(tabs2)
		length = lengthpart.split("\n")[0]
		hmmlengthsdict[name] = int(length)
	return hmmlengthsdict

def hmmscanparse(hmmscanoutputfile,hmmlengthsdict):
	domaindict = {}
	file = open(hmmscanoutputfile,"r")
	filetext = file.read()
	filetext = filetext.replace("\r","\n")
	outputs = filetext.split("Query:       ")[1:]

	for i in outputs:
		protname = i.split("\n")[0]
		protname = protname.split(" ")[0]
		domainresults = i.split("Domain annotation for each model:\n")[1]
		domainresults = domainresults.split("\n\nInternal pipeline statistics summary:")[0]
		domains = domainresults.split(">> ")
		domainlist = []

		#Find all domains
		for i in domains:
			domainname = i.split("\n")[0]
			domainname = domainname.split(" ")[0]
			domainresults = i.split("\n")[3:-2]

			for i in domainresults:
				tabs = i.split(" ")
				tabs2 = []
				for i in tabs:
					if i != "":
						tabs2.append(i)
				tabs = tabs2
				start = int(tabs[12])
				end = int(tabs[13])
				evalue = tabs[5]
				score = float(tabs[2])
				domainlist.append([domainname,start,end,evalue,score])
		# domainlist.sort(sortonsecondvalueoflist)
		domainlist.sort(key = lambda tup: tup[1])

		#Purify domain list to remove overlapping domains, only keeping those with the highest scores
		if len(domainlist) > 1:
			domainlist2 = [domainlist[0]]
			for i in domainlist[1:]:
				maxoverlap = 20
				if i[1] < (domainlist2[-1][2] - maxoverlap):
					if i[4] < domainlist2[-1][4]:
						pass
					elif i[4] > domainlist2[-1][4]:
						del domainlist2[-1]
						domainlist2.append(i)

				else:
					domainlist2.append(i)

			domainlist = domainlist2

		#Merge domain fragments which are really one domain
		if len(domainlist) > 1:
			domainlist2 = [domainlist[0]]
			for i in domainlist[1:]:
				alilength1 = int(domainlist2[-1][2]) - int(domainlist2[-1][1])
				alilength2 = int(i[2]) - int(i[1])
				domainlength = hmmlengthsdict[i[0]]

				if i[0] == domainlist2[-1][0] and (alilength1 < (0.75 * domainlength) or alilength2 < (0.75 * domainlength)) and (alilength1 + alilength2) < (1.5 * domainlength):
					name = i[0]
					start = domainlist2[-1][1]
					end = i[2]
					evalue = str(float(domainlist2[-1][3]) * float(i[3]))
					score = str(float(domainlist2[-1][4]) + float(i[4]))
					del domainlist2[-1]
					domainlist2.append([name,start,end,evalue,score])
				else:
					domainlist2.append(i)
			domainlist = domainlist2

		#Remove incomplete domains (covering less than 60% of total domain hmm length)
		if len(domainlist) > 1:
			domainlist2 = []
			for i in domainlist:
				alilength = int(i[2]) - int(i[1])
				domainlength = hmmlengthsdict[i[0]]
				if alilength > (0.6 * domainlength):
					domainlist2.append(i)
			domainlist = domainlist2

		#Save domainlist to domaindict
		domaindict[protname] = domainlist
	return domaindict



def sortonsecondvalueoflist(first,second):
	if int(first[1]) > int(second[1]):
		value = 1
	if int(first[1]) < int(second[1]):
		value = -1
	if int(first[1]) == int(second[1]):
		value = 0
	return value



##Core script
import sys
import os
from datetime import datetime

ARGS = sys.argv
if __name__ == '__main__':

	print('------------------------------------------------------------------------------------------------------------------------------')
	print('Script Start')
	print('------------------------------------------------------------------------------------------------------------------------------')
	#Takes folder named 'knownclusters' and finds all Genbank files in this folder
	topdir = os.getcwd()
	os.chdir('./all_gb')
	
	# gbkfiles = [filename for filename in os.listdir(".") if ".gbff" in filename]

	gbkfiles = [ARGS[1]]		

	print(str(datetime.now()) + '   ----    Looking for file -- ' + ARGS[1] + ' -- in directory ' + os.getcwd())
	sys.stdout.flush()

	#Read gbk file
	for iteration_number, gbkfile in enumerate(gbkfiles):
		print('-----------------------------------------------------------------------------------------------------')
		print(gbkfile + ' (' + str(iteration_number + 1) + ' of ' + str(len(gbkfiles)) + ') started at ' + str(datetime.now()))
		print('-----------------------------------------------------------------------------------------------------')
		sys.stdout.flush()
		
		out_path = topdir + '/all_pfam/'
		outfile1 = open(out_path + gbkfile.split(".gbff")[0] + ".pfs","w+")
		outfile2 = open(out_path + gbkfile.split(".gbff")[0] + ".pfd","w+")
		outfile1.close()
		outfile2.close()

		all_gbk_results = gbk2proteins(gbkfile)
		proteins_fasta_name = topdir + '/hmm_aux_files/' + gbkfile.split(".gbff")[0] + '_proteins.fasta'
		hmm_output_name = topdir + '/hmm_aux_files/' + gbkfile.split(".gbff")[0] + '_hmm_output.txt'
		tbl_output_name = topdir + '/hmm_aux_files/' + gbkfile.split(".gbff")[0] + '_hmm.txt'
		f = open(proteins_fasta_name, 'w+')
		f.close()

		accessiondict = {}
		seqdict = {}
		fullnamedict = {}
		strandsdict = {}
		z = 0
		locustags = []

		for proteins in all_gbk_results:
			# proteins = gbk2proteins(gbkfile)
			genomic_accnr = proteins[1]
			dnaseqlength = proteins[2]
			proteins = proteins[0]

			
			writefasta(proteins[0],proteins[1],proteins_fasta_name)
		

			# accessiondict = proteins[4]
			for key in proteins[4]:
				val = proteins[4][key]
				accessiondict[key] = val

			# seqdict = {}
			# fullnamedict = {}
			# strandsdict = {}
			z = 0
			# locustags = []

			for i in proteins[0]:
				name = i.split("|")[4]
				locustags.append(name)
				seq = proteins[1][z]
				seqdict[name] = seq
				strand = i.split("|")[3]
				strandsdict[name] = strand
				fullnamedict[name] = i
				z += 1


			#Output proteins fasta


		#Run PFAM search
		# hmmsearch = "hmmscan --cpu 4 -o hmm_output.txt --noali --cut_tc --tblout hmm.txt Pfam-A.hmm proteins.fasta"
		# hmmsearch = "/netapp/home/gkreder/Fischbach/hmmer/binaries/hmmscan --cpu 4 -o hmm_output.txt --noali --cut_tc --tblout hmm.txt /netapp/home/gkreder/Fischbach/hmmer/binaries/Pfam-A.hmm " + proteins_fasta_name
		hmmsearch = "/netapp/home/gkreder/Fischbach/hmmer/binaries/hmmscan --cpu 4 -o " + hmm_output_name +" --noali --cut_tc --tblout " + tbl_output_name + " /netapp/home/gkreder/Fischbach/hmmer/binaries/Pfam-A.hmm " + proteins_fasta_name
		print(str(datetime.now()) + '   --    Running Command hmmsearch with command line: ')
		print('\t\t\t' + hmmsearch)
		sys.stdout.flush()
		os.system(hmmsearch)

		# hmmlengthsdict = hmmlengths("PFAM-A.hmm")
		# hmmlengthsdict = hmmlengths("/netapp/home/gkreder/Fischbach/hmmer/binaries/Pfam-A.hmm")
		# domaindict = hmmscanparse("/netapp/home/gkreder/Fischbach/hmm_output.txt",hmmlengthsdict)
		print(str(datetime.now()) + '   --    Running hmmlengths()')
		sys.stdout.flush()
		hmmlengthsdict = hmmlengths("/netapp/home/gkreder/Fischbach/hmmer/binaries/Pfam-A.hmm")
		print(str(datetime.now()) + '   --    Running hmmparse()')
		sys.stdout.flush()
		domaindict = hmmscanparse(hmm_output_name, hmmlengthsdict)

		outfile1 = open(out_path + gbkfile.split(".gbff")[0] + ".pfs","a")
		outfile2 = open(out_path + gbkfile.split(".gbff")[0] + ".pfd","a")


		print(str(datetime.now()) + '   --    Creating Pfam ID dictionary')
		sys.stdout.flush()
		#Create PFAM ID dictionary
		pfamid_dict = {}
		# pfamdat = open(topdir+"/Pfam-A.hmm","r")
		# pfamdat = open("/netapp/home/gkreder/Fischbach/hmmer/binaries/Pfam-A.hmm.dat","r")
		pfamdat = open("/netapp/home/gkreder/Fischbach/hmmer/binaries/Pfam-A.hmm.dat","r")
		pfamdat = pfamdat.read()
		pfamdatparts = pfamdat.split("\n//\n")

		for i in pfamdatparts:
			lines = i.split("\n")
			for i in lines:
				if "#=GF ID" in i:
					pfamid = i.split("   ")[1]
					for j in lines:
						if "#=GF AC" in j:
							pfamacc = j.split("   ")[1]
							pfamid_dict[pfamid] = pfamacc.split(".")[0]

		# pkl.dump(pfamid_dict, open('pfamid_dict.pkl', 'wb'))

		#Parse locustags, positions, PFAM domains, strands, output to txt file; output string of PFAM domains to separate text file
		print(str(datetime.now()) + '   --    Parsing Locustags, Positions, and PFAM domains')
		sys.stdout.flush()
		for i in locustags:
			fullname = fullnamedict[i]
			print('\t\t' + str(i))
			print('\t\t\t\t' + fullname)
			sys.stdout.flush()

			# if domaindict.has_key(fullname):
			# print(fullname)
			if fullname in domaindict:
				print('\t\t\t\t\t Found in domaindict and writing to outfiles')
				sys.stdout.flush()
				for j in domaindict[fullname]:
					start = fullname.split("|")[2].split("-")[0]
					end = fullname.split("|")[2].split("-")[1]
					outfile1.write(pfamid_dict[j[0]] + " ")
					outfile2.write(gbkfile + "\t" + fullname.split("|")[4] + "\t" + start + "\t" + end + "\t" + fullname.split("|")[3] + "\t" + pfamid_dict[j[0]] + "\t" + j[0] + "\n")

		outfile1.close()
		outfile2.close()
		print(str(datetime.now()) + '   --    Wrote outfiles at locations')
		print('\t\t' + str(outfile1))
		print('\t\t' + str(outfile2))

	print('Finished at ' + str(datetime.now()))
	sys.stdout.flush()

