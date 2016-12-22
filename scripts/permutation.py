from Bio import AlignIO
from Bio.Nexus import Nexus
import numpy as np
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from tqdm import tqdm
from io import StringIO
import os
from pathlib2 import Path


def process_permutations(dataset_path, n_permut):
	for index in range(n_permut):
		name = "%s-permut-%s" % (os.path.basename(dataset_path).rstrip(".nex"), index)
		export_suffle_alns(dataset_path, name)


def export_suffle_alns(dataset_path, outfilename):
	
	dat = Nexus.Nexus()
	dat.read(dataset_path)
	aln = AlignIO.read(open(dataset_path), "nexus")

	aln_nex = []

	for key, value in dat.charsets.items():
		name  = key
		sites = value

		uce_aln = aln[:, min(sites) : max(sites) + 1]
		uce_perm = permut_aln_per_site(uce_aln)

		#nexus
		out_nex = StringIO()
		AlignIO.write(uce_perm, out_nex, "nexus")
		nex = Nexus.Nexus(out_nex.getvalue())
		aln_nex.append((name, nex))

	#output nexus
	concat_nex = Nexus.combine(aln_nex)
	concat_nex.write_nexus_data("%s.nex" % (outfilename))


def permut_aln_per_site(aln):
	
	aln_array   = np.array([list(rec) for rec in aln], order="F")
	aln_permut  = aln_array[: , np.random.permutation(aln_array.shape[1])]

	seqs_list = []
	for seqs in aln_permut:
		seq_str = "".join(seqs.tolist())
		new_seq = Seq(seq_str, IUPAC.unambiguous_dna)
		seqs_list.append(new_seq)

	multi_seqs = []
	for index, seqs in  enumerate(aln):
		multi_seqs.append(SeqRecord(seqs_list[index], seqs.id))

	return(MultipleSeqAlignment(multi_seqs))




