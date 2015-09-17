"""
Objective: 
- parse a FASTA file with sample collection dates
- ancestral sequence reconstruction at root (genomic sequence, S)
- estimate date of infection (TMRCA, t in units of days)
- estimate number of founder viruses (k)

It is assumed that the FASTA input represents a multiple sequence alignment.
"""

import argparse
from csv import DictReader
from datetime import date
import dateutil.parser as dateup
import tempfile
import os
import sys
import re
import subprocess
from pyphy import PyPhy
from beauti import Beauti
from Bio import Phylo, SeqIO
from StringIO import StringIO

class IronChef:
    def __init__(self,
                 csv=None, iso_format=True, origin=date(1970,1,1),
                 delimiter='_', field=-1,
                 ft2path=None, Rpath=None, java=None,
                 tmpfile='ironchef-tmp',
                 beast_xml_template='beast-template.xml'):
        self.csv = csv
        self.iso_format = iso_format
        self.origin = origin
        self.delimiter = delimiter
        self.field = field

        # paths to binaries
        self.ft2path = ft2path
        self.Rpath = Rpath
        self.java = java

        self.pyphy = PyPhy(os.getcwd(), 1)  # instance of HyPhy
        self.beauti = Beauti(beast_xml_template)

        # if given, parse dates from csv
        self.dates = {}
        if self.csv is not None:
            self.parse_date_csv()

        # store sequence records
        self.fasta = []
        self.seqs = {}  # grouped by collection date

        self.tmp = tempfile.gettempdir()
        self.tmpfile = os.path.join(self.tmp, tmpfile)
        self.test()


    def test(self):
        """
        Check whether expected binaries are accessible
        :return:
        """
        if not os.path.exists(self.ft2path):
            print 'ERROR: Failed to detect FastTree2 at', self.ft2path
            sys.exit()
        if not os.path.exists(self.Rpath):
            print 'ERROR: Failed to detect R at', self.ft2path
            sys.exit()


    def parse_date_csv (self):
        """
        Parse dates from CSV
        """
        reader = DictReader(open(self.csv, 'rU'))
        for row in reader:
            if self.iso_format:
                days = dateup.parse(row['date']).date() - self.origin
            else:
                days = int(row['date'])

            self.dates.update({row['header']: days})

    def parse_fasta (self):
        """
        :return:
        """
        for h, s in self.fasta:
            if self.csv is None:
                # parse date from header
                this_date = int(h.split(self.delimiter)[self.field])
            else:
                try:
                    this_date = self.dates[h]
                except:
                    print 'ERROR: sequence header', h, 'not found in dates parsed from CSV'
                    raise
            if this_date not in self.seqs:
                self.seqs.update({this_date: []})

            self.seqs[this_date].append(s)

    def read (self, handle):
        """
        Parse open file as FASTA.  Clean sequence labels.
        """
        self.fasta = []  # reset container
        sequence = ''
        for i in handle:
            if i[0] == '$': # skip h info
                continue
            elif i[0] == '>' or i[0] == '#':
                if len(sequence) > 0:
                    self.fasta.append([h, sequence])
                    sequence = ''   # reset containers
                h = i.strip('\n')[1:]
                h = re.sub('[-:|.]', '_', h)
            else:
                sequence += i.strip('\n').upper()
        self.fasta.append([h, sequence])
        self.parse_fasta()


    def plurality_consensus(self, column, alphabet='ACGT', resolve=False):
        """
        Plurality consensus - nucleotide with highest frequency.
        In case of tie, report mixtures.
        """
        mixture_dict = {'W':'AT', 'R':'AG', 'K':'GT', 'Y':'CT', 'S':'CG',
                'M':'AC', 'V':'AGC', 'H':'ATC', 'D':'ATG',
                'B':'TGC', 'N':'ATGC', '-':'ATGC'}
        ambig_dict = dict(("".join(sorted(v)), k) for k, v in mixture_dict.iteritems())
        freqs = {}

        for char in alphabet:
            freqs.update({char: 0})
        #freqs = {"A": 0, "T": 0, "C": 0, "G": 0, "-": 0}
        for char in column:
            if char in alphabet:
                freqs[char] += 1
            elif mixture_dict.has_key(char):
                # handled ambiguous nucleotides with equal weighting
                resolutions = mixture_dict[char]
                for char2 in resolutions:
                    freqs[char2] += 1./len(resolutions)
            else:
                # unrecognized nucleotide character
                pass

        base = max(freqs, key=lambda n: freqs[n])
        max_count = freqs[base]
        possib = filter(lambda n: freqs[n] == max_count, freqs)
        if len(possib) == 1:
            return possib[0]
        elif "-" in possib:
            if resolve:
                possib.remove("-")
                if len(possib) == 0:
                    return "-"
                elif len(possib) == 1:
                    return possib[0]
                else:
                    return ambig_dict["".join(sorted(possib))]
            else:
                # gap character overrides ties
                return "-"
        else:
            return ambig_dict["".join(sorted(possib))]

    def consensus(self, seqs, alphabet='ACGT', resolve=False):
        """
        Return plurality consensus of alignment.
        """
        # transpose the alignment
        n_columns = len(seqs[0])
        columns = []
        for c in range(n_columns):
            columns.append ( [ s[c] for s in seqs ] )

        consen = []
        for column in columns:
            consen.append(self.plurality_consensus(column, alphabet=alphabet, resolve=resolve))

        return "".join(consen)

    def earliest_sample (self):
        dates = self.seqs.keys()
        dates.sort()  # defaults to increasing order
        return self.seqs[dates[0]]


    def consensus_earliest (self):
        """
        Return the consensus of sequences from the earliest sample.
        :param fasta:
        :return:
        """
        if not self.seqs:
            # no sequences have been parsed
            return None

        sample = self.earliest_sample()  # list of sequences
        return self.consensus(sample)

    def consensus_all (self):
        """
        Return the consensus of all samples.
        :return:
        """
        all_seqs = []
        for date, sample in self.seqs.iteritems():
            all_seqs.extend(sample)
        return self.consensus(all_seqs)

    def output_fasta(self):
        """
        Write contents of self.fasta to temporary file
        :return:
        """
        with open(self.tmpfile, 'w') as f:
            for h, s in self.fasta:
                f.write('>%s\n%s\n' % (h, s))

    def output_seqs (self):
        """
        Write out parsed sequences to temporary file.
        :return:
        """
        fasta = []
        for date, sample in self.seqs.iteritems():
            for i, seq in enumerate(sample):
                fasta.append(['%d_%d' % (i, date), seq])
        with open(self.tmpfile, 'w') as f:
            for h, s in fasta:
                f.write('>%s\n%s\n' % (h, s))


    def call_fasttree2 (self, raw=True):
        """
        Call FastTree2 on FASTA file
        :param raw: if True, retain original sequence headers
        :return:
        """
        if raw:
            self.output_fasta()  # writes to self.tmpfile
        else:
            self.output_seqs()

        p = subprocess.Popen([self.ft2path, '-quiet', '-nosupport', '-nt', '-gtr'],
                             stdin=open(self.tmpfile, 'rU'),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        return stdout

    def call_rtt (self, tree):
        """
        Call an R script that implements Rosemary's rtt() function for re-rooting
        a tree based on tip dates.
        :param tree: Newick tree string
        :return: dictionary with two key-value pairs for rooted and dated trees
        """
        os.chdir('R/')
        p = subprocess.Popen([self.Rpath, 'rtt.r', tree], stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
        # clean kludge from R stdout
        #rooted_tree, dated_tree = map(lambda s: s.replace('[1] "', '').replace('NA;"', '0:0;'),
        #                              stdout.split('\n')[:2])
        #res = {'rooted': rooted_tree, 'dated': dated_tree}
        rooted_tree = stdout.replace('[1] ', '').strip('"\n')
        os.chdir('../')
        return rooted_tree

    def call_root2tip(self, tree):
        """
        Call jar file that implements a modified version of Andrew Rambaut's
        root-to-tip method (Path-O-Gen).
        :param tree: a Newick tree string
        :return: a dictionary that includes the time-scaled tree
        """
        # write tree to temporary file
        with open(self.tmpfile, 'w') as handle:
            handle.write(tree)

        out1 = os.path.join(self.tmp, 'ironchef-tmp.r2t.timetree')
        out2 = os.path.join(self.tmp, 'ironchef-tmp.r2t.csv')

        p = subprocess.check_call([self.java, '-jar', 'java/RLRootToTip.jar',
                              '-timetree', out1,
                              '-newick', self.tmpfile,
                              out2])

        # read outputs
        with open(out1, 'rU') as handle:
            timetree = Phylo.read(handle, 'nexus')
        with open(out2, 'rU') as handle:
            coef = handle.readlines()

        # convert NEXUS to Newick string
        output = StringIO()
        Phylo.write(timetree, output, 'newick')
        res = {'timetree': output.getvalue()}
        values = coef[1].strip('\n').split(',')
        for i, key in enumerate(coef[0].strip('\n').split(',')):
            res.update({key: values[i]})

        return res

    def call_hyphy_ancre(self, tree, is_codon):
        """

        :param tree: Newick tree string
        :param is_codon: if True, interpret alignment as codon sequences
        :return: [ancseq] is a dictionary of header/sequence pairs.
                   "Node0" keys the root node.
                 [lf] is a serialization of the likelihood function.
        """
        # make sure the tree labels match the sequence headers
        headers = [h for h, s in self.fasta]
        headers.sort()
        handle = StringIO(tree)
        phy = Phylo.read(handle, 'newick')
        tips = phy.get_terminals()
        tipnames = [tip.name for tip in tips]
        tipnames.sort()
        if headers != tipnames:
            print 'Warning: tree labels do not match FASTA in call_hyphy_ancre()'
            return None, None

        ancseq, lf = self.pyphy.ancre(self.fasta, tree, is_codon)
        return dict(ancseq), lf


    def call_beast(self, sample_size=10):
        """
        Use BEAST to sample trees from the posterior density under a
        strict molecular clock model.  If you want different settings,
        modify the template XML file.
        :return: a list of Newick tree strings
        """
        log, treelog = self.beauti.populate(fasta=self.fasta,
                                            stem=os.path.join(self.tmp, 'beast'))
        self.beauti.write(self.tmpfile)
        # this was tested on version 1.8.1
        p = subprocess.check_call([self.java, '-jar', 'java/beast.jar',
                                   '-beagle_off',
                                   '-overwrite',
                                   self.tmpfile])
        with open(treelog, 'rU') as f:
            trees, mctree = self.beauti.parse_treelog(f, sample_size=sample_size)

        return trees, mctree







def main():
    # command-line execution
    parser = argparse.ArgumentParser(description='Timing and ancestral reconstruction from dated'
                                                 'sequence alignments.')

    parser.add_argument('fasta', help='<input FASTA> sequence alignment')
    parser.add_argument('-csv', default=None, help='(input CSV) optional file containing sequence headers '
                                                   'and sample collection dates')
    parser.add_argument('-iso', action='store_true',
                        help='Specify date format as either ISO or number of days since some'
                             'time in the past.')
    parser.add_argument('-sep', choices=['_', '.', '-', ';', ':', ' '], default='_',
                        help='Sequence header field separator.')
    parser.add_argument('-pos', type=int, default=-1, help='Python-style index for date in sequence header.')
    parser.add_argument('-ft2', default='/usr/local/bin/fasttree2', help='Absolute path to FastTree2')
    parser.add_argument('-R', default='/usr/bin/Rscript', help='Absolute path to Rscript')
    parser.add_argument('-java', default='/usr/bin/java', help='Absolute path to Java interpreter')

    args = parser.parse_args()

    csv = None if args.csv is None else open(args.csv, 'rU')
    ichef = IronChef(csv=csv, iso_format=args.iso,
                     delimiter=args.sep, field=args.pos,
                     ft2path=args.ft2, Rpath=args.R, java=args.java)

    # load sequences
    with open(args.fasta, 'rU') as f:
        ichef.read(f)

    # construct a tree
    tree = ichef.call_fasttree2(raw=True)
    print tree

    # root the tree using tip dates
    rooted = ichef.call_rtt(tree)
    print rooted

    # re-estimate node heights
    res = ichef.call_root2tip(rooted)
    timetree = res['timetree']
    print timetree

    # get the consensus sequence
    conseq = ichef.consensus_earliest()
    print conseq

    # load sequences into HyPhy
    ancseq, lf = ichef.call_hyphy_ancre(tree=rooted, is_codon=False)
    print ancseq['Node0']

    # BEAST!
    ichef.call_beast()




if __name__ == "__main__":
    main()
