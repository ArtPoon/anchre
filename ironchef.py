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
import subprocess
from pyphy import PyPhy

class IronChef:
    def __init__(self, handle, csv=None, iso_format=True, origin=date(1970,1,1),
                 delimiter='_', field=-1, ft2path=None, Rpath=None, tmpfile='ironchef-tmp'):
        self.csv = csv
        self.iso_format = iso_format
        self.origin = origin
        self.delimiter = delimiter
        self.field = field

        # paths to binaries
        self.ft2path=ft2path
        self.Rpath=Rpath
        self.pyphy = PyPhy(os.getcwd(), 1)  # instance of HyPhy

        # if given, parse dates from csv
        self.dates = {}
        if self.csv is not None:
            self.parse_date_csv()

        # read sequence records from file handle
        fasta = self.iter_fasta(handle)
        self.seqs = {}
        self.parse_fasta(fasta)

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

    def parse_fasta (self, fasta):
        """
        :return:
        """
        for h, s in fasta:
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

    def iter_fasta (self, handle):
        """
        Parse open file as FASTA.  Returns a generator
        of handle, sequence tuples.
        """
        sequence = ''
        for i in handle:
            if i[0] == '$': # skip h info
                continue
            elif i[0] == '>' or i[0] == '#':
                if len(sequence) > 0:
                    yield h, sequence
                    sequence = ''   # reset containers
                h = i.strip('\n')[1:]
            else:
                sequence += i.strip('\n').upper()
        yield h, sequence

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

    def output_seqs (self, filename='ironchef-tmp.fa'):
        """
        Write out sequences to a temporary FASTA file.
        :return:
        """
        with open(filename, 'w') as f:
            for date, sample in self.seqs.iteritems():
                for i, seq in enumerate(sample):
                    f.write('>%d_%d\n%s\n' % (i, date, seq))
        return f.name

    def call_fasttree2 (self, handle):
        """
        Call FastTree2 on temporary FASTA
        :param filename:
        :return:
        """
        p = subprocess.Popen([self.ft2path, '-quiet', '-nt', '-gtr'],
                             stdin=handle,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        return stdout

    def call_root2tip (self, tree):
        """
        Call an R script that implements Rosemary's rtt() function for re-rooting
        a tree based on tip dates.
        :param tree: Newick tree string
        :return: dictionary with two key-value pairs for rooted and dated trees
        """
        p = subprocess.Popen([self.Rpath, 'rtt.r', tree], stdout=subprocess.PIPE)
        stdout, stderr = p.communicate()
        # clean kludge from R stdout
        rooted_tree, dated_tree = map(lambda s: s.replace('[1] "', '').replace('NA;"', '0:0;'), stdout.split('\n')[:2])
        res = {'rooted': rooted_tree, 'dated': dated_tree}
        return res

    def ancreml(self):
        """
        Ancestral sequence reconstruction at root by maximum likelihood.
        :return:
        """





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

    args = parser.parse_args()
    infile = open(args.fasta, 'rU')
    csvfile = None
    if args.csv is not None:
        csvfile = open(args.csv, 'rU')

    ichef = IronChef(handle=infile, csv=csvfile, iso_format=args.iso,
                     delimiter=args.sep, field=args.pos,
                     ft2path=args.ft2, Rpath=args.R)
    tmpfile = ichef.output_seqs()
    with open(tmpfile, 'rU') as f:
        tree = ichef.call_fasttree2(f)

    res = ichef.call_root2tip(tree)
    print res['dated']

    #conseq = ichef.consensus_earliest()


    infile.close()


if __name__ == "__main__":
    main()
