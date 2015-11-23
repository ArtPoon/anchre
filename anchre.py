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

class Anchre:
    def __init__(self,
                 csv=None, iso_format=True, origin=date(1970,1,1),
                 delimiter='_', field=-1,
                 ft2path=None, Rpath=None, java=None,
                 tmpfile='anchre-tmp',
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
        self.fasta = {}
        self.last_date = None

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


    def parse_date(self, date_str):
        """
        Convert a string representation of a sample collection date into
        an integer value (number of days since some time in the past).
        """
        if self.iso_format:
            try:
                days = (dateup.parse(date_str).date() - self.origin).days
            except ValueError:
                print 'ERROR: Failed to parse date', date_str
                raise
        else:
            # expressed as number of days since some time in the past (BEAST style)
            try:
                days = int(date_str)
            except:
                print 'ERROR: Expected integer value for sequence date, found', date_str
                raise
        return days

    def parse_dates_csv (self):
        """
        Parse dates from CSV.
        """
        reader = DictReader(open(self.csv, 'rU'))
        for row in reader:
            self.dates.update({row['header']: self.parse_date(row['date'])})


    def read (self, handle):
        """
        Parse open file as FASTA.  Clean sequence labels.
        """
        self.fasta = {}  # reset container
        h = None
        sequence = ''
        count = 0
        for line in handle:
            if line.startswith('$'):  # skip comments
                continue

            if line.startswith('>') or line.startswith('#'):
                if sequence:
                    days = self.get_date(h)
                    self.fasta.update({h: {'header': '%d_%d' % (count, days),
                                           'sequence': sequence,
                                           'days': days}})
                    sequence = ''   # reset container

                h = line.strip('>#\n')
                count += 1
            else:
                sequence += line.strip('\n').upper()

        # append last entry
        days = self.get_date(h)
        self.fasta.update({h: {'header': '%d_%d' % (count, days),
                               'sequence': sequence,
                               'days': days}})

        # determine most recent sample date
        all_dates = [v['days'] for v in self.fasta.itervalues()]
        assert len(set(all_dates)) > 0, 'ERROR: Only one sample date in data'
        all_dates.sort(reverse=True)
        self.last_date = all_dates[0]



    def get_date (self, h):
        """
        If dates were provided as a CSV input file, then return the
        date associated with the sequence header supplied as the first argument.
        Otherwise, parse date field from sequence header in FASTA object.
        """
        if self.csv:
            try:
                # this will always be days since X
                return self.dates[h]
            except:
                print 'ERROR: sequence header', h, 'not found in dates parsed from CSV'
                raise

        # otherwise, parse date from sequence headers
        date_field = h.split(self.delimiter)[self.field]
        return self.parse_date(date_field)


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
        # determine the earliest sample date
        dates = [v['days'] for v in self.fasta.itervalues()]
        dates.sort()  # defaults to increasing order
        earliest_date = dates[0]

        # retrieve all sequences with this date
        first_sample = [v['sequence'] for k, v in self.fasta.iteritems() if v['days'] == earliest_date]
        return first_sample


    def consensus_earliest (self):
        """
        Return the consensus of sequences from the earliest sample.
        :param fasta:
        :return:
        """
        if not self.fasta:
            # no sequences have been parsed
            return None

        sample = self.earliest_sample()  # list of sequences
        return self.consensus(sample)

    def consensus_all (self):
        """
        Return the consensus of all sequences.
        :return:
        """
        all_seqs = [v['sequence'] for v in self.fasta.itervalues()]
        return self.consensus(all_seqs)

    def output_fasta(self):
        """
        Write contents of self.fasta to temporary file
        :return:  Absolute path to temporary file
        """
        with open(self.tmpfile, 'w') as f:
            for i, (h, data) in enumerate(self.fasta.iteritems()):
                f.write('>%s\n%s\n' % (data['header'], data['sequence']))


    def call_fasttree2 (self, raw=False):
        """
        Call FastTree2 on FASTA file
        :param raw: if True, retain original sequence headers
        :return:
        """
        self.output_fasta()  # writes to self.tmpfile

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

        out1 = os.path.join(self.tmp, 'anchre.r2t.timetree')
        out2 = os.path.join(self.tmp, 'anchre.r2t.csv')

        p = subprocess.check_call([self.java, '-jar', 'java/RLRootToTip.jar',
                              '-timetree', out1,
                              '-newick', self.tmpfile,
                              out2], stdout=subprocess.PIPE)

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

    def call_hyphy_ancre(self, tree, model_spec='010010', is_codon=False):
        """
        Ancestral reconstruction with HyPhy
        :param tree: Newick tree string
        :param is_codon: if True, interpret alignment as codon sequences
        :return: [ancseq] is a dictionary of header/sequence pairs.
                   "Node0" keys the root node.
                 [lf] is a serialization of the likelihood function.
        """

        # cast Newick tree string as Phylo object to extract tip labels
        handle = StringIO(tree)
        phy = Phylo.read(handle, 'newick')
        tips = phy.get_terminals()
        tipnames = [tip.name for tip in tips]
        tipnames.sort()

        # make sure the tree labels match the sequence headers
        headers = [v['header'] for v in self.fasta.itervalues()]
        headers.sort()
        if headers != tipnames:
            print 'Warning: tree labels do not match FASTA in call_hyphy_ancre()'
            return None, None

        ancseq, lf = self.pyphy.ancre(fasta=self.fasta, newick=tree,
                                      model_spec=model_spec, is_codon=is_codon)
        return dict(ancseq), lf


    def call_beast(self, chain_length=1E6,
                   screen_step=1E5,
                   log_step=1E4,
                   treelog_step=1E4,
                   sample_size=100):
        """
        Use BEAST to sample trees from the posterior density under a
        strict molecular clock model.  If you want different settings,
        modify the template XML file.
        :return: a list of Newick tree strings
        """
        log, treelog = self.beauti.populate(fasta=self.fasta,
                                            stem=os.path.join(self.tmp, 'beast'),
                                            chain_length=chain_length,
                                            screen_step=screen_step,
                                            log_step=log_step,
                                            treelog_step=treelog_step)
        self.beauti.write(self.tmpfile)
        # this was tested on version 1.8.1
        p = subprocess.check_call([self.java, '-jar', 'java/beast.jar',
                                   '-beagle_off',
                                   '-overwrite',
                                   self.tmpfile], stdout=subprocess.PIPE)

        with open(log, 'rU') as f:
            traces = self.beauti.parse_log(f)

        with open(treelog, 'rU') as f:
            trees = self.beauti.parse_treelog(f, sample_size=sample_size)

        return traces, trees




def main():
    # command-line execution
    parser = argparse.ArgumentParser(description='Timing and ancestral reconstruction from dated'
                                                 'sequence alignments.')

    # positional arguments
    parser.add_argument('fasta', help='<input FASTA> sequence alignment')
    parser.add_argument('anc', help='<output FASTA> ancestral sequence reconstructions')
    parser.add_argument('mrca', help='<output CSV> time to MRCA')

    # optional arguments
    parser.add_argument('-csv', default=None, help='(input CSV) optional file containing sequence headers '
                                                   'and sample collection dates')
    parser.add_argument('-iso', action='store_true',
                        help='Specify date format as either ISO or number of days since some'
                             'time in the past.')
    parser.add_argument('-origin', help='ISO formatted date.  For expressing ISO dates as number of days '
                                        'since date of origin.', default=None)
    parser.add_argument('-sep', choices=['_', '.', '-', ';', ':', ' '], default='_',
                        help='Sequence header field separator; make sure to enclose in single quotes.')
    parser.add_argument('-pos', type=int, default=-1, help='Python-style index for date in sequence header.')
    parser.add_argument('-ft2', default='/usr/local/bin/fasttree2', help='Absolute path to FastTree2')
    parser.add_argument('-R', default='/usr/bin/Rscript', help='Absolute path to Rscript')
    parser.add_argument('-java', default='/usr/bin/java', help='Absolute path to Java interpreter')

    args = parser.parse_args()

    if args.sep != '_':
        print 'Using custom separator "%s"' % (args.sep, )

    if args.iso:
        print 'Expecting ISO formatted dates.'

    csv = None if args.csv is None else open(args.csv, 'rU')
    anchre = Anchre(csv=csv, iso_format=args.iso,
                     delimiter=args.sep, field=args.pos,
                     ft2path=args.ft2, Rpath=args.R, java=args.java,
                     origin = dateup.parse(args.origin).date() if args.origin else date(1970, 1, 1))

    # load sequences
    with open(args.fasta, 'rU') as f:
        anchre.read(f)

    # prepare outputs
    fanc = open(args.anc, 'w')
    fmrca = open(args.mrca, 'w')
    fmrca.write('method,value,lower.CI,upper.CI\n')

    # construct a tree
    print 'fasttree2'
    tree = anchre.call_fasttree2()  # writes sequences to temporary FASTA file
    #print tree

    # root the tree using tip dates
    print 'rtt'
    rooted = anchre.call_rtt(tree)
    #print rooted

    # re-estimate node heights
    print 'root2tip'
    res = anchre.call_root2tip(rooted)
    timetree = res['timetree']
    origin = res['x-intercept']
    fmrca.write('root2tip,%s,,\n' % (origin,))

    # get the consensus sequence
    conseq = anchre.consensus_earliest()
    fanc.write('>consensus_earliest\n%s\n' % (conseq,))
    conseq = anchre.consensus_all()
    fanc.write('>consensus_all\n%s\n' % (conseq,))

    # load sequences into HyPhy
    print 'ML anc res'
    ancseq, lf = anchre.call_hyphy_ancre(tree=rooted, model_spec='010020', is_codon=False)
    try:
        fanc.write('>ML_rooted\n%s\n' % (ancseq['Node0'], ))
    except:
        print ancseq
        raise

    ancseq, _ = anchre.call_hyphy_ancre(tree=timetree, model_spec='010020', is_codon=False)
    fanc.write('>ML_timetree\n%s\n' % (ancseq['Node0'], ))

    # BEAST!
    print 'BEAST'
    traces, trees = anchre.call_beast(sample_size=1000)
    #print traces['posterior']

    sample_origin = map(lambda t: anchre.last_date - t, traces['treeModel.rootHeight'][10:])
    sample_origin.sort()
    n = len(sample_origin)
    median = ((sample_origin[n/2-1] + sample_origin[n/2])/2.) if n % 2 == 0 else sample_origin[n/2]
    lower = sample_origin[int(0.025*n)]
    upper = sample_origin[int(0.975*n)]
    fmrca.write('BEAST,%f,%f,%f\n' % (median, lower, upper))

    print 'BEAST anc res'
    mctree = anchre.beauti.max_credible(trees)
    #print mctree

    # get additional ancestral sequence reconstructions
    ancseq, _ = anchre.call_hyphy_ancre(tree=mctree, model_spec='010020', is_codon=False)
    fanc.write('>ML_BEAST_mctree\n%s\n' % (ancseq['Node0'], ))
    for rep, tree in enumerate(trees):
        ancseq, _ = anchre.call_hyphy_ancre(tree=tree, model_spec='010020', is_codon=False)
        fanc.write('>ML_BEAST_tree%d\n%s\n' % (rep, ancseq['Node0']))

    fanc.close()
    fmrca.close()


if __name__ == "__main__":
    main()
