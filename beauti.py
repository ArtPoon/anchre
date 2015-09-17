from xml.etree.ElementTree import ElementTree as Tree
from xml.etree.ElementTree import Element as Node
import re
import subprocess
import math
from Bio import Phylo
from StringIO import StringIO


class Beauti:
    def __init__(self, template_file):
        self.template = Tree()
        _ = self.template.parse(template_file)

    def populate(self, fasta, stem, time_unit='days',
                 chain_length=None, screen_step=None, log_step=None, treelog_step=None):
        """
        Load sequences from FASTA object into BEAST XML template
        :param fasta: a Python list object containing sublists of header/sequence pairs
        :param stem: file path and prefix to write *.log and *.tree files
        :param time_unit: used by BEAST for annotation only (e.g., days, years)
        :param chain_length: optional setting for number of steps in MCMC chain
        :return: paths to BEAST log and tree log files
        """
        logfile = stem+'.log'
        treefile = stem+'.trees'

        # reset TAXA and ALIGNMENT blocks
        t_taxa = self.template.findall('taxa')[0]
        t_taxa._children = []
        t_aln = self.template.find('alignment')
        t_aln._children = []

        for h, s in fasta:
            date_val = float(h.split('_')[-1])
            date = Node('date', {'units': time_unit,
                                 'direction':'forwards',
                                 'value':str(date_val)})
            # TAXA
            taxon = Node('taxon', {'id': h})
            taxon.append(date)
            t_taxa.append(taxon)

            # SEQUENCE
            seqtag = Node('sequence', {})
            staxon = Node('taxon', {'idref': h})
            staxon.tail = '\n\t\t\t'+s.upper()  # mimic formatting in BEAST XML
            seqtag.append(staxon)
            t_aln.append(seqtag)

        # revise log settings
        t_mcmc = self.template.find('mcmc')
        if chain_length:
            t_mcmc.set('chainLength', str(int(chain_length)))  # number of MCMC steps

        for log in t_mcmc.findall('log'):
            if log.get('id') == 'fileLog':
                log.set('fileName', logfile)
                if log_step:
                    log.set('logEvery', str(int(log_step)))
            elif log.get('id') == 'screenLog':
                if screen_step:
                    log.set('logEvery', str(int(screen_step)))

        log_tree_element = t_mcmc.find('logTree')
        log_tree_element.set('fileName', treefile)
        if treelog_step:
            log_tree_element.set('logEvery', str(int(treelog_step)))

        return logfile, treefile

    def write(self, handle):
        self.template.write(handle)

    def parse_log(self, handle, mod=1):
        """
        Parse contents of a BEAST log file
        :param handle: an open file handle
        :return:
        """
        keys = None
        result = {}
        count = 0
        for line in handle:
            if line.startswith('#'):
                # ignore comment line
                continue
            if not result:
                # reached the first non-comment line
                keys = line.strip('\n').split('\t')
                result = dict([(key, []) for key in keys])
                continue
            if count % mod == 0:
                values = line.strip('\n').split('\t')
                for i, k in enumerate(keys):
                    v = values[i]
                    result[k].append(int(v) if k == 'state' else float(v))
            count += 1
        return result


    def parse_treelog(self, handle, sample_size):
        """
        Extract sampled trees from BEAST tree log and convert into
        Newick tree strings.
        :param handle:
        :return:
        """

        # use grep to figure out how many trees are in the file
        p = subprocess.Popen(['grep', '--count', '-e', '^tree', handle.name],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        ntrees = int(stdout.strip('\n'))

        # figure out a modulus to get target sample
        mod = max(ntrees / sample_size, 1)

        trBlock = False
        label_d = {}
        newicks = []
        tree_count = 0

        for line in handle:
            if 'Translate' in line:
                # entered taxa block
                trBlock = True
                continue

            if trBlock:
                if ';' in line:
                    # end of block
                    trBlock = False
                    continue
                index, label = line.strip(',\n').split()[-2:]
                label_d.update({index: label})
                continue

            # this should follow the Translate block
            if line.startswith('tree'):
                tree_count += 1
                if tree_count % mod == 0:
                    tree = line.split()[-1].strip(';')
                    # remove figtree annotations
                    tree1 = re.sub('\[[^]]+\]', '', tree)
                    # replace all indices with labels
                    tree2 = re.sub('([0-9]+):', lambda x: label_d[x.groups()[0]]+':', tree1)
                    newicks.append(tree2)
                continue

        return newicks


    def max_credible(self, newicks):
        """
        Return the maximum clade crediblity tree
        :param newicks: a list of Newick tree strings
        :return:
        """
        def get_tips(tree):
            return tuple(sorted([tip.name for tip in tree.root.get_terminals()]))

        def get_clades(tree):
            # recursive function to collect all monophyletic clades of tips
            result = [get_tips(tree)]
            for subtree in tree.root:
                if not subtree.is_terminal():
                    result.extend(get_clades(subtree))
            return result

        ntrees = len(newicks)
        counts = {}
        clades_cache = {}

        for newick in newicks:
            handle = StringIO(newick)
            tree = Phylo.read(handle, 'newick')

            clades = get_clades(tree)
            for clade in clades:
                if clade not in counts:
                    counts.update({clade: 0.})
                counts[clade] += 1./ntrees
            clades_cache.update({newick: clades})

        def credibility(clades, counts):
            return sum(math.log(counts[clade]) for clade in clades)

        return max(clades_cache.keys(), key=lambda t: credibility(clades_cache[t], counts))



def main():
    # derived from shankarappa test set
    test_fasta = [
        ['AF137634_p1c003_20_90',
         'GAAGAAGAGGTAGTAATTAGATCTGAAAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAAGGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAATAACAATACAAGAAAAAGTATACAGGTAGGACCAGGGAAAGCAATTTATACAACAGGAG'],
        ['AF137631_p1c003_17_90',
         'GAAGAAGAGGTAGTAATTAGATCTGAAAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAAGGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAATAACAATACAAGAAAAAGTATACAGGTAGGACCAGGGAAGGCAATTTATACAACAGGAG'],
        ['AF137725_p1p080_341_2400',
         'GAAGAAGAGGTAGTAATTAGATCTGAAAATTTCACGAACAATGCTAAGACCATAATAGTACAGCTGAATGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAGTAACAATACAAGAAGAAGTATAGCGGTAGGACCAGGGAGAGCATTTTATGCAACAGATA'],
        ['AF137724_p1p080_340_2400',
         'GAAAAAGAGGTAGTAATCAGATCTGAAAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAATGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAGTAACAATACAAGAAGAAGTATAGCGGTAGAACCAGGGAGAGCATTTTATGCAACAGATC'],
        ['AF137729_p1p087_206_2610',
         'GAAAAAGAGGTAATAATTAGATCTGAAAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAATGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAGTAACAATACAAGAAGAAGTATAGCGGTAGGACCAGGGAGAGCATTTTATGCAACAGATC'],
        ['AF137741_p1p094_347_2820',
         'GAAAAAGAGGTAGTAATCAGATCTGAAAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAATGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAGTAACAACACAAGAAGAAGTATAGCGGTAGGACCAGGGAGAGCATTTTATGCAACAGATC'],
        ['AF137759_p1p105_272_3150',
         'GAAAAAGAGGTAGTAATTAGATCTGAAAATTTCACGAACAATGCTAAAACCATAATAGTACAGCTAAATGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAGTAACAATACAAGAAGAAGTATAGCAGTAGGACCAGGGAGAGCATTTTATGCAACAGATA'],
        ['AF137758_p1p105_269_3150',
         'GAAGAAGAGGTAGTAATTAGATCTGAAAATTTCACGAACAATGCTAAAACCATAATAGTACAGCTGAATGAGTCTGTAGT'
         'AATTAATTGTACAAGACCCAGTAACAATACAAGAAGAAGTATAGCGGTAGGACCAGGAAGAGCATTTTATGCAACAGATC']
    ]

    # some tests
    beauti = Beauti('beast-template.xml')
    beauti.populate(fasta=test_fasta,
                    stem='tests/beast-test',
                    chain_length='100000',
                    screen_step='10000',  # 10 lines
                    log_step='100',
                    treelog_step='100')

    with open('tests/beast-test.xml', 'w') as f:
        beauti.write(f)

    p = subprocess.check_call(['java', '-jar', 'java/beast.jar',
                               '-beagle_off',
                               '-overwrite',
                               'tests/beast-test.xml'])

    with open('tests/beast-test.trees', 'rU') as f:
        # this log contains 100 trees
        trees = beauti.parse_treelog(f, 100)

    mctree = beauti.max_credible(trees)
    print mctree

if __name__ == '__main__':
    main()





