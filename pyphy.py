import os
import HyPhy

class PyPhy:
    def __init__(self, cwd='', nthreads=1,
                 dsid='my_ds', trid='my_tree', lfid='my_lf'):
        self.__instance = HyPhy._THyPhy(cwd, nthreads)
        self.version = self._get_version()

        self.dsid = dsid
        self.dsfid = self.dsid+'f'
        self.trid = trid
        self.lfid = lfid

    def _get_version(self):
        self.call('GetString(v, HYPHY_VERSION, 1);')
        return self.get('v', as_string=True)

    def call(self, msg, flush=False):
        """
        Executes batch language
        :param msg: HBL command
        :param flush: if False, preserves execution state of the system
        :return: any stdout from HBL execution
        """
        p = self.__instance.ExecuteBF(msg, False)
        return p.sData

    def get(self, expr, as_string=False):
        """
        Retrieve parameter value by name
        :param expr: HBL expression
        :return:
        """
        val = None
        self.call('val=%s;' % (expr,))
        self.call('fprintf(stdout, val);')
        if as_string:
            exec 'val="""%s"""' % (self.stdout(),)
        else:
            exec 'val=%s' % (self.stdout(),)  #FIXME: only works for scalars
        return val

    def stdout(self):
        return self.__instance.GetStdout().sData

    def stderr(self):
        return self.__instance.GetErrors().sData

    def warnings(self):
        return self.__instance.GetWarnings().sData

    def read_data(self, fasta, is_codon=False):
        """
        Read in FASTA as DataSet object in HyPhy instance
        :param fasta: Python list containing sublists of header/sequence pairs
        :return: Number of sequences read from string
        """
        input = ''
        for h, s in fasta:
            input += '>%s\n%s\n' % (h, s)
        self.call('DataSet %s = ReadFromString("%s");' % (self.dsid, input))
        if is_codon:
            # assumes universal code
            self.call('DataSetFilter %sf = CreateFilter(%s, 3, "", "", '
                      '"TAA,TGA,TAG");' % (self.dsid, self.dsid))
        else:
            self.call('DataSetFilter %sf = CreateFilter(%s, 1);' % (self.dsid, self.dsid))

        # retrieve the number of rows in alignment
        return self.get('%sf.species' % (self.dsid, ))

    def read_tree(self, newick):
        """
        Read Newick tree string as Tree object in HyPhy
        :param newick: Newick tree string
        :param trid: name for Tree object
        :return:
        """
        self.call('ACCEPT_ROOTED_TREES=1;')  # disable automatic unrooting
        self.call('ACCEPT_BRANCH_LENGTHS=1;')
        self.call('Tree %s = "%s";' % (self.trid, newick))

        return self.get('TipCount(%s)' % (self.trid,), as_string=True)

    def set_model(self, model_spec='010010', is_codon=False):
        """
        Define a nucleotide or codon substitution model in HyPhy.
        :param model_spec: PAUP* style model specification string
        :param is_codon: if True, use MG94
        :return:
        """

        self.call('modelDesc = "%s";' % (model_spec,))
        #print self.get('modelDesc', as_string=True)

        if is_codon:
            self.call('ExecuteAFile("%s/HBL/loadGeneticCode");' % (os.getcwd(),))
            self.call('HarvestFrequencies (observedFreq,%s,3,1,1);' % (self.dsfid,))
            #print self.get('observedFreq', as_string=True)
            self.call('ExecuteAFile("%s/HBL/MG94custom");' % (os.getcwd(),))
            #print self.get('MG94custom', as_string=True)
        else:
            self.call('HarvestFrequencies(vectorOfFrequencies,%s,1,1,0);' % (self.dsfid,))
            self.call('FREQUENCY_SENSITIVE=1;')
            self.call('ExecuteAFile("%s/HBL/custm4x4");' % (os.getcwd(),))


    def ancre(self, fasta, newick, model_spec='010010', is_codon=False):
        """
        Ancestral reconstruction using ReconstructAncestors() function in HyPhy
        :param dsfid:  name of DataSetFilter object
        :param trid:  name of Tree object
        :param lfid:  name for LikelihoodFunction object to build
        :return:
        """
        # FIXME: this isn't working (see hyphy issue #329)
        self.call('DATA_FILE_PRINT_FORMAT=0;')

        # set up analysis
        self.read_data(fasta, is_codon)
        self.set_model(model_spec=model_spec, is_codon=is_codon)
        self.read_tree(newick)

        if is_codon:
            # constrain branch lengths
            self.call('branchNames = BranchName(%s, -1);' % (self.trid,))
            self.call('branchLengths = BranchLength(%s, -1);' % (self.trid,))
            self.call('for(k=0; k<Columns(branchNames)-1; k=k+1) {'
                      'ExecuteCommands("%s."+branchNames[k]+".synRate:="+'
                      'branchLengths[k]/1000.+"/scalingB;");}')

        # optimize likelihood function
        self.call('LIKELIHOOD_FUNCTION_OUTPUT=7;')
        self.call('LikelihoodFunction %s = (%s, %s);' % (self.lfid, self.dsfid, self.trid))
        self.call('Optimize (res, %s);' % (self.lfid,))
        result = self.get('res', as_string=True)

        # ancestral reconstruction
        self.call('DataSet anc = ReconstructAncestors(%s);' % (self.lfid,))
        self.call('DataSetFilter ancf = CreateFilter(anc,1);')

        self.call('fprintf(stdout, ancf);')

        anc = {}
        header = None
        sequence = ''
        lines = self.stdout().split('\n')
        for line in lines:
            if line.startswith('#'):
                if sequence:
                    anc.update({header: sequence})
                    sequence = ''
                header = line.strip('#\n')
            else:
                sequence += line.upper()
        anc.update({header: sequence})

        # check for errors
        stderr = self.stderr()
        if stderr:
            print stderr

        return anc, result



def main():
    # test script
    p = PyPhy()

    from seqUtils import convert_fasta
    with open('test/shankarappa-p1.fa', 'rU') as f:
        fasta = convert_fasta(f)

    out = p.read_data(fasta, is_codon=True)
    print out

    t = '((A:1,B:2):0.5,C:0.2):0;'
    out = p.read_tree(t)
    print out

    p.set_model(is_codon=True)

if __name__ == '__main__':
    main()
