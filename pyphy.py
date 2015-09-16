import HyPhy

class PyPhy:
    def __init__(self, cwd='', nthreads=1):
        self.__instance = HyPhy._THyPhy(cwd, nthreads)

    def call(self, msg, flush=False):
        """
        Executes batch language
        :param msg: HBL command
        :param flush: if False, preserves execution state of the system
        :return: any stdout from HBL execution
        """
        p = self.__instance.ExecuteBF(msg, False)
        return p.sData

    def get(self, key, cast=None):
        response = self.__instance.AskFor(key)
        if not response:
            return None
        if self.__instance.CanCast(response, cast):
            result = self.__instance.CastResult(response, cast)
            if cast == HyPhy.THYPHY_TYPE_NUMBER:
                return result.castToNumber()
            elif cast == HyPhy.THYPHY_TYPE_MATRIX:
                return result.castToMatrix()
            elif cast == HyPhy.THYPHY_TYPE_STRING:
                return result.castToString()
            else:
                print 'Warning: Unreconigzed cast request in PyPhy.get(), ', cast
            return result
        return None

    def stdout(self):
        return self.__instance.GetStdout()

    def stderr(self):
        return self.__instance.GetErrors()

    def warnings(self):
        return self.__instance.GetWarnings()

    def read_data(self, fasta, dsid='my_ds'):
        """
        Read in FASTA
        :param fasta: Python list containing sublists of header/sequence pairs
        :return:
        """
        input = ''
        for h, s in fasta:
            input += '>%s %s ' % (h, s)
        self.call('DataSet %s = ReadFromString(%s);' % (dsid, fasta))
        nspecies = self.call('return %s.species;' % (dsid,))
        return nspecies
