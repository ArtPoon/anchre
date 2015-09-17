from xml.etree.ElementTree import ElementTree as Tree
from xml.etree.ElementTree import Element as Node

class Beauti:
    def __init__(self, template_file):
        self.template = Tree()
        self.t_root = self.template.parse(template_file)

    def populate(self, fasta, outfile, stem, time_unit='days'):
        """
        Load sequences from FASTA object into BEAST XML template
        :param fasta: a Python list object containing sublists of header/sequence pairs
        :return:
        """
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
            taxon = Node('taxon', {'id':id})
            taxon.append(date)
            t_taxa.append(taxon)

            # SEQUENCE
            seqtag = Node('sequence', {})
            staxon = Node('taxon', {'idref': h})
            staxon.tail = '\n\t\t\t'+s.upper()  # mimic formatting in BEAST XML
            seqtag.append(staxon)
            t_aln.append(seqtag)

        t_mcmc = self.template.find('mcmc')
        log_elements = t_mcmc.findall('log')
        log_tree_element = t_mcmc.find('logTree')

        for log in log_elements:
            if log.get('id') == 'fileLog':
                log.set('filename', stem+'.log')
                break

        log_tree_element.set('fileName', stem+'.trees')
        self.template.write(outfile)




