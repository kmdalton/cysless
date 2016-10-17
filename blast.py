###############################################################################
#                                                                             #
# Module for retrieving blast hits.                                           #
#                                                                             #
###############################################################################

import urllib2,re,subprocess,requests,datetime
from BeautifulSoup import BeautifulStoneSoup

AminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 
              'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y']

class ConnectivityError(Exception):
    """
    This is just a dummy class to help diagnose the source of an error. It may do more later.
    """
    pass

class blast_handle():
    """
    This class is used to make calls to the NCBI
    BLAST common url api (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo) 
    """
    def __init__(self, query):
        """
        Parameters
        ----------
        query : string
            The query for your blast search. Usually this will be a bare sequence, but users may also supply an 
            NCBI accession, GI, or FASTA
        """
        self.query  = query
        self.status = None
        self.rid    = None

    def request(self, **kw):
        """
        Parameters
        ----------
        **kwargs : optional
            Supply keyword arguments to define the options of the BLAST search. The available keyword options are 
            documented one the NCBI github (https://ncbi.github.io/blast-cloud/dev/api.html). The default arguments
            are reproduced here. Any properly formatted api option is supported.

                QUERY : string, default=self.query TODO: this needs to be seq basically always
                    Accession, GI, or FASTA.
                DATABASE : string, default='nr'
                    Blast database to use for this query
                PROGRAM BLAST : string, default='blastp'
                    One of blastn, megablast, blastp, blastx, tblastn, tblastx
                FORMAT_TYPE : string, default='XML'
                    Report type String  HTML, Text, XML, XML2, JSON2, or Tabular. HTML is the default.
                HITLIST_SIZE : string, default='3000'
                    Number of databases sequences to keep. The default is 3,000.
        Returns
        -------
        RID : string
            The request ID used allocated to your search by the BLAST server. 
        WAITTIME: int
            The number of seconds NCBI estimates your request will take. 
        """

        PutKwargs = ['AUTO_FORMAT', 'COMPOSITION_BASED_STATISTICS', 'DATABASE', 'DB_GENETIC_CODE', 'ENDPOINTS', 'ENTREZ_QUERY', 'EXPECT', 'FILTER', 'FORMAT_TYPE', 'GAPCOSTS', 'GENETIC_CODE', 'HITLIST_SIZE', 'I_THRESH', 'LAYOUT', 'LCASE_MASK', 'MATRIX_NAME', 'NUCL_PENALTY', 'NUCL_REWARD', 'OTHER_ADVANCED', 'PERC_IDENT', 'PHI_PATTERN', 'PROGRAM', 'QUERY', 'QUERY_FILE', 'QUERY_BELIEVE_DEFLINE', 'QUERY_FROM', 'QUERY_TO', 'SEARCHSP_EFF', 'SERVICE', 'THRESHOLD', 'UNGAPPED_ALIGNMENT', 'WORD_SIZE']

        #Default keyword arguments
        kw['QUERY'] = kw.get("QUERY", self.query)
        kw['FORMAT_TYPE'] = kw.get("FORMAT_TYPE", 'XML')
        kw['HITLIST_SIZE'] = kw.get('HITLIST_SIZE', 3000)
        kw['DATABASE'] = kw.get('DATABASE', 'nr')
        kw['PROGRAM'] = kw.get('PROGRAM', 'blastp')

        for kwarg in kw:
            if kwarg not in PutKwargs:
                raise TypeError("blast_handle.request got an unexpected keyword argument {}".format(kwarg))

        BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&"
        QueryString = '&'.join(['='.join((i, str(kw[i]))) for i in PutKwargs if i in kw])
        
        U  = urllib2.urlopen(BaseURL + QueryString)
        html = U.read()
        QBlastInfo = re.search(r"\<\!\-\-QBlastInfoBegin.+QBlastInfoEnd", html, re.DOTALL)
        QBlastInfo = QBlastInfo.group()
        RID        = QBlastInfo.split()[3]
        WAITTIME   = QBlastInfo.split()[6]
        try:
            WAITTIME = int(WAITTIME)
        except ValueError:
            print "Warning, invalid wait time returned by blast"

        return RID, WAITTIME

    def check_status(self):
        """
        Query the NCBI database to check the status of self.rid which is the ncbi server
        request id for the blast search initiated by self.blast_request()

        Returns
        -------
            bool
                returns True if the blast search is complete, False otherwise, and None if self.rid is unset
                also sets the instance variable self.status to the same value
        """
        if self.rid is None:
            self.status = None
            raise TypeError('The type of self.rid is None. Run self.blast_request() before checking status.')
        else:
        #Now we check to see if the server will return formatted results
        #This test relies on the idea that qblast only returns FORMAT_TYPE==XML
        #if the search is __complete__. Otherwise, we get an html waiting page
            text = self.ncbi_get(RID = self.rid, ALIGNMENTS='0', DESCRIPTIONS='0', FORMAT_TYPE='Text')
            if '<!DOCTYPE html' in text.split('\n')[0]:
                self.status = False
                return False
            else:
                self.status = True
                return True

    def fetch_result(self, **kw):
        """
        make an ncbi get request

        Parameters
        ----------
        **kwargs : optional
            Supply keyword arguments to define the options of the BLAST GET request. The available keyword options are 
            documented one the NCBI github (https://ncbi.github.io/blast-cloud/dev/api.html). The default values are:

        """
        kw['FORMAT_TYPE'] = kw.get('FORMAT_TYPE', 'XML')
        status = self.check_status()
        if status != True:
            raise TypeError("blast_handle.check_status() = {}. Status must be True if the results are ready to be downloaded from NCBI.".format(status))
        return self.ncbi_get(**kw)

    def ncbi_get(self, **kw):
        """
        make an arbitrary ncbi get request

        Parameters
        ----------
        **kwargs : optional
            Supply keyword arguments to define the options of the BLAST GET request. The available keyword options are 
            documented one the NCBI github (https://ncbi.github.io/blast-cloud/dev/api.html). The default values are:

        """
        BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&"
        GetKwargs = ['ALIGNMENTS', 'ALIGNMENT_VIEW', 'DESCRIPTIONS', 'ENTREZ_LINKS_NEW_WINDOW', 'EXPECT_LOW', 'EXPECT_HIGH', 'FORMAT_ENTREZ_QUERY', 'FORMAT_OBJECT', 'FORMAT_TYPE', 'NCBI_GI', 'RID', 'RESULTS_FILE', 'SERVICE', 'SHOW_OVERVIEW']

        for kwarg in kw:
            if kwarg not in GetKwargs:
                raise TypeError("blast_handle.fetch_result got an unexpected keyword argument {}".format(kwarg))
        kw['RID'] = kw.get('RID', self.rid)
        QueryString = '&'.join(['='.join((i, str(kw[i]))) for i in GetKwargs if i in kw])
        #print BaseURL + QueryString
        URL = urllib2.urlopen(BaseURL + QueryString)
        return URL.read()

    def __exit__(self):
        """
        delete a qblast RID from NCBI's servers

        Parameters
        ----------
        RID : the request ID you wish to delete
        """
        if self.rid != None:
            BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Delete&"
            URL = urllib2.urlopen(BaseURL + "RID=%s" %RID)
            URL.read()

class blast_results():
    """
    Turn the BLAST XML dump into a more useable object
    """
    def __init__(self, XML):
        """
        Parameters
        ----------
        XML : str
            The XML formatted BLAST results
        """
        self.soup = BeautifulStoneSoup(XML)
        query_length = int(self.soup.find("iteration_query-len").text)
        self.hits = sorted([blast_hit(str(hit), query_length) for hit in self.soup.findAll('hit')], key = lambda x: -x.score)
        self.uids = [i.accession for i in self.hits]
        #Add the aligned hit sequences to self.seqs

    def recommend_mutant(self, residues):
        """
        Parameters
        ----------
        residues : iterable
            The residue numbers (type int) that you wish to mutate
        """
        for hit in self.hits:
            if False not in [hit.qseq[i-1] != hit.hseq[i-1] and hit.hseq[i-1] in AminoAcids for i in residues]:
                return hit
        return None

    def fetch_full_sequences(self):
        """
        blaster.fetch_full_sequences()
        If blaster.uids is not None, download the full sequence of every uid and store it in blaster.seqs
        You should probably not use this unless you're desperate. As a caution, sometimes blast hits contain
        full genomes which might be awkard. 
        
        """
        if self.uids is not None:
            fasta= efetch(self.uids) #flat text fasta format is the default for efetch
            self.headers,self.seqs = zip(*(('>'+i.split('\n')[0], ''.join(i.split('\n')[1:])) for i in fasta.split('>')[1:]))

class blast_hit():
    """
    Create a simple object out of a blast hit's XML representation which exposes useful features. 
    """
    def __init__(self, XML, query_length=None):
        """
        Parameters
        ----------
        XML : str
            The XML corresponding to a single blast hit
        """
        self.soup  = BeautifulStoneSoup(XML)
        self.score = float(self.soup.hsp_score.text)
        self.accession = self.soup.hit_accession.text

        qseq = self.soup.hsp_qseq.text
        hseq = self.soup.hsp_hseq.text
        #Remove query sequence gaps from the alignment
        qseq,hseq = zip(*((i,j) for i,j in zip(qseq,hseq) if i!='-'))
        qseq,hseq = ''.join(qseq),''.join(hseq)

        #Number of gaps to pad the ends of the hit sequences
        lpad = int(self.soup.find("hsp_query-from").text) - 1 

        rpad = 0
        if query_length is not None:
            rpad = query_length - int(self.soup.find("hsp_query-to").text)
        
        self.qseq = '-'*lpad + qseq + '-'*rpad
        self.hseq = '-'*lpad + hseq + '-'*rpad

    def __str__(self):
       text = "Accession: {}\n{}\n".format(self.accession, self.soup.hit_def.text)
       text = text + format_alignment(self.qseq, self.hseq)
       return text

class efetch():
    def __init__(uids, **kw):
        self.kwargs = kw
        self.kwargs['db'] = kw.get('db', 'protein')
        self.kwargs['rettype'] = kw.get('rettype', 'fasta')
        self.kwargs['retmode'] = kw.get('retmode', 'text')
        self.kwargs['id'] = ','.join(uids)
        BaseURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        r = requests.post(BaseURL, data=kw)
        return r.text

class smith_waterman():
    """
    align two sequences using the smith-waterman algorithm via water of the EMBOSS suite

    Parameters
    ----------
        seq1: an amino acid sequence to be aligned with seq2
        seq2: an amino acid sequence
    Class Methods
    -------------
        smith_waterman.align(): make an external call to the emboss suite and align seq1&2 supplied at __init__

    Instance Variables
    ------------------
        smith_waterman.aln str is where we keep the text output from water
        smith_waterman.seq1 is the first sequence provided at instantiation
        smith_waterman.seq2 is the second sequence provided at instantiation
        smith_waterman.indentity  (float) is the sequence identity calculated by water
        smith_waterman.similarity (float) is the sequence similarity calculated by water
        smith_waterman.registered_seq2 is the second sequence with the appropriate gap structure such that it has a
            1-1 correspondence with seq1, residue by residue. this may include a series of gaps at the n and/or c-terminus.
    """
    def __init__(self, seq1, seq2, **kw):
        self.header1,self.header2 = kw.get('h1', ''),kw.get('h2','')
        self.seq1,self.seq2 = seq1,seq2
        self.aln,self.identity,self.similarity = None,None,None
        self.registered_seq2 = None
        self.align()

    def align(self):
        self.aln = subprocess.Popen(
                ["/bin/bash",
                    "-cs", 
                    "water <(echo %s) <(echo %s) stdout -gapopen=10 -gapextend=0.5 2>/dev/null" %(self.seq1, self.seq2)
                    ],
                stdout=subprocess.PIPE
            ).communicate()[0]
        self.similarity = float(re.findall(r'Similarity:.*\n', self.aln)[0].split()[-1].strip('%()'))
        self.identity   = float(re.findall(r'Identity:.*\n', self.aln)[0].split()[-1].strip('%()'))

        #This whole block just registers the second seq with the first so we can examine the snps later
        lines = (i for i in self.aln.split('\n') if len(i)>0 and i[0]!='#')
        S = ""
        for line1 in lines:
            lines.next()#The middle line in each group of three has no amino acids
            line2 = lines.next()
            if len(S) == 0:
                S = '-'*(int(line1.split()[0]) - 1) #Padding for the nterm of seq1
            s1,s2 = line1.split()[1],line2.split()[1]
            for aa1,aa2 in zip(s1,s2):
                if aa1 != '-':
                    S = S + aa2
        S = S.ljust(len(self.seq1), '-')
        self.registered_seq2 = S

def format_alignment(seq1, seq2, width=50):
    """
    Parameters
    ----------
    seq1 : str
        First aligned sequence
    seq2 : str
        Second aligned sequence
    width : int=50
        Optionaly specify the number of residues per line to print
    Returns
    -------
    text : str
        The multiline formatted alignment text
    """
    text = ''
    l1,l2,l3 = '','',''

    for i,(res1,res2) in enumerate(zip(seq1, seq2), 1):
        l1 += res1
        l3 += res2
        if res1 == res2:
            l2 += "|"
        else:
            l2 += ":"
        if i%width == 0 or i==len(seq1):
            l1 += " {}\n".format(i)
            l2 += "\n"
            l3 += " {}\n\n".format(i)
            text += l1 + l2 + l3
            l1, l2, l3 = '','',''

    return text

