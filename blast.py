###############################################################################
#                                                                             #
# Module for retrieving blast hits.                                           #
#                                                                             #
###############################################################################

import urllib2,re,subprocess,requests,datetime
from BeautifulSoup import BeautifulStoneSoup

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

        PutKwargs = ['AUTO_FORMAT', 'COMPOSITION_BASED_STATISTICS', 'DATABASE', 'DB_GENETIC_CODE', 'ENDPOINTS', 'ENTREZ_QUERY', 'EXPECT', 'FILTER', 'GAPCOSTS', 'GENETIC_CODE', 'HITLIST_SIZE', 'I_THRESH', 'LAYOUT', 'LCASE_MASK', 'MATRIX_NAME', 'NUCL_PENALTY', 'NUCL_REWARD', 'OTHER_ADVANCED', 'PERC_IDENT', 'PHI_PATTERN', 'PROGRAM', 'QUERY', 'QUERY_FILE', 'QUERY_BELIEVE_DEFLINE', 'QUERY_FROM', 'QUERY_TO', 'SEARCHSP_EFF', 'SERVICE', 'THRESHOLD', 'UNGAPPED_ALIGNMENT', 'WORD_SIZE']

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
            text = ncbiGet(self.rid, ALIGNMENTS='0', DESCRIPTIONS='0', FORMAT_TYPE='Text')
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
        BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&"
        GetKwargs = ['ALIGNMENTS', 'ALIGNMENT_VIEW', 'DESCRIPTIONS', 'ENTREZ_LINKS_NEW_WINDOW', 'EXPECT_LOW', 'EXPECT_HIGH', 'FORMAT_ENTREZ_QUERY', 'FORMAT_OBJECT', 'FORMAT_TYPE', 'NCBI_GI', 'RID', 'RESULTS_FILE', 'SERVICE', 'SHOW_OVERVIEW']

        status = self.check_status()
        if status != True:
            raise TypeError("blast_handle.check_status() = {}. Status must be True if the results are ready to be downloaded from NCBI.".format(status))

        for kwarg in kw:
            if kwarg not in GetKwargs:
                raise TypeError("blast_handle.fetch_result got an unexpected keyword argument {}".format(kwarg))
        kw['RID'] = str(RID)
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

class alignment():
    """
    """
    def __init__(self, XML):
        """
        """
        pass
    def 



class legacy_blast_handle():
    """
    blast_handle(seq, **kwargs)

    This class wraps a variety of calls the NCBI
    BLAST common url api (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo) 
    and
    NCBI eutils (https://www.ncbi.nlm.nih.gov/books/NBK25501/). 

    It can be used to perform a BLAST search and to retrieve the full sequences of the hits.
    """
    def __init__(self, seq):
        """Create blast handle with a query sequence

        Parameters
        ----------
        seq : string
            the query sequence for your blast search
        """
        self.seq     = seq
        self.rid     = None
        self.seqs    = None
        self.uids    = None
        self.soup    = None 
        self.status  = None
        self.numhits = None

    def uptime(self):
        """Check the uptime of a blast_handle instance. A server might do this to nuke old searches.
        This relies on the datetime library, and is not intended to be particularly portable

        Returns
        -------
        int
            the number of seconds since the blast_handle was instantiated
        """
        return (datetime.datetime.now() - self.starttime).seconds

    def blast_request(self, **kw):
        """
        Submit your query sequence to the BLAST server. This will set two instance variables:
            self.numhits -- the requested HITLIST_SIZE for the the search
            self.rid -- request ID from the BLAST server

        Parameters
        ----------
        **kwargs

        Keyword Arguments
        -----------------
        HITLIST_SIZE : int=3000
            max number of seqs to ask blast for
        DATABASE : str='nr'
            which blast database to search
        PROGRAM : str='blastp' 
            which blast program to use, default blastp
        """
        #Make sure default list of parameters is populated
        kw['HITLIST_SIZE'] = str(kw.get('HITLIST_SIZE', 20000))
        kw['DATABASE'] = kw.get('DATABASE', 'nr')
        kw['PROGRAM'] = kw.get('PROGRAM', 'blastp')
        try:
            self.rid = ncbiPut(self.seq, **kw)[0]
        except:
            raise ConnectivityError('Unable to reach blast server')
        self.numhits = int(kw['HITLIST_SIZE'])

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
            text = ncbiGet(self.rid, ALIGNMENTS='0', DESCRIPTIONS='0', FORMAT_TYPE='Text')
            if '<!DOCTYPE html' in text.split('\n')[0]:
                self.status = False
                return False
            else:
                self.status = True
                return True

    def get_hitlist(self):
        """
        Download the blast results and store them in self.soup as a BeautifulStoneSoup object. 
        Store the accession numbers for the hits in self.uids (list)
        This function requires that self.status == True
        
        """
        self.check_status()
        if self.status is None:
            pass
        elif self.status == False:
            pass
        elif self.status == True:
            text = ncbiGet(self.rid, DESCRIPTIONS='{}'.format(self.numhits), FORMAT_TYPE='XML')
            self.soup = sorted(BeautifulStoneSoup(text).findAll('hit'), key = lambda x: -float(x.hsp_score.text))
            self.uids = [i.hit_accession.text for i in self.soup]
            #Add the aligned hit sequences to self.seqs
            self.seqs = []
            for hit in self.soup:
                qseq = hit.hsp_qseq.text
                hseq = hit.hsp_hseq.text
                #Remove query sequence gaps from the alignment
                qseq,hseq = zip(*((i,j) for i,j in zip(qseq,hseq) if i!='-'))
                qseq,hseq = ''.join(qseq),''.join(hseq)

                #Number of gaps to pad the ends of the hit sequences
                lpad = int(hit.find("hsp_query-from").text) - 1 
                rpad = len(self.seq) - int(hit.find("hsp_query-to").text)
                self.seqs.append('-'*lpad + hseq + '-'*rpad)
        else:
            raise TypeError('The type of blaster.status must be True,False, or None. Current type is: {}'.format(type(self.status)))

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

    def recommend_mutant(self, residues):
        """
        blaster.recommend_mutant(residues)
        Parameters
        ----------
        residues: An iterable containing integer residue numbers to be mutated.
        Returns
        -------
        smith_waterman object or None: smith_waterman object of the most closely related sequence which is mutated at the requested positions or None if no suitable sequence exists
        """
        if self.seqs is None:
            return None
        else:
            WT = [self.seq[i-1] for i in residues]
            hits = []
            for hseq in self.seqs:
                MUT = [hseq[i-1] for i in residues]
                if '-' not in MUT:
                    if True not in [x == y for x,y in zip(WT,MUT)]:
                        hits.append(hseq)
            if len(hits) > 0.:
                return hits[0]
            else:
                return None


    def __exit__(self):
        if self.rid != None:
            ncbiDelete(self.rid)

def ncbiPut(seq, **kw):
    """
    webBlast(seq, **kwargs)
        Make a call to the ncbi webserver to run BLAST.

    Parameters
    ----------
    seq : Input protein or nucleic acid sequence

    kwargs : { AUTO_FORMAT, COMPOSITION_BASED_STATISTICS, DATABASE, DB_GENETIC_CODE, ENDPOINTS, ENTREZ_QUERY, EXPECT, FILTER, GAPCOSTS, GENETIC_CODE, HITLIST_SIZE, I_THRESH, LAYOUT, LCASE_MASK, MATRIX_NAME, NUCL_PENALTY, NUCL_REWARD, OTHER_ADVANCED, PERC_IDENT, PHI_PATTERN, PROGRAM, QUERY, QUERY_FILE, QUERY_BELIEVE_DEFLINE, QUERY_FROM, QUERY_TO, SEARCHSP_EFF, SERVICE, THRESHOLD, UNGAPPED_ALIGNMENT, WORD_SIZE }
    
        kwargs are verbatim those supported by the blast REST API with CMD=Put. The values will default to a vanilla blastp search which returns fasta formatted sequences with gids for headers. (http://www.ncbi.nlm.nih.gov/blast/Doc/node5.html)

    Returns
    -------
    Tuple (RID\string, WAITTIME\int) 
        RID: blast request ID (RID) which will be used to retrieve the results later with a Get request
        WAITTIME: blast returns an estimate of the amount of time in seconds the search will take
    """
    BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&"
    PutKwargs = ['AUTO_FORMAT',
                'COMPOSITION_BASED_STATISTICS',
                'DATABASE',
                'DB_GENETIC_CODE',
                'ENDPOINTS',
                'ENTREZ_QUERY',
                'EXPECT',
                'FILTER',
                'GAPCOSTS',
                'GENETIC_CODE',
                'HITLIST_SIZE',
                'I_THRESH',
                'LAYOUT',
                'LCASE_MASK',
                'MATRIX_NAME',
                'NUCL_PENALTY',
                'NUCL_REWARD',
                'OTHER_ADVANCED',
                'PERC_IDENT',
                'PHI_PATTERN',
                'PROGRAM',
                'QUERY',
                'QUERY_FILE',
                'QUERY_BELIEVE_DEFLINE',
                'QUERY_FROM',
                'QUERY_TO',
                'SEARCHSP_EFF',
                'SERVICE',
                'THRESHOLD',
                'UNGAPPED_ALIGNMENT',
                'WORD_SIZE'
                ]
    kw['QUERY'] = seq
    kw['HITLIST_SIZE'] = kw.get('HITLIST_SIZE', 20000)
    kw['DATABASE'] = kw.get('DATABASE', 'nr')
    kw['PROGRAM'] = kw.get('PROGRAM', 'blastp')

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

def ncbiGet(RID, **kw):
    """
    make an ncbi get request

    Parameters
    ----------
    RID : the request ID you wish to query blast about
    kwargs : { ALIGNMENTS, ALIGNMENT_VIEW, DESCRIPTIONS, ENTREZ_LINKS_NEW_WINDOW, EXPECT_LOW, EXPECT_HIGH, FORMAT_ENTREZ_QUERY, FORMAT_OBJECT, FORMAT_TYPE, NCBI_GI, RID, RESULTS_FILE, SERVICE, SHOW_OVERVIEW }
        kwargs are verbatim those supported by the blast REST API with CMD=Get. (http://www.ncbi.nlm.nih.gov/blast/Doc/node6.html)
    """
    BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&"
    GetKwargs = ['ALIGNMENTS',
                 'ALIGNMENT_VIEW',
                 'DESCRIPTIONS',
                 'ENTREZ_LINKS_NEW_WINDOW',
                 'EXPECT_LOW',
                 'EXPECT_HIGH',
                 'FORMAT_ENTREZ_QUERY',
                 'FORMAT_OBJECT',
                 'FORMAT_TYPE',
                 'NCBI_GI',
                 'RID',
                 'RESULTS_FILE',
                 'SERVICE',
                 'SHOW_OVERVIEW'
                ]
    kw['RID'] = str(RID)
    QueryString = '&'.join(['='.join((i, str(kw[i]))) for i in GetKwargs if i in kw])
    #print BaseURL + QueryString
    URL = urllib2.urlopen(BaseURL + QueryString)
    return URL.read()


def ncbiDelete(RID):
    """
    delete a qblast RID from NCBI's servers

    Parameters
    ----------
    RID : the request ID you wish to delete
    """
    BaseURL = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Delete&"
    URL = urllib2.urlopen(BaseURL + "RID=%s" %RID)
    URL.read()

def efetch(uids, **kw):
    kw['db'] = kw.get('db', 'protein')
    kw['rettype'] = kw.get('rettype', 'fasta')
    kw['retmode'] = kw.get('retmode', 'text')
    kw['id'] = ','.join(uids)
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
