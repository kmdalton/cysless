###############################################################################
#                                                                             #
# Module for retrieving blast hits.                                           #
#                                                                             #
###############################################################################

import urllib2,re,subprocess,requests,datetime
from time import sleep

class ConnectivityError(Exception):
    pass

class blaster():
    def __init__(self, seq = None):
        """
        var = blaster(seq) --> blaster object
        """
        if seq is None:
            self.seq = ''
        else:
            self.seq = seq
        self.starttime = datetime.datetime.now()
        self.rid = None
        self.status = None #None=>RID not set;False=>Search in progress;True=>Search Complete
        self.numhits= None
        self.headers= None
        self.seqs   = None
        self.uids   = None
        #self.blast_request()

    def full_analysis(self):
        """
        blaster.full_analysis()
        makes external blast call. downloads sequence hits and makes the pairwise smith-waterman 
        alignments between the initial seqs and the hits.
        --------
        Returns:
            self (blaster object)
        """
        self.blast_request()
        while self.check_status() is False:
            sleep(20)
        self.get_hitlist()
        self.make_pairwise_alignments()
        return self

    def uptime(self):
        """blaster.uptime() - return lifetime of a blaster object in seconds"""
        return (datetime.datetime.now() - self.starttime).seconds

    def blast_request(self, **kw):
        """
        blaster.blast_request(**kwargs)
            -------
            kwargs:
                HITLIST_SIZE (str|int): max number of seqs to ask blast for, default 3000
                DATABASE (str): which blast database to search, default nr
                PROGRAM (str): which blast program to use, default blastp
            --------
            returns:
                nothing. but it sets self.numhits, self.rid
        """
        #Make sure default list of parameters is populated
        kw['HITLIST_SIZE'] = str(kw.get('HITLIST_SIZE', 3000))
        kw['DATABASE'] = kw.get('DATABASE', 'nr')
        kw['PROGRAM'] = kw.get('PROGRAM', 'blastp')
        try:
            self.rid = ncbiPut(self.seq, **kw)[0]
        except:
            raise ConnectivityError('Unable to reach blast server')
        self.numhits = int(kw['HITLIST_SIZE'])

    def check_status(self):
        """
        blaster.check_status()
            query the ncbi database to check the status of self.rid which is the ncbi server
            request id for the blast search initiated by self.blast_request
        --------
        Returns:
            bool/None - returns True if the blast search is complete, False otherwise, and None if self.rid is unset
                        also sets the instance variable self.status to the same value
        """
        if self.rid is None:
            self.status = None
            return None
        else:
        #Now we check to see if the server will return formatted results
        #This test relies on the idea that qblast only returns FORMAT_TYPE==Text
        #if the search is __complete__. Otherwise, we get an html waiting page
            text = ncbiGet(self.rid, ALIGNMENTS='0', DESCRIPTIONS='0', FORMAT_TYPE='Text')
            if '<!DOCTYPE html' in text.split('\n')[0]:
                self.status = False
                return False
            else:
                self.status = True
                return True

    def get_hitlist(self):
        self.check_status()
        if self.status is None:
            return None,None
        elif self.status == False:
            return None,None
        elif self.status == True:
            text = ncbiGet(self.rid, ALIGNMENTS='0', DESCRIPTIONS='{}'.format(self.numhits), FORMAT_TYPE='Text')
            self.uids = [i.split('|')[1] for i in text.split('\n') if '|' in i and len(i.split('|')) == 3]
            self.uids = [i for i in self.uids if len(i)>3]#Shameful hack ... TODO Fix later
            self.fetch_sequences()
        else:
            raise TypeError('The type of blaster.status must be True,False, or None. Current type is: {}'.format(type(self.status)))

    def fetch_sequences(self):
        if self.uids is not None:
            fasta= efetch(self.uids) #flat text fasta format is the default for efetch
            self.headers,self.seqs = zip(*(('>'+i.split('\n')[0], ''.join(i.split('\n')[1:])) for i in fasta.split('>')[1:]))

    def make_pairwise_alignments(self):
        self.alignments = [smith_waterman(self.seq, i) for i in self.seqs]

    def recommend_mutant(self, residues):
        if self.alignments is None:
            return None
        else:
            WT = [self.seq[i-1] for i in residues]
            hits = []
            for a in self.alignments:
                MUT = [a.registered_seq2[i-1] for i in residues]
                if '-' not in MUT:
                    if True not in [x == y for x,y in zip(WT,MUT)]:
                        hits.append(a)
            hits.sort(key = lambda x: x.identity)
            if len(hits) > 0.:
                return hits[-1]
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

def blastFormatter(RID, **kw):
    outfmt   = kw.get('outfmt', '6 evalue sgi staxids sseq')
    outfmt   = "'%s'" %outfmt
    arguments = ["blast_formatter","-rid",str(RID),"-outfmt",outfmt]
    p = subprocess.Popen(' '.join(arguments), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    #p = subprocess.Popen(' '.join(arguments), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    lines = p.communicate()
    return lines
    #lines = lines.split("\n")
    #return lines[:-1] #The last line is empty

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
    def __init__(self, seq1, seq2):
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

