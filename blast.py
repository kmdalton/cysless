###############################################################################
#                                                                             #
# Module for retrieving blast hits.                                           #
#                                                                             #
###############################################################################

from time import sleep
import urllib2,re

class ConnectivityError(Exception):
    pass

class blaster():
    def __init__(self, seq = None):
        if seq is None:
            self.seq = ''
        else:
            self.seq = seq
        self.rid = None
        self.status = None #None=>RID not set;False=>Search in progress;True=>Search Complete
        self.numhits= None
        self.headers= None
        self.seqs   = None

    def blast_request(self, **kw):
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
        if self.rid is None:
            return None
        else:
        #Now we check to see if the server will return formatted results
        #This test relies on the idea that qblast only returns FORMAT_TYPE==Text
        #if the search is __complete__. Otherwise, we get an html waiting page
            text = ncbiGet(self.rid, ALIGNMENTS='0', DESCRIPTIONS='0', FORMAT_TYPE='Text')
            if '<!DOCTYPE html' in text.split('\n')[0]:
                return False
            else:
                return True

    def update_status(self):
        status = self.check_status()
        if status is None:
            self.status = None
        elif status == True:
            self.status = True
        else:
            self.status = False

    def get_hitlist(self):
        self.update_status()
        if self.status is None:
            return None,None
        elif self.status == False:
            return None,None
        elif self.status == True:
            text = ncbiGet(self.rid, ALIGNMENTS='0', DESCRIPTIONS='{}'.format(self.numhits), FORMAT_TYPE='Text')
            uids = [i.split('|')[1] for i in text.split('\n') if '|' in i and len(i.split('|')) == 3]
            fasta= efetch(uids) #flat text fasta format is the default for efetch
            headers,seqs = zip(*(('>'+i.split('\n')[0], ''.join(i.split('\n')[1:])) for i in fasta.split('>')))
            self.headers,self.seqs = headers,seqs
        else:
            raise TypeError('The type of blaster.status must be True,False, or None. Current type is: {}'.format(type(self.status)))

    def target_hit_pairwise_alignment():
        pass #TODO -- import alignment & registration functions from coupling_paper directory

    def fetch_seqs(self):
        hits = self.get_hitlist()
        uids = [i.split('|')[1] for i in hits]
        fasta= efetch(gids)

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
    kw  = '&'.join(['{}={}'.format(k,v)  for k,v in kw.items()])
    BaseURL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + kw
    data = ''
    for start in range(0, len(uids), 200):
        URL = urllib2.urlopen(BaseURL + "&id={}".format(','.join(uids[start:start+200])))
        data= data + URL.read()
    return data
    
