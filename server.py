import re,blast,pickle
from time import sleep
from multiprocessing import Pool,cpu_count
from tornado.ioloop import IOLoop
from tornado.web import RequestHandler, Application, url, asynchronous
from tornado.gen import coroutine

blast_polling_period = 60 #Number of seconds to wait between blast queries -- minimum sixty seconds according to the blast docs
blast_rid_lifetime = 24*60*60 #Cache results for 24 hours -- blast says it caches for approximately 36 hours fwiw
redis_url  = "localhost"
redis_port = 6379

def sanitize(seq):
    """sanitize(str): convert fasta or bare sequence to bare sequence with no whitespace. returns a string of upper case letters"""
    seq = u''.join([re.sub(r"^s+", "", i.strip()) for i in seq.split(u'\n') if i[0] != '>']) #In case of FASTA
    return seq.upper()

def is_sane(seq):
    """is_sane(str): are all characters in str.upper() amino acids, return True or False"""
    if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq.upper()) is None:
        return True
    else:
        return False

class MainHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    def get(self, **kw):
        header = kw.get('header', "Please enter your amino acid sequence below:")
        self.render('templates/frontpage.html', header = header)

    def post(self):
        seq = self.get_argument("usersequence")
        seq = sanitize(seq)
        if is_sane(seq):
            h = blast.blast_handle(seq)
            rid, waittime = h.request()
            self.db[rid] = "PROVISIONAL str(waittime)" #I think we need to maintain some state here to not be evil
            self.redirect("/blast/{}".format(rid))
        elif not is_sane(seq):
            self.get(header="Invalid sequence. Ensure all characters are amino acids")
        else:
            self.get()

class BlastHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    def get(self, rid):
        if rid in self.db:
            self.redirect("/sequence/{}".format(rid))
            handle = blast.blast_handle()
            handle.rid = rid
            if handle.check_status() != True:
                sleep(blast_polling_period)
                self.get(rid)
            elif handle.check_status == True:
                xml = handle.
        else:
            self.rdirect("/")

    def post(self)
        sleep(blast_polling_period)

class PreferenceHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    def get(self, **kw):
        sessionid = int(self.get_argument("sessionid"))
        if sessionid in self.db:
            seq = self.db[sessionid].seq
            if 'mutant' in self.request.arguments:
                mutants = self.request.arguments['mutant']
            else:
                mutants = []

            #Remove mutants slated for deletion
            deletes = [int(k) for k,v in self.request.arguments.items() if v[0].lower() == '-' and k.isdigit()]
            if len(deletes) == 1: #must be in {0,1}
                delete = deletes[0]
                mutants = mutants[:delete] + mutants[delete+1:]

            mutant_validity = [int(i) in range(1, len(seq)+1) if i.isdigit() else False for i in mutants]
            if False in mutant_validity:
                message = "<font color='firebrick'>*Invalid residue numbers detected</font>"
            else:
                message = kw.get('message', '')

            #If '+' was clicked, add an empty mutant; also never render an empty list
            if 'add' in self.request.arguments or len(mutants) == 0:
                mutants.append('')
                mutant_validity.append(True)

            self.render('templates/userprefs.html',
                    usersequence = seq,
                    mutants=mutants,
                    validity=mutant_validity,
                    sessionid=sessionid,
                    highlighted_markup=self.highlight(seq, mutants),
                    message = message, 
                    )
        else:
            self.render('templates/frontpage.html', header='Invalid session ID. Please enter a new sequence')

    def highlight(self, seq, mutants):
        #DON'T LOOK AT ME!!! I'M HIDEOUS AWWWWW!~
        mutants = sorted(list(set([int(i) for i in mutants if str(i).isdigit()])))[::-1]
        mutants = [i for i in mutants if i <= len(seq)]
        markup = ''
        currseq = seq
        for i in mutants:
            i = i-1 #zero indexing
            markup = u"<mark>" + currseq[i] + u"</mark>" + currseq[i+1:] + markup
            currseq = currseq[:i]
        markup = currseq + markup
        return markup

    def post(self, **kw):
        sessionid = int(self.get_argument("sessionid"))
        mutants = self.request.arguments['mutant']
        if self.db[sessionid].check_status():
            seq = self.db[sessionid].seq
            mutants = [int(i) for i in mutants if str(i).isdigit()]
            SW = self.db[sessionid].recommend_mutant(mutants) #There may be dragons in here
            if SW is not None:
                self.render('templates/results.html',
                        usersequence = seq,
                        mutants=mutants,
                        sessionid=sessionid,
                        source_seq=SW.registered_seq2,
                        source_header=SW.header2,
                        highlighted_markup=self.highlight(seq, mutants),
                        message = kw.get('message', ''),
                        )
            else:
                message = ''
                if len(mutants) > 1:
                    message = '<h3> No suitable match found. Try decreasing the number of simultaneous mutations. </h3>'
                else:
                    message = '<h3> No suitable match found. </h3>'
                self.render('templates/userprefs.html',
                        usersequence = seq,
                        mutants=mutants,
                        sessionid=sessionid,
                        validity = [i in range(1, len(seq)+1) for i in mutants],
                        highlighted_markup=self.highlight(seq, mutants),
                        message = kw.get('message', message),
                        )
        else:
            message = '<h3> Sorry, your BLAST results are not ready yet. Queries can take up to 15 minutes depending on load. Please wait several minutes and resubmit this form. </h3>'
            self.get(message = message)

RID_DB = redis.Redis(host=redis_url, port=redis_port, db=0)
RID_DB.set('debug', open('blast_results.xml').read())

application = Application([
    (r"/", MainHandler, {'DB': RID_DB}),
    (r"/blast/(.*)", BlastHandler, {'DB': RID_DB}),
    (r"/sequence/(.*)", PreferenceHandler, {'DB' : RID_DB}),
])

if __name__ == "__main__":
    #Don't be a jerk error
    if blast_polling_period < 60:
        raise ValueError("blast_polling_period set to {}. This must be at least sixty seconds to comply with the BLAST terms of service".format(blast_polling_period))

    application.listen(8888)
    IOLoop.instance().start()
