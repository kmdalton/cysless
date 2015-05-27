import re,blast,pickle
from multiprocessing import Pool,cpu_count
from tornado.ioloop import IOLoop
from tornado.web import RequestHandler, Application, url, asynchronous
from tornado.gen import coroutine

dummy_blaster_filename = 'dummy_blaster.pickle'


def blaster(seq):
    b = blast.blaster(seq)
    return b.full_analysis()

class infiniblaster(blast.blaster):
    def uptime(self):
        return 0

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
        NUMPROCS = cpu_count() - 1 | 1
        self.pool = Pool(NUMPROCS)

    def get(self, **kw):
        header = kw.get('header', "Please enter your amino acid sequence below:")
        self.render('templates/frontpage.html', header = header)

    def post(self):
        seq = self.get_argument("usersequence")
        seq = sanitize(seq)
        if is_sane(seq):
            k  = max(self.db.keys()) + 1
            self.db[k] = blast.blaster(seq) #This may cause bugs. The right thing would be to put a lock on or use a real db
            def callback(blaster_obj):
                self.db[int(k)] = blaster_obj
            self.pool.apply_async(blaster, (seq,), callback=callback)
            self.redirect('/sequenceprefs?sessionid=' + str(k))
        elif not is_sane(seq):
            self.get(header="Invalid sequences. Ensure all characters are amino acids")
        else:
            self.get()

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
        print mutants
        markup = ''
        currseq = seq
        #print "Highlighter is here! To save the markup!"
        for i in mutants:
            print currseq
            i = i-1 #zero indexing
            print i
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
            SW = self.db[sessionid].recommend_mutant(mutants)
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
            message = '<h3> Sorry, your BLAST results are not ready yet. Please wait a minute and resubmit. </h3>'
            self.get(message = message)


RID_DB = {0: pickle.load(open(dummy_blaster_filename))}
#print RID_DB

application = Application([
    (r"/", MainHandler, {'DB': RID_DB}),
    (r"/sequenceprefs", PreferenceHandler, {'DB' : RID_DB}),
])

if __name__ == "__main__":
    application.listen(8888)
    IOLoop.instance().start()
