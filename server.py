import re,blast
from multiprocessing import Pool,cpu_count
from tornado.ioloop import IOLoop
from tornado.web import RequestHandler, Application, url, asynchronous
from tornado.gen import coroutine



def blaster(seq):
    b = blast.blaster(seq)
    return b.full_analysis()

class infiniblaster():
    def __init__(self):
        pass
    def uptime():
        return 0

def sanitize(seq):
    """sanitize(str): convert fasta or bare sequence to bare sequence with no whitespace. returns a string"""
    seq = u''.join([re.sub(r"^s+", "", i.strip()) for i in seq.split(u'\n') if i[0] != '>']) #In case of FASTA
    return seq

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
            #self.pool.apply_async(blaster, (seq,), callback=callback)
            self.redirect('/sequenceprefs?sessionid=' + str(k))
        elif not is_sane(seq):
            self.get(header="Invalid sequences. Ensure all characters are amino acids")
        else:
            self.get()

class PreferenceHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    def get(self):
        sessionid = int(self.get_argument("sessionid"))
        seq = self.db[sessionid].seq
        if 'mutant' in self.request.arguments:
            mutants = self.sanitize_mutants(self.request.arguments['mutant'], seq)
        else:
            mutants = []
        deletes = [k for k,v in self.request.arguments.items() if v[0].lower() == 'delete']
        if len(deletes) == 1:
            delete = int(deletes[0])
            mutants = mutants[:delete] + mutants[delete+1:]

        if 'add' in self.request.arguments:
            mutants.append('')

        mutants = self.sanitize_mutants(mutants, seq)

        if len(mutants) == 0: #fuckall prevent returning an empty list
            mutants = ['']

        self.render('templates/userprefs.html',
                usersequence = seq,
                mutants=mutants,
                sessionid=sessionid,
                highlighted_markup=self.highlight(seq, mutants)
                )

    def sanitize_mutants(self, mutant_list, sequence):
        mutants = [i if i.isdigit() and int(i) <= len(sequence) and int(i) > 0 else '' for i in mutant_list]
        #remove duplicated DIGITS
        mutants = list(set([i for i in mutants if i.isdigit()])) + [i for i in mutants if not i.isdigit()]
        return mutants

    def highlight(self, seq, mutants):
        #DON'T LOOK AT ME!!! I'M HIDEOUS AWWWWW!~
        mutants = sorted(list(set([int(i) for i in mutants if i.isdigit()])))[::-1]
        mutants = [i for i in mutants if i <= len(seq)]
        print mutants
        markup = ''
        #from copy import deepcopy
        #currseq = deepcopy(seq)
        currseq = seq
        print "Highlighter is here! To save the markup!"
        for i in mutants:
            print currseq
            i = i-1 #zero indexing
            print i
            markup = u"<mark>" + currseq[i] + u"</mark>" + currseq[i+1:] + markup
            currseq = currseq[:i]
        markup = currseq + markup
        return markup

    def post(self):
        sessionid = int(self.get_argument("sessionid"))
        seq = self.db[sessionid].seq
        mutants = self.get_argument('mutant')
        self.write("sessionid:{}\nseq:{}\nmutants:{}\n".format(sessionid, seq, mutants))
        self.flush()

RID_DB = {0:infiniblaster()}

application = Application([
    (r"/", MainHandler, {'DB': RID_DB}),
    (r"/sequenceprefs", PreferenceHandler, {'DB' : RID_DB}),
])

if __name__ == "__main__":
    application.listen(8888)
    IOLoop.instance().start()
