import re,blast
from tornado.ioloop import IOLoop
from tornado.web import RequestHandler, Application, url, asynchronous

class infiniblaster():
    def __init__(self):
        pass
    def uptime():
        return 0

RID_DB = {0:infiniblaster()}

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

    def get(self, **kw):
        header = kw.get('header', "Please enter your amino acid sequence below:")
        self.render('templates/frontpage.html', header = header)

    def post(self):
        seq = self.get_argument("usersequence")
        seq = sanitize(seq)
        if is_sane(seq):
            k  = max(self.db.keys()) + 1
            self.db[k] = blast.blaster(seq) 
            self.blast(k) 
            self.redirect('/sequenceprefs?sessionid=' + str(k))
        elif not is_sane(seq):
            self.get(header="Invalid sequences. Ensure all characters are amino acids")
        else:
            self.get()

    @asynchronous #TODO: make this actually run asynchronously
    def blast(self, k):
        b = blast.blaster(self.db[k].seq)
        self.db[k] = b.full_analysis()

class PreferenceHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    def get(self):
        sessionid = int(self.get_argument("sessionid"))
#TODO: add + delete mutants fix form
        mutants = {k:v for k,v in self.request.arguments.items() if 'mutant' in k}
        deletes = {k:v for k,v in self.request.arguments.items() if 'delete' in k}
        seq = self.db[sessionid].seq
        self.render('templates/userprefs.html', usersequence = seq, mutants=[0])

application = Application([
    (r"/", MainHandler, {'DB': RID_DB}),
    #(r"/blastout", RIDHandler, {'DB' : RID_DB}),
    (r"/sequenceprefs", PreferenceHandler, {'DB' : RID_DB}),
])

if __name__ == "__main__":
    application.listen(8888)
    IOLoop.instance().start()
