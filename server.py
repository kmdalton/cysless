import re,blast
import redis
from uuid import uuid4
from bs4 import BeautifulSoup
from time import sleep,time
from tornado.ioloop import IOLoop
from tornado import gen
from tornado.web import RequestHandler, Application

blast_polling_period = 60 #Number of seconds to wait between blast queries -- minimum sixty seconds according to the blast docs
blast_rid_lifetime = 24*60*60 #Cache results for 24 hours -- blast says it caches for approximately 36 hours fwiw
redis_url  = "localhost"
redis_port = 6379
numhits = 3000 #Number of blast hits to ask for. During production this should be 20000

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
            rid, waittime = h.request(HITLIST_SIZE = numhits)
            value = """<status>provisional</status>
            <waittime>{}</waittime>
            <uptime>{}</uptime>
            <rid>{}</rid>
            <sequence>{}</sequence>""".format(waittime, time(), rid, seq)

            uid = uuid4()
            while uid in self.db:
                uid = uuid4()
            self.db.setex(uid, value, blast_rid_lifetime) #I think we need to maintain some state here to not be evil
            #self.("/blast/{}".format(uid))
            self.render("templates/waiting.html", uid=uid)
        elif not is_sane(seq):
            self.get(header="Invalid sequence. Ensure all characters are amino acids")
        else:
            self.get()

class BlastHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    @gen.coroutine
    def get(self, uid):
        if uid in self.db:
            soup = BeautifulSoup(self.db[uid], features='lxml')
            if soup.status is not None:
                rid = soup.rid.text
                handle = blast.blast_handle(soup.sequence.text)
                handle.rid = rid
                uptime = float(soup.uptime.text)
                yield gen.sleep(blast_polling_period)
                if time() - uptime > blast_polling_period:
                    print("Checking status for UID: {}".format(uid))
                    status = handle.check_status()
                    print("Status {} for UID: {}".format(status, uid))
                    if status == True:
                        self.db[uid] = handle.fetch_result() + "\n<sequence>{}</sequence>".format(soup.sequence.text)
                        self.redirect("/sequence/{}".format(uid))
                    elif status == False:
                        soup.uptime.string = str(time())
                        self.db[uid] = str(soup)
                        self.redirect("/blast/{}".format(uid))
                    else:
                        raise TypeError("blast_handle.check_status returned a variable of type: {}".format(type(status)))
                else:
                    self.get(uid)
            else:
                self.redirect('/sequence/{}'.format(uid))
        else:
            self.redirect('/')

class SequenceHandler(RequestHandler):
    def initialize(self, **kw):
        self.db = kw['DB']

    def get(self, uid):
        if uid in self.db:
            soup = BeautifulSoup(self.db[uid], features='lxml')
            if soup.status is None:
                sequence = soup.sequence.text
                self.render("templates/userprefs.html", sequence=sequence, uid=uid)
            else:
                self.redirect("/blast/{}".format(uid))
        else:
            self.redirect("/".format(uid))

    def post(self, uid):
        if uid in self.db:
            results = blast.blast_results(self.db[uid])
            if results.soup.status is None:
                try:
                    mutants = self.get_arguments('mutant')
                    mutants = list(map(int, mutants))
                    hit = results.recommend_mutant(mutants)
                    print(hit)
                    print(type(hit))
                    if hit is None:
                        self.render("templates/nohit.html")
                    self.render("templates/hit.html", hit=hit, mutants=mutants, uid=uid)
                except:
                    self.redirect("/sequence/{}".format(uid))
                
            else:
                self.redirect("/blast/{}".format(uid))
        else:
            self.redirect("/")

RID_DB = redis.Redis(host=redis_url, port=redis_port, db=0)
RID_DB.set('debug', open('blast_results.xml').read())

application = Application([
    (r"/", MainHandler, {'DB': RID_DB}),
    (r"/blast/(.*)", BlastHandler, {'DB': RID_DB}),
    (r"/sequence/(.*)", SequenceHandler, {'DB' : RID_DB}),
])

if __name__ == "__main__":
    #Don't be a jerk error
    if blast_polling_period < 60:
        raise ValueError("blast_polling_period set to {}. This must be at least sixty seconds to comply with the BLAST terms of service".format(blast_polling_period))

    application.listen(8889)
    IOLoop.instance().start()
