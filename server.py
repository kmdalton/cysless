import re
from tornado.ioloop import IOLoop
from tornado.web import RequestHandler, Application, url

def is_sane(seq):
    if re.search(r'[^ACDEFGHIKLMNPQRSTVWY]', seq.upper()) is None:
        return True
    else:
        return False

class MainHandler(RequestHandler):
    def get(self, **kw):
        header = kw.get('header', "Please enter your amino acid sequence below:")
        self.render('templates/frontpage.html', header = header)

    def post(self):
        seq = self.get_argument("usersequence")
        seq = seq.upper()
        if is_sane(seq):
            self.render('templates/userprefs.html', usersequence = seq)
        elif not is_sane(seq):
            self.get(header="Invalid sequences. Ensure all characters are amino acids")
        else:
            self.get()

class PreferenceHandler(RequestHandler):
    def get(self):
        seq = self.get_argument("usersequence")
        self.render('templates/userprefs.html', usersequence = seq)

application = Application([
    (r"/", MainHandler),
    (r"/sequenceprefs", PreferenceHandler),
])

if __name__ == "__main__":
    application.listen(8888)
    IOLoop.instance().start()
