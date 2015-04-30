from tornado import ioloop, web, template

loader = template.Loader("templates")

class MainHandler(tornado.web.RequestHandler):
    def get(self):
        

application = tornado.web.Application([
    (r"/", MainHandler),
])

if __name__ == "__main__":
    application.listen(8888)
    tornado.ioloop.IOLoop.instance().start()
