# REP-X
a python webapp for designing amino acid point mutants from homology



## This software has the following dependencies:
* python 3.6+
* beautifulsoup4 (https://www.crummy.com/software/BeautifulSoup/bs4/doc/)
* lxml is the preferred parser for beautifulsoup4 (https://lxml.de/)
* requests (http://docs.python-requests.org/)
* tornado  (http://www.tornadoweb.org/)
* redis  (https://redis.io/)
* water from the EMBOSS suite (http://emboss.sourceforge.net/)
* bash: the water interface is simply a shell call wrapped by the subprocess library. i don't expect this to work on a non-unix.

## Installation:
This has been tested on Anaconda Python 3.7 (https://www.anaconda.com/) on Ubuntu 18.04. Installing dependencies is easy with the conda package manager.
```bash
conda install requests tornado beautifulsoup4 lxml redis redis-py
```

On Ubuntu, EMBOSS is available through the package manager. 
```bash
sudo apt-get install emboss
```

## Running REP-X:
REP-x consists of two processes that need to run together. The backend is a Redis database server which is queried by the Python server. The Python servery is built with the Tornado web framework. In order to start REP-X, one must first start the Redis server by calling the redis-server command from the terminal. 
```bash
redis-server
```
This process needs to be running before the Python server is started. To start the REP-X server, open a new terminal and call
```bash
python server.py
```
. Then navigate to http://localhost:8889 in the browser to use the webapp. The port number can be changed at the bottom of the server.py file to suit your needs. 
