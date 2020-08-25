#Copyright (C) Systems & Technology Research
#Use of this software is subject to the restrictions in license.txt
#Distribution A: Approved for Public Release, Distribution Unlimited
import os
import sys
import traceback
import requests
import json
import time

from bs4 import BeautifulSoup
from scholarmetrics import hindex
from Bio import Entrez
from ratelimit import limits, sleep_and_retry
import threading
import _thread as thread
from tqdm import tqdm


tool = "biopython"
Entrez.tool = tool
Entrez.sleep_between_tries = 60
Entrez.max_tries = 5


def quit_function(fn_name):
    print('{0} took too long'.format(fn_name), file=sys.stderr)
    sys.stderr.flush() # Python 3 stderr is likely buffered.
    thread.interrupt_main() # raises KeyboardInterrupt


def exit_after(s):
    '''
    use as decorator to exit process if
    function takes longer than s seconds
    '''
    def outer(fn):
        def inner(*args, **kwargs):
            timer = threading.Timer(s, quit_function, args=[fn.__name__])
            timer.start()
            try:
                result = fn(*args, **kwargs)
            except KeyboardInterrupt as e:
                raise e
            finally:
                timer.cancel()
            return result
        return inner
    return outer


def get_links_ids(pmids):
    link_db = {}

    while True:
        try:
            links = Entrez.elink(dbfrom="pubmed", cmd='neighbor_score', id=pmids, linkname="pubmed_pmc_refs")

            record = Entrez.read(links)
            break
        except:
            print("Error... Waiting 2 Minutes before re-trying")

            print('-' * 60)
            traceback.print_exc(file=sys.stdout)
            print('-' * 60)

            time.sleep(120)

    for rr in record:
        id_list = {int(i) for i in rr[u'IdList']}
        link_set_db = rr[u'LinkSetDb']
        if len(link_set_db) > 0:
            for link in link_set_db[0][u'Link']:
                l_score = int(link[u'Score'])
                l_id = int(link[u'Id'])

                if l_score not in id_list:
                    print("ERROR: score doesn't match id set")

                if l_score not in link_db:
                    link_db[l_score] = []

                link_db[l_score].append(l_id)

    return link_db


@sleep_and_retry
@limits(calls=30, period=10)
def pull_url(author, email, apikey=None):
    if apikey is not None:
        req = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&datetype=pdat&mindate=1990' \
          f'&maxdate=2020&retmax=500&tool={tool}&term={author}&email={email}&apikey={apikey}'
    else:
        req = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&datetype=pdat&mindate=1990' \
              f'&maxdate=2020&retmax=500&tool={tool}&term={author}&email={email}'

    while True:
        try:
            return requests.get(req)
        except:
            print("Error... Waiting 1 Minute before re-trying")

            print('-' * 60)
            traceback.print_exc(file=sys.stdout)
            print('-' * 60)

            time.sleep(60)
            continue


def calculate_h_index(df_authors, db_path, has_key=False):
    # Collect three stats: (i) author name and his/her h-index, (ii) citation list of each pmid, and (iii) author pmids

    author_2_hindex = dict()
    author_2_hindex_return = dict()
    pmid_2_cite = dict()
    author_2_pmids = dict()
    h_index_db = dict()

    with open(db_path) as h_index_db_file:
        h_index_db = json.load(h_index_db_file)
    
    for name in tqdm(df_authors):

        if name[0] == ' ':
            author = name[1:] # + ' ' + surname
        else:
            author = name # + ' ' + surname

        author_pmids = []


        # First check if the name already exists in our db, if 0, then try again
        call_api_flag = True

        if author in h_index_db['h_indices']:
            h_index = h_index_db['h_indices'][author]
            # print(f"{author} found in current database ({db_path}) with a value of {h_index}")
            if h_index != 0:
                # print(f"{author} found in current database ({db_path}) with a value of {h_index}")
                author_2_hindex_return[author] = h_index
                call_api_flag = False
            else:
                # print(f"{author} found in current database ({db_path}) with a value of {h_index} ... retrying API call")
                
                # BELOW IS JUST FOR DEMO!!!
                #TODO
                author_2_hindex_return[author] = h_index
                call_api_flag = False
        else:
            print(f"{author} not found in current database ({db_path})")
        # This ensures that we are not checking short and very common names which takes forever to collect information
        if len(author) > 5 and call_api_flag:
            if has_key:
                page = pull_url(author, Entrez.email, Entrez.api_key)
            else:
                page = pull_url(author, Entrez.email)

            soup = BeautifulSoup(page.content, 'xml')

            ids = soup.find_all('Id', {})

            for id_ in ids:
                author_pmids.append(id_.get_text())
            author_2_pmids[author] = author_pmids
            citations = []

            retrieved = get_links_ids([int(pmid) for pmid in author_pmids if int(pmid) not in pmid_2_cite.keys()])

            for pmid in retrieved.keys():
                link_list = []
                if pmid in retrieved:
                    link_list = retrieved[pmid]
                pmid_2_cite[pmid] = link_list
                citations.append(len(link_list))

            author_2_hindex[author] = int(hindex(citations))
            author_2_hindex_return[author] = author_2_hindex[author]
        elif call_api_flag:
            author_2_hindex[author] = -1
            author_2_hindex_return[author] = -1

    h_index_db['h_indices'].update(author_2_hindex)
    h_index_db['pmids'].update(author_2_pmids)
    h_index_db['citations'].update(pmid_2_cite)

    #TODO: make changes to path to overwrite current file
    with open(os.path.join(os.path.dirname(__file__), f'./data/author_h_indexes.json'), 'w') as output:
        output.write(json.dumps({
            'h_indices': h_index_db['h_indices'],
            'pmids': h_index_db['pmids'],
            'citations': h_index_db['citations']
        }))
    
    return author_2_hindex_return
    

# START HERE: Public Function
'''
Args:
email (string) -> valid email address for api calls
author_list (list) -> valid list of author names in the format "name surname"
db_path (string) -> optional path to current h_index database

Return:
author_2_hindex_return (dict) -> key: author name, value: h_index for authors requested in author_list
'''
def find_h_index(email, author_list, api_key=None, db_path=os.path.join(os.path.dirname(__file__), './data/author_h_indexes.json')):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        return calculate_h_index(author_list, db_path, True)
    else:
        return calculate_h_index(author_list, db_path)
