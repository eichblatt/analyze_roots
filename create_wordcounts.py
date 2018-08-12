# -*- coding: utf-8 -*- 
import os,io,re,sys,pdb,itertools,argparse,pickle
import pandas as pd
import numpy as np
import atomigi
from lxml import etree
from datetime import datetime
from collections import OrderedDict
from matplotlib import pyplot as plt
from sklearn.cluster.bicluster import SpectralBiclustering

xmlparser=etree.XMLParser(dtd_validation=True,no_network=False)

reload(sys)  
sys.setdefaultencoding('utf-8')

def flattern(A):
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flattern(i))
        else: rt.append(i)
    return rt

def get_text(f):
# p=etree.parse(f,xmlparser)
 t=[]
 aa=etree.iterparse(f)
 for action,elem in aa:
    if re.match(r'.*p$',elem.tag): t.append(elem.text)
 return t

def canonize(wl):
  pattern=re.compile('[\W_]+',re.UNICODE)    
  result=[]
  if wl is None:return result
  for word in wl:
    w=word.lower()
    if re.match(r'.*\'',w): w=w+'o'
    w=pattern.sub('',w)
    if len(w)>0: result.append(w)
  return result

def get_word_df(textpath): 
  wordlist=[]
  for filename in os.listdir(textpath):
    print 'reading file',filename
    try: wlist=get_text(os.path.join(textpath,filename))
    except: pass
    wlist=[canonize(w) for w in wlist]
    wordlist.extend([item for sublist in wlist for item in sublist])  
  word_df=pd.DataFrame({'word':wordlist})
  return word_df

def get_wcdf():
  ''' get the word count dataframe. This function needs to be sped up. Probably best to cache the table.'''
  wcdfpath=os.path.join(args.basepath,args.wordcount_name)
  if ((os.path.exists(wcdfpath)) and not args.force): return pickle.load(open(wcdfpath,'rb'))

  print '{}: leganta tekstojn...'.format(datetime.now().strftime('%H:%M:%S'))
  word_df=get_word_df(os.path.join(args.textpath,'tekstoj'))
  print '{}: kreanta vortkvanton...'.format(datetime.now().strftime('%H:%M:%S'))

  word_counts=word_df['word'].value_counts()
  wcdf=pd.DataFrame(word_counts)
  wcdf=wcdf.rename(columns={'word':'N'})
  wcdf=wcdf[wcdf.N>=args.min_count]
  atomigitaj=[atomigi.vorto(x) for x in wcdf.index]  # this is the slow part. 
  aaa=pd.DataFrame([(x.vorto,x.radik,x.prefs,x.sufs,x.fins,x.prefs+x.sufs+x.fins,x.valid,x.__repr__()) for x in atomigitaj],columns=['word','root','prefs','sufs','fins','mods','valid','repr'])
  wcdf=pd.concat([wcdf,aaa],axis=1).sort_values(by='N',ascending=False)
  wcdf['root_count']=wcdf['N'].groupby(wcdf['root']).transform('sum')
  wcdf['root_freq']=wcdf['N']/wcdf['root_count']
  wcdf.loc[wcdf.index==wcdf.root,'mods']=[u'NENIU']
  wcdf['modstring']=[''.join(x) for x in wcdf['mods']]
  pickle.dump(wcdf,open(args.wcdfpath,'wb'))
  return wcdf

def main(args):
  for k in sorted(args.__dict__.keys()): print '{0:20} : {1}'.format(k,args.__dict__[k])
  pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])

  wcdf0=get_wcdf()

  return(wcdf0)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-basepath", dest='basepath', default='/home/steve/esperanto/analyze_roots',help="top folder of corpus")
    argparser.add_argument("-textpath", dest='textpath', default='/home/steve/esperanto/tekstaro_de_esperanto_2017-07-25',help="top folder of corpus")
    argparser.add_argument("-wordcount_name", dest='wordcount_name', default='word_counts.pkl',help="word count dataframe path")
    argparser.add_argument("-force", dest='force', default=False,type=bool,help="create word count dataframe, even if it already exists")
    argparser.add_argument("-min_count", dest='min_count', default=3,type=int,help="minimum number of appearances for a word to count")
    argparser.add_argument("-topn", dest='topn', default=None,type=int,help="number of roots to consider (ranked by frequency)")

    args = argparser.parse_args()
    import pdb; pdb.set_trace()
    wcdf=main(args)
