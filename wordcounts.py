# -*- coding: utf-8 -*- 
import os,io,re,sys,pdb,itertools,argparse,pickle
import pandas as pd
import numpy as np
import atomigi
from lxml import etree
from datetime import datetime
from collections import OrderedDict

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
  for word in wl.split():
    w=word.lower()
    if re.match(r'.*\'',w): w=w+'o'
    w=pattern.sub('',w)
    if len(w)>0: result.append(w)
  return result

def get_word_df(textpath): 
  wordlist=[]
  for filename in os.listdir(textpath):
    print 'reading file',filename
    try: 
       wlist=get_text(os.path.join(textpath,filename))
       wlist=[canonize(w) for w in wlist]
       wordlist.extend([item for sublist in wlist for item in sublist])  
    except: 
       print "exception reading file"
       pass
  word_df=pd.DataFrame({'word':wordlist})
  return word_df

def add_flags(wcdf,neniu_thresh=0.3):
  root_wordcount=pd.DataFrame(wcdf.groupby(['root'])['root'].agg('count').sort_values(ascending=False))
  root_wordcount=root_wordcount.rename(columns={'root':'Nwords'})
  atomvortoj=pd.DataFrame(wcdf.loc[(wcdf.modstring=='NENIU')&(wcdf.root_freq>neniu_thresh)].root)
  atomvortoj['atomvorto']=True
  atomvortoj=atomvortoj.set_index('root')
  root_wordcount=pd.merge(root_wordcount,atomvortoj,left_index=True,right_index=True,how='left').fillna(False)
  wcdf=pd.merge(wcdf,root_wordcount,left_on='root',right_index=True,how='left')
  return wcdf


def get_wcdf(wcdfpath,args):
  ''' get the word count dataframe. This function needs to be sped up. Probably best to cache the table.'''
  if ((os.path.exists(wcdfpath)) and not args.force): return pickle.load(open(wcdfpath,'rb'))

  print '{}: leganta tekstojn...'.format(datetime.now().strftime('%H:%M:%S'))
  word_df=get_word_df(os.path.join(args.textpath,'tekstoj'))
  print '{}: kreanta vortkvanton...'.format(datetime.now().strftime('%H:%M:%S'))

  word_counts=word_df['word'].value_counts()
  aaa=pd.DataFrame(word_counts)
  aaa=aaa.rename(columns={'word':'N'})
  aaa=aaa[aaa.N>=args.min_count]
  atomigitaj=[atomigi.vorto(x) for x in aaa.index]  # this is the slow part. 
  wcdf=pd.DataFrame([(x.vorto,x.radik,x.prefs,x.sufs,x.fins,x.prefs+x.sufs+x.fins,x.valid,x.__repr__()) for x in atomigitaj],columns=['word','root','prefs','sufs','fins','mods','valid','repr'])
  wcdf['N']=aaa.values
  wcdf['root_count']=wcdf['N'].groupby(wcdf['root']).transform('sum')
  wcdf['root_freq']=wcdf['N']/wcdf['root_count']
  wcdf.loc[wcdf.index==wcdf.root,'mods']=[u'NENIU']
  wcdf['modstring']=[''.join(x) for x in wcdf['mods']]
  wcdf=add_flags(wcdf)
  pickle.dump(wcdf,open(wcdfpath,'wb'))
  return wcdf

def main(args):
  for k in sorted(args.__dict__.keys()): print '{0:20} : {1}'.format(k,args.__dict__[k])
  pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])

  wcdfpath=os.path.join(args.basepath,args.wordcount_name)
  wcdf0=get_wcdf(wcdfpath)

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
