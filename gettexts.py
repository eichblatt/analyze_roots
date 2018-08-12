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
  if os.path.exists(wcdfpath): return pickle.load(open(wcdfpath,'rb'))

  print '{}: leganta tekstojn...'.format(datetime.now().strftime('%H:%M:%S'))
  word_df=get_word_df(os.path.join(args.basepath,args.textpath,'tekstoj'))
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

def get_modstring_rack(wcdf):
  www=wcdf[['root','root_freq','modstring']] 
  g=pd.DataFrame(www.groupby(['root','modstring'])['root_freq'].apply(sum)).reset_index()
  g=g.rename(columns={'modstring':'mod'})
  rack=pd.DataFrame(list (itertools.product(np.unique(wcdf.root),np.unique(wcdf.modstring))),columns=['root','mod'])
  rack=pd.merge(rack,g,on=['root','mod'],how='left')
  return rack 

def get_rack(wcdf):
  ''' create a rack with roots by prefix/suffix endings '''
  bb=wcdf.apply(lambda x: pd.Series(x['mods']),axis=1).stack().reset_index(level=1,drop=True)
  www=pd.merge(pd.DataFrame(wcdf[['root','root_freq']]),pd.DataFrame(bb,columns=['mod']),left_index=True,right_index=True)
  www['mod_freq']=www.groupby(['root','mod'])['root_freq'].transform('sum')
#  www['mod_freq']=www['root_freq'].groupby(www[['root','mod']]).transform('sum')
  g=pd.DataFrame(www.groupby(['root','mod'])['root_freq'].apply(sum)).reset_index()
#  rack=pd.DataFrame(list (itertools.product(atomigi.allroots,atomigi.ksufs+atomigi.kprefs+atomigi.fins)),columns=['root','mod'])
  rack=pd.DataFrame(list (itertools.product(np.unique(wcdf.root),atomigi.ksufs+atomigi.kprefs+atomigi.fins+[u'NENIU'])),columns=['root','mod'])
  rack=pd.merge(rack,g,on=['root','mod'],how='left')
  return rack

def pivot_rack(rack,colname='mod'):
  prack=rack.pivot(index='root',columns=colname,values='root_freq')
  return prack

def get_corr(prack):
  pcorr=prack.fillna(0).T.corr()
  bcorr=pd.isnull(prack).T.corr()
  return (pcorr,bcorr)

def roots_compare(wcdf,r1,r2):
  r1=atomigi.addhats(r1)
  r2=atomigi.addhats(r2)
  wc1=wcdf[wcdf.root==r1][['root','root_freq','mods','modstring']]
  wc2=wcdf[wcdf.root==r2][['root','root_freq','mods','modstring']]
  wc12=pd.merge(wc1,wc2,on='modstring',how='outer',suffixes=['_1','_2']) 
  overlap=wc12[np.logical_and(pd.notnull(wc12['root_freq_1']),pd.notnull(wc12['root_freq_2']))]
  differs=wc12[np.logical_or(pd.isnull(wc12['root_freq_1']),pd.isnull(wc12['root_freq_2']))]
  return (overlap,differs)

def roots_similarity(wcdf,r1,r2):
  (overlap,differs)=roots_compare(wcdf,r1,r2)
  overlap_sum=overlap[['root_freq_1','root_freq_2']].min(axis=1).sum()
  differ_sum=0.5*(differs['root_freq_1'].sum()+differs['root_freq_2'].sum())
  return (overlap_sum,differ_sum)

def get_sim_matrix(wcdf,corr,thresh=0.5):
  ''' create a table with the similarity computed for each pair of roots. Very slow '''
  roots=corr.index
  simdf=pd.DataFrame()
  for r in roots:
    others=corr[corr[r]>thresh].index
    others=others[np.where(others>r)] 
    print "working on",r,"with",len(others),"others."
    sim=pd.DataFrame([(r,o,roots_similarity(wcdf,r,o)) for o in others],columns=['root1','root2','sim_diff'])
    simdf=pd.concat([simdf,sim])
  simdf['similarity']=simdf['sim_diff'].apply(lambda x:x[0]-x[1])
  simdf=simdf.sort_values(by='similarity')
  return simdf
  
def spectral_cluster(dataframe,n_clusters=(30,30),show_plots=False):
  model = SpectralBiclustering(n_clusters=n_clusters, method='log',random_state=0)
  data=dataframe.fillna(0.0).values
  model.fit(data)

  fit_data = data[np.argsort(model.row_labels_)]
  fit_data = fit_data[:, np.argsort(model.column_labels_)]

  if show_plots:
    plt.matshow(fit_data, cmap=plt.cm.Blues)
    plt.title("After biclustering; rearranged to show biclusters")
    plt.matshow(np.outer(np.sort(model.row_labels_) + 1, np.sort(model.column_labels_) + 1), cmap=plt.cm.Blues)
    plt.title("Checkerboard structure of rearranged data")
  
  return model

def group_roots(wcdf,n,min_root_count=50,min_n=10):
  return np.unique(wcdf.loc[np.logical_and(wcdf.groupnum==n,np.logical_and(wcdf.root_count>min_root_count,wcdf.N>min_n))].root)
  
def cluster_example(wcdf,prack):
  ''' toy function to demonstrate the case of finding a few obvious groups '''
  newdata=prack.loc[[u'ruĝ','blu','verd','hav','est','ir','hund','kat','kunikl','kaj','sed','la']][['NENIU','as','os','is',u'iĝi','anta','inta','onta','a','aj','o','oj']]
  clust_example=spectral_cluster(newdata,(4,4))
  newdata['groupnum']=clust_example.row_labels_
  wc_example=pd.merge(wcdf,newdata.reset_index()[['root','groupnum']],on='root',how='left')
  print {n:group_roots(wc_example,n) for n in range(4)}
  return newdata

def main(args):
  for k in sorted(args.__dict__.keys()): print '{0:20} : {1}'.format(k,args.__dict__[k])
  pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])

  wcdf0=get_wcdf()

  rootlist=list(OrderedDict.fromkeys(wcdf0[wcdf0.valid].sort_values(by='root_count',ascending=False)['root'].values))[:args.topn]
  wcdf=wcdf0[wcdf0['root'].isin(rootlist)]
  
  print '{}: enrakigante...'.format(datetime.now().strftime('%H:%M:%S'))
  rack=get_rack(wcdf)
  
  print '{}: komputante korelacion...'.format(datetime.now().strftime('%H:%M:%S'))
  prack=pivot_rack(rack)
  clusters=spectral_cluster(prack,(args.ngroups,50))

  prack['groupnum']=clusters.row_labels_

  wcdf=pd.merge(wcdf,prack[['root','groupnum']],on='root',how='left')

  groups={n:group_roots(wcdf,n) for n in range(args.ngroups)}
  #(pcorr,bcorr)=get_corr(prack)
  #sim_matrix=get_sim_matrix(wcdf,bcorr)
  return(wcdf,rack,sim_matrix,groups)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-basepath", dest='basepath', default='/home/steve/esperanto',help="top folder of corpus")
    argparser.add_argument("-textpath", dest='textpath', default='/home/steve/esperanto/tekstaro_de_esperanto_2017-07-25',help="top folder of corpus")
    argparser.add_argument("-wordcount_name", dest='wordcount_name', default='word_counts.pkl',help="word count dataframe path")
    argparser.add_argument("-min_count", dest='min_count', default=3,type=int,help="minimum number of appearances for a word to count")
    argparser.add_argument("-topn", dest='topn', default=None,type=int,help="number of roots to consider (ranked by frequency)")
    argparser.add_argument("-ngroups", dest='ngroups', default=30,type=int,help="minimum number of appearances for a word to count")

    args = argparser.parse_args()
    import pdb; pdb.set_trace()
    (pcorr,bcorr,wcdf,rack)=main(args)
