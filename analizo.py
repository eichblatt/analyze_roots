# -*- coding: utf-8 -*- 
import os,io,re,sys,pdb,itertools,argparse,pickle
import pandas as pd
import numpy as np
import atomigi
import wordcounts
import math
from datetime import datetime
from collections import OrderedDict
from matplotlib import pyplot as plt
from sklearn.cluster.bicluster import SpectralBiclustering

#pd.options.display.max_columns=30
#pd.set_option('display.width', pd.util.terminal.get_terminal_size()[0])

pd.set_option('display.height', 30)
pd.set_option('display.max_rows', 30)
pd.set_option('display.max_columns', 30)
pd.set_option('display.width', 1000)

#reload(sys)  
#sys.setdefaultencoding('utf-8')

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
  www['mod_count']=www.groupby(['mod'])['root'].transform('count')
  www['mod_freq']=www.groupby(['root','mod'])['root_freq'].transform('sum')
  www=www.loc[www.mod_count>args.min_mod_count]
  g=pd.DataFrame(www.groupby(['root','mod'])['root_freq'].apply(sum)).reset_index()
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

def build_sim_matrix(wcdf,groups):
  result=pd.DataFrame()
  for group in groups:
    smat=get_sim_matrix(wcdf.loc[wcdf.groupnum==group],None)
    smat['groupnum']=group
    result=pd.concat((result,smat))
  return result

def get_sim_matrix(wcdf,corr=None,thresh=0.5):
  ''' create a table with the similarity computed for each pair of roots. Very slow '''
  simdf=pd.DataFrame()
  roots=np.unique(wcdf['root'].values)
  for r in roots:
    if not corr is None:
      others=corr[corr[r]>thresh].index
      others=others[np.where(others>r)] 
    else:
      others=roots[np.where(roots>r)]
    print("working on",r,"with",len(others),"others.")
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
  
def describe_all(pstack,prack,wcdf,depth=5,n_examples=5,n_digits=4):
  g1=pd.DataFrame({x:pstack.sort_values(by=x,ascending=False).index[:depth].values for x in pstack.columns}).T
  g1=g1.rename(columns={x:str(x)+'_X' for x in g1.columns}) 
  g2=pd.DataFrame({x:pstack.sort_values(by=x,ascending=False).reset_index()[x][:depth].values for x in pstack.columns}).T
  g2=(g2*math.pow(10,n_digits)).fillna(0).astype('int')
  g2=g2.rename(columns={x:str(x)+'_n' for x in g2.columns}) 
  g=pd.concat((g1,g2),axis=1).sort_index(axis=1)
  g.index=g.index.astype('int')
  g=pd.merge(g,pd.DataFrame(prack.fillna(0.0).groupby('groupnum')['NENIU'].count()).rename(columns={'NENIU':'N'}),left_index=True,right_index=True)
  g['ekz']=[list(wcdf.loc[wcdf.groupnum==x].groupby('root').first().sort_values(by='root_count',ascending=False)[:n_examples].index) for x in g.index if not pd.isnull(x)]
  g=g.sort_values(by='N',ascending=False)
  colwidth=pd.get_option('max_colwidth')
  pd.set_option('max_colwidth',200)
  print(g)

  pd.set_option('max_colwidth',colwidth)
  return g

def describe_group(wcdf,groupnum):
  ppp=wcdf.pivot_table(index=['groupnum','root'],columns=['modstring'],values='root_freq')
  pmeans=ppp.reset_index().fillna(0.0).groupby('groupnum').mean()
  pstack=pd.DataFrame(pmeans.unstack()).reset_index().pivot(index='modstring',columns='groupnum')[0].sort_values(by=groupnum,ascending=False)
  pstack=pstack.replace({0.00:None})
  print(pstack)
  return pstack

def atomwords(wcdf,neniu_thresh=0.3):
  root_wordcount=pd.DataFrame(wcdf.groupby(['root'])['root'].agg('count').sort_values(ascending=False))
  root_wordcount=root_wordcount.rename(columns={'root':'Nwords'})
  atomvortoj=pd.DataFrame(wcdf.loc[(wcdf.modstring=='NENIU')&(wcdf.root_freq>neniu_thresh)].root)
  atomvortoj['atomvorto']=True
  atomvortoj=atomvortoj.set_index('root')
  root_wordcount=pd.merge(root_wordcount,atomvortoj,left_index=True,right_index=True,how='left').fillna(False)
  wcdf=pd.merge(wcdf,root_wordcount,left_on='root',right_index=True,how='left')
  return wcdf

def examine_group(wcdf,groupnum,topn=5,min_freq=0.05,means_only=False):
  part=wcdf.loc[np.in1d(wcdf.root,groups[groupnum][:topn])].sort_values(['root','root_freq'],ascending=False)[['root','word','mods','modstring','groupnum','root_count','root_freq']].loc[wcdf.root_freq>min_freq]
  #print part
  p=part.pivot(index='root',columns='modstring',values='root_freq')
  if not means_only: print(p)
  pmeans=p.fillna(0.).mean().sort_values(ascending=False)
  print(pmeans.loc[pmeans>min_freq])
  return part

def cluster_example(wcdf,prack):
  ''' toy function to demonstrate the case of finding a few obvious groups '''
  newdata=prack.loc[[u'ruĝ','blu','verd','hav','est','ir','hund','kat','kunikl','kaj','sed','la']][['NENIU','as','os','is',u'iĝi','anta','inta','onta','a','aj','o','oj']]
  clust_example=spectral_cluster(newdata,(4,4))
  newdata['groupnum']=clust_example.row_labels_
  wc_example=pd.merge(wcdf,newdata.reset_index()[['root','groupnum']],on='root',how='left')
  print ({n:group_roots(wc_example,n) for n in range(4)})
  return newdata

def get_wcdict(wcdf):
  wcd={}
  wcd['narrow']=wcdf.loc[(wcdf.atomvorto==False)&(wcdf.Nwords<=7)]
  wcd['broad']=wcdf.loc[(wcdf.atomvorto==False)&(wcdf.Nwords>7)]
  wcd['atom']=wcdf.loc[(wcdf.atomvorto==True)&(wcdf.root_count>100)]
  return wcd

def get_sim(wcd):
  try:
    sim=pd.read_pickle(os.path.join(args.basepath,'sim_matrix.pkl'))
  except:
    print('failed to load similarity matrix')
    sim={k:build_sim_matrix(wcd[k],np.unique(wcd[k].groupnum)) for k in wcd.keys()}
  return sim
  
def analyze(wcd,rackd,k):
  prack=pivot_rack(rackd[k])
  clusters=spectral_cluster(prack,(args.n_groups,args.n_modgroups))
  prack['groupnum']=clusters.row_labels_

  wcdf=wcdf.reset_index()
  wcdf=pd.merge(wcdf,prack.reset_index()[['root','groupnum']],on='root',how='left')

  groups={n:group_roots(wcdf,n) for n in range(args.n_groups)}
  pstack=describe_group(wcdf,0)
  description=describe_all(pstack,prack,wcdf,n_digits=2)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("-basepath", dest='basepath', default=os.path.join(os.getenv('HOME','C:/Users/steve'),'esperanto/analyze_roots'),help="top folder of corpus")
    argparser.add_argument("-textpath", dest='textpath', default=os.path.join(os.getenv('HOME','C:/Users/steve'),'esperanto/tekstaro_de_esperanto_2017-07-25'),help="top folder of corpus")
    argparser.add_argument("-wordcount_name", dest='wordcount_name', default='word_counts.pkl',help="word count dataframe path")
    argparser.add_argument("-force", dest='force', default=False,type=bool,help="create word count dataframe, even if it already exists")
    argparser.add_argument("-min_count", dest='min_count', default=3,type=int,help="minimum number of appearances for a word to count")
    argparser.add_argument("-min_mod_count", dest='min_mod_count', default=5,type=int,help="minimum number of appearances for a modification to count")
    argparser.add_argument("-topn", dest='topn', default=None,type=int,help="number of roots to consider (ranked by frequency)")
    argparser.add_argument("-n_groups", dest='n_groups', default=30,type=int,help="number of root clusters")
    argparser.add_argument("-n_modgroups", dest='n_modgroups', default=50,type=int,help="number of modstring clusters")

    args = argparser.parse_args()
    #import pdb; pdb.set_trace()
    for k in sorted(args.__dict__.keys()): print ('{0:20} : {1}'.format(k,args.__dict__[k]))
    pstackd={}; description={}
    wcdfpath=os.path.join(args.basepath,'word_counts.pkl')
    wcdf0=wordcounts.get_wcdf(wcdfpath,args)
    rootlist=list(OrderedDict.fromkeys(wcdf0[wcdf0.valid].sort_values(by='root_count',ascending=False)['root'].values))[:args.topn]
    wcdf=wcdf0[wcdf0['root'].isin(rootlist)]
    wcd=get_wcdict(wcdf)  
  
    rackd={k:get_modstring_rack(wcd[k]) for k in wcd.keys()}
    prackd={k:pivot_rack(rackd[k]) for k in rackd.keys()}
    clusters={k:spectral_cluster(prackd[k],(args.n_groups,args.n_modgroups)) for k in prackd.keys()}
    
    for k in clusters.keys(): 
      prackd[k]['groupnum']=clusters[k].row_labels_
      wcd[k]=wcd[k].reset_index()
      wcd[k]=pd.merge(wcd[k],prackd[k].reset_index()[['root','groupnum']],on='root',how='left')
    
      groups={n:group_roots(wcd[k],n) for n in range(args.n_groups)}
      pstackd[k]=describe_group(wcd[k],0)
      description[k]=describe_all(pstackd[k],prackd[k],wcd[k],n_digits=2)
    
    sim=get_sim(wcd)
    sim=pd.merge(sim,wcdf.groupby('root')[['Nwords','atomvorto']].agg('first'),how='left',right_index=True,left_on='root1')
    sim=pd.merge(sim,wcdf.groupby('root')[['Nwords','atomvorto']].agg('first'),how='left',right_index=True,left_on='root2')
    most_similar=sim.loc[(sim.atomvorto_x==False)&(sim.atomvorto_y==False)&(sim.Nwords_x>9)&(sim.Nwords_y>9)]
    
    most_similar=most_similar.loc[most_similar.similarity>.8]
  
