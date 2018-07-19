# -*- coding: utf-8 -*- 
import re,string,itertools,math,sys
import pandas as pd

#reload(sys)  
#sys.setdefaultencoding('utf-8')

rootfile="/home/steve/esperanto/analyze_roots/roots.txt" 
#rootfile="/home/seichblatt/python/roots.txt" 
hatdict={'au':'aŭ','eu':'eŭ','ou':'oŭ','cx':'ĉ','ch':'ĉ','gx':'ĝ','gh':'ĝ','hx':'ĥ','hh':'ĥ','jx':'ĵ','jh':'ĵ','sx':'ŝ','sh':'ŝ'}
atomvortoj=['adiaŭ','almenaŭ','ambaŭ','ankaŭ','ankoraŭ','anstataŭ','antaŭ','apenaŭ','aŭ','baldaŭ','hieraŭ','hodiaŭ','kontraŭ','kvazaŭ','malgraŭ','preskaŭ','ĉirkaŭ','apud','ekster','preter','tamen','kvankam','ili','oni','unu']

tkap=['neni','ti','ĉi','ki','i']
tvost=['u','o','am','el','a','om','es','al','e']
tvort=[k+v for k in tkap for v in tvost]
atomvortoj=atomvortoj+tvort
suf_dict={'ig':1,'iĝ':1,'il':1,'ul':1,'uj':1,'in':1,'aĵ':1,'ebl':1,'eg':1,'ej':1,'et':1,'it':1,'ant':1,\
          'ad':2,'an':2,'ar':2,'ec':2,'em':2,'ind':2,'ism':2,'ist':2,'ont':2,'int':2,'at':2,'um':2,\
          'aĉ':3,'end':3,'obl':3,'on':3,'op':3,'ot':3,'er':3,'id':3,'estr':3,'end':3}
ksufs=suf_dict.keys()

pref_dict={'mal':1,'re':1,'ne':1,'for':1,'sub':1,'en':1,'ĉio':1,'al':1,\
           'ĉe':2,'ek':2,'el':2,'super':2,'inter':2,'dis':2,'sam':2,'post':2,'antaŭ':2,'sen':2,'kun':2,\
           'ĉef':3,'pra':3,'ge':3,'bo':3,'eks':3,'fi':3}
kprefs=pref_dict.keys()
 
fins=['a','e','i','o','u','j','n','as','is','os','us']

def addhats(x):
    if not (re.match(r'.*x.*',x) or re.match(r'.*[c,g,h,j,s]h.*',x) or re.match(r'.*[a,e,o]u.*',x)): return x
    for (k,v) in hatdict.items():
      x=re.sub(k,v,x)
    return x

def readroots(rootfile):
    with open(rootfile,"r") as f:
      allroots=f.read().splitlines()
    f.close()
    allroots=[string.lower(addhats(x)) for x in allroots]
    allroots.sort()
    allroots=[key for key,_ in itertools.groupby(allroots)]
    return allroots

def vowelcount(v):
    vc=0
    for c in v:
      if c in 'aeiou': vc+=1
    return vc

allroots=readroots(rootfile)

class vorto:
  """ Vorto, kaj metodoj por analizi ĝin """

  chiuradik=allroots
 
  def __init__(mem,vorto):
    mem.sufs=[]; mem.prefs=[]; mem.radik=''; mem.fins=[]; mem.valid=False
    mem._vorto=vorto
    mem.vorto=addhats(string.lower(vorto))
    mem.ebloj=pd.DataFrame()
    mem.disigi()
    
  def __str__(mem):
    stringout=''
    stringout+='vorto '+str(mem.vorto)
    stringout+='\nprefs: '
    for s in mem.prefs: stringout+=str(s)
    stringout+='\nroot: '+str(mem.radik)
    stringout+='\nsufs: '
    for s in mem.sufs: stringout+=str(s)
    stringout+='\nfins: '
    for s in mem.fins: stringout+=str(s)
    stringout+='\nebloj: '
    stringout+='\n'+mem.ebloj.to_string()
    stringout+='\nvalid:',str(mem.valid)
    return stringout
    
  def __repr__(mem):
    pword=''
    for p in mem.prefs: pword+=p+','
    pword+=mem.radik if len(mem.prefs)==0 else mem.radik.title()
    for p in mem.sufs:pword+=','+p
    for p in mem.fins:pword+=','+p
    return pword

  def vowelpos(mem,v):
    vp=[];i=0
    for c in v:
      if c in 'aeiou': vp.append(i)
      i=i+1
    return vp
  
  def matchpref(mem,v):
    for p in kprefs:
      if re.match(p,v): return p
    return None
    
  def headtail(mem,v):
    lv=mem.vowelpos(v)[-1]
    tail=v[lv:]; head=v[:lv]
    return [head,tail]
   
  def tiru_radik(mem,sufs,prefs):
    allsufs=''.join(sufs); allprefs=''.join(prefs); allfins=''.join(mem.fins)
    if len(allfins)+len(allsufs)+len(allprefs)<=len(mem.vorto):
      radik=mem.vorto[len(allprefs):len(mem.vorto)-(len(allsufs)+len(allfins))]
    else: 
      radik=None
    return {'radik':radik,'prefs':prefs,'sufs':sufs,'fins':mem.fins,'valid':False}

  def ebloj_prob(mem):
    mem.ebloj['prob']=1.0/len(vorto.chiuradik)
    prob_list=[]
    for i,ebl in mem.ebloj.iterrows():
       prob=ebl['prob']
       if not ebl['valid']:prob=prob*math.pow(len(vorto.chiuradik),-1*vowelcount(ebl['radik']))
       for j in ebl.sufs:
          prob=prob*(.8/len(ksufs))/math.pow(2,(suf_dict[j.encode('utf8')]-1))
       for j in ebl.prefs:
          prob=prob*(.3/len(kprefs))/math.pow(2,(pref_dict[j.encode('utf8')]-1))
       prob_list.append(prob) 
    mem.ebloj['prob']=prob_list
    return

  def tiru_ebloj(mem):
    ebloj=[]
    for isuf in range(1+len(mem.sufs)):
      for ipref in range(1+len(mem.prefs)):
        eblo=mem.tiru_radik(mem.sufs[::-1][:isuf][::-1],mem.prefs[:ipref])
        er=eblo['radik']
        if ((er is None) or (vowelcount(er)==0)): continue
        eblo['valid']=er in vorto.chiuradik
        ebloj.append(eblo)
    mem.ebloj=pd.DataFrame(ebloj)
    mem.ebloj_prob()
    return(mem.ebloj)

  def disigi(mem):
    if mem.vorto=='malpli': mem.prefs=['mal']; mem.radik='pli'
    if ((mem.vorto in atomvortoj) or (vowelcount(mem.vorto)<=1)): mem.radik=mem.vorto
    if len(mem.radik)==0:
      mem.finajhoj()  
      if len(mem.radik)==0:
        mem.sufiksoj()  # chiuj eblaj sufiksoj
        mem.prefiksoj() # chiuj eblaj prefiksoj
    ebloj=mem.tiru_ebloj()
    #import pdb; pdb.set_trace()
    if len(ebloj)>0: 
      eblo=ebloj.loc[ebloj['prob']==max(ebloj['prob'])]
      (mem.radik,mem.prefs,mem.sufs,mem.fins)=ebloj.loc[ebloj['prob'].idxmax()][['radik','prefs','sufs','fins']]
    mem.valid=mem.radik in vorto.chiuradik
    return
 
  def prefiksoj(mem):
    mem.prefs=[]
    peco=mem.vorto
    if vowelcount(peco)<2: return
    pmatch=mem.matchpref(peco)
    while not pmatch is None:
      mem.prefs.append(pmatch)
      peco=peco[len(pmatch):]
      pmatch=mem.matchpref(peco)
    return mem.prefs 

  def finajhoj(mem):
    mem.fins=[]
    if vowelcount(mem.vorto)<2: return
    (head,tail)=mem.headtail(mem.vorto)
    if tail[-2:]=='jn':
      if mem.vorto[:-2] in atomvortoj: 
         mem.radik=mem.vorto[:-2]
         mem.fins = ['j','n'];return
      for c in tail: mem.fins.append(c); 
      return
    if tail[-1] in ['n','j']:
      if mem.vorto[:-1] in atomvortoj: 
         mem.radik=mem.vorto[:-1]
         mem.fins=[tail[-1]]; return
      for c in tail: mem.fins.append(c); 
      return
    if tail[-1]=='s':
      mem.fins.append(tail)
    if len(tail)==1:
      mem.fins.append(tail)
    return mem.fins

  def sufiksoj(mem):
    mem.sufs=[]
    if len(mem.fins)==0:return
    peco=mem.vorto[:len(mem.vorto)-len(''.join(mem.fins))]  # vorto sen finajhoj
    smatch=True
    while smatch:
      if vowelcount(peco)<2: break
      (head,tail)=mem.headtail(peco)
      if tail in ksufs:
        peco=head
        mem.sufs.append(tail)
        smatch=True
        continue
      smatch=False
    mem.sufs=mem.sufs[::-1]
    return mem.sufs
"""

vlist=[vorto(x) for x in ['peco','tamen','atomighi','baldau','esperantistojn','enretigi']]
problemlist=[vorto(x) for x in ['ekvidis','oni','chiu','ilin','nenion','malpli','ĉefruntan','antaŭsciigi','altamoŝta']

word='aferojn'
word='malaperas'
#word='aperas'

#word='sed'
#word='tamen'
word='nemalkontentajn'  
word='Germanio'
word='Germanujo'
res=decompose(word)

vorto='tiam'
vorto='ĉevalo'
decompose('esperantistojn')
frazo='la esperantistoj devus ŝati tion eblon se ĝi funkcios korekte tamen neniu ajn ŝatos ĝin alie'
[decompose(x) for x in string.split(frazo)]
"""
