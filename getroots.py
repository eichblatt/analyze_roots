import os,io
from lxml import etree
parser=etree.XMLParser(dtd_validation=True)


basepath='/home/steve/esperanto/revo/xml'
radlist=[]

def get_radical(f):
 p=etree.parse(f,parser)
 for elem in p.iter('*'):
    if elem.tag=='rad': return elem.text

for filename in os.listdir(basepath):
  radlist.append(get_radical(os.path.join(basepath,filename)))

radlist.sort()
with io.open("/home/steve/esperanto/roots.txt","w",encoding='utf8') as f:
  [f.write(unicode(r+"\n")) for r in radlist]
  
