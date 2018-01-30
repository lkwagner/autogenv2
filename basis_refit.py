import numpy as np
import scipy.optimize as optimize
import scipy
from math import factorial
from xml.etree.ElementTree import ElementTree,SubElement,Element
def generate_norm(e,n):
  return np.sqrt(2* (4*e)**n * np.sqrt(2*e/np.pi))
      

n_map={'S':1,'P':2,'D':3,'F':5 } 

def evaluate(exp,coeff,el,x):
  f=np.zeros(x.shape)
  for e,c in zip(exp,coeff):
    n=generate_norm(e,n_map[el] )
    f+=n*c*np.exp( - e*x**2)
  return f
  
def fit_exp(basis,expnew=(.2,.4,.8,1.6,3.2,6.4,12.8,25.6) ):
  p0=[1.]*len(expnew)

  xdata=np.logspace(-5,1,1000)
  ydata=evaluate(basis['exp'],basis['coeff'],basis['el'],xdata)
  #print(ydata)
  pf,pcov=optimize.curve_fit(lambda x,p1,p2,p3,p4,p5,p6,p7,p8: \
                             evaluate(expnew,[p1,p2,p3,p4,p5,p6,p7,p8],basis['el'],x),
                             xdata, ydata,p0)
  err=evaluate(expnew,pf,basis['el'],xdata)-ydata
  
  
  return expnew,pf,np.sqrt(np.sum(err**2))
                  


if __name__=="__main__":
  tree=ElementTree()
  tree.parse("BFD_Library.xml")
  root=tree.getroot()
  for child in root:
    print(child.tag,child.attrib['symbol'])
    for child2 in child:
      print(child2.tag)
      if child2.tag=='Basis-set':
        for child3 in child2:
          print(child3.tag)
          if child3.tag=='Contraction':
            n=int(child3.attrib['nterms'])

            basis={}
            basis['el']=child3.attrib['Angular_momentum'].upper()
            if n>1:
              basis['exp']=[]
              basis['coeff']=[]
              for contract in child3:
                basis['exp'].append(float(contract.attrib['Exp']))
                basis['coeff'].append(float(contract.attrib['Coeff']))
              if np.min(basis['exp']) >=0.2:
                expnew=basis['exp']
                pf=basis['coeff']
                err=0.
              elif n==1:
                expnew=[]
                pf=[]
                err=0
              else:
                expnew,pf,err=fit_exp(basis)
              for contract in child3.findall('Basis-term'):
                print(contract)
                child3.remove(contract)
              child3.attrib['nterms']=str(len(expnew))
              for e,c in zip(expnew,pf):
                child3.append(Element("Basis-term",
                  attrib={'Exp':str(e),'Coeff':str(c)}))

  tree.write("BFD_PBC.xml")             
          

