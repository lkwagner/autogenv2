

class PropertiesWriter:
  def __init__(self):
    self.kmesh=[2,2,2]
    self.cryapi=False
    self.gmesh=16


  def properties_input(self):
    outlines=['NEWK']
    if self.cif!=None:
      outlines+=[ "0 %i"%self.gmesh,
                  " ".join(map(str,self.kmesh)),
    outlines+=["1 1"]
    if self.cryapi:
      outlines+=["CRYAPI_OUT"]
    else:
      outlines+=["67 999"]
    outlines+=["END"]
    return "\n".join(outlines)
    
  
if __name__=="__main__":
  pwriter=PropertiesWriter()
  print(pwriter.properties_input())
