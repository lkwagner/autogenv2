

class PropertiesWriter:
  # Should we require a CrystalWriter object?
  def __init__(self,crystal_writer,cryapi=False):
    self.kmesh=crystal_writer.kmesh
    self.cryapi=cryapi
    self.gmesh=crystal_writer.kmesh
    self.boundary=crystal_writer.boundary

  def properties_input(self):
    outlines=['NEWK']
    if self.boundary=='3d':
      outlines+=[ "0 %i"%self.gmesh,
                  " ".join(map(str,self.kmesh))]
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
