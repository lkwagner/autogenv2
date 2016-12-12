
class CrystalManager:
  def __init__(self,writer,crys_reader,crys_runner,prop_reader,prop_runner):
    self.writer=writer
    self.creader=crys_reader
    self.crunner=crys_runner
    self.preader=prop_reader
    self.prunner=prop_runner

  ####################################
  def update_status(self):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. """ 
    write_status,read_status,run_status,status = 'unknown','unknown','unknown','unknown'
    crysbase='crys.in' # (BB) Not autogen.d12: trying to be consistent about input file names.
    propbase='prop.in'

    print("### This will print something like job ID")

    # Generate input files.
    with open(crysbase,'w') as f:
      f.write(self.writer.crystal_input())
    with open(propbase,'w') as f:
      f.write(self.writer.properties_input())

    self.crunner.run(crysbase)
    stat=self.creader.check_status(crysbase+".o")
    print("Crystal status",stat)
    self.creader.collect(crysbase+".o")

    self.prunner.run(propbase)
    stat=self.preader.check_status(crysbase+".o")
    print("Crystal properties status",stat)

    return 'ok'

  ####################################
  def to_json(self):
    raise NotImplementedError
