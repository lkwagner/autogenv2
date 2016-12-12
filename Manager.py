
# Skeleton of manager.
class CrystalManager:
  def __init__(self,writer,reader,runner):
    self.writer=writer
    self.reader=reader
    self.runner=runner

  ####################################
  def update_status(self):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. """ 
    write_status,read_status,run_status,status = 'unknown','unknown','unknown','unknown'
    crysbase='autogen.d12'

    print("### This will print something like job ID")

    with open(crysbase,'w') as f:
      f.write(self.writer.crystal_input())
    self.runner.run()
    
    stat=self.reader.check_status(crysbase+".o")
    print("Status",stat)
    self.reader.collect()

    return 'ok'

  ####################################
  def to_json(self):
    raise NotImplementedError
