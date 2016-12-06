
# Skeleton of manager.
class Manager:
  def __init__(self,writer,reader,runner):
    self.writer=writer
    self.reader=reader
    self.runner=runner

  ####################################
  def update_status(self,job_record):
    """ Check write status, then if it's running. If not running check if
    finished. If not finished, attempt to run. """ 
    write_status,read_status,run_status,status = 'unknown','unknown','unknown','unknown'

    print("### This will print something like job ID")

    write_status=self.writer.check_status()
    if write_status=='failed':
      status=='failed'
    elif write_status=='not_started':
      inpfn=self.writer.write()
      write_status=self.writer.check_status()
    print("Writer: %s"%write_status)

    run_status=self.runner.check_status()
    print("Runner: %s"%run_status)

    if run_status!='running':
      read_status=self.reader.check_status()
      print("Reader: %s"%read_status)
      if read_status=='not_started':
        run_status=self.runner.run(job_record)
        print("Runner: %s"%run_status)
        status = run_status
    if write_status=='failed' or read_status=='failed' or run_status=='failed':
      status='failed'
    if write_status=='ok' and read_status=='ok' and run_status=='ok':
      status='ok'

    print("Manager: %s"%status)
    return status

  ####################################
  def to_json(self):
    raise NotImplementedError
