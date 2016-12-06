
# Skeleton of manager.
class Manager:
  def __init__(writer,reader,runner):
    self.writer=writer
    self.reader=reader
    self.runner=runner

  ####################################
  def update_status(self):
    status = 'unknown'
    print("This will print something like job ID")

    write_status = self.writer.check_status()
    print("Writer: %s"%write_status)
    if write_status=='not started':
      status=writer.write()
    elif write_status=='failed':
      status='failed'
    if status=='failed':
      print("Manager: %s"%status)
      return 'failed'

    read_status=self.reader.check_status()
    print("Reader: %s"%read_status)
    if read_status=='not started':
      run_status=self.runner.run()
      print("Runner: %s"%read_status)
    if read_status=='failed' or run_status=='failed':
      status='failed'
    if status=='failed':
      print("Manager: %s"%status)
      return 'failed'

  ####################################
  def to_json(self):
    raise NotImplementedError
