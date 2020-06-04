class Task:
    """Handler of making/running qsubs for given run"""
    STATUS = 'status'
    CURRENT = 'current_job'
    MEM = 'scale_memory_by'

    SINGLE = 'single'
    PAIRED = 'paired'

    def __init__(self, sample_id, run_ids, jobs, read_format=None, mem_scale=1, download_paths=None):
        self.sample_id = sample_id
        self.run_ids = run_ids
        self.jobs = jobs
        self.status = 0
        self.mem_scale = mem_scale
        self.download_paths = download_paths
        if read_format is None:  # todo, autodetect on None, or somewhere?
            self.read_format = Task.SINGLE
        else:
            self.read_format = read_format
        self.trimmed_paired_extras = ['_1P', '_2P', '_1U', '_2U']
        self.untrimmed_paired_extras = ['_1', '_2']

    def get_status_from_directories(self):
        pass

    def get_status_from_log(self, log_dict):
        self.status = log_dict[self.sample_id][Task.STATUS]
        self.mem_scale = log_dict[self.sample_id][Task.MEM]

    def check_latest(self):
        pass

    def log_dict(self):
        return {self.sample_id: {Task.STATUS: self.status,
                                 Task.CURRENT: self.jobs[self.status].name,
                                 Task.MEM: self.mem_scale}}

    def run_next_qsub(self):
        pass

    def write_next_script(self):
        self.jobs[self.status].write_script(self)

    def write_next_qsub(self, n_total):
        self.jobs[self.status].write_qsub(n_total=n_total)

    def reset(self):
        pass
