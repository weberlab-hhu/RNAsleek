import os
import jobs  # used via eval
import csv
from configobj import ConfigObj
from tasks import Task
import pandas as pd
import subprocess


class Project:
    """Project wide settings and master control"""
    def __init__(self, directory, job_config=None, scale_memory_by=1):
        self.directory = directory
        self.jobs_single = self.define_jobs(job_config, extra='Single')
        self.jobs_paired = self.define_jobs(job_config, extra='Paired')
        self.mem_scale = scale_memory_by
        self.tasks = []

    def setup_directories(self):
        d = self.directory + '/'
        to_make = [d, d + 'logs', d + 'qsubs', d + 'scripts']
        for job in self.jobs_single + self.jobs_paired:
            dirs = job.output_dirs()
            if dirs is not None:
                to_make += [d + targ for targ in dirs]
        for fullpath in to_make:
            if not os.path.exists(fullpath):
                os.mkdir(fullpath)

    def define_jobs(self, job_config=None, extra=None):
        ok_formats = [Task.SINGLE, Task.PAIRED]
        if extra.lower() not in ok_formats:
            raise ValueError("attempt to define jobs with read format not in ".format(ok_formats))

        directory = self.directory  # used via eval
        if job_config is None:
            targ_jobs = ['Fetch', 'Trimming', 'Fastqc', 'Hisat', 'Coverage']
            job_list = []
            for job_name in targ_jobs:
                job_list.append(eval('jobs.{}Job{}(directory)'.format(job_name, extra)))
        else:
            config = ConfigObj(job_config)
            print('adding jobs from {}:'.format(job_config))
            print(config.sections)
            job_list = []
            for section in config.sections:
                if section not in ["All", "Filters"]:
                    new_job = eval('jobs.{}Job{}(directory)'.format(section, extra))
                    new_job.configure(config, section, extra.lower())
                    job_list.append(new_job)
        return job_list

    def define_tasks(self, filtered_run_info=None):
        if filtered_run_info is None:
            filtered_run_info = []
        for sample in filtered_run_info:
            read_format = sample['library_layout']
            if read_format == Task.SINGLE:
                job_list = self.jobs_single
            elif read_format == Task.PAIRED:
                job_list = self.jobs_paired
            else:
                raise ValueError("read_format not in [{}, {}]".format(Task.SINGLE, Task.PAIRED))
            self.tasks.append(Task(sample_id=sample['sample'], jobs=job_list,
                                   run_ids=sample['runs'], read_format=read_format,
                                   download_paths=sample["download_paths"]))

    def check_output(self):
        allsamples = [x.sample_id for x in self.tasks]
        output_report = self.directory + '/output_report.txt'
        with open(output_report, 'w') as f:
            for task in self.tasks:
                for job in task.jobs:
                    for err_type, err in job.check_output(task, allsamples=allsamples):
                        f.write('{}|{}|{}|{}\n'.format(task.sample_id, job.name, err_type, err))

    def write_scripts(self):
        for i in range(len(self.jobs_single)):
            for t in range(len(self.tasks)):
                self.tasks[t].status = i
                self.tasks[t].write_next_script()

    def write_qsubs(self):
        n_total = len(self.tasks)
        # only called for first task, because it sets up job arrays to run all. Todo, refactor to be task independent?
        for i in range(len(self.jobs_single)):
            self.tasks[0].status = i
            self.tasks[0].write_next_qsub(n_total=n_total)
        # and the list of tasks from which to start the job array
        sample_ids = [x.sample_id + '\n' for x in self.tasks]
        sample_id_file = self.tasks[0].jobs[0].sample_id_file()
        with open(sample_id_file, 'w') as f:
            f.writelines(sample_ids)

    def clean_up(self):
        pass

    def prep_multiqc(self):
        """add symlinks so that SRS numbers are used as widely as possible and only latest log files used"""
        # todo, remove or force symbolic links?
        outpath = '{}/multiqc'.format(self.directory)
        fastqc_trimmed_path = '{}/multiqc/fastqc/'.format(self.directory)
        fastqc_untrimmed_path = '{}/multiqc_untrimmed/'.format(self.directory)

        allsamples = [x.sample_id for x in self.tasks]
        for path in [outpath, fastqc_trimmed_path, fastqc_untrimmed_path]:
            if not os.path.exists(path):
                os.makedirs(path)
        for task in self.tasks:
            for job in task.jobs:
                for rtype in ['out', 'err', 'log']:
                    old = job.get_report_path(sample_id=task.sample_id,
                                              all_samples=allsamples,
                                              report_type=rtype)
                    new = '{}/{}.{}.{}'.format(outpath, task.sample_id, job.name, rtype)

                    if job.name == "hisat" and rtype == "err":
                        # because multiQC seems to give up if it sees to many warnings before hisat output
                        self.copy_tail(old, new, 30)
                    else:
                        os.symlink(os.path.relpath(old, start=outpath), new)

        picard = self.directory + '/picard'
        if os.path.exists(picard):
            print('picard relative path: {}'.format(os.path.relpath(picard, start='multiqc/')))
            os.symlink(os.path.relpath(picard, start='multiqc/'), 'multiqc/picard')
            # as multiQC doesn't export summarized coverage automatically for whatever reason, do so here
            infiles = ['{}/{}'.format(picard, x) for x in os.listdir(picard)]
            outfile = '{}/coverage.tsv'.format(self.directory)
            summarize_hist(infiles, outfile)

        trimmed = self.directory + '/fastqc/trimmed'
        untrimmed = self.directory + '/fastqc/untrimmed'
        if os.path.exists('fastqc/trimmed'):
            for item in os.listdir(trimmed):
                if "U_" not in item:
                    os.symlink(os.path.relpath('{}/{}'.format(trimmed, item),
                                               start=fastqc_trimmed_path),
                               '{}/{}'.format(fastqc_trimmed_path, item))

        if os.path.exists(untrimmed):
            print('from: {}\n to: {}'.format(fastqc_untrimmed_path + '/fastqc',
                                             os.path.relpath(untrimmed, start=fastqc_untrimmed_path)))
            os.symlink(os.path.relpath(untrimmed, start=fastqc_untrimmed_path), fastqc_untrimmed_path + '/fastqc')

    @staticmethod
    def copy_tail(source, destination, n):
        with open(source) as f:
            lines = f.readlines()
        lines = lines[-n:]
        with open(destination, 'w') as f:
            f.writelines(lines)


def summarize_hist(infiles, outfile):
    """make tabular summary of picard coverage data"""
    # AUTHOR: Daniela Dey (with minor modification)
    hist = dict()
    rel_position = []
    for f in infiles:
        with open(f, 'r') as hist_file:
            all_values = []
            all_pos = []
            lines = hist_file.readlines()
            filename = f.split('/')[-1].rstrip('.stats')
            for l in lines[11:112]:
                vals = l.split()[1]
                pos = l.split()[0]
                all_pos.append(pos)
                all_values.append(vals)
            if all_values:  # ignore files where no coverage could be calculated
                hist[filename] = all_values
                if not rel_position:
                    rel_position = all_pos  # within loop and if just in case the last file has no coverage
    hist['relative_position'] = rel_position
    print([(key, len(hist[key])) for key in hist])
    df = pd.DataFrame(hist)
    df.to_csv(outfile, sep='\t', index=False)


# Section: sra run info parsing
def parse_runinfo(filein):
    """takes SraRunInfo.csv ncbi format to list of sample-based dictionaries"""
    index_errs_so_far = 0
    with open(filein) as handlein:
        f = csv.reader(handlein)
        next(f)  # skip title line
        group_by_sample = {}
        for csvline in f:
            try:
                line = parse_runinfo_line(csvline)

                # setup dictionary[sample_id] = [1st_sample_line, 2nd_sample_line, ...]
                try:
                    group_by_sample[line['sample']].append(line)
                except KeyError:
                    group_by_sample[line['sample']] = [line]
            except IndexError as e:  # allow exactly one IndexError (for empty last line)
                if not index_errs_so_far:
                    index_errs_so_far += 1
                else:
                    raise IndexError(e)
    out = []
    for sample_name in sorted(group_by_sample.keys()):
        line_list = group_by_sample[sample_name]
        runs = [x['run'] for x in line_list]
        download_paths = [x["download_path"] for x in line_list]
        out.append(format_sample(line_list[0], runs, download_paths))
    return out


def format_sample(line, runs, download_paths):
    pairing_key = {'SINGLE': Task.SINGLE, 'PAIRED': Task.PAIRED}
    out = line
    out.pop('run', None)  # delete run
    out.pop('download_path', None)
    out['runs'] = runs  # replace with accumulated runs
    out['download_paths'] = download_paths
    out['library_layout'] = pairing_key[out['library_layout']]
    return out


def parse_runinfo_line(line):
    """get important bits for a line in SraRunInfo.csv ncbi format"""
    # todo, check column names
    name_cols = {"sample": 24, "sample_name": 29, "run": 0, "library_strategy": 12, "library_layout": 15,
                 "platform": 18, "study": 20, "taxid": 27, "scientific_name": 28, "consent": 44, "download_path": 9}
    out = {}
    for key in name_cols:
        out[key] = line[name_cols[key]]
    return out


def filter_runinfo(parsed_runinfo, config_file):

    #              scientific_name='Arabidopsis thaliana', taxid='3702', platform='ILLUMINA',
    #               library_strategy='RNA-Seq', consent='public'):
    targets = config_based_filters(config_file)
    out = parsed_runinfo
    for key in targets:
        out = [x for x in out if x[key] == targets[key] or x[key] is None]
    return out


def config_based_filters(config_file):
    config = ConfigObj(config_file)
    targets = {
        'scientific_name': None,
        'taxid': None,
        'platform': "ILLUMINA",
        'library_strategy': "RNA-Seq",
        'consent': "public"
    }
    for key in targets:
        try:
            targets[key] = config['Filters'][key]
        except KeyError:
            pass
    return targets
