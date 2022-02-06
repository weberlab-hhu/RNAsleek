import os
import re
import copy
from tasks import Task


class FileTooSmallError(Exception):
    pass


class FileStdErrError(Exception):
    pass


class FileStdErrKilled(Exception):
    pass


class Job:
    """Parent class for all jobs with generics setup"""

    TO_OE_TEXT = "2> logs/$PBS_JOBNAME.e${PBS_JOBID%\\[*}.$PBS_ARRAY_INDEX 1> logs/$PBS_JOBNAME.o${PBS_JOBID%\\[*}.$PBS_ARRAY_INDEX"

    def __init__(self, directory, is_array=False):
        self.directory = directory
        self.name = None
        self.threads = 1
        self.time = "01:55:00"
        self.mb = 1000
        self.project = "HelixerOpt"
        self.more_resources = ""
        self.modules = []
        self.user_verbatim = ''
        self.is_array = is_array

    @property
    def extra_loading_verbatim(self):
        return None

    def scale_memory(self, val):
        pass

    def scale_time(self, val):
        pass

    def configure(self, config, job_key, read_format):
        try:
            self.configure_section(config["All"])
        except KeyError:
            pass
        self.configure_section(config[job_key])  # the core job section has to be there
        try:
            self.configure_section(config[job_key][read_format])
        except KeyError:
            pass

    def configure_section(self, config_section):
        config_section = copy.deepcopy(config_section)
        # recognize numbers from config file as numbers (I'm sure there's a better way...)
        for intkey in ["scale_memory", "scale_time", "threads"]:
            try:
                config_section[intkey] = int(config_section[intkey])
            except KeyError:
                pass
        # special case handling
        if "scale_memory" in config_section:
            val = config_section.pop("scale_memory")
            self.scale_memory(val)

        if "scale_time" in config_section:
            val = config_section.pop("scale_time")
            self.scale_time(val)

        # set attributes from rest of non-sub section keys
        for key in config_section:
            # skip any sub-sections
            if isinstance(config_section[key], dict):
                pass
            else:
                if key in self.__dict__:
                    setattr(self, key, config_section[key])
                else:
                    raise ValueError("cannot parse config file for job: {}, unrecognized key: {}".format(self.name,
                                                                                                         key))

    def get_report_path(self, sample_id, all_samples, report_type='out'):
        report_types = {
            'out': self.is_out_path,
            'err': self.is_err_path,
            'log': self.is_log_path
        }
        task_id = self.guess_task_id(sample_id, all_samples)
        is_type_fn = report_types[report_type]
        allfiles = os.listdir(self.directory)
        targfiles = [x for x in allfiles if is_type_fn(x, task_id)]
        targfiles = ['{}/{}'.format(self.directory, x) for x in targfiles]
        if len(targfiles) == 1:
            out = targfiles[0]
        elif len(targfiles) > 1:
            log_times = [os.path.getmtime(x) for x in targfiles]
            pre_out = sorted(zip(targfiles, log_times), key=lambda x: x[1])
            print("time sorted??")
            print(pre_out)
            out = pre_out[-1][0]
        else:
            raise ValueError("Zero matches for {} file {}".format(report_type, targfiles))
        return out

    def is_out_path(self, a_path, task_id):
        return a_path.startswith('{}.o'.format(self.name)) and a_path.endswith('.{}'.format(task_id))

    def is_err_path(self, a_path, task_id):
        return a_path.startswith('{}.e'.format(self.name)) and a_path.endswith('.{}'.format(task_id))

    def is_log_path(self, a_path, task_id):
        right_start = a_path.startswith('{}.'.format(self.name))
        right_task = '[{}].hpc'.format(task_id) in a_path
        right_end = a_path.endswith('.log')
        return right_start and right_task and right_end

    @staticmethod
    def guess_task_id(sample_id, all_samples):
        return all_samples.index(sample_id) + 1  # +1 bc switch python -> bash numbering

    def job_array_info(self, n_total):
        if self.is_array:
            return f"#PBS -J1-{n_total}"
        else:
            return ""

    def preface(self, n_total):
        pfx = """#!/bin/bash
#PBS -l select=1:ncpus={threads}:mem={mb}mb{more}
#PBS -l walltime={time}
#PBS -A "{project}"
#PBS -N {name}
{job_array_info}
#PBS -r y
#PBS -e /dev/null
#PBS -o /dev/null

source $HOME/.bashrc

cd $PBS_O_WORKDIR

## Log-File setup
export LOGFILE=$PBS_O_WORKDIR/logs/$PBS_JOBNAME"."$PBS_JOBID".log"
echo "$PBS_JOBID ($PBS_JOBNAME) @ `hostname` at `date` in "$RUNDIR" START" > $LOGFILE
echo "`date +"%d.%m.%Y-%T"`" >> $LOGFILE 

## Job specific bits
""".format(threads=self.threads,
           mb=self.mb,
           time=self.time,
           project=self.project,
           name=self.name,
           job_array_info=self.job_array_info(n_total),
           more=self.more_resources)

        pfx += self.load_module_text()
        pfx += '\n## main\n'
        return pfx

    def load_module_text(self):
        if isinstance(self.modules, str) or not isinstance(self.modules, list):
            raise ValueError("job modules must be a non-string list")
        text = '## load modules & similar\n'
        for mod in self.modules:
            text += 'module load {}\n'.format(mod)
        if self.extra_loading_verbatim is not None:
            text += self.extra_loading_verbatim + '\n'
        text += '\n'
        return text

    @staticmethod
    def epilog():
        sfx = """
## and more logging

qstat -f $PBS_JOBID >> $LOGFILE  

echo "$PBS_JOBID ($PBS_JOBNAME) @ `hostname` at `date` in "$RUNDIR" END" >> $LOGFILE
echo "`date +"%d.%m.%Y-%T"`" >> $LOGFILE
"""
        return sfx

    def start_main_text(self):
        out = """sample_id=`sed $PBS_ARRAY_INDEX"q;d" {sample_id_file}`
bash scripts/{name}$sample_id".sh {redirect}"
""".format(sample_id_file=self.sample_id_file(True), name=self.name, redirect=self.TO_OE_TEXT)
        return out

    def sample_id_file(self, relative=False):
        if relative:
            direc = ''
        else:
            direc = '{}/'.format(self.directory)
        return '{}sample_ids.txt'.format(direc)

    def main_text(self, run_task):
        raise NotImplementedError("no main_text method implemented for parent Job class")

    def verbatimable(self, **kwargs):
        raise NotImplementedError("no verbatimable method implemented for parent Job Class")

    def write_qsub(self, n_total):
        qsub_text = self.preface(n_total=n_total) + self.start_main_text() + self.epilog()
        qsub_file = self.qsub_file()
        with open(qsub_file, 'w') as f:
            f.write(qsub_text)

    def qsub_file(self):
        qsub_file = '{dir}/qsubs/{name}.qsub'.format(dir=self.directory,
                                                     name=self.name)
        return qsub_file

    def script_file(self, sample_id):
        script_file = '{dir}/scripts/{name}{sample_id}.sh'.format(dir=self.directory,
                                                                  name=self.name,
                                                                  sample_id=sample_id)
        return script_file

    def write_script(self, run_task):
        script_text = self.main_text(run_task)
        script_file = self.script_file(run_task.sample_id)
        with open(script_file, 'w') as f:
            f.write(script_text)

    def output_dirs(self):
        return []

    def expected_output(self, run_task):
        raise NotImplementedError('parent Job has no expected_output')

    @property
    def expected_output_size(self):
        # default to something larger than an empty text file
        return 7000

    def check_output(self, run_task, allsamples):
        sample_id = run_task.sample_id
        try:
            report_path = self.get_report_path(sample_id, allsamples, 'err')
        except ValueError as e:
            yield ('ValueError', e)
            report_path = None

        for f in self.expected_output(run_task):
            if not os.path.exists(f):
                yield ('FileNotFoundError', 'expected output file {} not found, err path {}'.format(f, report_path))
            elif os.path.getsize(f) <= self.expected_output_size:
                yield ('FileTooSmallError', 'expected output file {} to be larger than {}, err path {}'.format(
                    f, self.expected_output_size, report_path)
                )
        if report_path is not None:
            with open(report_path) as f:
                try:
                    stderr_txt = f.read()
                except OSError as e:
                    print(report_path)
                    raise e
                stderr_txt = stderr_txt.lower()
                if 'error' in stderr_txt:
                    yield ('FileStdErrError', 'stderr for job: {}, and sample: {} at {} contains "error"'.format(
                        self.name, sample_id, report_path))
                if 'kill' in stderr_txt:
                    yield ('FileStdErrKilled', 'stderr for job: {}, and sample: {}  at {} contains "kill"'.format(
                        self.name, sample_id, report_path))
        return

    def clean_up(self, run_task):
        # some jobs will want to delete their output if the next job has completed successfully
        pass


# user customizable job (useless w/o config of course)
class CustomJob(Job):
    def __init__(self, directory):
        super(CustomJob, self).__init__(directory)
        self.name = 'custom'

    def expected_output(self, run_task):
        print("Warn: no expected output for CustomJob with name {}".format(self.name))  # todo, how can user specify?
        return []

    def output_dirs(self):
        print("Warn: no known output_dirs for CustomJob with name {}".format(self.name))  # todo, how can user specify?
        return []

    def main_text(self, run_task):
        text = "sample_id={} ;\nrun_ids=({}) ;\n".format(run_task.sample_id, ' '.join(run_task.run_ids))
        text += self.verbatimable()
        text += '\n'
        return text

    def verbatimable(self, **kwargs):
        return self.user_verbatim


class CustomJobSingle(CustomJob):
    pass


class CustomJobPaired(CustomJob):
    pass


class WgetJob(Job):
    def __init__(self, directory):
        super(WgetJob, self).__init__(directory)
        self.name = "wget"
        self.time = "23:59:00"
        self.mb = 400

    def main_text(self, run_task):
        text = "cd $PBS_O_WORKDIR/fastqs\n\n"
        for i, run_id in enumerate(run_task.run_ids):
            download_path = run_task.download_paths[i]
            sra = download_path.split('/')[-1]
            text += "{}\n".format(self.verbatimable(run_id, download_path, sra, run_task))

        text += "\ncd $PBS_O_WORKDIR\n"
        return text

    def verbatimable(self, run_id, download_path, sra_file, run_task):
        out = """wget {} 
fastq-dump ./{} --gzip --split-3 {}
""".format(download_path, sra_file, self.user_verbatim)
        # link to expected name
        if sra_file != run_id:
            if run_task.read_format == Task.SINGLE:
                out += "\nln -s {}.fastq.gz {}.fastq.gz".format(sra_file, run_id)
            else:  # paired
                out += "\n" + "\n".join(
                    ["ln -s {}{}.fastq.gz {}{}.fastq.gz".format(sra_file, x, run_id, x)
                     for x in run_task.untrimmed_paired_extras])
        # todo ln -s
        return out + '\n'

    def output_dirs(self):
        return ['fastqs']

# todo, this is 1:1 with fetch... maybe make more generic?
class WgetJobSingle(WgetJob):
    def expected_output(self, run_task):
        out = ['{}/fastqs/{}.fastq.gz'.format(self.directory, x) for x in run_task.run_ids]
        return out


class WgetJobPaired(WgetJob):
    def expected_output(self, run_task):
        out = []
        for run_id in run_task.run_ids:
            out += ['{}/fastqs/{}{}.fastq.gz'.format(self.directory, run_id, x)
                    for x in run_task.untrimmed_paired_extras]
        return out

# fetching section
class FetchJob(Job):
    def __init__(self, directory):
        super(FetchJob, self).__init__(directory)
        self.name = 'fetch'
        self.time = "23:59:00"
        self.mb = 200

    def main_text(self, run_task):
        text = "cd $PBS_O_WORKDIR/fastqs\n\n"
        for run_id in run_task.run_ids:
            text += "{}\n".format(self.verbatimable(run_id))

        text += "\ncd $PBS_O_WORKDIR\n"
        return text

    def verbatimable(self, run_id):
        return "fastq-dump {} --gzip --split-3 {}".format(run_id, self.user_verbatim)

    def output_dirs(self):
        return ['fastqs']

    # todo, move this elsewhere
    def get_read_format(self, run_task):
        s_files = self.expected_output_single(run_task)
        p_files = self.expected_output_paired(run_task)
        if all([os.path.exists(x) for x in s_files]):
            out = Task.SINGLE
        elif all([os.path.exists(x) for x in p_files]):
            out = Task.PAIRED
        else:
            raise FileNotFoundError("""neither single: {}, 
            nor paired: {}
            output could be found for {}""".format(s_files, p_files, run_task.sample_id))
        return out

    def expected_output(self, run_task):
        raise NotImplementedError('parent FetchJob has no expected_output')


class FetchJobSingle(FetchJob):
    def expected_output(self, run_task):
        out = ['{}/fastqs/{}.fastq.gz'.format(self.directory, x) for x in run_task.run_ids]
        return out


class FetchJobPaired(FetchJob):
    def expected_output(self, run_task):
        out = []
        for run_id in run_task.run_ids:
            out += ['{}/fastqs/{}{}.fastq.gz'.format(self.directory, run_id, x)
                    for x in run_task.untrimmed_paired_extras]
        return out


# trimming section
class TrimmingJob(Job):
    def __init__(self, directory):
        super(TrimmingJob, self).__init__(directory)
        self.name = "trimmomatic"
        self.time = "06:05:00"
        self.mb = 2400
        self.threads = 2
        self.install_directory = "$HOME/extra_programs/Trimmomatic-0.36/"
        self.MAXINFO = "36:0.7"

    @property
    def extra_loading_verbatim(self):
        return "export trimmopath={}\n".format(self.install_directory)

    def output_dirs(self):
        return ['trimmed']

    def main_text(self, run_task):
        text = "\n"
        for run_id in run_task.run_ids:
            text += "{}\n".format(self.verbatimable(run_id))
        return text

    def verbatimable(self, run_id):
        raise NotImplementedError("no verbatimable method defined for parent TrimmingJob")

    def expected_output(self, run_task):
        raise NotImplementedError('parent TrimmingJob has no expected_output')


class TrimmingJobSingle(TrimmingJob):
    def __init__(self, directory):
        super(TrimmingJobSingle, self).__init__(directory)
        self.ILLUMINACLIP = "TruSeq3-SE.fa:2:30:10"

    def verbatimable(self, run_id):
        text = """java -jar -Xmx2048m $trimmopath/trimmomatic-0.36.jar SE -threads {threads} fastqs/{run_id}.fastq.gz \
trimmed/{run_id}.fastq.gz ILLUMINACLIP:$trimmopath/adapters/{illumina_clip} MAXINFO:{max_info} {verbatim}\
""".format(run_id=run_id, threads=self.threads, verbatim=self.user_verbatim,
           illumina_clip=self.ILLUMINACLIP, max_info=self.MAXINFO)
        return text

    def expected_output(self, run_task):
        return ['{}/trimmed/{}.fastq.gz'.format(self.directory, x) for x in run_task.run_ids]


class TrimmingJobPaired(TrimmingJob):
    def __init__(self, directory):
        super(TrimmingJobPaired, self).__init__(directory)
        self.ILLUMINACLIP = "TruSeq3-PE-2.fa:3:30:10:1:true"  # todo, is this a good default, for coverage?

    def main_text(self, run_task):
        # raise error if naming assumptions aren't accurate
        if run_task.untrimmed_paired_extras != ['_1', '_2']:
            raise ValueError("mismatch between paired naming Task and TrimmingJob")
        return super(TrimmingJobPaired, self).main_text(run_task)

    def verbatimable(self, run_id):
        text = """java -jar -Xmx2048m $trimmopath/trimmomatic-0.36.jar PE -threads {threads} -basein \
fastqs/{run_id}_1.fastq.gz -baseout trimmed/{run_id}.fastq.gz ILLUMINACLIP:$trimmopath/adapters/{illumina_clip} \
MAXINFO:{max_info} {verbatim}""".format(run_id=run_id, threads=self.threads, illumina_clip=self.ILLUMINACLIP,
                                        max_info=self.MAXINFO, verbatim=self.user_verbatim)
        return text

    def expected_output(self, run_task):
        endings = ['{}.fastq.gz'.format(x) for x in run_task.trimmed_paired_extras]
        return ['{}/trimmed/{}{}'.format(self.directory, ri, end) for end in endings for ri in run_task.run_ids]


class FastqcJob(Job):
    def __init__(self, directory):
        super(FastqcJob, self).__init__(directory)
        self.name = "fastqc"
        self.time = "02:59:00"
        self.mb = 2000
        self.threads = 2

    def main_from_pfx_list(self, trimmed_pfx, untrimmed_pfx):
        text = ''
        qc_cmds = ['fastqc fastqs/{} &'.format(self.verbatimable(x)) for x in untrimmed_pfx]
        qc_cmds += ['fastqc trimmed/{} &'.format(self.verbatimable(x)) for x in trimmed_pfx]

        mv_cmds = ['mv fastqs/{0}_fastqc.zip fastqs/{0}_fastqc.html fastqc/untrimmed/'.format(x) for x in untrimmed_pfx]
        mv_cmds += ['mv trimmed/{0}_fastqc.zip trimmed/{0}_fastqc.html fastqc/trimmed/'.format(x) for x in trimmed_pfx]
        while qc_cmds:
            for i in range(self.threads):
                try:
                    text += qc_cmds.pop() + '\n'
                except IndexError:
                    pass
            text += 'wait\n\n'
        for cmd in mv_cmds:
            text += cmd + '\n'
        return text

    def main_text(self, run_task):
        # logic is in main_from_pfx_list, but the pfx_lists are required from subclasses
        raise NotImplementedError("main_text method not implemented for parent FastqcJob class")

    def expected_output(self, run_task):
        raise NotImplementedError('parent FastqcJob has no expected_output')

    def verbatimable(self, base_file):
        return '{}.fastq.gz {}'.format(base_file, self.user_verbatim)

    def output_dirs(self):
        return ['fastqc', 'fastqc/trimmed', 'fastqc/untrimmed']


class FastqcJobSingle(FastqcJob):
    def main_text(self, run_task):
        files_trimmed = run_task.run_ids
        files_untrimmed = run_task.run_ids
        text = self.main_from_pfx_list(files_trimmed, files_untrimmed)
        return text

    def expected_output(self, run_task):
        endings = ['.zip', '.html']
        out = []
        for direc in ['untrimmed', 'trimmed']:
            for run_id in run_task.run_ids:
                out += ['{}/fastqc/{}/{}_fastqc{}'.format(self.directory, direc, run_id, x) for x in endings]
        return out


class FastqcJobPaired(FastqcJobSingle):
    def main_text(self, run_task):
        files_trimmed = [ri + x for x in run_task.trimmed_paired_extras for ri in run_task.run_ids]
        files_untrimmed = [ri + x for x in run_task.untrimmed_paired_extras for ri in run_task.run_ids]
        text = self.main_from_pfx_list(files_trimmed, files_untrimmed)
        return text

    def expected_output(self, run_task):
        singles = super().expected_output(run_task)  # todo, double inheriting from single is weird, move to FastqcJob?
        out = []
        # todo, clean up, make legible and all that
        untrimmed = [x for x in singles if '/untrimmed/' in x]
        trimmed = [x for x in singles if '/trimmed/' in x]
        midbits = run_task.trimmed_paired_extras
        for single in trimmed:
            out += [re.sub('_fastqc', x + '_fastqc', single) for x in midbits]
        midbits = run_task.untrimmed_paired_extras
        for single in untrimmed:
            out += [re.sub('_fastqc', x + '_fastqc', single) for x in midbits]

        return out

class MappingJob(Job):
    def __init__(self, directory):
        super(MappingJob, self).__init__(directory)
        self.name = "mapping"
        self.time = "06:55:00"
        self.mb = 1000
        self.threads = 4

    def output_dirs(self):
        return ['mapped']

    def verbatimable(self, run_task):
        raise NotImplementedError('Parent Mapping class has no "verbatimable" defined')

    def expected_output(self, run_task):
        return ['{}/mapped/{}.{}'.format(self.directory, run_task.sample_id, x) for x in ['sam', 'bam']]

    # todo, make this into some sort of subordinate class so that it too can have user_verbatim, etc..
    def shared_text(self, run_task):
        text = """
samtools view mapped/{sample_id}.sam -b |samtools sort -T mapped/tmp{sample_id} -@{threads} -o mapped/{sample_id}.bam
""".format(sample_id=run_task.sample_id, threads=self.threads)
        return text

    def main_text(self, run_task):
        text = self.verbatimable(run_task)
        text += self.shared_text(run_task)
        return text


class HisatJob(MappingJob):
    def __init__(self, directory):
        super(HisatJob, self).__init__(directory)
        self.name = 'hisat'
        self.mb = 4000
        self.modules = ['SamTools/1.6']
        self.sp = 'Athaliana'

    def verbatimable(self, run_task):
        raise NotImplementedError('Parent Hisat class has no "verbatimable" defined')


class HisatJobSingle(HisatJob):

    def verbatimable(self, run_task):
        unpaired = ','.join(['trimmed/{}.fastq.gz'.format(x) for x in run_task.run_ids])
        text = """hisat2 -x ../genomes/{sp}/{sp} -U {unpaired} -S mapped/{sample_id}.sam -p{threads} {verbatim}
""".format(sp=self.sp, unpaired=unpaired, sample_id=run_task.sample_id, threads=self.threads,
           verbatim=self.user_verbatim)
        return text


class HisatJobPaired(HisatJob):
    def verbatimable(self, run_task):
        if run_task.trimmed_paired_extras != ['_1P', '_2P', '_1U', '_2U']:
            raise ValueError("miss match between expected trimming endings and Task")
        forward = ','.join(['trimmed/{}_1P.fastq.gz'.format(x) for x in run_task.run_ids])
        reverse = ','.join(['trimmed/{}_2P.fastq.gz'.format(x) for x in run_task.run_ids])
        unpaired = ','.join(['trimmed/{0}_1U.fastq.gz,trimmed/{0}_2U.fastq.gz'.format(x) for x in run_task.run_ids])
        text = """hisat2 -x ../genomes/{sp}/{sp} -1 {forward} -2 {reverse} -U {unpaired} -S mapped/{sample_id}.sam \
-p{threads} {verbatim}
""".format(sp=self.sp, sample_id=run_task.sample_id, forward=forward, unpaired=unpaired, reverse=reverse,
           threads=self.threads, verbatim=self.user_verbatim)
        return text


class BWAJob(MappingJob):
    def __init__(self, directory):
        super(MappingJob, self).__init__(directory)
        self.name = "bwa"
        self.time = "06:55:00"
        self.mb = 2000
        self.modules = ['SamTools/1.6', 'bwa-mem2/2.1']
        self.sp = 'Athaliana'  # todo, what is this doing hard coded?


class BWAJobPaired(BWAJob):
    def verbatimable(self, run_task):
        if run_task.trimmed_paired_extras != ['_1P', '_2P', '_1U', '_2U']:
            raise ValueError("miss match between expected trimming endings and Task")

        forward = ','.join(['trimmed/{}_1P.fastq.gz'.format(x) for x in run_task.run_ids])
        reverse = ','.join(['trimmed/{}_2P.fastq.gz'.format(x) for x in run_task.run_ids])
        text = f"bwa-mem2 mem ../genomes/{self.sp}/{self.sp} {forward} {reverse} {self.user_verbatim}"
        return text


class BWAJobSingle(BWAJob):
    def verbatimable(self, run_task):
        unpaired = ','.join(['trimmed/{}.fastq.gz'.format(x) for x in run_task.run_ids])
        text = f"bwa-mem2 mem ../genomes/{self.sp}/{self.sp} {unpaired} {self.user_verbatim}"
        return text


class CoverageJob(Job):
    def __init__(self, directory):
        super(CoverageJob, self).__init__(directory)
        self.name = "coverage"
        self.time = "01:55:00"
        self.mb = 2500
        self.threads = 1
        self.modules = ['Python/3.4.5']
        self.sp = 'Athaliana'

    @property
    def extra_loading_verbatim(self):
        return 'source ~/virtualenvs/venv3/bin/activate'

    def main_text(self, run_task):
        return "\n{}".format(self.verbatimable(run_task.sample_id))

    def verbatimable(self, sample_id):
        text = """python ~/repos/naivlix1/cov_vs_pos.py -o coverage/{sample_id} -f ../genomes/{sp}/{sp}.fa \
-b mapped/{sample_id}.bam --seperate_x --deterministic=puma -F raw --skip_x {verbatim}
""".format(sp=self.sp, sample_id=sample_id, verbatim=self.user_verbatim)
        return text

    def expected_output(self, run_task):
        sets = ['val', 'test', 'train']
        xy = ['x', 'y']
        endings = ['{}_{}.csv'.format(aset, c) for aset in sets for c in xy]
        return ['coverage/{}.{}'.format(run_task.sample_id, x) for x in endings]

    def output_dirs(self):
        return ['coverage']


class CoverageJobSingle(CoverageJob):
    pass


class CoverageJobPaired(CoverageJob):
    pass


class CollectRNAseqMetricsJob(Job):
    def __init__(self, directory):
        super(CollectRNAseqMetricsJob, self).__init__(directory)
        self.name = 'collect_rnaseq_metrics'
        self.mb = 1200
        self.modules = ['Java/1.8.0']
        self.sp = 'Athaliana'

    def output_dirs(self):
        return ['picard']

    def expected_output(self, run_task):
        return ['{}/picard/{}.stats'.format(self.directory, run_task.sample_id)]

    def verbatimable(self, run_task, strand="NONE"):
        text = """java -Xmx1G -jar $HOME/extra_programs/picard.jar CollectRnaSeqMetrics \
    REF_FLAT=../genomes/{sp}/{sp}.refflat3.gz I=mapped/{srs}.bam O=picard/{srs}.stats \
    {verbatim}""".format(sp=self.sp,
                         srs=run_task.sample_id,
                         verbatim=self.user_verbatim)
        return text

    def main_text(self, run_task):
        return self.verbatimable(run_task)

    @property
    def expected_output_size(self):
        # most stats files seem to be 2700 and some
        return 2500


class CollectRNAseqMetricsJobPaired(CollectRNAseqMetricsJob):
    pass


class CollectRNAseqMetricsJobSingle(CollectRNAseqMetricsJob):
    pass


class GenomePrepJob(Job):
    """single job for genome annotation, not combined with Tasks"""
    def __init__(self, directory, species):
        super(GenomePrepJob, self).__init__(directory)
        self.name = 'genome_prep'
        self.mb = 14000
        self.modules = ['Cufflinks', 'bwa-mem2/2.1']  # will have to doublecheck if using the avail bwa makes sense
        self.sp = species
        self.time = "23:55:00"

    def start_main_text(self):
        text = f"""
cd ..
## refflat for Picard Tools
species={self.sp}
basename=genomes/$species/$species

# remove gz file if exists (doesn't overwrite otherwise)
[ -e $basename.refflat3.gz ] && rm $basename.refflat3.gz

gffread $basename.gff3 -T -o $basename.gtf
gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons $basename.gtf $basename.refflat2
paste <(cut -f 12 $basename.refflat2) <(cut -f 1-10 $basename.refflat2) > $basename.refflat3
gzip $basename.refflat3

## hisat2 index
hisat2-build genomes/$species/$species.fa genomes/$species/$species

## bwa mem index
bwa-mem2 index -p genomes/$species/$species genomes/$species/$species.fa 
"""
        return text
