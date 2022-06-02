import controller
import click
import os
import shutil
import jobs


@click.group()
@click.option('-d', '--directory', help='project directory', required=True)
@click.option('-s', '--sra-run-info', help='info on target data in the format of SraRunInfo.csv')
@click.option('-c', '--run-config', required=True, help='config for commands in pipeline, filters, etc...')
@click.pass_context
def cli(ctx, **kwargs):
    ctx.ensure_object(dict)
    ctx.obj.update(kwargs)


@cli.command()
@click.pass_context
def setup(ctx):
    project = get_project(ctx)
    project.setup_directories()
    project.write_scripts()
    project.write_qsubs()


@cli.command()
@click.pass_context
def check(ctx):
    project = get_project(ctx)
    project.check_output()


@cli.command()
@click.pass_context
def multiqc(ctx):
    project = get_project(ctx)
    project.prep_multiqc()


@cli.command(help="run this only if the genome is not yet setup for target species")
@click.option('-f', '--fasta', required=True, help="genome sequence as fasta file")
@click.option('-g', '--gff3', required=True, help='reference structural annotation as gff3 file')
@click.option('-s', '--species', required=True, help='short species name (for naming files/must later match config)')
@click.pass_context
def genome(ctx, fasta, gff3, species):
    proj_directory = os.path.abspath(ctx.obj['directory'])
    print(proj_directory, type(proj_directory))
    if proj_directory.endswith('/'):
        proj_directory = proj_directory[:-1]
    # genomes directory is sister to project directory, and the species is therein
    genome_sp_directory = os.path.join('/', *proj_directory.split('/')[:-1], 'genomes', species)

    if not os.path.exists(genome_sp_directory):
        os.makedirs(genome_sp_directory)
    shutil.copy(fasta, os.path.join(genome_sp_directory, f'{species}.fa'))
    shutil.copy(gff3, os.path.join(genome_sp_directory, f'{species}.gff3'))

    j = jobs.GenomePrepJob(proj_directory, species)
    j.write_script()
    j.write_qsub(1)


def get_project(ctx):
    directory = ctx.obj['directory']
    run_file = ctx.obj['sra_run_info']
    config_file = ctx.obj['run_config']
    assert run_file is not None, "missing option. --sra-run-info required for this subcommand"
    samples = controller.parse_runinfo(run_file)
    filtered_samples = controller.filter_runinfo(samples, config_file)
    print(len(filtered_samples), "samples passed filter")
    project = controller.Project(directory, job_config=config_file)  # todo, cleanup old read format remnants
    project.define_tasks(filtered_samples)
    return project


if __name__ == "__main__":
    cli()
