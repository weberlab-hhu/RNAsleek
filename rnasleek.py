import controller
import click


@click.group()
@click.option('-d', '--directory', help='project directory', required=True)
@click.option('-s', '--sra-run-info', help='info on target data in the format of SraRunInfo.csv')
@click.option('-c', '--run-config', help='config for commands in pipeline, filters, etc...')
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


@cli.command()
def genome():
    pass


def get_project(ctx):
    directory = ctx.obj['directory']
    run_file = ctx.obj['sra_run_info']
    config_file = ctx.obj['run_config']
    assert run_file is not None, "missing option. --sra-run-info required for this subcommand"
    assert config_file is not None, "missing option. --run-config required for this subcommand"
    samples = controller.parse_runinfo(run_file)
    filtered_samples = controller.filter_runinfo(samples, config_file)
    print(len(filtered_samples), "samples passed filter")
    project = controller.Project(directory, job_config=config_file)  # todo, cleanup old read format remnants
    project.define_tasks(filtered_samples)
    return project


if __name__ == "__main__":
    cli()
