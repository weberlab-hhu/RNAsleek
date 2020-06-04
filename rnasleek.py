import controller
import argparse


def main(directory, run_file, config_file, check_output=False, prep_multiqc=False):
    samples = controller.parse_runinfo(run_file)
    print(controller.config_based_filters(config_file))
    filtered_samples = controller.filter_runinfo(samples, config_file)
    print(filtered_samples)
    project = controller.Project(directory, job_config=config_file)  # todo, cleanup old read format remnants
    project.define_tasks(filtered_samples)
    if check_output:
        project.check_output()
    else:
        project.setup_directories()
        project.write_scripts()
        project.write_qsubs()
    if prep_multiqc:
        project.prep_multiqc()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', help='project working directory')
    parser.add_argument('run_file', help='file with SRR numbers, one per line')
    parser.add_argument('-c', '--config', help='config file (see included example.ini)', default=None)
    parser.add_argument('--check_output', action='store_true', help='set to check output instead of writing qsubs')
    parser.add_argument('--prep_multiqc', action='store_true', help="sets up symlinks for multiqc")
    args = parser.parse_args()
    main(directory=args.directory, run_file=args.run_file, config_file=args.config, check_output=args.check_output,
         prep_multiqc=args.prep_multiqc)
