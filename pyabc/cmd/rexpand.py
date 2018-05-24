import click

from pyabc.crystal.derive import cells_nonredundant, ConfigurationGenerator

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@cli.command('cell', short_help="Expanding primitive cell to specific range of volumes.")
@click.argument('pcell', metavar='<primitive_cell_file>')
@click.option('--volume', '-v', multiple=True, nargs=2, type=int, metavar='<min> <max>',
              help="Expand primitive cell to supercell of volume <min> to <max>, set <min> as -1 for creating specific <max> volume expanded supercell")
@click.option('--verbose', '-vvv', is_flag=False, metavar='',
              help="Will print verbose messages.")
def cell():
    # if
    click.echo("expanding cell")


@cli.command('conf', short_help="Generating configurations in grid cells or supercells.")
def conf():
    click.echo("expanding configuration")
