import click
import time
import threading
import sys

from pyabc.crystal.derive import cells_nonredundant, ConfigurationGenerator
from pyabc.io.vasp import read_vasp, write_vasp

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@cli.command('cell', short_help="Expanding primitive cell to specific range of volumes.")
@click.argument('pcell_filename', metavar='<primitive_cell_file>',
                type=click.Path(exists=True, resolve_path=True, readable=True, file_okay=True))
@click.option('--comment', '-c', type=str, default='cell',
              help="identifier (first word) of output files, Defualt='cell'.")
@click.option('--volume', '-v', nargs=2, type=int, metavar='<min> <max>',
              help="Expand primitive cell to supercell of volume <min> to <max>, set <min> as -1 for creating only <max> volume expanded supercells")
@click.option('--symprec', '-s', type=float, default=1e-5,
              help="Symmetry precision to decide the symmetry of primitive cell. Default=1e-5")
@click.option('--comprec', '-p', type=float, default=1e-5,
              help="Compare precision to judging if supercell is redundant. Defalut=1e-5")
@click.option('--verbose', '-vvv', is_flag=True, metavar='',
              help="Will print verbose messages.")
def cell(pcell_filename, comment, volume, symprec, comprec, verbose):
    """
    <primitive_cell_file>  Primitive cell structure file, now only vasp POSCAR version5 supported.
    """
    pcell = read_vasp(pcell_filename)
    (min_v, max_v) = volume
    if min_v == -1:
        click.echo("Expanding primitive to volume {:d}".format(max_v))
        _export_supercell(pcell, comment, max_v, symprec, comprec, verbose)
    else:
        for v in range(min_v, max_v + 1):
            click.echo("Expanding primitive to volume {:d}".format(v))
            _export_supercell(pcell, comment, v, symprec, comprec, verbose)


def _export_supercell(pcell, comment, v, symprec, comprec, verbose):
    spinner = Spinner()
    spinner.start()
    cells = cells_nonredundant(
        pcell, v, symprec=symprec, comprec=comprec)
    for idx, c in enumerate(cells):
        if verbose:
            print("    " + "No.{:d}: Processing".format(idx))
        filename = '{:s}_v{:d}_id{:d}'.format(comment, v, idx)
        write_vasp(c, filename)
    spinner.stop()


@cli.command('conf', short_help="Generating configurations in grid cells or supercells.")
@click.argument('cell_filename', metavar="<parent_cell_filename>")
@click.option('--comment', '-c', type=str, default='confs')
@click.option('--pmode', '-mp', type=click.Choice(['varv', 'svc', 'sc']),
              help="[varv|svc|sc] represent ['variable volume cells'|'specific volume cells'|'specific cell'] respectively. Deciding what kinds of parent cell to be used to getting configurations")
@click.option('--cmode', '-mc', type=click.Choice(['vc', 'cc']),
              help="[vc|cc] for 'variable concentration' and 'certain concentration' respectively.")
@click.option('--volume', '-v', nargs=2, type=int, metavar='<min> <max>',
              help="Expand primitive cell to supercell and to generate configurations of volume <min> to <max>, set <min> as -1 for creating only <max> volume expanded supercells. ONLY USED WHEN --pmode=[varv|svc]")
@click.option('--element', '-e', type=str, metavar='<symbol of element>',
              help="Symbol of element of original sites")
@click.option('--substitute', '-s', type=str, multiple=True, metavar='<symbol of element>',
              help="Symbol of element to be disorderd substituting, 'Vac' for empty position aka vacancy, multiple optione supported for multielement alloy")
@click.option('--symprec', '-s', type=float, default=1e-5,
              help="Symmetry precision to decide the symmetry of cells. Default=1e-5")
@click.option('--comprec', '-p', type=float, default=1e-5,
              help="Compare precision to judging if supercell is redundant. Defalut=1e-5")
@click.option('--verbose', '-vvv', is_flag=True, metavar='',
              help="Will print verbose messages.")
def conf(cell_filename, comment, pmode, cmode, volume, element, substitute, symprec, comprec, verbose):
    """
    <parent_cell_file> is the parent cell to generating configurations by sites disorder.\n
    The non-primitive cell can only used as argument when '--pmode=sc'.\n
    Command line tool only provide configurations generator for elements disorder, for more flexible usage such as specific site disorder, please see document http:// , or use python library directly.
    """
    cell = read_vasp(cell_filename)
    cg = ConfigurationGenerator(cell, symprec)
    if pmode == 'varv' and cmode == 'vc':
        click.echo("Expanding and generating configurations: (may take much time)")
        spinner = Spinner()
        spinner.start()
        (min_v, max_v) = volume
        sites =
        confs = cg.cons_max_volume(sites, max_v, min_volume=min_v, symprec)
        for idx, c in enumerate(confs):
            filename = '{:s}_id{:d}'.format(comment, idx)
            write_vasp(c, filename)
        spinner.stop()
    elif pmode == 'svc' and cmode == 'vc':
        print(pmode, cmode)
    elif pmode == 'sc' and cmode == 'vc':
        print(pmode, cmode)
    elif pmode == 'sc' and cmode == 'cc':
        print(pmode, cmode)
    else:
        print("error")


class Spinner:
    busy = False
    delay = 0.1

    @staticmethod
    def spinning_cursor():
        while 1:
            for cursor in '|/-\\':
                yield cursor

    def __init__(self, delay=None):
        self.spinner_generator = self.spinning_cursor()
        if delay and float(delay):
            self.delay = delay

    def spinner_task(self):
        while self.busy:
            sys.stdout.write(next(self.spinner_generator))
            sys.stdout.flush()
            time.sleep(self.delay)
            sys.stdout.write('\b')
            sys.stdout.flush()

    def start(self):
        self.busy = True
        threading.Thread(target=self.spinner_task).start()

    def stop(self):
        self.busy = False
        time.sleep(self.delay)
