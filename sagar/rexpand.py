# -*- coding: utf-8 -*-
import click
import time
import threading
import sys

from sagar.crystal.derive import cells_nonredundant, ConfigurationGenerator
from sagar.io.vasp import read_vasp, write_vasp
from sagar.crystal.structure import symbol2number as s2n

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


@cli.command('cell', short_help="Expanding primitive cell to specific range of volumes.")
@click.argument('pcell_filename', metavar='<primitive_cell_file>',
                type=click.Path(exists=True, resolve_path=True, readable=True, file_okay=True))
@click.option('--comment', '-c', type=str, default='cell',
              help="identifier (first word) of output files, Defualt='cell'.")
@click.option('--dimension', '-d', type=int, default=3,
              help="Dimension of the system, 2 for slab. Defalut=3 for crystal")
@click.option('--volume', '-v', nargs=2, type=int, metavar='<min> <max>',
              help="Expand primitive cell to supercell of volume <min> to <max>, set <min> as -1 for creating only <max> volume expanded supercells")
@click.option('--symprec', '-s', type=float, default=1e-5,
              help="Symmetry precision to decide the symmetry of primitive cell. Default=1e-5")
@click.option('--comprec', '-p', type=float, default=1e-5,
              help="Compare precision to judging if supercell is redundant. Defalut=1e-5")
@click.option('--verbose', '-vvv', is_flag=True, metavar='',
              help="Will print verbose messages.")
def cell(pcell_filename, comment, dimension, volume, symprec, comprec, verbose):
    """
    <primitive_cell_file>  Primitive cell structure file, now only vasp POSCAR version5 supported.
    """
    pcell = read_vasp(pcell_filename)
    (min_v, max_v) = volume
    if min_v == -1:
        click.echo("Expanding primitive to volume {:d}".format(max_v))
        _export_supercell(pcell, comment, dimension, max_v, symprec, comprec, verbose)
    else:
        for v in range(min_v, max_v + 1):
            click.echo("Expanding primitive to volume {:d}".format(v))
            _export_supercell(pcell, comment, dimension, v, symprec, comprec, verbose)


def _export_supercell(pcell, comment, dimension, v, symprec, comprec, verbose):
    spinner = Spinner()
    # spinner.start()
    cells = cells_nonredundant(
        pcell, v, dimension, symprec=symprec, comprec=comprec)
    for idx, c in enumerate(cells):
        if verbose:
            print("    " + "No.{:d}: Processing".format(idx))
        filename = '{:s}_v{:d}_id{:d}'.format(comment, v, idx)
        write_vasp(c, filename)
    # spinner.stop()


@cli.command('conf', short_help="Generating configurations in grid cells or supercells.")
@click.argument('cell_filename', metavar="<parent_cell_filename>")
@click.option('--comment', '-c', type=str, default='confs')
@click.option('--pmode', '-mp', type=click.Choice(['varv', 'svc', 'sc']), default='sc',
              help="[varv|svc|sc] represent ['variable volume cells'|'specific volume cells'|'specific cell'] respectively. Deciding what kinds of parent cell to be used to getting configurations")
@click.option('--cmode', '-mc', type=click.Choice(['vc', 'cc']), default='cc',
              help="[vc|cc] for 'variable concentration' and 'certain concentration' respectively.")
@click.option('--dimension', '-d', type=int, default=3,
              help="Dimension of the system, 2 for slab. Defalut=3 for crystal")
@click.option('--volume', '-v', nargs=2, type=int, metavar='<min> <max>',
              help="Expand primitive cell to supercell and to generate configurations of volume <min> to <max>, set <min> as -1 for creating only <max> volume expanded supercells. ONLY USED WHEN --pmode=[varv|svc]")
@click.option('--element', '-e', type=str, metavar='<symbol of element>',
              help="Symbol of element of original sites")
@click.option('--substitutes', '-s', type=str, multiple=True, metavar='<symbol of element>',
              help="Symbol of element to be disorderd substituting, 'Vac' for empty position aka vacancy, multiple optione supported for multielement alloy")
@click.option('--number', '-n', type=int, multiple=True,
              help="number of substitutes element, only used when --cmode=cc")
@click.option('--symprec', type=float, default=1e-5,
              help="Symmetry precision to decide the symmetry of cells. Default=1e-5")
@click.option('--comprec', type=float, default=1e-5,
              help="Compare precision to judging if supercell is redundant. Defalut=1e-5")
@click.option('--verbose', '-vvv', is_flag=True, metavar='',
              help="Will print verbose messages.")
def conf(cell_filename, comment, pmode, cmode, dimension, volume, element, substitutes, number,  symprec, comprec, verbose):
    """
    <parent_cell_file> is the parent cell to generating configurations by sites disorder.\n
    The non-primitive cell can only used as argument when '--pmode=sc'.\n
    Command line tool only provide configurations generator for elements disorder, for more flexible usage such as specific site disorder, please see document http:// , or use python library directly.
    """
    cell = read_vasp(cell_filename)
    cg = ConfigurationGenerator(cell, symprec)
    if pmode == 'varv' and cmode == 'vc':
        click.secho("Expanding and generating configurations: ")
        click.secho(
            "(may take much time)", blink=True, bold=True, bg='magenta', fg='white')
        spinner = Spinner()
        spinner.start()
        (min_v, max_v) = volume
        if min_v == -1:
            min_v = 1
        sites = _get_sites(list(cell.atoms), element, substitutes)
        confs = cg.cons_max_volume(
            sites, max_v, min_volume=min_v, dimension=dimension, symprec=symprec)
        for idx, c in enumerate(confs):
            c = c.get_primitive_cell()
            filename = '{:s}_id{:d}'.format(comment, idx)
            write_vasp(c, filename)
        spinner.stop()
        click.secho("DONE", bold=True, bg='green', fg='white')
    elif pmode == 'svc' and cmode == 'vc':
        click.secho("Expanding and generating configurations: ")
        click.secho(
            "(may take much time)", blink=True, bold=True, bg='magenta', fg='white')
        spinner = Spinner()
        spinner.start()
        (min_v, max_v) = volume
        sites = _get_sites(list(cell.atoms), element, substitutes)
        confs = cg.cons_specific_volume(
            sites, volume=max_v, e_num=None, dimension=dimension, symprec=symprec)
        f_deg = open('deg.txt', 'a')
        for idx, (c, d) in enumerate(confs):
            filename = '{:s}_id{:d}'.format(comment, idx)
            write_vasp(c, filename)
            deg_line = filename + '{:10d}'.format(d) + '\n'
            f_deg.write(deg_line)
        f_deg.close()

        spinner.stop()
        click.secho("DONE", bold=True, bg='green', fg='white')
    elif pmode == 'svc' and cmode == 'cc':
        click.secho("Expanding and generating configurations: ")
        click.secho(
            "(may take much time)", blink=True, bold=True, bg='magenta', fg='white')
        spinner = Spinner()
        spinner.start()
        (min_v, max_v) = volume
        l_atoms = cell.atoms.tolist()
        sites = _get_sites(l_atoms, element, substitutes)
        # number to enum
        ele_n = s2n(element)
        e_total = l_atoms.count(ele_n) * max_v
        e_n = e_total - sum(number)    # 第一个元素的数量
        e_num = [e_n] + list(number)    # 各个元素的数量
        confs = cg.cons_specific_volume(
            sites, volume=max_v, e_num=e_num, dimension=dimension, symprec=symprec)
        f_deg = open('deg.txt', 'a')
        for idx, (c, d) in enumerate(confs):
            filename = '{:s}_id{:d}'.format(comment, idx)
            write_vasp(c, filename)
            deg_line = filename + '{:10d}'.format(d) + '\n'
            f_deg.write(deg_line)
        f_deg.close()

        spinner.stop()
        click.secho("DONE", bold=True, bg='green', fg='white')
    elif pmode == 'sc' and cmode == 'vc':
        click.secho("Generating configurations: ")
        click.secho(
            "(may take much time)", blink=True, bold=True, bg='magenta', fg='white')
        spinner = Spinner()
        spinner.start()
        l_atoms = cell.atoms.tolist()
        sites = _get_sites(l_atoms, element, substitutes)
        confs = cg.cons_specific_cell(sites, None, symprec=symprec)
        f_deg = open('deg.txt', 'a')
        for idx, (c, d) in enumerate(confs):
            filename = '{:s}_id{:d}'.format(comment, idx)
            write_vasp(c, filename)
            # import pdb; pdb.set_trace()
            deg_line = filename + '{:10d}'.format(d) + '\n'
            f_deg.write(deg_line)
        f_deg.close()

        spinner.stop()
        click.secho("DONE", bold=True, bg='green', fg='white')
    elif pmode == 'sc' and cmode == 'cc':
        click.secho("Generating configurations: ")
        click.secho(
            "(may take much time)", blink=True, bold=True, bg='magenta', fg='white')
        spinner = Spinner()
        spinner.start()
        l_atoms = cell.atoms.tolist()
        sites = _get_sites(l_atoms, element, substitutes)
        # number to enum
        ele_n = s2n(element)
        e_total = l_atoms.count(ele_n)
        e_n = e_total - sum(number)    # 第一个元素的数量
        e_num = [e_n] + list(number)    # 各个元素的数量
        confs = cg.cons_specific_cell(sites, e_num, symprec=symprec)
        f_deg = open('deg.txt', 'a')
        # TODO f.close()
        for idx, (c, d) in enumerate(confs):
            filename = '{:s}_id{:d}'.format(comment, idx)
            write_vasp(c, filename)
            deg_line = filename + '{:10d}'.format(d) + '\n'
            f_deg.write(deg_line)
        f_deg.close()

        spinner.stop()
        click.secho("DONE", bold=True, bg='green', fg='white')
    else:
        click.secho("ERROR: --pmode={:s} --cmode={:s} not supported.".format(
            pmode, cmode), bold=True, bg='red', fg='white')


def _get_sites(l_atoms, ele, l_sub):
    ele_n = s2n(ele)
    l_sub_n = [s2n(sub_n) for sub_n in l_sub]
    sites = []
    for a in l_atoms:
        if a == ele_n:
            sites.append(tuple([a] + l_sub_n))
        else:
            sites.append(tuple([a]))
    return sites


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
