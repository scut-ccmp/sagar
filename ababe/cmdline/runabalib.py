# coding: utf-8
# Distributed under the terms of the MIT License.
import ababe
import sys, ast
import click
import yaml
from ababe.cmdline.apps import *

def run():
    try:
        ababe.cmdline.runabalib.exec_from_cmdline()
    except KeyboardInterrupt:
        pass
    except EOFError:
        pass

@click.group()
@click.version_option(version='0.1.0')
def exec_from_cmdline():
    pass

#######################################################################
## In the module all inputs from cmd-line is checked
## And then call the modules in directory
## ababe/cmdline/command
#######################################################################

@exec_from_cmdline.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--comment', '-c', default=None)
@click.option('--volumn', '-v', type=int, default=1)
@click.option('--ld/--no-ld', '-L/', default=False)
@click.option('--outmode', '-o', type=click.Choice(['vasp', 'yaml']), default='yaml')
def suplat(input, comment, volumn, ld, outmode):
    infile = click.format_filename(input)

    appsuperlattice = superlattice.App(infile, comment, volumn, ld, outmode)
    appsuperlattice.run()

@exec_from_cmdline.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--comment', '-c', default=None)
@click.option('--element', '-e', default=None)
@click.option('--speckle', '-s', default=None)
@click.option('--number-speckle', '-n', 'nspeckle', type=int, default=None)
@click.option('--dist-restrict', '-r', 'trs', nargs=2, type=click.Tuple([str, float]), multiple=True)
@click.option('--refined/--no-refined', default=True)
@click.option('--outmode', '-o', type=click.Choice(['vasp', 'yaml']), default='vasp')
@click.option('--move-supercell/--no-move-supercell', '-S/-N', 'mpr', default=True,
              help='Whether move no primitive structures')
def ocumaker(input, comment, element, speckle, nspeckle, trs, refined, outmode, mpr):
    infile = click.format_filename(input)

    appoccupymaker = occupymaker.App(infile, comment, element, speckle,
                                     nspeckle, trs, refined, outmode, mpr)
    appoccupymaker.run()

@exec_from_cmdline.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--comment', '-c', default=None)
@click.option('--element', '-e', default=None)
@click.option('--speckle', '-s', default=None)
@click.option('--number-speckle', '-n', 'nspeckle', type=int, default=None)
@click.option('--dist-restrict', '-r', 'trs', nargs=2, type=click.Tuple([str, float]), multiple=True)
@click.option('--refined/--no-refined', default=True)
@click.option('--outmode', '-o', type=click.Choice(['vasp', 'yaml']), default='vasp')
@click.option('--move-supercell/--no-move-supercell', '-S/-N', 'mpr', default=True,
              help='Whether move no primitive structures')
def ocubiter(input, comment, element, speckle, nspeckle, trs, refined, outmode, mpr):
    infile=click.format_filename(input)

    appoccupybiter = occupybiter.App(infile, comment, element, speckle,
                                     nspeckle, trs, refined, outmode, mpr)
    appoccupybiter.run()

@exec_from_cmdline.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--comment', default=None)
@click.option('--element', default=None)
@click.option('--speckle', nargs=2)
@click.option('--number-speckle', '-n', 'nspeckle', type=int, nargs=2)
@click.option('--zoom', type=float, default=None)
@click.option('--dist-restrict', '-r', 'trs', nargs=2, type=click.Tuple([str, float]), multiple=True)
@click.option('--refined/--no-refined', default=True)
@click.option('--outmode', type=click.Choice(['vasp', 'yaml']), default='vasp')
@click.option('--move-supercell/--no-move-supercell', '-S/-N', 'mpr', default=True,
              help='Whether move no primitive structures')
def ocutenter(input, comment, element, speckle, nspeckle, zoom, trs, refined, outmode, mpr):
    infile=click.format_filename(input)
    y = yaml.load(open(infile, "r"))

    appoccupytenter = occupytenter.App(y, comment, element, speckle,
                                       nspeckle, zoom, trs, refined, outmode, mpr)
    appoccupytenter.run()

@exec_from_cmdline.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--comment', '-c', default=None)
@click.option('--element', '-e', default=None)
@click.option('--speckle', '-s', nargs=2)
@click.option('--dist-restrict', '-r', 'trs', nargs=2, type=click.Tuple([str, float]), multiple=True)
@click.option('--refined/--no-refined', default=True)
@click.option('--outmode', '-o', type=click.Choice(['vasp', 'yaml']), default='vasp')
def ocumakert(input, comment, element, speckle, trs, refined, outmode):
    infile=click.format_filename(input)

    appoccupymakert = occupymakert.App(infile, comment, element, speckle,
                                       trs, refined, outmode)
    appoccupymakert.run()

@exec_from_cmdline.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--comment', default=None)
@click.option('--exch', nargs=2)
@click.option('--number-exchange', '-n', 'nexch', type=int)
@click.option('--dist-restrict', '-r', 'trs', nargs=2, type=click.Tuple([str, float]), multiple=True)
@click.option('--refined/--no-refined', default=True)
@click.option('--outmode', type=click.Choice(['vasp', 'yaml']), default='vasp')
@click.option('--move-supercell/--no-move-supercell', '-S/-N', 'mpr', default=False,
              help='Whether move no primitive structures')
def exch(input, comment, exch, nexch, trs, refined, outmode, mpr):
    infile=click.format_filename(input)

    appoccupyexch = occupyexch.App(infile, comment, exch, nexch,
                                   trs, refined, outmode, mpr)
    appoccupyexch.run()
