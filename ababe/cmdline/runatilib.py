# coding: utf-8
# Distributed under the terms of the MIT License.
import ababe
import click
import sys, ast
from ababe.cmdline.apps import *

def run():
    try:
        ababe.cmdline.runatilib.exec_from_cmdline()
    except KeyboardInterrupt:
        pass
    except EOFError:
        pass

@click.group()
@click.version_option(version='0.1.0')
def exec_from_cmdline():
    pass

@exec_from_cmdline.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--center-element', '-c', 'cenele', required=True)
@click.option('--radius', '-r', type=float, default=0)
@click.option('--element-remove', '-e', 'ele', required=True)
@click.option('--refined/--no-refined', default=True)
def atclear(files, cenele, radius, ele, refined):
    appatomclarifier = atomclarifier.App(files, cenele, radius, ele, refined)
    appatomclarifier.run()

@exec_from_cmdline.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--radius', '-r', type=float, default=0)
@click.option('--outmode', type=click.Choice(['vasp', 'yaml', 'stdio']), default='vasp')
@click.option('--verbose/--no-verbose', '-v/ ', default=False)
def perturb(files, radius, outmode, verbose):
    from ababe.io.io import GeneralIO
    from ababe.stru.scaffold import ModifiedCell
    from ababe.stru.element import Specie
    import os
    for infile in files:
        basefname = os.path.basename(infile)

        # read
        gcell = GeneralIO.from_file(infile)
        # process
        mcell = ModifiedCell.from_gcell(gcell)
        mcell.perturb(radius)

        # write
        out = GeneralIO(mcell.to_gcell())

        if outmode == 'stdio':
            out.write_file(fname=None, fmt='vasp')
        else:
            if verbose: print("PROCESSING: {:}".format(infile))
            ofname = "{:}_PURB.{:}".format(basefname.split('.')[0], outmode)
            out.write_file(ofname)

@exec_from_cmdline.command()
@click.argument('file', type=click.Path(exists=True))
@click.option('--scale', '-s', nargs=3, type=int)
@click.option('--outmode', type=click.Choice(['vasp', 'yaml', 'stdio']), default='stdio')
def supcell(file, scale, outmode):
    from ababe.io.io import GeneralIO
    import os
    import numpy as np

    basefname = os.path.basename(file)

    gcell = GeneralIO.from_file(file)

    scale_matrix = np.diag(np.array(scale))
    sc = gcell.supercell(scale_matrix)

    out = GeneralIO(sc)

    print("PROCESSING: {:}".format(file))
    if outmode == 'stdio':
        out.write_file(fname=None, fmt='vasp')
    else:
        ofname = "{:}_SUPC.{:}".format(basefname.split('.')[0], outmode)
        out.write_file(ofname)

@exec_from_cmdline.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--outmode', type=click.Choice(['vasp', 'yaml', 'stdio']), default='stdio')
@click.option('--verbose/--no-verbose', '-v/ ', default=False)
def pcell(files, outmode, verbose):
    from ababe.io.io import GeneralIO
    import os
    for infile in files:
        basefname = os.path.basename(infile)

        # read
        gcell = GeneralIO.from_file(infile)
        # process
        pcell = gcell.get_refined_pcell()

        # write
        out = GeneralIO(pcell)

        if outmode == 'stdio':
            out.write_file(fname=None, fmt='vasp')
        else:
            if verbose: print("PROCESSING: {:}".format(infile))
            ofname = "{:}_PRIMC.{:}".format(basefname.split('.')[0], outmode)
            out.write_file(ofname)

@exec_from_cmdline.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--outmode', type=click.Choice(['vasp', 'yaml', 'stdio']), default='stdio')
@click.option('--verbose/--no-verbose', '-v/ ', default=False)
def concell(files, outmode, verbose):
    from ababe.io.io import GeneralIO
    import os
    for infile in files:
        basefname = os.path.basename(infile)

        gcell = GeneralIO.from_file(infile)
        pcell = gcell.get_refined_cell()

        out = GeneralIO(pcell)

        if outmode == 'stdio':
            out.write_file(fname=None, fmt='vasp')
        else:
            if verbose: print("PROCESSING: {:}".format(infile))
            ofname = "{:}_CONC.{:}".format(basefname.split('.')[0], outmode)
            out.write_file(ofname)

@exec_from_cmdline.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--zlength', '-z', type=float, default=15)
@click.option('--outmode', type=click.Choice(['vasp', 'yaml', 'stdio']), default='vasp')
@click.option('--verbose/--no-verbose', '-v/ ', default=False)
def d2cell(files, zlength, outmode, verbose):
    from ababe.io.io import GeneralIO
    from ababe.stru.scaffold import ModifiedCell
    from ababe.stru.element import Specie
    import os
    for infile in files:
        basefname = os.path.basename(infile)

        # read
        gcell = GeneralIO.from_file(infile)
        # process
        mcell = ModifiedCell.from_gcell(gcell)
        mcell.d2_at_Z(zlength)

        # write
        out = GeneralIO(mcell.to_gcell())

        if outmode == 'stdio':
            out.write_file(fname=None, fmt='vasp')
        else:
            if verbose: print("PROCESSING: {:}".format(infile))
            ofname = "{:}_d2C.{:}".format(basefname.split('.')[0], outmode)
            out.write_file(ofname)

@exec_from_cmdline.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--sym', '-s', type=float, default=1e-3)
@click.option('--verbose/--no-verbose', '-v/ ', default=False)
def pspg(files, sym, verbose):
    from ababe.io.io import GeneralIO
    import os
    for infile in files:
        basefname = os.path.basename(infile)

        # read
        gcell = GeneralIO.from_file(infile)

        # print spacegroup
        print(gcell.get_spacegroup(sym))
