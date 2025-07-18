from ase.optimize import FIRE
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase import units
from ase.md import MDLogger
from ase.io import write
from sevenn.calculator import SevenNetCalculator
from pymatgen.io.ase import AseAtomsAdaptor
from siman.core.structure import Structure as STS
from bvlain import Lain 
import numpy as np 
from siman.geo import supercell
import os
import sys


def siman2ase(sts):
    stp = sts.convert2pymatgen()
    sta = AseAtomsAdaptor.get_atoms(stp)
    return sta

def ase2siman(sta):
    stp = AseAtomsAdaptor.get_structure(sta)
    st = STS().update_from_pymatgen(stp)
    return st

def make_inst(st, at):
    st_sc = supercell(st, [10, 10, 10])
    sta = siman2ase(st_sc)
    calc = Lain()

    calc.read_atoms(sta)

    _ = calc.bvse_distribution(mobile_ion='Na+')
    index = np.argmin(calc.distribution)
    coord = calc.mesh_[index]%1
    print(coord)
    # int = Atoms(at, positions=[list(coord)])
    # sta = sta + int
    st = ase2siman(sta)
    st = st.add_atom(xr = coord, element = at)
    return st


    
def Lnvt_7net(sta, init_T, targ_T, steps, write_poscar = 50):

    formula = sta.get_chemical_formula()
    logfile = formula + '_' + 'Lnvt_7net.log'
    poscar_file = './' + formula + '_' + 'Lnvt_7net_poscar'
    os.mkdir(poscar_file)
    MaxwellBoltzmannDistribution(sta, temperature_K=init_T)
    sta.calc = SevenNetCalculator("7net-0")

    timestep = 1 * units.fs            # Шаг по времени
    friction = 0.02                    # Коэффициент трения (нормально ~0.01–0.1)


    # Создание интегратора (Langevin = NVT с термостатом)
    dyn = Langevin(sta, 
                  timestep=timestep, 
                  temperature=targ_T * units.kB,
                  friction=friction,
                  trajectory=formula + f'_{init_T}.traj',
                  loginterval=20  # Записывать данные каждые 20 шагов
                 )

    logfile = open(logfile + ".log", "w")
    dyn.attach(MDLogger(dyn, sta, sys.stdout, header=True, stress=False, peratom=False), interval=10)

    def write_poscar_step():
        step = dyn.nsteps  # текущий номер шага
        filename = f'POSCAR_mlip_{step:05d}'
        write(poscar_file + '/' + filename, sta, format='vasp')

    # Подключаем запись на каждый шаг
    if write_poscar:
        dyn.attach(write_poscar_step, interval=write_poscar)

    # Запускаем MD
    dyn.run(steps)