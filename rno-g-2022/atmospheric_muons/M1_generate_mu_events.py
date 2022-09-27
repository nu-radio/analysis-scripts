from NuRadioReco.utilities import units
from NuRadioMC.EvtGen.generator import generate_surface_muons
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description='Generate events for muon simulations.')
parser.add_argument('--n_events', type=int, default=10000, help='Choose the number of events for this file')
parser.add_argument('--Emin', type=float, default=1e18, help='min energy of muon energy to simulate in eV')
parser.add_argument('--Emax', type=float, default=1e19, help='max energy of muon energy to simulate in eV')
parser.add_argument('--cos_theta_max', type=float, default=0.0, help='cos zenith bin max (0.8...0')
parser.add_argument('--cos_theta_min', type=float, default=0.1, help='cos zenith bin max (1...0.2')
parser.add_argument('--job', type=int, default=0, help='number of part to simulate')
parser.add_argument('--base_dir', default="/lustre/fs22/group/radio/lpyras/muon_sim/M1_simulation_input",
                    help='"base" directory for the simulations')
args = parser.parse_args()
print("Running generation with arguments", args)

if not os.path.exists(args.base_dir):
    os.makedirs(args.base_dir)

volume = {
    'fiducial_rmin': 0 * units.km,
    'fiducial_rmax': 8 * units.km,
    'fiducial_zmin': -3 * units.km,
    'fiducial_zmax': 0 * units.km}

# volume['full_rmax'] = 6 * units.km
# volume['full_rmin'] = 0 * units.km
# volume['full_zmax'] = 0 * units.km
# volume['full_zmin'] = -4 * units.km

# See https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/INSTALL.md

filename = f'input_{args.Emin:.1e}_{args.Emax:.1e}_cos_theta_{args.cos_theta_min:.1f}_{args.cos_theta_max:.1f}_job_{args.job}.hdf5'
file_path = os.path.join(args.base_dir, filename)


generate_surface_muons(file_path, args.n_events, args.Emin, args.Emax,
                       volume,
                       thetamin=np.arccos(args.cos_theta_min), thetamax=np.arccos(args.cos_theta_max),
                       phimin=0. * units.rad, phimax=2 * np.pi * units.rad,
                       start_event_id=1,
                       plus_minus='mix',
                       n_events_per_file=None,
                       spectrum='log_uniform',
                       start_file_id=0,
                       config_file='Greenland',
                       proposal_kwargs={},
                       log_level=None,
                       max_n_events_batch=args.n_events,
                       seed=None)
