import numpy as np

def stepped_path(edges, bins):
    
    """
    Create a stepped path suitable for histogramming
    :param edges: bin edges
    :param bins: bin contents
    """
    if len(edges) != len(bins) + 1:
        raise ValueError("edges must be 1 element longer than bins")

    x = np.zeros((2 * len(edges)))
    y = np.zeros((2 * len(edges)))

    x[0::2], x[1::2] = edges, edges
    y[1:-1:2], y[2::2] = bins, bins
    return x,y

def load_gen2opticalehe_aeffs(flavor):
    f = np.load("data/Gen2_EHE_effective_area_{}.npz".format(flavor))
    cos_theta = f['cos_theta_bins']
    energies = f['energy_bins']

    areas = f['area_in_sqm']
    areas = np.asarray(areas)
    areas_vs_energy_zenith = np.sum(areas, axis=1)
    areas_vs_energy = np.sum(areas_vs_energy_zenith, axis=0) / len(cos_theta) * np.pi
    
    return energies, areas_vs_energy

def get_gen2opticalehe():
    
    energies, e_aeff = load_gen2opticalehe_aeffs('e')
    energies, mu_aeff = load_gen2opticalehe_aeffs('mu')
    energies, tau_aeff = load_gen2opticalehe_aeffs('tau')

    aeff = (e_aeff + mu_aeff + tau_aeff)/3.
    return energies, aeff

def get_radio_review_array():
    n_deep = 144*1
    n_shallow = (169+144)*1

    data_i = np.genfromtxt(f'data/tabulated_veff_aeff_review_hybrid.csv', delimiter=',', skip_header=6, names=['logE', 'dveff', 'daeff', 'sveff', 'saeff', 'overlap_frac', 'deep_only', 'shallow_only'])
    energies = np.power(10.,data_i['logE'])

    average_total_veff = ((data_i['dveff']*n_deep + data_i['sveff']*n_shallow)) * (1/(1+data_i['overlap_frac']))
    average_total_aeff = ((data_i['daeff']*n_deep + data_i['saeff']*n_shallow)) * (1/(1+data_i['overlap_frac']))
    return np.asarray(energies), np.asarray(average_total_veff), np.asarray(average_total_aeff)