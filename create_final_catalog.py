# Script to combine the aperture and timeline fluxes into a final catalog.
# Which one to use will be based on the photometric classification that I gave
# each source in each band.
# If the class is 'P' then the timeline flux will be used.
# If the class is 'E', 'C', or 'U' then the aperture flux will be used.
# All aperture fluxes will be corrected if a point source aperture was used
# using the appropriate aperture corrections.

import numpy as np
import pandas as pd

# Import the final aperture fluxes
ap_final = pd.read_csv('results_06-24-2015/final_aperture_fluxes_06-24-2015.csv',
                        index_col=[0,1])
                        
# Import the timeline fitting fluxes
tf_final = pd.read_table('point_source_timeline_fitting_results.txt', delimiter='\t',
                          index_col=0)

# Import the classifications
phot_class = pd.read_csv('bat_spire_phot_class.csv', index_col=0)

# Import the methods for the flux calculation
methods = pd.read_csv('results_06-24-2015/bat_spire_flux_methods.csv', index_col=0)

# Aperture corrections for point source apertures
AP_CORR = {'PSW': 1.2697, 'PMW': 1.2271, 'PLW': 1.2194}

# Calibration error for the timeline fitting fluxes
CAL_ERR_TF = 0.055

# Initialize the DataFrame to hold all of the final results
flux_final = pd.DataFrame(columns=['flux_PSW', 'err_PSW', 'method_PSW',
                                   'flux_PMW', 'err_PMW', 'method_PMW',
                                   'flux_PLW', 'err_PLW', 'method_PLW'],
                                   index=phot_class.index)

# Loop over all of the sources and filters
for s in phot_class.index.values:
    for f in ['PSW', 'PMW', 'PLW']:
        
        pclass = phot_class.loc[s, f]
        m = methods.loc[s, f]
        
        if m == 'TF':
            if f == 'PSW':
                flux = tf_final.loc[s, 'flux_250']/1000.
                flux_err = tf_final.loc[s, 'flux_err_250']/1000.
                flux_err = np.sqrt(flux_err**2 + (flux*CAL_ERR_TF)**2)
                snr = flux/flux_err
            elif f == 'PMW':
                flux = tf_final.loc[s, 'flux_350']/1000.
                flux_err = tf_final.loc[s, 'flux_err_350']/1000.
                flux_err = np.sqrt(flux_err**2 + (flux*CAL_ERR_TF)**2)
                snr = flux/flux_err
            elif f == 'PLW':
                flux = tf_final.loc[s, 'flux_500']/1000.
                flux_err = tf_final.loc[s, 'flux_err_500']/1000.
                flux_err = np.sqrt(flux_err**2 + (flux*CAL_ERR_TF)**2)
                snr = flux/flux_err
            method = 'TF'
            
        elif m == 'AP':
            phot_data = ap_final.loc[s, f]
            if phot_data['type'] == 'extended':
                flux = phot_data['bkgsub_flux']
                flux_err = phot_data['total_err']
                snr = flux/flux_err
            else:
                flux = phot_data['bkgsub_flux']*AP_CORR[f]
                flux_err = phot_data['total_err']*AP_CORR[f]
                snr = flux/flux_err
            method = 'AP'
        
        if snr >= 5:
            flux_final.loc[s, 'flux_'+f] = flux
            flux_final.loc[s, 'err_'+f] = flux_err
            flux_final.loc[s, 'method_'+f] = method
        else:
            flux_final.loc[s, 'flux_'+f] = 0.0
            flux_final.loc[s, 'err_'+f] = 5*flux_err
            flux_final.loc[s, 'method_'+f] = method

