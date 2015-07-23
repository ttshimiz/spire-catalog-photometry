# Script to write the file containing the photometry classifications for all of 
# BAT sources in the SPIRE images.
# P = point source
# E = extended
# U = undetected
# C = cirrus
# Sources are given a P if the source is below a certain FWHM and reduced chi-square limit
# from the HIPE timeline fitting routine. Otherwise they are considered Extended (E).
# Undetected and Cirrus classifications come from a visual classification of the images
# contained in the file bat_spire_visual_class.csv

import pandas as pd

# FWHM and reduced chi-square cutoffs
fwhm_250_cut = 21.0
fwhm_350_cut = 28.0
fwhm_500_cut = 40.0
rchi2_250_cut = 0.0022420542605145711
rchi2_350_cut = 0.0023140627264933229
rchi2_500_cut = 0.0033905236792526215

# Upload the visual classifications
viz_class = pd.read_csv('bat_spire_visual_class.csv', index_col=0)

# Upload the timeline fitting results
tl_results = pd.read_table('point_source_timeline_fitting_results.txt', delimiter='\t',
                           index_col=0)

# Create a new dataframe to hold the classification for each source
phot_class = pd.DataFrame(columns=['PSW', 'PMW', 'PLW'], index=tl_results.index)

# Find all of the point sources
ind_point = ((tl_results['fwhm_250'] < fwhm_250_cut) & (tl_results['rchi2_250'] < rchi2_250_cut) &
            (tl_results['fwhm_350'] < fwhm_350_cut) & (tl_results['rchi2_350'] < rchi2_350_cut) &
            (tl_results['fwhm_500'] < fwhm_500_cut) & (tl_results['rchi2_500'] < rchi2_500_cut))

# Cirrus dominated sources
ind_cirrus = (viz_class['250 Classification'] == 'Cirrus') | (viz_class['250 Classification'] == 'Cirrus (1)')

# Undetected sources
ind_undetect_250 = viz_class['250 Classification'] == 'Undetected'
ind_undetect_350 = viz_class['350 Classification'] == 'Undetected'
ind_undetect_500 = viz_class['500 Classification'] == 'Undetected'

phot_class.loc[ind_point, ['PSW', 'PMW', 'PLW']] = 'P'
phot_class.loc[ind_cirrus, ['PSW', 'PMW', 'PLW']] = 'C'
phot_class.loc[ind_undetect_250, 'PSW'] = 'U'
phot_class.loc[ind_undetect_350, 'PMW'] = 'U'
phot_class.loc[ind_undetect_500, 'PLW'] = 'U'
phot_class[phot_class.isnull()] = 'E'

