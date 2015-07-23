#Script to convert the Planck 857 and 545 GHz fluxes to SPIRE 350 and 500 um values based on an
#assumed SED.
#The corrections are from the file SpireHfiColourCorrTab_v2.4.fits which gives corrections based on a 
#model graybody SED with varying temperature and beta fixed at 1.8.
#It also provides values for 545/857 flux ratio given the model SED.
#This script will find the nearest 545/857 in the table then use the corrections to determine SPIRE
#350 and 500 um values.
#Then it will compare that to the actual SPIRE values to try and see how close the two surveys are.

from pylab import *
from astropy.table import Table
#Upload the data
execfile('/Users/ttshimiz/Dropbox/Research/Thesis/scripts/upload_bat_ir_database.py')

data_dir = '/Users/ttshimiz/Github/bat-data/'
bat_herschel = pd.read_csv(data_dir+'bat_herschel.csv', index_col=0)
h250 = bat_herschel['PSW']
h350 = bat_herschel['PMW']
h500 = bat_herschel['PLW']

#Upload the table of corrections
table = Table.read('/Users/ttshimiz/Dropbox/Research/Thesis/SPIRE_photometry/SpireHfiColourCorrTab_v2.4.fits')
ratio_model = table['ratio545_857']
k545toPLW = table['k545toPLW']
k857toPMW = table['k857toPMW']


#Use sources which have both a 545 and 857 GHz Planck value
ind = (planck857 != 0) & (planck545 != 0)

planck857_select = planck857[ind]
planck545_select = planck545[ind]
h350_select = h350[ind]
h500_select = h500[ind]
ratio_src = planck545_select/planck857_select

correct_545 = zeros(len(ratio_src))
correct_857 = zeros(len(ratio_src))

for i in range(len(ratio_src)):

	index_correct = argmin(abs(ratio_model - ratio_src[i]))
	correct_545[i] = k545toPLW[i]
	correct_857[i] = k857toPMW[i]
	
h350_from_planck = planck857_select*correct_857
h500_from_planck = planck545_select*correct_545

flux_ratio_350 = h350_select/h350_from_planck
flux_ratio_500 = h500_select/h500_from_planck