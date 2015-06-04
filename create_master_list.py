from glob import glob
import os

f = open('bat_agn_spire_obsids.txt', 'r')
lines = f.readlines()
names = [x.split('\t')[0] for x in lines]
f.close()
#names = ['CenA', 'IC4329A']

mlist = open('master_list_bat', 'w')

for name in names:
    stru_files = glob('/ricci5nb/tshimizu/input_scans/*'+name+'*.xdr')
    mlist.write(name+'_scanamorphos\n')
    mlist.write('/ricci5nb/tshimizu/input_scans/\n')
    mlist.write('/spire\n')
    mlist.write('/nogains\n')
    mlist.write('/one_plane_fits\n')
    if len(stru_files) <= 2:
    	mlist.write('/minimap\n')
    mlist.write("orient='astro'\n")

    for sf in stru_files:
        mlist.write(os.path.basename(sf)+'\n')

    if name != names[-1]:
        mlist.write('#########################################################\n')
    else:
        mlist.write('#########################################################')

mlist.close()
    
