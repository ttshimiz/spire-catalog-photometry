# 
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2015 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
# 
###########################################################################
###            SPIRE Scan Map 2-Pass Pipeline Reprocessing Script       ###
###########################################################################
#  Purpose:  A simplified version of the SPIRE 2-Pass Pipeline 
#             This is for data reprocessing using the latest SPIRE calibration products.
#             
# 
#            The results are 
#            - a reprocessed observation
#            -  3 FITS files calibrated for point sources (Jy/beam) for the PSW, PMW, PLW arrays,
#            -  3 FITS files calibrated for extended emission (MJy/sr) for the PSW, PMW, PLW arrays,
#            -  Optionally 3 FITS files for SSO moving objects (Jy/beam) for the PSW, PMW, PLW arrays,
#            All individual FITS files contain the final image map, error map, coverage map.   
#            Creation of Extended Emission maps and SSO maps are controlled by the keywords 
#            createExtdMaps    default = True
#            createSSOMaps     default = False
#
#    
#  Usage:    The user needs to specify the options in the simple user input 
#            section at the beginning of the script;
#            - Observation ID (obsid)  
#            - Data pool name  (if accessing a data pool on disk)
#            - Output directory for final fits files  
#
# Note:  it is possible to save entire observation back to a pool by
#     uncommenting the saveProducts command at the end of the script
#
#  Assumptions: A data pool has already been created on disk (else access from HSA). 
#               For Extended Emission maps, the Planck maps need to be accesible
#               Planck map location is set by adding the "spire.spg.hfi.545map" 
#               and "spire.spg.hfi.857map" with the relevant path-to-file, in user.props
#
#  Updated: 19/02/2015
#  Updated: 02/06/2015: Taro Shimizu
#                       Modified to re-process all of the BAT AGN using the same pixel
#                       sizes as those used in Scanamorphos
###########################################################################


# Load in all of the OBSIDs for the BAT AGN
f_obsid = open('/Users/ttshimiz/Github/spire-catalog-photometry/bat_agn_spire_obsids.txt', 'r')
lines = f_obsid.readlines()[0]
lines = lines.split('\r')
obsids = [int(x.split('\t')[-1]) for x in lines]
#obsids = [1342245154]
sources = [x.split('\t')[0] for x in lines]
#sources = ['MCG -01-24-012']
myDataPool = "spire-bat-agn"
outDir     = "/Users/ttshimiz/Research/Thesis/data/Herschel/SPIRE_naive/"

# Calibration Context and Calibration Files 
# Read the latest calibration tree relevant to HCSS v13 from the local disc:
cal = spireCal(pool="spire_cal_13_1")
#cal = spireCal(calTree="spire_cal_13_1", saveTree=True)

# Global processing environment options
saveMaps          = True       # Save all final maps to outDir path
tempStorage       = True       # Use Temporary Storage on disk if obs is very large
includeTurnaround = True       # Include the scan line turnarounds (recommended)
createExtdMaps    = True       # Created extended emission maps in MJy/sr (requires the Planck maps)
createSSOMaps     = False      # Create moving object frame maps for SSO

# Deglitching Parameters
l2DeglitchRepeat = 100
kappa = 5.0
kappa2 = 5.0

# Destriping Parameters
offsetFunction = "perScan"
polyDegree = 0
withMedianCorrected = True
nThreads = 2
jumpThresh =-1.0
jumpIter = 100
brightSourceThresh = 1.5
roi = 0

# Mapmaking Parameters
pswSize = 4.5      # Recommended map pixel size for PSW
pmwSize = 6.25     # Recommended map pixel size for PMW
plwSize = 9.0     # Recommended map pixel size for PLW
minVel  = 5.0      # Recommended min scan velocity to be included in map

# Extended Emission map creation Planck ZeroPointCorrection Parameters
# The Planck Maps are required in order to produce the maps calibrated for extended emission
# change Paths to point to Planck maps on YOUR local disk
hfi545Gain = zeroPointCorrection["hfi545Gain"]  # Recommended gain for Planck HFI 545GHz channel
hfi857Gain = zeroPointCorrection["hfi857Gain"]  # Recommended gain for Planck HFI 857Hz channel
hfiFwhm    = zeroPointCorrection["hfiFwhm"]     # Recommended Planck HFI FWHM
#*************************************************************************
# Run the loop for all of the BAT AGN.
for i in range(len(obsids)):
#for i in [0]:
    myObsid = obsids[i]
    src = sources[i]
	
#*************************************************************************
##  Load in an observation context from your data pool into HIPE:
    obs=getObservation(myObsid,useHsa=True,instrument="SPIRE")           # from the HSA
    #obs=getObservation(myObsid,poolName=myDataPool,instrument="SPIRE")   # from a pool
    print
    print "Processing %s observation %i (0x%X)"%(src, myObsid, myObsid)


#*************************************************************************
# Check that the data are really SPIRE data
    if obs.instrument != "SPIRE": 
        raise BadDataException("This ObservationContext cannot be processed with this pipeline: it contains "+obs.instrument+" data, not SPIRE data")

#*************************************************************************
# RUN THE modified 2 Pass PIPELINE
    spirePhotPipeline(obs, cal=cal,  \
                           tempStorage=tempStorage, includeTurnaround=includeTurnaround, \
                           createExtdMaps=createExtdMaps, createSSOMaps=createSSOMaps,   \
                           l2DeglitchRepeat=l2DeglitchRepeat, kappa=kappa, kappa2=kappa2,\
                           offsetFunction=offsetFunction, polyDegree=polyDegree,         \
                           withMedianCorrected=withMedianCorrected, nThreads=nThreads,   \
                           jumpThresh=jumpThresh, jumpIter=jumpIter,                     \
                           brightSourceThresh=brightSourceThresh,                        \
                           pswSize=pswSize, pmwSize=pmwSize, plwSize=plwSize,            \
                           minVel=minVel,                                                \
                           hfi545Gain=hfi545Gain, hfi857Gain=hfi857Gain, hfiFwhm=hfiFwhm)

#*************************************************************************


#*************************************************************************
# Save Maps to output directory
    if saveMaps:
        simpleFitsWriter(obs.level2.refs["psrcPSW"].product, "%s%s_psrcPSW.fits"%(outDir, src))
        simpleFitsWriter(obs.level2.refs["psrcPMW"].product, "%s%s_psrcPMW.fits"%(outDir, src))
        simpleFitsWriter(obs.level2.refs["psrcPLW"].product, "%s%s_psrcPLW.fits"%(outDir, myObsid))
        print "Map saved as FITS files to %s"%(outDir)
        #
        if createExtdMaps:
            simpleFitsWriter(obs.level2.refs["extdPSW"].product, "%s%s_extdPSW.fits"%(outDir, src))
            simpleFitsWriter(obs.level2.refs["extdPMW"].product, "%s%s_extdPMW.fits"%(outDir, src))
            simpleFitsWriter(obs.level2.refs["extdPLW"].product, "%s%s_extdPLW.fits"%(outDir, src))
		#
        if createSSOMaps:
            simpleFitsWriter(obs.level2.refs["ssoPSW"].product, "%sssoPSW_%i.fits"%(outDir, src))
            simpleFitsWriter(obs.level2.refs["ssoPMW"].product, "%sssoPMW_%i.fits"%(outDir, src))
            simpleFitsWriter(obs.level2.refs["ssoPLW"].product, "%sssoPLW_%i.fits"%(outDir, src))
	#

	#Save the entire observation to your Local Store pool
	# Finally we can save the new reprocessed observation back to your local pool on your hard disk
	# Uncomment the next line and choose a poolName, either the existing one or a new one
	#
        saveProduct(product=obs, pool=myDataPool)

	#
        print
        print "Completed the processing of %s, OBSID= %i"%(src, myObsid)

#### End of the script ####