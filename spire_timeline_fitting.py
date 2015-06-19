#Script to run the timeline fitter on the BAT sources

github_dir = '/ricci9nb/tshimizu/Github/'
# Upload file with OBSIDs
file = open(github_dir+'spire-catalog-photometry/bat_agn_spire_obsids.txt', 'r')
lines_obsid = file.read().splitlines()
file.close()

# Results file
fn_results = open(github_dir+'spire-catalog-photometry/point_source_timeline_fitting_results.txt', 'w')

# Upload coordinates for BAT AGN
file_bat = open(github_dir+'bat-data/bat_info.csv', 'r')
lines_bat = file_bat.read().splitlines()[1:]
file_bat.close()

# Upload the visual classifications for each source
file_class = open(github_dir+'spire-catalog-photometry/bat_spire_visual_class.csv', 'r')
lines_class = file_class.read().splitlines()[1:]
file_class.close()

for i in range(len(lines_obsid)):
	
	parse = lines_bat[i].split(',')
	name = parse[0]
	ra = float(parse[3])
	dec = float(parse[4])
	cl = lines_class[i].split(',')[1]

	obsid = int((lines_obsid[i].split()[1]))
	print name
	obs = getObservation(obsid, useHsa=True)
	
	# Use different settings for sources in Cirrus dominated regions
	if ((cl == 'Cirrus') | (cl == 'Cirrus (1)')):
		allowTiltBack = True
		fitEllipticalGauss2d = False
		useBackInFit = False
		allowVaryBack = False
	else:
		allowTiltBack = False
		fitEllipticalGauss2D = False
		useBackInFit = True
		allowVaryBack = True

	timeline_250 = sourceExtractorTimeline(input=obs.refs["level1"].product,
					       array='PSW', rPeak=22.0,
					       inputSourceList=[ra, dec],
					       allowTiltBack=allowTiltBack,
					       fitEllipticalGauss2d=fitEllipticalGauss2D,
					       useBackInFit=useBackInFit,
					       allowVaryBackground=allowVaryBack)
	timeline_350 = sourceExtractorTimeline(input=obs.refs["level1"].product,
					       array='PMW', rPeak=30.0,
					       inputSourceList=[ra, dec],
					       allowTiltBack=allowTiltBack,
					       fitEllipticalGauss2d=fitEllipticalGauss2D,
					       useBackInFit=useBackInFit,
					       allowVaryBackground=allowVaryBack)
	timeline_500 = sourceExtractorTimeline(input=obs.refs["level1"].product,
					       array='PLW', rPeak=42.0,
					       inputSourceList=[ra, dec],
					       allowTiltBack=allowTiltBack,
					       fitEllipticalGauss2d=fitEllipticalGauss2D,
					       useBackInFit=useBackInFit,
					       allowVaryBackground=allowVaryBack)
	
	flux_250 = timeline_250['sources']['flux'].data[0]
	flux_err_250 = timeline_250['sources']['fluxPlusErr'].data[0]
	fwhm_250 = timeline_250['sources']['sigma'].data[0]*3600*2.3548
	fwhm_err_250 = timeline_250['sources']['sigmaErr'].data[0]*3600*2.3548
	ra_250 = timeline_250['sources']['ra'].data[0]
	dec_250 = timeline_250['sources']['dec'].data[0]
	rchi2_250 = timeline_250['sources']['reducedChiSquare'].data[0]

	flux_350 = timeline_350['sources']['flux'].data[0]
	flux_err_350 = timeline_350['sources']['fluxPlusErr'].data[0]
	fwhm_350 = timeline_350['sources']['sigma'].data[0]*3600*2.3548
	fwhm_err_350 = timeline_350['sources']['sigmaErr'].data[0]*3600*2.3548
	ra_350 = timeline_350['sources']['ra'].data[0]
	dec_350 = timeline_350['sources']['dec'].data[0]
	rchi2_350 = timeline_350['sources']['reducedChiSquare'].data[0]

	flux_500 = timeline_500['sources']['flux'].data[0]
	flux_err_500 = timeline_500['sources']['fluxPlusErr'].data[0]
	fwhm_500 = timeline_500['sources']['sigma'].data[0]*3600*2.3548
	fwhm_err_500 = timeline_500['sources']['sigmaErr'].data[0]*3600*2.3548
	ra_500 = timeline_500['sources']['ra'].data[0]
	dec_500 = timeline_500['sources']['dec'].data[0]
	rchi2_500 = timeline_500['sources']['reducedChiSquare'].data[0]
	
	fn_results.write(name+'\t'+str(ra_250)+'\t'+str(dec_250)+'\t'+str(flux_250)+'\t'+str(flux_err_250)+'\t'+str(fwhm_250)+'\t'+str(fwhm_err_250)+'\t'+str(rchi2_250)+'\t'+
                         str(ra_350)+'\t'+str(dec_350)+'\t'+str(flux_350)+'\t'+str(flux_err_350)+'\t'+str(fwhm_350)+'\t'+str(fwhm_err_350)+'\t'+str(rchi2_350)+'\t'+
                         str(ra_500)+'\t'+str(dec_500)+'\t'+str(flux_500)+'\t'+str(flux_err_500)+'\t'+str(fwhm_500)+'\t'+str(fwhm_err_500)+'\t'+str(rchi2_500)+'\n')
fn_results.close()