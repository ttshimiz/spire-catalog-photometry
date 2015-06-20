# Standard scientific imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.special import hyp2f1
from mmm import mmm

# Astropy imports
import astropy.units as u
import astropy.coordinates as coord
import astropy.io.fits as pyf
from astropy import wcs
from astropy.stats import sigma_clipped_stats, gaussian_sigma_to_fwhm

# Photutils imports
from photutils import detect_sources, segment_properties, aperture_photometry, remove_segments
from photutils import CircularAperture, EllipticalAperture
from photutils import SkyCircularAperture, SkyEllipticalAperture
from photutils import CircularAnnulus, EllipticalAnnulus
from photutils import SkyCircularAnnulus, SkyEllipticalAnnulus
from photutils.geometry import circular_overlap_grid, elliptical_overlap_grid
from photutils.aperture_funcs import get_phot_extents

# Plotting imports
import aplpy

# Saving modules
import datetime as dt
import pickle

# Global variables
PIX_SIZES = {'PSW':4.50, 'PMW': 6.25, 'PLW': 9.00}
FWHM = {'PSW':17.6, 'PMW': 24.0, 'PLW': 35.0}
BEAM_AREAS = {'PSW': 469.7, 'PMW': 831.7, 'PLW': 1793.5}
PS_SRC_APS = {'PSW': 22., 'PMW': 30., 'PLW': 42.}
PS_BKG_AP = {'r_in': 60., 'r_out': 90.}
AP_CORR = {'PSW': 1.2697, 'PMW': 1.2271, 'PLW': 1.2194}
CALIBRATION_ERR = 0.095
WAVES = {'PSW': 250, 'PMW': 350, 'PLW': 500}

# Directories
bat_dir = '/Users/ttshimiz/Github/bat-data/'
image_dir = '/Users/ttshimiz/Dropbox/Herschel_Images/SPIRE_reprocessed/'

# Function to estimate the global background
def estimate_bkg(image, clip=3.0, snr=2.0, npix=5, tol=1e-6, max_iter=10):

    # Create the NaN mask
    nanmask = np.isnan(image)
    
    # First estimate of the global background from sigma-clipping
    im_mean0, im_med0, im_std0 = sigma_clipped_stats(image, sigma=clip, mask=nanmask)
    
    # Create segmentation image using threshold = im_med + snr*im_std
    # Greater than npix pixels have to be connected to be a source
    thresh = im_med0 + snr*im_std0
    segm_img = detect_sources(image, thresh, npixels=npix)
    
    # Create new mask that masks all the detected sources
    mask = segm_img.astype(np.bool) | nanmask
    
    # Re-calculate background and rms
    im_mean, im_med, im_std = sigma_clipped_stats(image, sigma=clip, mask=mask)
    
    # Calculate percent change in the background estimate
    perc_change = np.abs(im_med0 - im_med)/ im_med
    
    if perc_change > tol:
    	im_med0 = im_med
        for i in range(max_iter):
 
    	    thresh = im_med + snr*im_std
            segm_img = detect_sources(image, thresh, npixels=npix)
            mask = segm_img.astype(np.bool) | nanmask
            im_mean, im_med, im_std = sigma_clipped_stats(image, sigma=clip, mask=mask)
            perc_change = np.abs(im_med0 - im_med)/ im_med

            if perc_change < tol:
            	break
            else:
            	im_med0 = im_med
            	
    return im_med, im_std


# Function to find the BAT counterpart in the segmentation image
def find_bat_source(coord_bat, props, rdist):

    ra_src = [p.ra_icrs_centroid for p in props]
    dec_src = [p.dec_icrs_centroid for p in props]

    coord_src = coord.SkyCoord(ra_src, dec_src, frame='icrs')
    dist = coord_bat.separation(coord_src).to(u.arcsec)
    ind_closest = np.argmin(dist)
    cdist = np.min(dist).value

    if cdist < rdist:
        return ind_closest
    else:
        return None

def create_spire_aperture(prop, filter, extent=3., coord_bat=None, wcs=None):
    
    if prop is None:
        position = wcs.all_world2pix([[coord_bat.ra.value, coord_bat.dec.value]], 0)
        ap = CircularAperture(position, PS_SRC_APS[filter]/PIX_SIZES[filter])
        type = 'undetected'
    else:
        position = [prop.xcentroid.value, prop.ycentroid.value]
        a = prop.semimajor_axis_sigma.value * extent
        b = prop.semiminor_axis_sigma.value * extent
        theta = prop.orientation.value
    
        if a < PS_SRC_APS[filter]/PIX_SIZES[filter]:
            ap = CircularAperture(position, PS_SRC_APS[filter]/PIX_SIZES[filter])
            type = 'point'
        else:
            ap = EllipticalAperture(position, a, b, theta=theta)
            type = 'extended'
    
    return ap, type

    
def calc_bkg_rms(ap, image, type, filter, mask=None, min_ap=6):

    rpsrc = PS_SRC_APS[filter]/PIX_SIZES[filter]

    if (type == 'point') | (type == 'undetected'):
        aback = bback =  PS_BKG_AP['r_in']/PIX_SIZES[filter] + rpsrc
        ap_theta = 0
    elif type == 'extended':
        aback = ap.a + rpsrc
        bback = ap.b + rpsrc
        ap_theta = ap.theta

    ecirc = ellip_circumference(aback, bback)
    diam = 2*rpsrc
    
    # Estimate the number of background apertures that can fit around the source
    # aperture.
    naps = np.int(np.round(ecirc/diam))
    
    # Use a minimum of 6 apertures
    naps = np.max([naps, min_ap])
    
    theta_back = np.linspace(0, 2*np.pi, naps, endpoint=False)
    
    # Get the x, y positions of the background apertures
    x, y = ellip_point(ap.positions[0], aback, bback, ap_theta, theta_back)

    # Create the background apertures and calculate flux within each
    bkg_aps = CircularAperture(np.vstack([x,y]).T, rpsrc)
#    flux_bkg = aperture_photometry(image, bkg_aps, mask=mask)
    flux_bkg, area_bkg = calc_masked_aperture(bkg_aps, image.data, method='sum', mask=mask)
    flux_bkg_adj = flux_bkg/area_bkg * bkg_aps.area()
	
	# Use sigma-clipping to determine the RMS of the background
	# Scale to the area of the source aperture
    me, md, sd = sigma_clipped_stats(flux_bkg_adj, sigma=3)
    bkg_rms = sd/bkg_aps.area()*ap.area()
    
    return bkg_rms, bkg_aps


def prep_image(im, filt):

    im.data *= im.header['PFOV']**2/BEAM_AREAS[filt]
    im.header['BUNIT'] = 'Jy/pixel'
    
    return im

# Function to calculate the circumference of an ellipse
def ellip_circumference(a, b):
    t = ((a-b)/(a+b))**2
    return np.pi*(a+b)*hyp2f1(-0.5, -0.5, 1, t)


# Function to determine x, y, position of an ellipse
def ellip_point(pos, a, b, theta, alpha):
    x = a*np.cos(alpha)*np.cos(theta) - b*np.sin(alpha)*np.sin(theta) + pos[0]
    y = a*np.cos(alpha)*np.sin(theta) + b*np.sin(alpha)*np.cos(theta) + pos[1]    
    return x, y


def calc_masked_aperture(ap, image, method='mmm', mask=None):

    positions = ap.positions
    extents = np.zeros((len(positions), 4), dtype=int)
    
    if isinstance(ap, EllipticalAnnulus):
        radius = ap.a_out
    elif isinstance(ap, CircularAnnulus):
        radius = ap.r_out
    elif isinstance(ap, CircularAperture):
        radius = ap.r
    
    extents[:, 0] = positions[:, 0] - radius + 0.5
    extents[:, 1] = positions[:, 0] + radius + 1.5
    extents[:, 2] = positions[:, 1] - radius + 0.5
    extents[:, 3] = positions[:, 1] + radius + 1.5
    
    ood_filter, extent, phot_extent = get_phot_extents(image, positions,
                                                       extents)
    
    x_min, x_max, y_min, y_max = extent
    x_pmin, x_pmax, y_pmin, y_pmax = phot_extent
    
    bkg = np.zeros(len(positions))
    area = np.zeros(len(positions))
    
    for i in range(len(bkg)):
        if isinstance(ap, EllipticalAnnulus):
            fraction = elliptical_overlap_grid(x_pmin[i], x_pmax[i],
                                               y_pmin[i], y_pmax[i],
                                               x_max[i] - x_min[i],
                                               y_max[i] - y_min[i],
                                               ap.a_out, ap.b_out, ap.theta, 0,
                                               1)
            b_in = ap.a_in * ap.b_out / ap.a_out
            fraction -= elliptical_overlap_grid(x_pmin[i], x_pmax[i],
                                                y_pmin[i], y_pmax[i],
                                                x_max[i] - x_min[i],
                                                y_max[i] - y_min[i],
                                                ap.a_in, b_in, ap.theta,
                                                0, 1)
        elif isinstance(ap, CircularAnnulus):
            fraction = circular_overlap_grid(x_pmin[i], x_pmax[i],
                                             y_pmin[i], y_pmax[i],
                                             x_max[i] - x_min[i],
                                             y_max[i] - y_min[i],
                                             ap.r_out, 0, 1)
            
            fraction -= circular_overlap_grid(x_pmin[i], x_pmax[i],
                                                y_pmin[i], y_pmax[i],
                                                x_max[i] - x_min[i],
                                                y_max[i] - y_min[i],
                                                ap.r_in, 0, 1)
        elif isinstance(ap, CircularAperture):
            fraction = circular_overlap_grid(x_pmin[i], x_pmax[i],
                                             y_pmin[i], y_pmax[i],
                                             x_max[i] - x_min[i],
                                             y_max[i] - y_min[i],
                                             ap.r, 0, 1)  
        
        pixel_data = image[y_min[i]:y_max[i], x_min[i]:x_max[i]] * fraction
        if mask is not None:
            pixel_data[mask[y_min[i]:y_max[i], x_min[i]:x_max[i]]] = 0
        
        good_pixels = pixel_data[pixel_data != 0].flatten()
        
        if method == 'mmm':
            skymod, skysigma, skew = mmm(good_pixels)
            bkg[i] = skymod
        elif method == 'sum':
            bkg[i] = np.sum(good_pixels)
        area[i] = len(good_pixels)

    return bkg, area

     
def spire_aperture_photometry(ap, image, error, type, filter, mask=None,
                              ellip_ann_fac=1.5):
	
    image.data = np.array(image.data, dtype=np.float)
    error = np.array(error, dtype=np.float)
	
    # Calculate the flux in the global background subtracted image using the input aperture
    flux_tbl = aperture_photometry(image, ap, error=error, mask=mask)
    
    # Create a background annulus aperture
    if (type == 'point') | (type == 'undetected'):
        
        bkg_ap = CircularAnnulus(ap.positions, r_in=PS_BKG_AP['r_in']/PIX_SIZES[filter], r_out=PS_BKG_AP['r_out']/PIX_SIZES[filter])
        
    elif type == 'extended':
    
        bkg_ap = EllipticalAnnulus(ap.positions, a_in=ap.a, a_out=ap.a*ellip_ann_fac,
                                   b_out=ap.b*ellip_ann_fac, theta=ap.theta)
    
    bkg_flux, bkg_area = calc_masked_aperture(bkg_ap, image.data, method='mmm', mask=mask)
    
    # Use a series of apertures around the source aperture to calculate the rms of the
    # background
    bkg_rms, bkg_circles = calc_bkg_rms(ap, image, type, filter, mask=mask)
    
    w = wcs.WCS(image.header)
    
    # Final fluxes and errors
    raw_flux = flux_tbl['aperture_sum'].data[0]
    bkg_flux = bkg_flux[0]*ap.area()
    bkgsub_flux = raw_flux - bkg_flux
    ap_err = flux_tbl['aperture_sum_err'].data[0]
    cal_err = bkgsub_flux * CALIBRATION_ERR
    total_err = np.sqrt(ap_err**2 + cal_err**2 + bkg_rms**2)
    
    result = {'raw_flux': raw_flux,
              'bkg_flux': bkg_flux,
              'bkgsub_flux': bkgsub_flux,
              'aperture_err': ap_err,
              'bkg_rms': bkg_rms,
              'calib_err': cal_err,
              'total_err': total_err}
    apertures = {'source': ap,
                 'background_annulus': bkg_ap,
                 'background_circles': bkg_circles}

    return result, apertures


def plot_aperture_photometry(image, apertures, filter, type, global_bkg=0, stretch='arcsinh',
                             figure=None, subplot=(1,1,1), colorscale='cubehelix',
                             pixel_max=None, title=None, plot_bat_loc=None):

    w = wcs.WCS(image.header)
    fig = aplpy.FITSFigure(image, figure=figure, subplot=subplot)
    if pixel_max is None:
        fig.show_colorscale(vmin=global_bkg, cmap=colorscale, stretch=stretch)
    else:
        fig.show_colorscale(vmin=global_bkg, vmax=pixel_max, cmap=colorscale,
                            stretch=stretch)

    src_ap = apertures['source']
    bkg_ann = apertures['background_annulus']
    bkg_circs = apertures['background_circles']
    
    # Recenter the image to the center of the source aperture
    x_c = src_ap.positions[0][0]
    y_c = src_ap.positions[0][1]
    ra_c, dec_c = fig.pixel2world(x_c+1, y_c+1)
    if type == 'extended':
    	radius = (bkg_ann.a_out*PIX_SIZES[filter] + PS_SRC_APS[filter])/3600.
    else:
        radius = (PS_BKG_AP['r_out'] + PS_SRC_APS[filter]*1.5)/3600.
    fig.recenter(ra_c, dec_c, radius=radius)
    
    # Plot the apertures
    # Need to adjust the pixel positions by adding 1
    src_ap.positions = src_ap.positions + 1
    bkg_ann.positions = bkg_ann.positions + 1
    bkg_circs.positions = bkg_circs.positions + 1
    
    src_ap.plot(color='w', lw=1.5)
    bkg_ann.plot(color='w', lw=1.5, ls='dashed')
    bkg_circs.plot(color='w', lw=1.5, ls='dashed')
    
    src_ap.positions = src_ap.positions - 1
    bkg_ann.positions = bkg_ann.positions - 1
    bkg_circs.positions = bkg_circs.positions - 1
    
    # Add a beam
    fig.add_beam(major=FWHM[filter]*u.arcsec, minor=FWHM[filter]*u.arcsec, angle=0)
    fig.beam.set_facecolor('white')
    fig.beam.set_hatch('/')
    
    # Add a scalebar
    fig.add_scalebar(30*u.arcsec)
    fig.scalebar.set_color('white')
    fig.scalebar.set_linewidth(3) 
    fig.scalebar.set_label('30"')
    
    # Add a colorbar
    fig.add_colorbar()
    fig.colorbar.set_axis_label_text('Jy/pixel')
    
    # Add title
    if title is not None:
        fig.add_label(0.05, 0.97, title, color='white', relative=True, size=20,
                      horizontalalignment='left', verticalalignment='top')
    
    # Add a marker for the BAT position
    if plot_bat_loc is not None:
        fig.show_markers(plot_bat_loc[0], plot_bat_loc[1], marker='+', color='r', s=40)
    
    return fig
    
 
def run_bat_sources(source, filter, save_results=False, outdir=None, plot=False, results_file=None,
                    ap_file=None, save_fig=False):
 
    bat_info = pd.read_csv(bat_dir+'bat_info.csv', index_col=0)
    
    if np.all(source != 'all'):
        if np.isscalar(source):
            srcs = [source]
        else:
            srcs = source
    else:
        srcs = bat_info.index.values
    
    if np.isscalar(filter):
        filter = [filter]
    
    # Setup the DataFrame to store the photometry results
    index = pd.MultiIndex.from_product([srcs, filter], names=['Name', 'Filter'])
    src_df = pd.DataFrame(columns=['raw_flux', 'bkg_flux', 'bkgsub_flux', 'total_err', 'aperture_err', 'bkg_rms', 'calib_err', 'type'], index=index)
    
    # Setup a dictionary to hold the apertures used
    src_aps = {}
    
    for s in srcs:
        for f in filter:
            
            # Load the data
            hdu_image = pyf.open(image_dir+s+'_scanamorphos_spire'+str(WAVES[f])+'_signal.fits')[0]
            hdu_err = pyf.open(image_dir+s+'_scanamorphos_spire'+str(WAVES[f])+'_error.fits')[0]
            
            # Convert to Jy/pixel
            hdu_image = prep_image(hdu_image, f)
            hdu_err = prep_image(hdu_err, f)
            
            # Calculate the global background
            im = hdu_image.data
            im_med, im_std = estimate_bkg(im)
            
            # Use a 1.5 sigma threshold to detect the BAT source
            thresh = im_med + 1.5*im_std
            segm_img = detect_sources(im, thresh, npixels=5)
            props = segment_properties(im-im_med, segm_img, wcs=wcs.WCS(hdu_image.header))
            
            # Find the BAT source in the properties list
            ra_bat = bat_info.loc[s, 'RA_(J2000)']
            dec_bat = bat_info.loc[s, 'DEC_(J2000)']
            coord_bat = coord.SkyCoord(ra=ra_bat, dec=dec_bat, frame='fk5')
            ind_bat = find_bat_source(coord_bat, props, FWHM[f]*2)
            
            if ind_bat is None:
                
                ap, type = create_spire_aperture(None, f, coord_bat=coord_bat, wcs=wcs.WCS(hdu_image.header)) 
            
            else:
            
                ap, type = create_spire_aperture(props[ind_bat], f, coord_bat=coord_bat)
            
            # Create new mask based on 3-sigma threshold and remove the bat source if detected
            thresh2 = im_med + 3*im_std
            segm_img2 = detect_sources(im, thresh2, npixels=5)
            props2 = segment_properties(im-im_med, segm_img2, wcs=wcs.WCS(hdu_image.header))
            nanmask = np.isnan(im)
            nanerr = np.isnan(hdu_err.data) | np.isinf(hdu_err.data)
            
            if len(props2) != 0:
                ind_bat2 = find_bat_source(coord_bat, props2, FWHM[f]*2)
            else:
                ind_bat2 = None
    
            if ind_bat2 is None:
                mask = segm_img2.astype(np.bool) | nanmask | nanerr
            else:
                si = remove_segments(segm_img2, ind_bat2 + 1)
                mask = si.astype(np.bool) | nanmask | nanerr
            
            results, apertures = spire_aperture_photometry(ap, hdu_image, hdu_err.data, type, f, mask=mask)
            results['type'] = type
            src_df.loc[s, f] = pd.Series(results)
            src_aps[s+'_'+f] = apertures  
            
            if plot:

                cbat = [coord_bat.ra.deg, coord_bat.dec.deg]
                
                if (type == 'undetected'):
                    pmax = None
                else:
                    pmax = props[ind_bat].max_value + im_med
                
                if (pmax < im_med):
                    pmax = None

                fig = plot_aperture_photometry(hdu_image, apertures, f, type, global_bkg=im_med,
                                                   pixel_max=pmax, title=s+' ['+f+']', plot_bat_loc=cbat)
           
                if save_fig:
                    if outdir is None:
                        fig.save(s+'_'+f+'.png')
                    else:
                        fig.save(outdir+s+'_'+f+'.png')
                    fig.close()
            
    if save_results:
        if (outdir is None) and (results_file is None):
            src_df.to_csv('photometry_results_'+str(dt.datetime.today().isoformat()))
            f_ap = open('apertures_'+str(dt.datetime.today().isoformat()), 'wb')
        elif (outdir is not None) and (results_file is None):
            src_df.to_csv(outdir+'photometry_results_'+str(dt.datetime.today().isoformat()))
            f_ap = open(outdir+'apertures_'+str(dt.datetime.today().isoformat()), 'wb')
        else:
            src_df.to_csv(outdir+results_file)
            f_ap = open(outdir+ap_file, 'wb')
        
        pickle.dump(src_aps, f_ap)
        f_ap.close()
    
    if plot and not save_fig:
        return src_df, src_aps, fig  
    else:
        return src_df, src_aps             