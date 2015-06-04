PRO create_scanam_scan_files
; Load in the names of the BAT AGN
    readcol, 'bat_agn_spire_obsids.txt', names, obsids, DELIMITER=string(9b), FORMAT='A,UL'
;    names = ['CenA', 'IC4329A']
    n_src = n_elements(names)

; Directories
    d_in = '/ricci5nb/tshimizu/L1_HIPE13/'
    d_out = '/ricci5nb/tshimizu/input_scans/'

    for i = 0, n_src-1 do begin

;   Determine the number of scans
        name=names[i]
        scanlegs = FILE_SEARCH('/ricci5nb/tshimizu/L1_HIPE13/'+name+'*')
        scans = strarr(n_elements(scanlegs))

        for j = 0, n_elements(scanlegs)-1 do begin
            split = strsplit(scanlegs[j], '_', /extract)
            scans[j] = split[5]
        endfor

        nscans = n_elements(uniq(scans))

        for k = 0, nscans - 1 do convert_hcssfits_spire, dir_in=d_in, dir_out=d_out, root=name+'_processed_subscans_level1_scan0'+chain(k+1)

    endfor
end
