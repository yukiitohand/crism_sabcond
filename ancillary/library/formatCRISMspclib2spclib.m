function [spclib_crism] = formatCRISMspclib2spclib(CRISMspclib)
%  [spclib_jesslib] = formatCRISMspclib2spclib(CRISMspclib)
%    convert the format of CRISMspclib read by "readCRISMspclib" to the
%    normal one depth spectral library format
%  Input Parameters
%    CRISMspclib: output of readCRISMspclib
%  Output Parameters
%    spclib_crism: struct, having fields
%         'product_id'
%         'product_name'
%         'name'            
%         'spc_name'   
%         'wavlength'
%         'reflectance' 
%         'wavelength_units'
%         'subfolder'
%         'filename' 
%         'LineIdx'

spclib_crism = [];
flds= fieldnames(CRISMspclib);
for i=1:length(flds)
    fld = flds{i};
    CRISMspclibi = CRISMspclib.(fld);
    L = CRISMspclibi.hdr.lines;
    spclibi = struct(...
        'product_id'      , repmat({''},[L,1]),...
        'product_name'    , repmat({''},[L,1]),...
        'name'            , repmat({''},[L,1]),...
        'spc_name'    , repmat({''},[L,1]),...
        'wavelength'       , repmat({''},[L,1]),...
        'reflectance'     , repmat({''},[L,1]),...
        'wavelength_units', repmat({''},[L,1]),...
        'subfolder'       , repmat({''},[L,1]),...
        'filename'        , repmat({''},[L,1]),...
        'LineIdx'         , repmat({''},[L,1])...
        );
    % exception: 47&48 in synth
    if strcmpi(fld,'synth')
        CRISMspclibi.hdr.spectra_names{47} = 'NA-JAROSITE BKR1JB440';
        CRISMspclibi.hdr.spectra_names{48} = 'NA-JAROSITE C1JB440';
    end
    [name_info] = get_spectrumIDfromspectra_namesCRISMspclib(CRISMspclibi.hdr.spectra_names);
    for j=1:L
        spclibi(j).product_id       = name_info(j).product_id;
        spclibi(j).product_name       = name_info(j).product_name;
        spclibi(j).name             = name_info(j).name;
        spclibi(j).spc_name         = CRISMspclibi.hdr.spectra_names{j};
        spclibi(j).wavelength       = CRISMspclibi.hdr.wavelength;
        spclibi(j).reflectance      = CRISMspclibi.spc(j,:);
        spclibi(j).wavelength_units = CRISMspclibi.hdr.wavelength_units;
        spclibi(j).subfolder        = CRISMspclibi.subfolder;
        spclibi(j).filename         = fld;
        spclibi(j).LineIdx          = j;
    end
    
    spclib_crism = merge_struct(spclib_crism,spclibi);
    
end

end