function [spclib_relab_actv] = load_spclib_relab(opt_relab)
    % relab (some that are not in the crism spectral library)
switch opt_relab
    case 1
        % opal
        % hydrated silica and gypsum: referred in the CRISM Type
        % Spectral Library
        load('spclib_relab2016Dec.mat','spclib_relab');
        opal_smpIDList = {'JB-JLB-874','JB-JLB-A33','JM-TGS-066',...
            'NV-RBS-022','OP-MCG-001','OP-MCG-002','OP-MCG-003',...
            'OP-MCG-004','OP-MCG-005','OP-MCG-006','OP-MCG-007',...
            'OP-MCG-008','OP-MCG-009','OP-MCG-010','OP-MCG-011',...
            'OP-MCG-012','SI-DWM-008-A','SI-DWM-008-B'};
        opal = searchby('sampleID',opal_smpIDList,spclib_relab,...
                        'COMP_FUNC','strcmpi');
        gypsum = searchby('spectrumID','LASF41A',spclib_relab,...
                          'COMP_FUNC','strcmpi');
        hyd_silica = searchby('spectrumID','BKR1JC329',spclib_relab);
        spclib_relab_actv = [opal gypsum hyd_silica];

    case 2
        % only opal
        load('spclib_relab2016Dec.mat','spclib_relab');
        opal_smpIDList = {'JB-JLB-874','JB-JLB-A33','JM-TGS-066',...
            'NV-RBS-022','OP-MCG-001','OP-MCG-002','OP-MCG-003',...
            'OP-MCG-004','OP-MCG-005','OP-MCG-006','OP-MCG-007',...
            'OP-MCG-008','OP-MCG-009','OP-MCG-010','OP-MCG-011',...
            'OP-MCG-012','SI-DWM-008-A','SI-DWM-008-B'};
        opal = searchby('sampleID',opal_smpIDList,spclib_relab,...
                        'COMP_FUNC','strcmpi');
        spclib_relab_actv = [opal];

    case 3
        % opal
        % hydrated silica and gypsum: referred in the CRISM Type
        % Spectral Library
        load('spclib_relab2016Dec.mat','spclib_relab');
        opal_smpIDList = {'JB-JLB-874','JB-JLB-A33','JM-TGS-066',...
            'NV-RBS-022','OP-MCG-001','OP-MCG-002','OP-MCG-003',...
            'OP-MCG-004','OP-MCG-005','OP-MCG-006','OP-MCG-007',...
            'OP-MCG-008','OP-MCG-009','OP-MCG-010','OP-MCG-011',...
            'OP-MCG-012','SI-DWM-008-A','SI-DWM-008-B'};
        szomolnokite_smpIDList =...
            {'JB-JLB-622-A','JB-JLB-824','JB-JLB-B66','JB-JLB-B67-A',...
             'JB-JLB-B67-B','LH-JFM-026'};
        melanterite_smpIDList = {'LH-JFM-041','SF-EAC-044-A'};
        opal = searchby('sampleID',opal_smpIDList,spclib_relab,...
                        'COMP_FUNC','strcmpi');
        gypsum = searchby('spectrumID','LASF41A',spclib_relab,...
                          'COMP_FUNC','strcmpi');
        szomolnokite = searchby('sampleID',szomolnokite_smpIDList,spclib_relab);
        melanterite = searchby('sampleID',melanterite_smpIDList,spclib_relab);
        hyd_silica = searchby('spectrumID','BKR1JC329',spclib_relab);
        spclib_relab_actv = [opal gypsum hyd_silica szomolnokite,melanterite];

    case 4
        load('spclib_relab2016Dec.mat','spclib_relab');
        idx = and([spclib_relab.wavelength_strt]<1000,[spclib_relab.wavelength_end]>2800);
        spclib_relab_actv = spclib_relab(idx);

    case 5
        load('spclib_relab2016Dec.mat','spclib_relab');
        opal_smpIDList = {'JB-JLB-874','JB-JLB-A33','JM-TGS-066',...
            'NV-RBS-022','OP-MCG-001','OP-MCG-002','OP-MCG-003',...
            'OP-MCG-004','OP-MCG-005','OP-MCG-006','OP-MCG-007',...
            'OP-MCG-008','OP-MCG-009','OP-MCG-010','OP-MCG-011',...
            'OP-MCG-012','SI-DWM-008-A','SI-DWM-008-B'};
        opal = searchby('sampleID',opal_smpIDList,spclib_relab,...
                        'COMP_FUNC','strcmpi');
        gypsum = searchby('spectrumID','LASF41A',spclib_relab,...
                          'COMP_FUNC','strcmpi');
        hyd_silica = searchby('spectrumID','BKR1JC329',spclib_relab);
        load('spclib_relab_jess_exclusive.mat','spclib_relab_jess_exclusive');

        spclib_relab_actv = [opal gypsum hyd_silica spclib_relab_jess_exclusive];
    case 0
        spclib_relab_actv = []; % empty is added.
    case 6
        load('spclib_relab2018Dec31.mat','spclib_relab');
        wv_strt = cell2mat([spclib_relab.wavelength_strt]);
        wv_end = cell2mat([spclib_relab.wavelength_end]);
        idx = and(wv_strt<900,wv_end>2800);
        generalType1_list = {'mineral','Rock','Rocks','RockCoating','MineralPwdr','Coating','Regolith','Sediment','Soil','Ash'};
        spclib_relab_wv = spclib_relab(idx);
        spclib_relab_actv = searchby('generalType1',generalType1_list,spclib_relab_wv);        
    otherwise
        error('opt_relab %d is not defined',opt_relab);
end

if ~isempty(spclib_relab_actv)
    [spclib_relab_actv.spclib] = deal('RELAB spectral library');
    iList = num2cell(1:length(spclib_relab_actv));
    [spclib_relab_actv.cumindex] = iList{:};
end
    
end