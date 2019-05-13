[obs_dirname_list] = read_obs_dirname_list('obs_dirname_list_ffc_OlymMons');


for i=1:length(obs_dirname_list)
    obs_dirname = obs_dirname_list{i};
    obs_id = obs_dirname(4:11);
    crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L','download_trrif',2,...
        'download_trrra',2,'download_edrscdf',2);
end

%% create ffc info
L = length(obs_dirname_list);
ffc_info = struct('dirname',cell(L,1),...
    'sclk_start_01',cell(L,1),'sclk_stop_01',cell(L,1),...
    'existSC01',cell(L,1),'existTRR3IF01',cell(L,1),'existTRR3RA01',cell(L,1),...
    'sclk_start_03',cell(L,1),'sclk_stop_03',cell(L,1),...
    'existSC03',cell(L,1),'existTRR3IF03',cell(L,1),'existTRR3RA03',cell(L,1));

for i=1:length(ffc_info)
    dirname = obs_dirname_list{i};
    ffc_info(i).dirname = dirname;
    obs_id = dirname(4:11);
    crism_obs = CRISMObservation(obs_id,'SENSOR_ID','L');
    
    % scene 01 ------------------------------------------------------------
    [EDRSCdata1] = get_scene_CRISMdata_FFC(crism_obs.info.basenameSC,'',1);
    if ~isempty(EDRSCdata1)
        ffc_info(i).existSC01 = true;
        ffc_info(i).sclk_start_01 = EDRSCdata1.get_sclk_start();
        ffc_info(i).sclk_stop_01 = EDRSCdata1.get_sclk_stop();
    else
        ffc_info(i).existSC01 = false;
    end
    [TRRIFdata1] = get_scene_CRISMdata_FFC(crism_obs.info.basenameIF,'',1);
    ffc_info(i).existTRR3IF01 = ~isempty(TRRIFdata1);
    [TRRRAdata1] = get_scene_CRISMdata_FFC(crism_obs.info.basenameRA,'',1);
    ffc_info(i).existTRR3RA01 = ~isempty(TRRRAdata1);
    
    % scene 03 ------------------------------------------------------------
    [EDRSCdata3] = get_scene_CRISMdata_FFC(crism_obs.info.basenameSC,'',3);
    if ~isempty(EDRSCdata3)
        ffc_info(i).existSC03 = true;
        ffc_info(i).sclk_start_03 = EDRSCdata3.get_sclk_start();
        ffc_info(i).sclk_stop_03 = EDRSCdata3.get_sclk_stop();
    else
        ffc_info(i).existSC03 = false;
    end
    [TRRIFdata3] = get_scene_CRISMdata_FFC(crism_obs.info.basenameIF,'',3);
    ffc_info(i).existTRR3IF03 = ~isempty(TRRIFdata3);
    [TRRRAdata3] = get_scene_CRISMdata_FFC(crism_obs.info.basenameRA,'',3);
    ffc_info(i).existTRR3RA03 = ~isempty(TRRRAdata3);
    
end