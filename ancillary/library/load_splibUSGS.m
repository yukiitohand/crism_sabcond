function [splib_usgs_actv] = load_splibUSGS(opt_splibUSGS)
    % usgs spectral library
switch opt_splibUSGS
    case 1
        % samples referred in the CRISM Type Spectral Library
        load('splib06a.mat','splib06a');
        % irecno:  1930,  'Analcime GDS1 Zeolite        W1R1Ba AREF'
        % irecno:  3422,  'Bassanite GDS145 (syn)       W1R1Ba AREF'
        % irecno:  8863,  'Halite HS433.3B              W1R1Ba AREF'
        % irecno:  7546,  'Epidote GDS26.a 75-200um     W1R1Bb AREF'
        % irecno:  7560,  'Epidote GDS26.b <75um        W1R1Bb AREF'
        % irecno: 26841,  'H2O-Ice GDS136 77K           W1R1Ba AREF'
        % irecno:  9282,  'Hematite GDS27               W1R1Ba AREF'
        % irecno: 13344,  'Margarite GDS106             W1R1Bc AREF'
        % irecno: 21499,  'Talc GDS23 74-250um fr       W1R1Ba AREF'
        irecnoList = [1930 3422 8863 7546 7560 26841 9282 13344 21499];
        splib_usgs_actv = searchby('irecno',irecnoList,splib06a);
        [splib_usgs_actv.spclib] = deal('USGS splib06a');
        [splib_usgs_actv.cumindex] = deal(1:length(splib_usgs_actv));
    case 2
        splib_usgs_actv = [];
    case 3
        load('splib07b.mat','splib07b');
        % irecno:  1607,  'Analcime GDS1 Zeolite         BECKa AREF'
        % irecno:  2879,  'Bassanite GDS145 (syn)        BECKa AREF'
        % irecno:  7518,  'Halite HS433.3B               BECKa AREF'
        % irecno:  6222,  'Epidote GDS26.a 75-200um      BECKb AREF'
        % irecno:  6269,  'Epidote GDS26.b <75um         BECKb AREF'
        % irecno: 22531,  'H2O-Ice GDS136 77K            BECKa AREF'
        % irecno:  7928,  'Hematite GDS27                BECKa AREF'
        % irecno: 11092,  'Margarite GDS106              BECKc AREF'
        % irecno: 17452,  'Talc GDS23 74-250um fr        BECKa AREF'

        % irecno: 21801,  'Opalized Tuff CU00-15E       ASDFRb AREF'
        irecnoList = [1607 2879 7518 6222 6269 22531 7928 11092 17452,...
            21801];
        splib_usgs_actv1 = searchby('irecno',irecnoList,splib07b);
        splib_usgs_actv2 = searchby('documentation_format','MIXTURE',splib07b);
        splib_usgs_actv2_1 = [];
        for i=1:length(splib_usgs_actv2)
            if ~isempty(regexp(splib_usgs_actv2(i).ititl,'BECK[a-c]{1} [A-Z]{1}REF'))
                splib_usgs_actv2_1 = [splib_usgs_actv2_1 splib_usgs_actv2(i)];
            end
        end
        splib_usgs_actv = [splib_usgs_actv1 splib_usgs_actv2_1];
        [splib_usgs_actv.spclib] = deal('USGS splib07b');
    case 4
        load('s07CRSMj.mat', 's07CRSMj');
        [splib_usgs_actv,j] = searchby('documentation_format',{'MINERAL','MIXTURE'},s07CRSMj);
        splib_usgs_actv1 = [];
        for i=1:length(splib_usgs_actv)
            if ~isempty(regexp(splib_usgs_actv(i).ititl,'BECK[a-c]{1} [A-Z]{1}REF'))
                splib_usgs_actv1 = [splib_usgs_actv1 splib_usgs_actv(i)];
            end
        end
        splib_usgs_actv = splib_usgs_actv1;
        [splib_usgs_actv.spclib] = deal('USGS s07CRSMj');
    case 5
        load('splib07b.mat','splib07b');
        idx_actv = cellfun(@(x) ~isempty(regexpi(x,'BECK[a-c]{1} [A-Z]{1}REF','ONCE')),{splib07b.ititl});
        splib_usgs_actv = splib07b(idx_actv);
        [splib_usgs_actv.spclib] = deal('USGS splib07b');
    case 6
        % no ice...
        % samples referred in the CRISM Type Spectral Library
        load('splib06a.mat','splib06a');
        % irecno:  1930,  'Analcime GDS1 Zeolite        W1R1Ba AREF'
        % irecno:  3422,  'Bassanite GDS145 (syn)       W1R1Ba AREF'
        % irecno:  8863,  'Halite HS433.3B              W1R1Ba AREF'
        % irecno:  7546,  'Epidote GDS26.a 75-200um     W1R1Bb AREF'
        % irecno:  7560,  'Epidote GDS26.b <75um        W1R1Bb AREF'
        % %  irecno: 26841,  'H2O-Ice GDS136 77K           W1R1Ba AREF'
        % irecno:  9282,  'Hematite GDS27               W1R1Ba AREF'
        % irecno: 13344,  'Margarite GDS106             W1R1Bc AREF'
        % irecno: 21499,  'Talc GDS23 74-250um fr       W1R1Ba AREF'
        irecnoList = [1930 3422 8863 7546 7560 9282 13344 21499];
        splib_usgs_actv = searchby('irecno',irecnoList,splib06a);
        [splib_usgs_actv.spclib] = deal('USGS splib06a');
        [splib_usgs_actv.cumindex] = deal(1:length(splib_usgs_actv));
    otherwise
        error('opt_splibUSGS %d is not defined',opt_splibUSGS);
end

if ~isempty(splib_usgs_actv)
    iList = num2cell(1:length(splib_usgs_actv));
    [splib_usgs_actv.cumindex] = iList{:};
end

end