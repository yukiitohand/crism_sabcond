classdef DATA_sabcond < handle
    % DATA_sabcond class
    %   
    
    properties
        info = [];
        Xall = [];
        Bg = [];
        AB = [];
        Ym = [];
        ancillaries = [];
    end
    
    methods
        function [obj] = DATA_sabcond(TRRIFdata,varargin)
            [ obj.info ] = const_sabcond_info( TRRIFdata,varargin{:} );
            
        end
        function [Xall,infoA,infoAicelib] = constAbundance3D(obj)
            if isempty(obj.ancillaries)
                load(obj.info.fname_supple,'ancillaries');
                obj.ancillaries = ancillaries;
                clear ancillaries
            end
            Alibcachefname = constAlibcachefname(obj.info.libprefix_Alib,obj.info.basenameWA,...
                                      obj.info.optInterpid,obj.info.bands_opt,300);
            Alibcachefilepath = joinPath(obj.info.Alibdir,Alibcachefname);
            load(Alibcachefilepath,'Aall','valid_idx','logAallNrmedvalid',...
                 'logAcntrmvd','infoAall','setups');
             
            if ~isempty(obj.info.opticelib)
                Aicelibcachefname = constAicelibcachefname(obj.info.basenameWA,...
                               obj.info.opticelib,obj.info.bands_opt,300,obj.info.cntRmvl_ice);
                Aicelibcachefilepath = joinPath(obj.info.Alibdir,Aicelibcachefname);
                iceA_data = load(Aicelibcachefilepath,'Aicelib','infoAicelib','bands','valid_idx');
                Aicelib = iceA_data.Aicelib;
                infoAicelib = iceA_data.infoAicelib;
                valid_idx_ice = iceA_data.valid_idx;
            else
                Aicelib = [];
                infoAicelib = [];
                valid_idx_ice = [];
            end
             
            
            Xall = cat(3,obj.ancillaries.X);
            Xall = permute(Xall,[2 3 1]);
            infoA = infoAall(valid_idx);
            infoAicelib = infoAicelib(valid_idx_ice);
        end
        
        function [hsi] = readimg(obj)
            hsi = HSI(obj.info.basename,obj.info.save_dir);
            hsi.wa = obj.info.wa;
        end
        
        function [hsi] = readAB(obj)
            hsi = HSI([obj.info.basename '_AB'],obj.info.save_dir);
            hsi.wa = obj.info.wa;
        end
        
        function [hsi] = readBg(obj)
            hsi = HSI([obj.info.basename '_Bg'],obj.info.save_dir);
            hsi.wa = obj.info.wa;
        end
        
        function [hsi] = readCR(obj)
            hsi = envireadx([obj.info.save_dir obj.info.basename '_CR.img']);
            hsi.wa = obj.info.wa(obj.info.bands,:);
        end
        
        function [ancillaries] = load_ancillaries(obj)
            load(obj.info.fname_supple,'ancillaries');
            Alibcachefname = constAlibcachefname(obj.info.libprefix,obj.info.basenameWA,...
                                      obj.info.optInterpid,obj.info.bands_opt,200);
            Alibcachefilepath = joinPath(obj.info.Alibdir,Alibcachefname);
            load(Alibcachefilepath,'Aall','valid_idx','logAallNrmedvalid',...
                 'logAcntrmvd','infoAall');
            obj.info.infoAall = infoAall;
            valid_idxall = false(size(Aall,2),length(ancillaries));
            for c=1:length(ancillaries)
                if ~isempty(ancillaries(c).X)
                    Alibcachefname = constAlibcachefname(obj.info.libprefix,obj.info.basenameWA,...
                                      obj.info.optInterpid,obj.info.bands_opt,c);
                    Alibcachefilepath = joinPath(obj.info.Alibdir,Alibcachefname);
                    load(Alibcachefilepath,'Aall','valid_idx','logAallNrmedvalid',...
                                                    'logAcntrmvd','infoAall');
                    valid_idxall(:,c) = valid_idx';
                    if obj.info.cntRmvl
                        ancillaries(c).A = logAcntrmvd;
                    else
                        ancillaries(c).A = logAallNrmedvalid;
                    end
                end
            end
            obj.ancillaries = ancillaries;
            obj.info.valid_idxall = valid_idxall;
            
            
        end
    end
    
end