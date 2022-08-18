classdef SABCONDdataset < dynamicprops
    % SABCONDdataset
    %   Container for the output of the sabcond correction. 
    
    properties
        basename_cor
        dirpath
        prop_source
        source_trr3
        cor
        nr
        nr_ds
        nr_ds_rp
        ori
        AB
        Bg
        Ice
    end
    
    methods
        function obj = SABCONDdataset(basename,dirpath,varargin)
            % Constructor
            
            suffix = '';
            if (rem(length(varargin),2)==1)
                error('Optional parameters should always go by pairs');
            else
                for i=1:2:(length(varargin)-1)
                    switch upper(varargin{i})
                        % 
                        case 'SUFFIX'
                            suffix = varargin{i+1};
                        otherwise
                            error('Unrecognized option: %s',varargin{i});
                    end
                end
            end
            if isempty(suffix)
                basename_cor = basename;
            else
                basename_cor = [basename '_' suffix];
            end
            
            obj.basename_cor = basename_cor;
            obj.dirpath      = dirpath;
            
            % Stock different kinds of processed images
            % 
            obj.appendCAT('cor',basename_cor,dirpath,'');
            obj.appendCAT('nr',basename_cor,dirpath,'nr');
            obj.appendCAT('ori',basename_cor,dirpath,'ori');
            obj.appendCAT('AB',basename_cor,dirpath,'AB');
            obj.appendCAT('Bg',basename_cor,dirpath,'Bg');
            obj.appendCAT('mdl',basename_cor,dirpath,'mdl');
            obj.appendCAT('Ice',basename_cor,dirpath,'Ice');
            obj.appendCAT('nr_ds',basename_cor,dirpath,'nr_ds');
            obj.appendCAT('mdl_ds',basename_cor,dirpath,'mdl_ds');
            obj.appendCAT('nr_ds_rp',basename_cor,dirpath,'nr_ds_rp');
            
            % Decode "basename" to get information about the source image. 
            obj.prop_source = crism_getProp_basenameOBSERVATION(basename);
            prop_source_trr3 = obj.prop_source;
            prop_source_trr3.version = 3;
            basename_source_trr3 = crism_get_basenameOBS_fromProp(prop_source_trr3);
            
            TRR3IFdata = CRISMdata(basename_source_trr3,'');
            if isempty(TRR3IFdata.lbl)
                TRR3IFdata.download(2);
                TRR3IFdata = CRISMdata(basename_source_trr3,'');
            end
            TRR3IFdata.readWAi();
            
            if ~isempty(obj.cor)
                obj.cor.wa = TRR3IFdata.wa;
            end
            if ~isempty(obj.nr)
                obj.nr.wa  = TRR3IFdata.wa;
            end
            if ~isempty(obj.ori)
                obj.ori.wa = TRR3IFdata.wa;
            end
            if ~isempty(obj.AB)
                obj.AB.wa  = TRR3IFdata.wa;
            end
            if ~isempty(obj.Bg)
                obj.Bg.wa  = TRR3IFdata.wa;
            end
            if ~isempty(obj.Ice)
                obj.Ice.wa = TRR3IFdata.wa;
            end
            
            obj.source_trr3 = TRR3IFdata;

        end
        
        function appendCAT(obj,propName,bname,dirpath,suffix_a)
            if isempty(suffix_a)
                basename_a = bname;
            else
                basename_a = [bname '_' suffix_a];
            end
            if exist(joinPath(dirpath,[basename_a '.img']),'file') && exist(joinPath(dirpath,[basename_a '.hdr']),'file')
                if ~isprop(obj,propName)
                    addprop(obj,propName);
                end
                obj.(propName) = CRISMdataCAT(basename_a,dirpath);
            end
            
        end
        
    end
end

