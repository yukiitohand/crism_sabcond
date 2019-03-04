function [suffix] = const_suffix_v2(mt,cntRmvl,libprefix,optInterpid,bands_opt,nIter,additional_suffix)

suffix = sprintf('_atcr_%s_%d_%s_%d_%d_%d',...
               mt,cntRmvl,libprefix,optInterpid,bands_opt,nIter);

if ~isempty(additional_suffix)
    suffix = [suffix '_' additional_suffix];
end

end