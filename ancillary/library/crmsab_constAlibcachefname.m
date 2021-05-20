function Alibcachefname = constAlibcachefname(libprefix,wabasename,optInterpid,bands_opt,c)

Alibcachefname = sprintf('Alib%s_%s_interp%d_b%d_c%03d.mat',...
                             libprefix,wabasename,optInterpid,bands_opt,c);

end