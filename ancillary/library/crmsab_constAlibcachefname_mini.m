function Alibcachefname = crmsab_constAlibcachefname_mini(libprefix,wabasename,optInterpid,bands_opt,c)

Alibcachefname = sprintf('Alib%s_%s_interp%d_b%d_mini_c%03d.mat',...
                             libprefix,wabasename,optInterpid,bands_opt,c);

end