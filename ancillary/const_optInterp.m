function optInterp = const_optInterp(optInterpid)

switch optInterpid
    case 1
        optInterp(1).method = 'interpCRISMspc';
        optInterp(1).retainRatio = 0.1;
        optInterp(2).method = 'interpCRISMspc';
        optInterp(2).retainRatio = 0.1;
        optInterp(3).method = 'interp1';
        optInterp(3).retainRatio = [];
        optInterp(4).method = 'interp1';
        optInterp(4).retainRatio = [];
    otherwise
        error('optIntepid %d is not defined',optInterpid);
end

end