function [masterbase] = const_masterbase(libname,opt,wabasename,method,retainRatio)
switch method
    case 'interp1'
        masterbase = sprintf('%s%d_%s_%s',libname,opt,wabasename,method);
    case 'interpCRISMspc'
        masterbase = sprintf('%s%d_%s_r%d',libname,opt,wabasename,retainRatio*10);
    otherwise
        error('method %s is not valid',method);
end
end