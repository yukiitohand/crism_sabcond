function [obs_dirname_list] = read_obs_dirname_list(fpath)
% [obs_dirname_list] = read_obs_dirname_list(fpath)
% read observation text files storing dir names 
%  INPUTS
%    fpath: path to the file
%  OUTPUTS
%    obs_dirname_list: cell array of the observation dirnames in the form
%    of
%       cccnnnnnn
%        ccc = class type of the obervation
%        nnnnnnn = observation id 

obs_dirname_list = [];

fid = fopen(fpath,'r');
while ~feof(fid)
    tline = fgetl(fid);
    obs_dirname_list = [obs_dirname_list {tline}];
end

fclose(fid);

end