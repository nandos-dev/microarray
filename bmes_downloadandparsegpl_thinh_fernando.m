function gpldata=bmes_downloadandparsegpl_thinh_fernando(gplid)
% * Download e.g.,
% http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?form=text&acc=GPL96&view=full
% * geosoftread()

%TODO: This function takes time! You should spend some effort to "cache"
%its results.

if ~exist('gplid','var'); gplid='GPL96'; end

url = sprintf('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?form=text&acc=%s&view=full',gplid);
file = [tempdir '/' gplid '.txt'];
websave(file, url);
gpldata = geosoftread( file );
