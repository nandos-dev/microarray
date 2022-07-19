function gsedata=bmes_downloadandparsegse_thinh_fernando(gseid)
% by AhmetSacan.

% * Download e.g.,
% https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5847/matrix/GSE5847_series_matrix.txt.gz
% * gunzip()
% * geoseriesread()

%TODO: This function takes time! You should spend some effort to "cache"
%its results.

if ~exist('gseid','var'); gseid='GSE5847'; end

gseid_temp = floor(str2double(gseid(4:end))/1000);

url = sprintf('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%dnnn/%s/matrix/%s_series_matrix.txt.gz', ...
	gseid_temp,gseid,gseid);
gzfile = [tempdir '/' sprintf('%s.txt.gz',gseid)];
fprintf('Downloading %s ...\n',url);
websave(gzfile, url);
file = gunzip( gzfile );
file = file{1};

fprintf('Reading %s ...\n',file);
gsedata = geoseriesread( file );
