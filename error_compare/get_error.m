% ---
% Gets data from files, given file name
% ---

function [e2, N] = get_error(ncfilename,res)

    Qname = strcat('Q',res{1});
    xname = strcat('x',res{1});
    yname = strcat('y',res{1});
    
    Q = nc_varget(ncfilename, Qname);
    x = nc_varget(ncfilename, xname);
    y = nc_varget(ncfilename, yname);
    t = nc_varget(ncfilename, 'time');

    tmpic = squeeze(Q(1,:,:));
    tmpfin = squeeze(Q(end,:,:));
        
    % Am I computing these norms correctly?
    nx = size(tmpic,1);
    ny = size(tmpic,2);
    
    e2 = sqrt( sum( (tmpfin(:)-tmpic(:)).^2 )/(nx*ny) );
    N = nx;
    
end