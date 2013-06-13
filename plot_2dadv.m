% Animated GIF plotting function for operator splitting advection
% By: Devin Light
% ----

function [Q,x,y,t] = plot_2dadv(methname,ncfilename,res,file,stat)
    
    Qname = strcat('Q',res{1});
    xname = strcat('x',res{1});
    yname = strcat('y',res{1});
    
    Q = nc_varget(ncfilename, Qname);
    x = nc_varget(ncfilename, xname);
    y = nc_varget(ncfilename, yname);
    t = nc_varget(ncfilename, 'time');
    
    nnodes = 5;
    nt = size(Q,1);
    nx = size(Q,2);
    nelem = nx/nnodes;
    elemdx = 1/nelem;

    scrsz = get(0,'ScreenSize');

    % Make initial frame
    figure('Position',[1 scrsz(4)/2 3*scrsz(4)/4 3*scrsz(4)/4])
    
    % Get figure size
    pos = get(gcf, 'Position');
    width = pos(3); height = pos(4);
    
    %set(gca, 'nextplot', 'replacechildren'); 
    %caxis manual; % allow subsequent plots to use the same color limits
    %caxis([-1 1]); % set the color axis scaling to your min and max color limits

    % Preallocate data (for storing frame data)
    len = length(t);
    mov = zeros(height, width, 1, len, 'uint8');

    % Loop through by changing XData and YData
    if(stat == 1)
    for id = 1:length(t)
            
        tmp = squeeze(Q(id,:,:));
        contourf(x,y,tmp)
        colorbar('location','EastOutside')
        ftitle = {strcat(methname, sprintf(', Time: %0.2f sec',t(id)));['nx=' num2str(size(tmp,1)) ' ny=' num2str(size(tmp,2))]};
        title(ftitle);
        axis square;
        % Get frame as an image
        f = getframe(gcf);

        % Create a colormap for the first frame. For the rest of the frames,
        % use the same colormap
        if id == 1
            [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
            set(gca, 'nextplot', 'replacechildren'); 
            caxis manual; % allow subsequent plots to use the same color limits
            caxis([-1 1]); % set the color axis scaling to your min and max color limits
        else
            mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
        end
       
    end
    name = strcat(file, res{1},'.gif');
% Create animated GIF
imwrite(mov, map, name, 'DelayTime', 0.1, 'LoopCount', inf);
    else
        tmp = squeeze(Q(end,:,:));
        contourf(x,y,tmp);
        colorbar('location','EastOutside')
        ftitle = {strcat(methname, sprintf(', Time: %0.2f sec',t(end)));['nx=' num2str(size(tmp,1)) ' ny=' num2str(size(tmp,2))]};
        title(ftitle);
        axis square;
        caxis([-1 1])
    end
    

end