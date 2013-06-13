% Plot Advection Tests using plot_2dadv.m
% By: Devin Light
% ------

%clear all;
%close all;
clc;

tests = {'adv_sine', ... % Uniform adv of sine^4
         'def_cosinebell', ... % LeVeque deformation test cosinebell
         'def_smooth_cosinebell', ... % Smoother version of LeVeque test
         'fdef_sqwave', ... % LeVeque deformation test square wave
         'hadv_cosinebell', ... % Horizontal advection of steep cosinebell
         };
res = {'1','2','3','4'};

which_test = tests(2);
which_res = res(2);


ncfilename = strcat('weno2d_' ,which_test{1}, '.nc');
nmax = 8;
ushoot = zeros(nmax,1);
oshoot = zeros(nmax,1);
%% Read in data, animate gifs
stat = -1;

for n=1:nmax
    if n==1
        methname = 'PPM, No Limiting';
        nc = ['pfctnon/' ncfilename];
        file = ['figures/ppm/ppm_' which_test{1}];
        pltn = 1;
        [Q1,x1,y1,t1] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q1(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==2
        methname = 'PPMDG, No Limiting';
        nc = ['ppmdghy/', ncfilename];
        file = ['figures/ppmdg/ppmdg_' which_test{1}];
        pltn = 4;
        [Q2,x2,y2,t2] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q2(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==3
        methname = 'DG, nolimiting';
        nc = ['dgnolim/', ncfilename];
        file = ['figures/dg/dg_' which_test{1}];
        pltn = 4;
        [Q3,x3,y3,t3] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q3(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==4
        methname = 'PPMDG, FCT, Positive';
        nc = ['ppmdgfc/', ncfilename];
        file = ['figures/ppmdgfct/ppmdgfct_' which_test{1}];
        pltn = 4;
        [Q4,x4,y4,t4] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q4(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==5
        methname = 'PPM, FCT, Positive';
        nc = ['pfctpos/' ncfilename];
        file = ['figures/ppmfct/ppmfct_' which_test{1}];
        pltn = 1;
        [Q5,x5,y5,t5] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q5(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==6
        methname = 'PPMDG, PMOD, Positive';
        nc = ['ppmdgpm/' ncfilename];
        file = ['figures/ppmdgpm/ppmdgpm_' which_test{1}];
        [Q6,x6,y6,t6] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q6(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==7
        methname = 'PPM, PMOD, Positive';
        nc = ['ppmdpos/' ncfilename];
        file = ['figures/ppmd/ppmd_' which_test{1}];
        [Q7,x7,y7,t7] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q7(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    elseif n==8
        methname = 'PPM, FCT, Sel, Pos';
        nc = ['pfctpse/' ncfilename];
        file = ['figures/ppmselp/ppmselp_' which_test{1}];
        [Q8,x8,y8,t8] = plot_2dadv(methname,nc,which_res,file,stat);
        tmp = squeeze(Q8(end,:,:));
        ushoot(n) = abs(min(tmp(:)));
        oshoot(n) = abs(max(tmp(:)));
    end
   
    
end

%% Make combined animation
% ---

%{
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2])
subplot(1,2,1)
%caxis manual; % allow subsequent plots to use the same color limits
%caxis([-0.1 1]); % set the color axis scaling to your min and max color limits

subplot(1,2,2)
set(gca, 'nextplot', 'replacechildren'); 
%caxis manual; % allow subsequent plots to use the same color limits
%caxis([-0.1 1]); % set the color axis scaling to your min and max color limits
pos = get(gcf, 'Position');
width = pos(3); height = pos(4);
mov = zeros(height, width, 1, length(t1), 'uint8');

for id = 1:length(t1)
        subplot(1,2,1)
        tmp1 = squeeze(Q1(id,:,:));
        contourf(x1,y1,tmp1)
        colorbar('location','EastOutside')
        ftitle = sprintf('PPM, No Limit - Time: %0.2f sec',t1(id));
        title(ftitle);
        
        subplot(1,2,2)
        tmp2 = squeeze(Q2(id,:,:));
        contourf(x2,y2,tmp2)
        colorbar('location','EastOutside')
        ftitle = sprintf('PPMDG - Time: %0.2f sec',t2(id));
        title(ftitle);

        % Get frame as an image
        f = getframe(gcf);

        % Create a colormap for the first frame. For the rest of the frames,
        % use the same colormap
        if id == 1
            [mov(:,:,1,id), map] = rgb2ind(f.cdata, 256, 'nodither');
        else
            mov(:,:,1,id) = rgb2ind(f.cdata, map, 'nodither');
        end
end

% Create animated GIF
imwrite(mov, map, 'figures/combo.gif', 'DelayTime', 0, 'LoopCount', inf);
%}

%% Make 6-panel figure
% ---
height = 717;
width = 630;
fig = figure('Position',[1 height width height]);
h = subplot(3,2,1);
nstep = size(t5,1);
tmp = squeeze(Q5(nstep/2,:,:));
contourf(x5,y5,tmp);
ftitle = {strcat('PPM FV, Pos limit, ',sprintf('Time: %0.2f sec',t5(nstep/2)));['nx=' num2str(size(tmp,1)) ' ny=' num2str(size(tmp,2))]};
title(ftitle);
axis square;
caxis([-1 1])

h = subplot(3,2,2);
tmp = squeeze(Q5(end,:,:));
contourf(x5,y5,tmp);
colorbar('location','EastOutside');
% caxis([0 1]);
ftitle = {strcat('PPM FV, Pos limit, ',sprintf('Time: %0.2f sec',t5(end)));['nx=' num2str(size(tmp,1)) ' ny=' num2str(size(tmp,2))]};
title(ftitle);
axis square;
caxis([-1 1])

h = subplot(3,2,3);
nstep = size(t3,1);
tmp = squeeze(Q3(nstep/2,:,:));
contourf(x3,y3,tmp);
ftitle = {strcat('DG, Unlimit, ',sprintf('Time: %0.2f sec',t3(nstep/2)));['x nodes=' num2str(size(tmp,1)) ' y nodes=' num2str(size(tmp,2))]};
title(ftitle);
axis square;
caxis([-1 1])

h = subplot(3,2,4);
tmp = squeeze(Q3(end,:,:));
contourf(x3,y3,tmp);
colorbar('location','EastOutside');
caxis([-1 1]);
ftitle = {strcat('DG, Unlimit, ',sprintf('Time: %0.2f sec',t3(end)));['x nodes=' num2str(size(tmp,1)) ' y nodes=' num2str(size(tmp,2))]};
title(ftitle);
axis square;
caxis([-1 1])


h = subplot(3,2,5);
nstep = size(t4,1);
tmp = squeeze(Q4(nstep/2,:,:));
contourf(x4,y4,tmp);
ftitle = {strcat('PPM/DG, Pos limit, , ',sprintf('Time: %0.2f sec',t4(nstep/2)));['x nodes=' num2str(size(tmp,1)) ' ny=' num2str(size(tmp,2))]};
title(ftitle);
axis square;
caxis([-1 1])

h = subplot(3,2,6);
tmp = squeeze(Q4(end,:,:));
contourf(x4,y4,tmp);
colorbar('location','EastOutside');
ftitle = {strcat('PPM/DG, Pos limit, , ',sprintf('Time: %0.2f sec',t4(end)));['x nodes=' num2str(size(tmp,1)) ' ny=' num2str(size(tmp,2))]};
title(ftitle);
axis square;
caxis([-1 1]);

print(fig,'-dpsc2','./figures/test_file')

%}

%% Compute largest magnitude undershoot over all time steps
tmpmin = 0;
for tt = 1:length(t4)
    tmp = squeeze(Q4(tt,:,:));
    min1 = min(tmp(:));
    tmpmin = min(min1,tmpmin);
end
Q40 = squeeze(Q4(1,:,:));
maxunder = abs(tmpmin)-abs(min(Q40(:)))


