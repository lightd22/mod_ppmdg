% ----------
%  Error comparison between DG/DG and PPM/PPM Strang splits
%  By : Devin Light
% ----------

clear all;
close all;
clc;

res = {'1','2','3','4'};
meth = {'ppm','N4','N5','N6','N7'};
nmax = 5;
nres = 4;

errors = zeros(nres,1);
nx = zeros(nres,1);


for n=1:nmax
    if n==1
        which_meth = meth(n);
        ncfilename = strcat('weno2d_adv_sine_',which_meth{1},'.nc');
        for i=1:nres
            [fooe,foon] = get_error(ncfilename,res(i));
            errors(i) = fooe;
            nx(i) = foon;
        end
        semilogy(nx,errors,'--bd','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',8);
        hold on
    
    elseif n==2
        which_meth = meth(n);
        ncfilename = strcat('weno2d_adv_sine_',which_meth{1},'.nc');
        for i=1:nres
            [fooe,foon] = get_error(ncfilename,res(i));
            errors(i) = fooe;
            nx(i) = foon;
        end
        holdr = errors;
        semilogy(nx,errors,'--rd','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',8);
    elseif n==3
        which_meth = meth(n);
        ncfilename = strcat('weno2d_adv_sine_',which_meth{1},'.nc');
        for i=1:nres
            [fooe,foon] = get_error(ncfilename,res(i));
            errors(i) = fooe;
            nx(i) = foon;
        end
        semilogy(nx,errors,'--gs','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',8);
    elseif n==4
        which_meth = meth(n);
        ncfilename = strcat('weno2d_adv_sine_',which_meth{1},'.nc');
        for i=1:nres
            [fooe,foon] = get_error(ncfilename,res(i));
            errors(i) = fooe;
            nx(i) = foon;
        end
        semilogy(nx,errors,'--co','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','c',...
                'MarkerSize',8);
    elseif n==5
        which_meth = meth(n);
        ncfilename = strcat('weno2d_adv_sine_',which_meth{1},'.nc');
        for i=1:nres
            [fooe,foon] = get_error(ncfilename,res(i));
            errors(i) = fooe;
            nx(i) = foon;
        end
        semilogy(nx,errors,'--mx','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','m',...
                'MarkerSize',8);
                
    
    end
    
end
hold off
title('L2 errors assoc. w/ sine adv');
xlabel('nx')
ylabel('error')
legend('PPM', 'DGN4','DGN5','DGN6','DGN7')
