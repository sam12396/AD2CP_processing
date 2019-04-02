%% Plotting of t,s,d and velocity data calculated from AD2CP
% This script will make two plots. (1) Comparing t,s,d over the duration of
% the glider's mission. (2) Comparing Shear Magnitude, N2, and Richardson
% number over the same duration. The shear magnitude and Richardson number
% are bin x segment, N2 is bin x profile, and tsd is averaged into the depth
% bins. I have already made most of these plots individually in different 
% programs but putting it all together here will help with formatting.

% Sam Coakley
% 3/6/19

%% Notes:
% 3/10/19: Going to add a prints for before and after Jan 31 to highlight an event
% in ru31 data
% 
% 3/20/19: Added water column depth in black on a seperate axis

%% Starters
clear all;close all;
% Define the deployment of interest
deployment='ru31_467_2016_01_13';
fig_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Figures/' deployment '/'];
data_path=['/Users/samcoa/Documents/MATLAB/ADCP_Wave/Data/' deployment '/'];

shear_data=load([data_path deployment '_ocean_shear.mat']);
g_data=load([data_path deployment '_dgroup.mat']);
ri_data=load([data_path deployment '_rich_num_grid.mat']);
vel_data=load([data_path deployment '_ocean_velo_grid.mat']);

%% Grab the data from the structures
% Shear
s_maggrid=shear_data.s_maggrid;
usgrid=shear_data.usgrid;
vsgrid=shear_data.vsgrid;
grid_bin=shear_data.grid_bin;
dz=grid_bin(2)-grid_bin(1);
% Create plotting bin var that uses bin edges
edge_grid_bin=grid_bin-dz./2; 

% Rich num
N2_grid=ri_data.N2_grid1m;
N2_time=ri_data.profile_time;
ri_grid=ri_data.ri_grid;

% Velocity
ugrid=vel_data.ugrid;
vgrid=vel_data.vgrid;

%% Load glider data
% Put the data to a grid with the same depth interval as the AD2CP data
[g_time, g_bin, temp]=g_data.dgroup.toGrid2d('drv_sea_water_temperature');
[~, ~, sal]          =g_data.dgroup.toGrid2d('drv_sea_water_salinity');
[~, ~, dens]         =g_data.dgroup.toGrid2d('drv_sea_water_density');
segs=g_data.dgroup.startDatenums;
sege=g_data.dgroup.endDatenums;
%% Calculate backscatter
beta=toArray(g_data.dgroup,'sensors','sci_flbbcd_bb_units');
bb_sal=toArray(g_data.dgroup,'sensors','drv_sea_water_salinity');
lambda=toArray(g_data.dgroup,'sensors','sci_flbbcd_bb_ref');
lambda=lambda(find(~isnan(lambda),1),3);
bb = beta2bb(beta, real(bb_sal), lambda);

% Put backscatter to the glider data grid
bb_grid=nan(length(g_bin),length(segs));
for kk=1:length(segs)
    seg_ind= bb(:,1)>=segs(kk) & bb(:,1)<sege(kk);
    seg_bb=bb(prof_ind,:);
    
    for jk=1:length(g_bin)-1
        ds_ind= seg_bb(:,2)>=g_bin(jk) & seg_bb(:,2)<g_bin(jk+1);
        if ~isempty(ds_ind)
            bb_grid(jk,kk)=nanmean(seg_bb(ds_ind,3));
        end
    end
end
clear bb bb_sal beta

% Glider data and AD2CP data do not need to be on the same grid b/c the
% glider data has better resolution so why not use it
g_dz=g_bin(2)-g_bin(1);

%% Get data for plotting depth of the water column
% Read in the water depth data from the glider
[col_depth]=toArray(g_data.dgroup,'sensors',{'m_water_depth'});
% Sometimes the m_water_depth records negative values for water column
% depth so we throw those out
ind= col_depth(:,3) > 1;
col_depth=col_depth(ind,:);

% Calculate the segment mean of the depth of the water column and put it on
% a grid
seg_mean_depth=nan(length(edge_grid_bin),length(segs));
for ii=1:length(segs)
    % Find all water depth measurements in the segment and take the mean
    seg_ind= col_depth(:,1)>= segs(ii) & col_depth(:,1)<sege(ii);
    seg_mean=nanmean(col_depth(seg_ind,3));
    
    [min_dif, dep_ind]=nanmin(abs(seg_mean-edge_grid_bin));
    seg_mean_depth(dep_ind:end,ii)=1;
end
% Add a few more rows on so that the depth extends to the max(g_bin) the max depth recorded by the glider 
f_grid_bin=[edge_grid_bin; (edge_grid_bin(end)+dz:dz:max(g_bin))'];
seg_mean_depth=[seg_mean_depth; ones(length((edge_grid_bin(end)+dz:dz:max(g_bin))),length(segs))];

% Fill a grid with water column depth data by finding nearest m_water_depth to
% scientific measurements
g_maxdep=nan(length(g_bin),length(g_time));
for ii=1:length(g_time)
    % Find the minimum difference between any timestamp of col_depth and
    % this stamp of g_time
    [min_dif, t_ind]=nanmin(abs(g_time(ii)-col_depth(:,1)));
    if min_dif*24*3600<60 % If the difference in time stamps is less than a minute 
        % Find the column depth at that time
        water_dep=col_depth(t_ind,3);
        
        % Find what bins to fill with bottom
        d_ind= ceil(water_dep)==g_bin;
        bin=find(d_ind==1,1);
        
        if sum(d_ind)~=0
            % Fill that bin and all bins below it with values that won't
            % interfere with the data
            g_maxdep(bin:end,ii)=1;
        end
    elseif ii~=1
        % If there is no good depth reading for the g_time stamp then take
        % the one before it and hope it doesn't look awful
        g_maxdep(:,ii)=g_maxdep(:,ii-1);
    end  
end

%% Set some plotting parameters
siz_title = 15;
siz_text = 12;
siz_marker= 8;
plot_width=8;
plot_height=1.5;
Ylim=[-floor(max(g_bin)) 0];
Xlim=[segs(1) sege(end)];

% We want to view the variables before and after the pivot time
pivot_time=[2016, 1 ,31];
pivot_datenum=datenum(pivot_time);
pivot_str=datestr(pivot_datenum,'yyyy_mm_dd');
start_str=datestr(segs(1),'yyyy_mm_dd');
end_str  =datestr(sege(end),'yyyy_mm_dd');

%% Plotting glider t,s,d
std_threshold=1;
% Make colorbar limits for each variable
tm=nanmean(temp(:));
ts=nanstd(temp(:));
tclim=[(tm-std_threshold*ts) (tm+std_threshold*ts)];
sm=nanmean(sal(:));
ss=nanstd(sal(:));
sclim=[(sm-std_threshold*ss) (sm+std_threshold*ss)];
dm=nanmean(dens(:));
ds=nanstd(dens(:));
dclim=[(dm-std_threshold*ds) (dm+std_threshold*ds)];

% Define colormaps for each variable
tcmap=cmocean('thermal');
scmap=cmocean('haline');
dcmap=cmocean('dense');

% Plot temperature
figure(1)
tax=axes;
% Plot the variable on one axis
pcolorjw(g_time,-g_bin,temp)
% Plot water column depth on a second axes
dep_ax=axes;
pcolorjw(g_time,-g_bin,-g_maxdep);
dep_ax.Colormap=[0,0,0];
% Merge the axes
linkaxes([dep_ax,tax]);
% Hide the depth axes
dep_ax.Visible='off';
% Turn off ticks for the depth axis
dep_ax.XTick=[];
dep_ax.YTick=[];
% Align the axes in the frame
set([tax,dep_ax],'Position',[.08 .2 .78 .63]);
% Make a colorbar and align it within the frame
tcb = colorbar(tax,'Position',[0.87 .2 .03 .63]);
% Turn on shading and the grid in the variable axis
shading(tax,'interp'); grid(tax,'on');
% Set texts sizes
set(tax,'fontsize',siz_text)
% Set colorbar limits and the colormap
tax.CLim=tclim;
tax.Colormap=[tcmap];
% Turn in datetick in the variable axis
datetick(tax,'x',6)
% Add a label on the y and set limits for the x and y in the variable axis
ylabel(tax,'Depth [m]'); ylim(tax,Ylim);
xlim(tax,Xlim);
% Set a color bar label in the variable axis
tcb_label=get(tcb,'Label');
set(tcb_label,'String','°C');
% Make a title in the variable axis
title(tax,[deployment ' Temperature ' num2str(g_dz) 'm Resolution'],'Interpreter','none','fontsize',15)
% Set the printing position and dimensions
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
% Show what the printed figure will look like
wysiwyg

% Print the figure
print([fig_path 'full_temp_' num2str(g_dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    % Take prints of the variable around the pivot time
    xlim([segs(1) pivot_datenum]) %Why does this make my images taller???
    print([fig_path start_str '_' pivot_str '_temp_' num2str(g_dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_temp_' num2str(g_dz) 'm_cs'],'-dpng');
end
close

% Plot salinity
% subplot(3,2,3)
figure(2)
sax=axes;
pcolorjw(g_time,-g_bin,sal)
% Plot water column depth
dep_ax=axes;
pcolorjw(g_time,-g_bin,-g_maxdep);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,sax]);
dep_ax.Visible='off'; % Hide the top axes
dep_ax.XTick=[];
dep_ax.YTick=[];
set([sax,dep_ax],'Position',[.08 .2 .78 .63]);
scb = colorbar(sax,'Position',[0.87 .2 .03 .63]);
shading(sax,'interp'); grid(sax,'on');
set(sax,'fontsize',siz_text)
sax.CLim=sclim;
sax.Colormap=[scmap];
datetick(sax,'x',6)
ylabel(sax,'Depth [m]'); ylim(sax,Ylim);
xlim(sax,Xlim);
scb_label=get(scb,'Label');
set(scb_label,'String','PSU');
title(sax,[deployment ' Salinity ' num2str(g_dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_sal_' num2str(g_dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_sal_' num2str(g_dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_sal_' num2str(g_dz) 'm_cs'],'-dpng');
end
close

% Plot density
%subplot(3,2,5)
figure(3)
dax=axes;
pcolorjw(g_time,-g_bin,dens)
% Plot water column depth
dep_ax=axes;
pcolorjw(g_time,-g_bin,-g_maxdep);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,dax]);
dep_ax.Visible='off'; % Hide the top axes
dep_ax.XTick=[];
dep_ax.YTick=[];
set([dax,dep_ax],'Position',[.08 .2 .78 .63]);
dcb = colorbar(dax,'Position',[0.87 .2 .03 .63]);
shading(dax,'interp'); grid(dax,'on');
set(dax,'fontsize',siz_text)
dax.CLim=dclim;
dax.Colormap=[dcmap];
datetick(dax,'x',6)
ylabel(dax,'Depth [m]'); ylim(dax,Ylim);
xlim(dax,Xlim);
dcb_label=get(dcb,'Label');
set(dcb_label,'String','kg m^-^3');
title(dax,[deployment ' Density ' num2str(g_dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_dens_' num2str(g_dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_dens_' num2str(g_dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_dens_' num2str(g_dz) 'm_cs'],'-dpng');
end
close
%% Plotting shear, N2, and rich num
% Shear magnitude and rich num are on the same grid but N2 is on a bin x
% profile grid.
% Make colorbar limits for each variable
n2m=nanmean(N2_grid(:));
n2s=nanstd(N2_grid(:));
n2clim=[0 (n2m+std_threshold*n2s)];
smm=nanmean(s_maggrid(:));
sms=nanstd(s_maggrid(:));
smclim=[0 (smm+std_threshold*sms)];
% Rm=nanmean(ri_grid(:));
% Rs=nanstd(ri_grid(:));
% Rclim=[0 (Rm+std_threshold*Rs)];

% Define colormaps for each variable
n2cmap=cmocean('speed');
smcmap=cmocean('speed');
Rcmap=cmocean('balance');


% Plot N2
figure(4)
n2ax=axes;
pcolorjw(N2_time,-g_bin(1:size(N2_grid,1)),N2_grid)
% Plot water column depth
dep_ax=axes;
pcolorjw(g_time,-g_bin,-g_maxdep);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,n2ax]);
dep_ax.Visible='off'; % Hide the top axes
dep_ax.XTick=[];
dep_ax.YTick=[];
set([n2ax,dep_ax],'Position',[.08 .2 .78 .63]);
n2cb = colorbar(n2ax,'Position',[0.87 .2 .03 .63]);
shading(n2ax,'interp'); grid(n2ax,'on');
set(n2ax,'fontsize',siz_text)
n2ax.CLim=n2clim;
n2ax.Colormap=[n2cmap];
datetick(n2ax,'x',6)
ylabel(n2ax,'Depth [m]'); ylim(n2ax,Ylim);
xlim(n2ax,Xlim);
n2cb_label=get(n2cb,'Label');
set(n2cb_label,'String','s^-^1');
title(n2ax,[deployment ' N2 ' num2str(g_dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_N2_' num2str(g_dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_N2_' num2str(g_dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_N2_' num2str(g_dz) 'm_cs'],'-dpng');
end
close

%subplot(3,2,4)
figure(5)
smax=axes;
pcolorjw(segs,-edge_grid_bin,s_maggrid)
dep_ax=axes;
pcolorjw(segs,-f_grid_bin,seg_mean_depth);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,smax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([smax,dep_ax],'Position',[.08 .2 .78 .63]);
smcb = colorbar(smax,'Position',[0.87 .2 .03 .63]);
shading(smax,'interp'); grid(smax,'on');
set(smax,'fontsize',siz_text)
smax.CLim=smclim;
smax.Colormap=[smcmap];
datetick(smax,'x',6)
ylabel(smax,'Depth [m]'); ylim(smax,Ylim);
xlim(smax,Xlim);
smcb_label=get(smcb,'Label');
set(smcb_label,'String','s^-^1');
title(smax,[deployment ' dU/dz ' num2str(dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_shear_mag_' num2str(dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_shear_mag_' num2str(dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_shear_mag_' num2str(dz) 'm_cs'],'-dpng');
end
close

% Plot rich num
% Need an even number of contours
contour_start=-3.5;
contour_end=3.5;
contour_step=1;
num_c=length(contour_start:contour_step:contour_end);

% Adjust the colormap. Need to make the amount of colors=the # of contour
% steps. Want transition from red to blue to be sharp and exactly at 0
csize=length(Rcmap);
% Want num_c/2 colors in red and num_c/2 in blue
cinds=floor(linspace(1,ceil(csize/4),num_c/2));
cinds=[cinds, floor(linspace(ceil(3*csize/4),csize,num_c/2))];
nmap=Rcmap(cinds,:);
% Make proper tick position and tick labels
ticks=round(linspace(contour_start,contour_end,num_c+1),2);

%subplot(3,2,6)
figure(6)
Rax=axes;
contourf(segs,-edge_grid_bin,real(log10(ri_grid)),[contour_start:contour_step:contour_end]);
dep_ax=axes;
pcolorjw(segs,-f_grid_bin,seg_mean_depth);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,Rax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([Rax,dep_ax],'Position',[.08 .2 .78 .63]);
Rcb = colorbar(Rax,'Position',[0.87 .2 .03 .63]);
shading(Rax,'interp'); grid(Rax,'on');
set(Rax,'fontsize',siz_text)
Rax.CLim=[contour_start contour_end];
Rax.Colormap=nmap;
datetick(Rax,'x',6)
ylabel(Rax,'Depth [m]'); ylim(Rax,Ylim);
xlim(Rax,Xlim);
Rcb_label=get(Rcb,'Label');
set(Rcb_label,'String','[no dim]');
title(Rax,[deployment ' log10(Ri) ' num2str(dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_rich_num_' num2str(dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_rich_num_' num2str(dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_rich_num_' num2str(dz) 'm_cs'],'-dpng');
end
close

%% Plotting component velocity and component shear
std_threshold=3;
% Make colorbar limits for each variable and center at 0
% Velocity
um=nanmean(ugrid(:));
us=nanstd(ugrid(:));
uclim=[-(um+std_threshold*us) (um+std_threshold*us)];
vm=nanmean(vgrid(:));
vs=nanstd(vgrid(:));
vclim=[-(vm+std_threshold*vs) (vm+std_threshold*vs)];
% Shear
usm=nanmean(usgrid(:));
uss=nanstd(usgrid(:));
usclim=[-(usm+std_threshold*uss) (usm+std_threshold*uss)];
vsm=nanmean(vsgrid(:));
vss=nanstd(vsgrid(:));
vsclim=[-(vsm+std_threshold*vss) (vsm+std_threshold*vss)];

% Define colormaps
ucmap=cmocean('delta');
vcmap=cmocean('delta');
uscmap=cmocean('delta');
vscmap=cmocean('delta');


% Plot U velocity
figure(7)
uax=axes;
pcolorjw(segs,-edge_grid_bin,ugrid)
dep_ax=axes;
pcolorjw(segs,-f_grid_bin,seg_mean_depth);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,uax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([uax,dep_ax],'Position',[.08 .2 .78 .63]);
ucb = colorbar(uax,'Position',[0.87 .2 .03 .63]);
shading(uax,'interp'); grid(uax,'on');
set(uax,'fontsize',siz_text)
uax.CLim=uclim;
uax.Colormap=[ucmap];
datetick(uax,'x',6)
ylabel(uax,'Depth [m]'); ylim(uax,Ylim);
xlim(uax,Xlim);
ucb_label=get(ucb,'Label');
set(ucb_label,'String','m s^-^1');
title(uax,[deployment ' East(+)/West(-) Velocity ' num2str(dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_u_vel_' num2str(dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_u_vel_' num2str(dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_u_vel_' num2str(dz) 'm_cs'],'-dpng');
end
close

% Plot V velocity
figure(8)
vax=axes;
pcolorjw(segs,-edge_grid_bin,vgrid)
dep_ax=axes;
pcolorjw(segs,-f_grid_bin,seg_mean_depth);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,vax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([vax,dep_ax],'Position',[.08 .2 .78 .63]);
vcb = colorbar(vax,'Position',[0.87 .2 .03 .63]);
shading(vax,'interp'); grid(vax,'on');
set(vax,'fontsize',siz_text)
vax.CLim=vclim;
vax.Colormap=[vcmap];
datetick(vax,'x',6)
ylabel(vax,'Depth [m]'); ylim(vax,Ylim);
xlim(vax,Xlim);
vcb_label=get(vcb,'Label');
set(vcb_label,'String','m s^-^1');
title(vax,[deployment ' North(+)/South(-) Velocity ' num2str(dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_v_vel_' num2str(dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_v_vel_' num2str(dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_v_vel_' num2str(dz) 'm_cs'],'-dpng');
end
close

% Plot U shear
figure(9)
usax=axes;
pcolorjw(segs,-edge_grid_bin,usgrid)
dep_ax=axes;
pcolorjw(segs,-f_grid_bin,seg_mean_depth);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,usax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([usax,dep_ax],'Position',[.08 .2 .78 .63]);
uscb = colorbar(usax,'Position',[0.87 .2 .03 .63]);
shading(usax,'interp'); grid(usax,'on');
set(usax,'fontsize',siz_text)
usax.CLim=usclim;
usax.Colormap=[uscmap];
datetick(usax,'x',6)
ylabel(usax,'Depth [m]'); ylim(usax,Ylim);
xlim(usax,Xlim);
uscb_label=get(uscb,'Label');
set(uscb_label,'String','s^-^1');
title(usax,[deployment ' East(+)/West(-) du/dz ' num2str(dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_u_shear_' num2str(dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_u_shear_' num2str(dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_u_shear_' num2str(dz) 'm_cs'],'-dpng');
end
close

% Plot V Shear
figure(10)
vsax=axes;
pcolorjw(segs,-edge_grid_bin,vsgrid)
dep_ax=axes;
pcolorjw(segs,-f_grid_bin,seg_mean_depth);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,vsax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([vsax,dep_ax],'Position',[.08 .2 .78 .63]);
vscb = colorbar(vsax,'Position',[0.87 .2 .03 .63]);
shading(vsax,'interp'); grid(vsax,'on');
set(vsax,'fontsize',siz_text)
vsax.CLim=vsclim;
vsax.Colormap=[vscmap];
datetick(vsax,'x',6)
ylabel(vsax,'Depth [m]'); ylim(vsax,Ylim);
xlim(vsax,Xlim);
vscb_label=get(vscb,'Label');
set(vscb_label,'String','s^-^1');
title(vsax,[deployment ' North(+)/South(-) dv/dz ' num2str(dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_v_shear_' num2str(dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_v_shear_' num2str(dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_v_shear_' num2str(dz) 'm_cs'],'-dpng');
end
close

%% Plot backscatter
std_threshold=2;
% Make colorbar limits for each variable
bbm=nanmean(bb_grid(:));
bbs=nanstd(bb_grid(:));
bbclim=[0 (bbm+std_threshold*bbs)];

% Define colormap
bbcmap=cmocean('turbid');

figure(11)
bbax=axes;
pcolorjw(segs,-g_bin,bb_grid)
dep_ax=axes;
pcolorjw(g_time,-g_bin,-g_maxdep);
dep_ax.Colormap=[0,0,0];
linkaxes([dep_ax,bbax]);
dep_ax.Visible='off';
dep_ax.XTick=[];
dep_ax.YTick=[];
set([bbax,dep_ax],'Position',[.08 .2 .78 .63]);
bbcb = colorbar(bbax,'Position',[0.87 .2 .03 .63]);
shading(bbax,'interp'); grid(bbax,'on');
set(bbax,'fontsize',siz_text)
bbax.CLim=bbclim;
bbax.Colormap=[bbcmap];
datetick(bbax,'x',6)
ylabel(bbax,'Depth [m]'); ylim(bbax,Ylim);
xlim(bbax,Xlim);
bbcb_label=get(bbcb,'Label');
set(bbcb_label,'String','no dim');
title(bbax,[deployment ' Total Backscattering Coeff ' num2str(g_dz) 'm Resolution'],'Interpreter','none','fontsize',siz_title)
set(gcf,'PaperPosition',[0 0 plot_width plot_height])
wysiwyg

print([fig_path 'full_bb_' num2str(g_dz) 'm_cs'],'-dpng');
if exist('pivot_time','var')
    xlim([segs(1) pivot_datenum])
    print([fig_path start_str '_' pivot_str '_bb_' num2str(g_dz) 'm_cs'],'-dpng');
    xlim([pivot_datenum sege(end)])
    print([fig_path pivot_str '_' end_str '_bb_' num2str(g_dz) 'm_cs'],'-dpng');
end
close