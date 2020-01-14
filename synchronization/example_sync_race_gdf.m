clearvars; clc;

subject = 'F1';
depthlevel  = 2;

%% General gdf and race paths
folder      = '2019_cybathlon_whiteam';
experiment  = 'mi_cybathlon';
gdfpath     = ['/home/ltonin/Desktop/' folder '/' subject '_' experiment '/'];

racpath     = '/home/ltonin/Desktop/2019_cybathlon_whiteam/F1_mi_cybathlon/races/logs/';
racinclude  = {'race'};

%% Get Race files
RaceFiles = whi_util_getfile(racpath, '.mat', 'include', racinclude, 'level', depthlevel);
NumRaces  = length(RaceFiles);

%% Find corresponding gdf
nfound = 0;
for rId = 1:NumRaces
    cfile = RaceFiles{rId};
    [cpath, cname, cext] = fileparts(cfile);
    [~, cfolder] = fileparts(cpath);
    whi_util_bdisp(['[io] - Loading ' num2str(rId) '/' num2str(NumRaces) ' from ''' cfolder '/'':']);
    disp(['     - Race file: ' [cname cext]]);
    
    % Find associated gdf and corresponding delay
    [foundgdf, founddelay] = whi_sync_findgdf(cfile, gdfpath);
    
    if(isnan(foundgdf))
        continue;
    end
    
    disp(['     - GDF file found: ' foundgdf]);
    disp(['     - Delay: ' num2str(founddelay) ' s']);
    nfound = nfound + 1;

    % Plotting [for debugging]
%     close all;
%     whi_sync_plot(foundgdf, cfile, founddelay);
%     whi_fig_set_position(gcf, 'Top');

    % Synchronize events
    mevent = whi_sync_merge(foundgdf, cfile, founddelay);
end



