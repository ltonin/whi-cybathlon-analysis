clearvars; clc;

subject = 'F1';
depthlevel  = 2;

%% General gdf and race paths
rootpath    = '/mnt/data/Research/';
folder      = 'cybathlon';
experiment  = 'mi_cybathlon';
gdfpath     = [rootpath '/' folder '/' subject '_' experiment '/'];

racpath     = [rootpath '/' folder '/' subject '_' experiment '/races/logs/extracted/'];
racinclude  = {'race', '2019'};

savepath = 'analysis/events/';
util_mkdir(pwd, savepath);

%% Get Race files
RaceFiles = whi_util_getfile(racpath, '.mat', 'include', racinclude, 'level', depthlevel);
NumRaces  = length(RaceFiles);

%% Find corresponding gdf

pfoundgdf = '';
nfound    = 1;
NumEvtFiles = 0;
for rId = 1:NumRaces
    
    cfile = RaceFiles{rId};
    
    [cpath, cname, cext] = fileparts(cfile);
    [~, cfolder] = fileparts(cpath);
    whi_util_bdisp(['[io] - Loading ' num2str(rId) '/' num2str(NumRaces) ' from ''' cfolder '/'':']);
    disp(['     - Race file: ' [cname cext]]);
    
    % Find associated gdf and corresponding delay
   % [cfoundgdf, founddelay] = whi_sync_findgdf(cfile, gdfpath, 'GDFCommands', [1670 1672]); % 2020 Races
    [cfoundgdf, founddelay] = whi_sync_findgdf(cfile, gdfpath); % 2019 Races
    
    if(isnan(cfoundgdf))
        continue;
    end
    
    
    [evtfilename, nfound] = generate_event_filename(cfoundgdf, pfoundgdf, nfound);
    

    disp(['     - GDF file found: ' cfoundgdf]);
    disp(['     - Delay: ' num2str(founddelay) ' s']);
    disp(['     - Event File: ' evtfilename ]);
    NumEvtFiles = NumEvtFiles + 1;

    
    [~, gdfH] = sload(cfoundgdf);
    nsamples = gdfH.FILE.POS;
    
    % Synchronize events
    h.EVENT = whi_sync_merge(cfoundgdf, cfile, founddelay);
    
    pos = h.EVENT.RAC.POS(h.EVENT.RAC.TYP == 800);
    posId = find(h.EVENT.RAC.TYP == 800);
    
   
    if (pos + h.EVENT.RAC.DUR(posId) > nsamples) 
       keyboard;
    end
    
    pfoundgdf = cfoundgdf;
    
    % Save event MAT file
    save([savepath evtfilename], 'h');
end

function [filename, nItem] = generate_event_filename(cGDF, pGDF, nItem)

    if (strcmp(pGDF, cGDF))
        nItem = nItem + 1;
    else
        nItem = 1;
    end
    
    [~, filename] = fileparts(cGDF);
    filename     = [filename '.events.' num2str(nItem) '.mat'];

end


