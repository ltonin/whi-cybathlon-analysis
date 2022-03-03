clearvars; clc;

includepat  = {'raceLog_', '2020'};
depthlevel  = 2;
datapath = '/mnt/data/Research/cybathlon/F1_mi_cybathlon/races/logs/';
%datapath = '/media/stefano/74A0406FA04039BE/cybathlon/F1_mi_cybathlon/races/logs/';
files  = whi_util_getfile(datapath, '.txt', 'include', includepat, 'level', depthlevel);
nfiles = length(files);
TotNumRaces = 0;

for fId = 1:nfiles
    cfilepath = files{fId};
    
    whi_util_bdisp(['[io] - Extracting log from file ' num2str(fId) '/' num2str(nfiles) ':']);
    disp(['     - ' cfilepath]);
    
    % Extracting log events
    fEVENT = whi_racelog_extract2(cfilepath);

    % Splitting per race
    Races    = unique(fEVENT.RAC);
    NumRaces = length(Races);
    TotNumRaces = TotNumRaces + NumRaces;
    
    % Getting current date and start timestamp (from filename)
    currDateStr = regexp(cfilepath, '.*\/(\d+)\/raceLog_.*', 'tokens');
    currDateStr = char(currDateStr{:});
    currDate = str2double(currDateStr);
    currTimeFileStr = regexp(cfilepath, 'raceLog_\d+-(\d*-\d*-\d*)', 'tokens');
    currTimeFileStr = char(currTimeFileStr{:});
    currdt = datetime([currDateStr ' ' currTimeFileStr], 'InputFormat', 'yyyyMMdd HH-mm-ss');

    
    % Saving file for each race
    for rId = 1:NumRaces

        cindex = fEVENT.RAC == Races(rId);

        EVENT.TYP = fEVENT.TYP(cindex);
        EVENT.DUR = fEVENT.DUR(cindex);
        EVENT.PLY = fEVENT.PLY(cindex);
        EVENT.RAC = fEVENT.RAC(cindex);

        T = fEVENT.T(cindex);

        startTime = T(1);
        startdt = currdt + seconds(startTime);
       
        startdt_str = datestr(startdt, 'yyyymmdd.HHMMss');

        T = T - T(1);
        nfilename = [datapath '/extracted/' startdt_str '.cybathlon.race.mat'];

        disp(['[out] - Saving race ' num2str(Races(rId)) ' in file: ' nfilename]);
        save(nfilename, 'T', 'EVENT');
    end
end

