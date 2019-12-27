clearvars; clc;

includepat  = {'raceLog_'};
depthlevel  = 2;
datapath = '/home/ltonin/Desktop/2019_cybathlon_whiteam/F1_mi_cybathlon/races/logs/';
files  = whi_util_getfile(datapath, '.txt', 'include', includepat, 'level', depthlevel);
nfiles = length(files);
TotNumRaces = 0;

for fId = 1:nfiles
    cfilepath = files{fId};
    util_bdisp(['[io] - Extracting log from file ' num2str(fId) '/' num2str(nfiles) ':']);
    disp(['     - ' cfilepath]);
    
    % Extracting log events
    [fT, fE, fEVENT] = whi_racelog_extract(cfilepath);

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

        cPOS = fEVENT.POS(cindex);
        cTYP = fEVENT.TYP(cindex);
        cPLY = fEVENT.PLY(cindex);

        T = fT(cPOS(1):cPOS(end));
        E = fE(cPOS(1):cPOS(end));

        EVENT.POS = cPOS;
        EVENT.TYP = cTYP;
        EVENT.PLY = cPLY;

        startTime = T(1);
        startdt = currdt + seconds(startTime);

        startdt_str = datestr(startdt, 'yyyymmdd.HHMMss');

        nfilename = [fileparts(cfilepath) '/' startdt_str '.cybathlon.race.mat'];

        disp(['[out] - Saving race ' num2str(Races(rId)) ' in file: ' nfilename]);
        save(nfilename, 'T', 'E', 'EVENT');
    end
end

