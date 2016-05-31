clear;

% Load SAXS data
filename = 'exp-data/PEG_film.csv';
% filename = 'exp-data/PEG_doped.csv';
% filename = 'exp-data/PEG_undoped.csv';
if strcmp(filename,'exp-data/PEG_doped.csv')
    n = 48;  % total number of conditions
elseif strcmp(filename,'exp-data/PEG_undoped.csv')
    n = 9;
elseif strcmp(filename,'exp-data/PEG_film.csv')
    n = 5;
end

% read header
headerstr = strcat(repmat('%[^,],,,,',1,n-1),'%[^,\r\n]');
fid      = fopen(filename,'r');
header = textscan(fid, headerstr, 1);
fclose(fid);

% convert header to array
head = cell(1,n);
for ii = 1:n
    head{ii} = cell2mat(header{ii});
end

% load scattering wavevector and intensity
s = csvread(filename,1);

NAME = 'ODPA-AP6F-ODPA-PEG';
MW = 1500;
WPV = [10,20,30,40,50];
TV = 40;

figure;hold
leg = cell(1,length(WPV));
for ii = 1:length(WPV)
    WP = WPV(ii);
    col = (ii-1)/(length(WPV)-1);
    
    if strcmp(filename,'exp-data/PEG_doped.csv')
        plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%) DMAc EAN',sprintf(' %d',T),'C');
    elseif strcmp(filename,'exp-data/PEG_undoped.csv')
        plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%) DMAc');
    elseif strcmp(filename,'exp-data/PEG_film.csv')
        plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%)');
    end
    plotsaxs(s,head,plotname,[col 0 1-col]);
    leg{ii} = strcat(sprintf('WP = %d',WP),'%');
end

set(gca,'xscale','log');set(gca,'yscale','log');
xlim([.1,25]);ylim([1e3,1e6])
legend(leg)