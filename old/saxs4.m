clear;

% Load SAXS data
filename = 'exp-data/PI.csv';
% filename = 'exp-data/PEG_film.csv';
% filename = 'exp-data/PEG_doped.csv';
% filename = 'exp-data/PEG_undoped.csv';
if strcmp(filename,'exp-data/PEG_doped.csv')
    n = 48;  % total number of conditions
elseif strcmp(filename,'exp-data/PEG_undoped.csv')
    n = 9;
elseif strcmp(filename,'exp-data/PEG_film.csv')
    n = 5;
elseif strcmp(filename,'exp-data/PI.csv')
    n = 4;
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

% NAME = 'ODPA-AP6F-ODPA-PEG';
NAME = 'ODPA-AP6F';
MW = 1500;
WP = 50;
TV = [40,50,60,70,80];

figure;hold

if strcmp(filename,'exp-data/PEG_doped.csv')
    plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%) DMAc EAN',sprintf(' %d',T),'C');
elseif strcmp(filename,'exp-data/PEG_undoped.csv')
    plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%) DMAc');
elseif strcmp(filename,'exp-data/PEG_film.csv')
    plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%)');
elseif strcmp(filename,'exp-data/PI.csv')
    plotname = strcat(NAME);
end
plotsaxs(s,head,plotname,'k');

set(gca,'xscale','log');set(gca,'yscale','log');
xlim([.1,3]);ylim([1e3,1e6])