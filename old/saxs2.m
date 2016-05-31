clear;

% Load SAXS data
% filename = 'exp-data/PEG_film.csv';
% filename = 'exp-data/PEG_doped.csv';
filename = 'exp-data/PEG_undoped.csv';
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
MWV = [990,1500,2000,6000];
WP = 50;
% TV = [40,50,60,70,80];

figure;hold
leg = cell(1,length(MWV));
for ii = 1:length(MWV)
    MW  = MWV(ii);
    col = (ii-1)/(length(MWV)-1);
    
    if strcmp(filename,'exp-data/PEG_doped.csv')
        plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%) DMAc EAN',sprintf(' %d',T),'C');
    elseif strcmp(filename,'exp-data/PEG_undoped.csv')
        plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%) DMAc');
    elseif strcmp(filename,'exp-data/PEG_film.csv')
        plotname = strcat(NAME,sprintf('%d (%d',MW,WP),'%)');
    end
    plotsaxs(s,head,plotname,[col 0 1-col]);
    leg{ii} = strcat(sprintf('MW = %d',MW),'');
end

set(gca,'xscale','log');set(gca,'yscale','log');
xlim([.1,20]);ylim([1e3,1e6])
legend(leg)