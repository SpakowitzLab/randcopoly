function [qp,sp] = plotsaxs(s,head,plotname,col)
% example use:
% plotsaxs(s,head,'ODPA-AP6F-ODPA-PEG1500 (50%) DMAc EAN 40C','k');
% plotsaxs(s,head,'ODPA-AP6F-ODPA-PEG1500 (50%) DMAc EAN RT','k');
% plotsaxs(s,head,'ODPA-AP6F-ODPA-PEG1500 (50%) DMAc','k');

% find desired condition
x = strfind(head,plotname);
y = cellfun(@(x) ~isempty(x),x);
name = find(y == 1);
name = name(1); % find the first occurance (look out for '... step2', second experiment label)

indq = 4*(name-1)+1;
inds = 4*(name-1)+2;

% filter zero values
ind = find(s(:,indq)~=0);

% make a plot
qp = s(ind,indq);
sp = s(ind,inds);
plot(s(ind,indq),s(ind,inds),'color',col)