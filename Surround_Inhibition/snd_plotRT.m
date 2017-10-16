%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping;
fs  = filesep;

Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;
clear D 

subs        = spm_select('List', [Fanalysis fs 'MEEG'], '^A.*.mat');
subs        = cellstr(subs(:,1:9));

load([Fanalysis fs 'Fourier_summaries']);

timstps = size(S(1).ssvep,1);
ssvep   = [];
rctm    = [];
cptrace = [];


% Extract variables of interest from summary variable
%--------------------------------------------------------------------------
for s = 1:length(S)
    if size(S(s).ssvep,2) ~= size(S(s).events), disp(s); end
    ssvep = [ssvep; S(s).ssvep'];
    rctm  = [rctm; [S(s).events.rctm]'];
    cptrace = [cptrace; squeeze(S(s).cptrace)];
end

% Smooth bins of reaction times
%--------------------------------------------------------------------------
clear ssvsort avgbin avgsteps srtd sortid ssmean
[srtd sortid] = sort(rctm);
avgbin      = 80;
avgsteps    = 1:avgbin:size(rctm,1)-1;
ssvsort     = ssvep(sortid,:);
cpsort      = cptrace(sortid,:);

for a = 1:length(avgsteps)-1
    ssmean(a,:)     = mean(ssvsort(avgsteps(a):avgsteps(a) + avgbin,:),1);
    cpmean(a,:)     = mean(cpsort(avgsteps(a):avgsteps(a)+avgbin,:),1);
end

% Plotting routines
%==========================================================================
% Plot smoothed SSVEP 
%--------------------------------------------------------------------------
subplot(2,2,1)
    imagesc(S(1).timesteps, 1:length(ssvep), ssmean, [1 2]), hold on
    set(gca, 'Ydir', 'normal');
    plot(sort(rctm)/1000, 1:length(ssvep), 'Linewidth', 2);
    xlim([-.3 2]);
    
% Plot smoothed ERP decision component
%--------------------------------------------------------------------------
subplot(2,2,3)
    imagesc(cpmean);
    
% Plot time trace over SSVEP
%--------------------------------------------------------------------------
sortss = ssvep(sortid,:);
third = floor(size(sortss,1) / 3);

low     = sortss(1:third,:);                mlow = median(low,1);
med     = sortss([third+1]:[2*third],:);    mmed = median(med,1);
hig     = sortss([2*third+1]:[3*third],:);  mhig = median(hig,1);

subplot(2,2,2)
    plot(S(1).timesteps, mlow), hold on
    plot(S(1).timesteps, mmed),
    plot(S(1).timesteps, mhig)

    legend({'Low', 'Med', 'Hig'});

% Plot time trace over ERP decision component
%--------------------------------------------------------------------------
low     = cpsort(1:third,:);                mlow = mean(low,1);
med     = cpsort([third+1]:[2*third],:);    mmed = mean(med,1);
hig     = cpsort([2*third+1]:[3*third],:);  mhig = mean(hig,1);    
    
subplot(2,2,4)
    plot(mlow), hold on;
    plot(mmed)
    plot(mhig)
    
    legend({'Low', 'Med', 'Hig'});