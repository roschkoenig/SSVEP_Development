% This code applies principal component analysis to the Allen Brain
% developmental transcriptome dataset

clear all

%% Housekeeping
%==========================================================================
D   = snd_housekeeping;
fs  = filesep;

Fscripts = D.Fscripts;
Fdata    = D.Fdata;
Fgenex   = D.Fgenex;

% Load data array
%==========================================================================
genex = csvread([Fgenex fs 'expression_matrix.csv']);
[col ID name age gend strID strct strctN] = textread([Fgenex fs 'columns_metadata.csv'], '%s %s %s %q %s %s %q %q', 'delimiter', ',');
[row geneID ensembl geneS entrez]         = textread([Fgenex fs 'rows_metadata.csv'], '%s %s %q %q %s', 'delimiter', ',');

%% Select regions defined during time of interest (i.e. 4 years+)
%==========================================================================
% Identify all subjects aged four or older
%--------------------------------------------------------------------------
min_age = find(ismember(age, '4 yrs'));  
fourplus    = min_age(1):length(age);
agelbls     = unique(age(fourplus));

% Generate relevant vectors containing ages and structures
%--------------------------------------------------------------------------
for a = 1:length(agelbls)
    ageno(a) = str2double(agelbls{a}(1:end-3));
end
[sd sg] = sort(ageno);
agelbls = agelbls(sg);

struct_list = strct(min_age);

for s = 1:length(struct_list)
for a = 1:length(agelbls)
    ai = find(ismember(age, agelbls{a}));
    si = find(ismember(strct, struct_list{s}));
    if ~isempty(intersect(ai,si))
        datacol         = intersect(ai, si);
        gxdata(a,s,:)   = genex(:,datacol);
    end
    clear ai si
end
end

gxdata(gxdata == 0) = NaN;
gxdata              = snd_expression_interpolate(gxdata);
mxage               = squeeze(nanmean(gxdata,1)); 


[coef score latent tsquared explained] = pca(mxage);

% Reconstruct age changes in the PC across all samples
%=========================================================================
pcage = [];
for a = 1:size(gxdata,1)
    recon = squeeze(gxdata(a,:,:)) / coef';
    for r = 1:size(recon,2)
        recon(:,r) = recon(:,r) - mean(recon(:,r));
    end
    pcage(a,:,:) = recon;
end
 
%% Plotting summary of the PCA decomposition
%==========================================================================
subcortical = [3 9 13];
allo = [2 5];

subplot(2,2,1)
    imagesc(abs((coef)))
    % Labels 
    title('Log-scaled absolute coefficient values');
    xlabel('Components'); ylabel('Genes');
    set(gca, 'YTick', [1:length(struct_list)], 'YTickLabel', struct_list);
    % Settings
    axis square

subplot(2,2,3) 
    pareto(explained);
    title('Scree plot');
    
subplot(2,2,2)
    scatter(score(:,1), score(:,2), 'b', 'filled'), hold on
    scatter(score(allo,1), score(allo,2), 'm', 'filled'),
    scatter(score(subcortical,1), score(subcortical,2), 'r', 'filled'), hold off
    axis square
    xlabel('Component 1');  ylabel('Component 2');
    legend({'Neocortex', 'Allocortex', 'Subcortical'});
   	title('Examples of components');
    
subplot(2,2,4)
    scatter(score(:,6), score(:,7), 'b', 'filled'); hold on
    scatter(score(allo,6), score(allo,7), 'm', 'filled'),
	scatter(score(subcortical,6), score(subcortical,7), 'r', 'filled'), hold off
    axis square
    xlabel('Component 6 ');  ylabel('Component 7');

    

%% Plot example traces across development
%==========================================================================
v1c = squeeze(pcage(:,14,:));
vfc = squeeze(pcage(:,15,:));
ageso = sort(ageno);

% Colors
greens  = cbrewer('seq', 'Greens', 25);  
greens  = greens(end-size(vfc,1):end-1 , :);
reds    = cbrewer('seq', 'Reds', 25);
reds    = reds(end-size(vfc,1):end-1 , :);

subplot(2,1,1)
    % Plotting
    scatter(vfc(:,1), vfc(:,2), 50, greens, 'filled'), hold on
    scatter(v1c(:,1), v1c(:,2), 50, reds, 'filled')

    % Labels
    legend({'Ventral Frontal Ctx', 'Primary Visual Ctx'});
    xlabel('PC1', 'Fontweight', 'bold');
    ylabel('PC2', 'Fontweight', 'bold');
    
    % Settings
    axis square
    set(gcf, 'color', 'w');
    
subplot(2,1,2)
    % Plotting
    ageall = 4:40;
    vfcint = interp1(ageso, vfc(:,2), ageall);
    plot(ageso, vfc(:,2), 'o', 'color', greens(end-4,:)), hold on
    plot(ageall, vfcint, 'color', greens(end-4,:), 'Linewidth', 2);
    
    v1cint = interp1(ageso, v1c(:,2), ageall, 'linear');
    plot(ageso, v1c(:,2), 'o', 'color', reds(end-4,:));
    plot(ageall, v1cint, 'color', reds(end-4,:), 'Linewidth', 2);
    
    % Labels
    xlabel('Age in years', 'Fontweight', 'bold');
    ylabel('PC2', 'Fontweight', 'bold');
    
%% Generate linear interpolations of PCA and save
%==========================================================================
for r = 1:size(pcage, 2)
for c = 1:size(pcage, 3)
    pcint(:,r,c) = interp1(ageso, squeeze(pcage(:,r,c)), ageall);
end
end

GX.orig = pcage;
GX.intp = pcint;
GX.ageo = ageso;
GX.agei = ageall;
GX.expl = explained;

cumul   = cumsum(explained);
above   = find(cumul>95);
thr     = above(1);

for c = 1:thr
    weights = coef(:,c);
    [sd si] = sort(weights);
    sid     = flip(si);
    ranked(:,c) = geneS(sid(1:2000));
end

save([Fgenex fs 'gx_pca.mat'], 'GX');
cell2csv([Fgenex fs 'Ranked_Genes.csv'], ranked);


