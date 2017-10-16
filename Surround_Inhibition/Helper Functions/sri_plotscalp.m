function sri_plotscalp(volfile, varargin)
% Plotting the results of an statistical parametric map
% Usage: 
%--------------------------------------------------------------------------
% sri_plotscalp(volfile): 
% Plotting Tmap collapsed averaged by dimension
% 
% sri_plotscalp(volfile, threshold):
% Define the threshold that will be considered significant
%
% sri_plotscalp(volfile, threshold, averagetype);
% Define single frequency plane within which results will be plotted
% averagetype:  'average' (string)      - collapse as average over the dimension
%               f (double)              - extract single frequency plane

switch nargin
    case 1;     Tmin = 4;               plotting = 'average';
    case 2;     Tmin = varargin{1};     plotting = 'average';
    case 3;     Tmin = varargin{1};     plotting = varargin{2};
end

V               = spm_vol(volfile);
[Y XYZ]         = spm_read_vols(V);
frequencies     = unique(XYZ(3,:));

Tmap = Y * 0;
Tmap(Y > Tmin) = Y(Y > Tmin);


for r = 1:size(Tmap,1);
for c = 1:size(Tmap,2)
    rl  = size(Tmap,1);     cl = size(Tmap,2);
    r2  = floor(rl/2);      c2 = floor(rl/2);

    if sqrt((r - r2)^2 + (c - c2)^2) > r2
        mask(r,c) = 0; 
    else mask(r,c) = 1;
    end
end
end

% Plot single average collapsed onto scalp
%==========================================================================

if ischar(plotting) && strcmp('average', plotting)

    for r = 1:size(Tmap,1)
    for c = 1:size(Tmap,2)
        Smap(r,c) = mean(find(Tmap(r,c,:)));
        if isnan(Smap(r,c)), Smap(r,c) = 0; end
    end
    end
    
    
% Plot 3D average collapsed onto individual dimensions
%==========================================================================

elseif ischar(plotting) && strcmp('average3', plotting)

    for r = 1:size(Tmap,1)
    for c = 1:size(Tmap,2)
        Smap(r,c) = mean(find(Tmap(r,c,:)));
        if isnan(Smap(r,c)), Smap(r,c) = 0; end
    end
    end
    
  	for r = 1:size(Tmap,1)
    for f = 1:size(Tmap,3)
        LRmap(r,f) = mean(find(Tmap(r,:,f)));
        if isnan(LRmap(r,f)), LRmap(r,f) = 0; end
    end
 	end
    
	for c = 1:size(Tmap,2)
    for f = 1:size(Tmap,3)
        APmap(c,f) = mean(find(Tmap(:,c,f)));
        if isnan(APmap(c,f)), APmap(c,f) = 0; end
    end
    end
    disp('yeah');
    
    subplot(3,3,[2 3]);
    imagesc(LRmap')
    set(gca, 'Ydir', 'normal');
    
    subplot(3,3,[4 7]);
    imagesc(APmap);
    
    subplot(3,3,[5 6 8 9]);
    plotsmap(Smap, Tmin);


% Plot single plane defined by frequency
%==========================================================================
elseif isnumeric(plotting)
    ftoplot = plotting;
    fid     = nearest(frequencies, ftoplot);
    
    for r = 1:size(Tmap,1)
    for c = 1:size(Tmap,2)
        Smap(r,c) = Tmap(r,c,fid);
        if isnan(Smap(r,c)), Smap(r,c) = 0; end
    end
    end
    
    Smap = Smap .* mask;
    plotsmap(Smap, Tmin);

end
end

function plotsmap(Smap, Tmin)

    imagesc(Smap'); axis square; hold on
    colorbar;
    set(gca, 'Ydir', 'normal');
    set(gca, 'clim', [0 Tmin * 3]);
    xlabel('Left to Right');
    ylabel('Posterior to Anterior');
    
    
%     t = linspace(0,2*pi);
%     plot(r2*cos(t)+r2+.5,r2*sin(t)+r2+.5)

end

