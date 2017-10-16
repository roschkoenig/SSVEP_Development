function cxdata = snd_expression_interpolate(gxdata)

for g = 1:size(gxdata,3)    % each gene
for r = 1:size(gxdata,2)    % each region
    v = squeeze(gxdata(:,r,g));
    empt = find(isnan(v));
    if ~isempty(empt)
        if empt(1) == 1 
        elseif empt(end) == length(v)
            for e = 1:length(empt)
                v(empt(e)) = v(empt(e)-1);
            end
        else
            for e = 1:length(empt)
                v(empt(e)) = mean([v(empt(e)+1) v(empt(e)-1)]);
            end
        end
    end
    cxdata(:,r,g) = v;
end
end

if ~isempty(find(isnan(cxdata)))
    error('I have not fixed NaNs, sorry')
end
   