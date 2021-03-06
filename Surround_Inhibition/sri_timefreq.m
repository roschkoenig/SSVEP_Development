%% Housekeeping
%==========================================================================
clear all
D   = snd_housekeeping('surround');
fs  = filesep;

Fscripts    = D.Fscripts;
Fdata       = D.Fdata;
Fgenex      = D.Fgenex;
Fanalysis   = D.Fanalysis;
Fmeeg       = D.Fmeeg;

clear D

% Load list of subjects with MEEG files
%--------------------------------------------------------------------------
meegs       = cellstr(spm_select('FPList', Fmeeg, '^ra.*.mat'));
subs        = spm_select('List', [Fanalysis fs 'MEEG'], '^A.*.mat');
dotpos      = find(subs(1,:) == '.');
subs        = cellstr(subs(:,1:dotpos(end)-1));

%% Loop through individual MEEGs
%==========================================================================
load(meegs{1});
chanxy  = spm_eeg_project3D(D.sensors.eeg, 'EEG');
chanlbl = D.sensors.eeg.label;

for m = 1:length(meegs)
    load(meegs{m});
    for c = 1:length(D.channels)
        chid    = find(strcmp(D.channels(c).label, chanlbl));
        D.channels(c).X_plot2D = chanxy(1,chid);
        D.channels(c).Y_plot2D = chanxy(2,chid);
    end 
    save(meegs{m}, 'D');
    clear D
end


for m = 1:length(meegs)
    clear job 
    job{1}.spm.meeg.tf.tf.D = meegs(m);
    job{1}.spm.meeg.tf.tf.channels{1}.all = 'all';
    job{1}.spm.meeg.tf.tf.frequencies = [1:0.5:30];
    job{1}.spm.meeg.tf.tf.timewin = [0 2400];
    job{1}.spm.meeg.tf.tf.method.mtmfft.taper = 'sine';
    job{1}.spm.meeg.tf.tf.method.mtmfft.freqres = 0.5;
    job{1}.spm.meeg.tf.tf.phase = 0;
    job{1}.spm.meeg.tf.tf.prefix = '';
    
    spm_jobman('run', job);
end

% Average across repeated trials
%--------------------------------------------------------------------------
tfs = cellstr(spm_select('FPList', Fmeeg, '^tf.*\.mat'));
for t = 1:length(tfs)
    clear job
    job{1}.spm.meeg.averaging.average.D = tfs(t);
    job{1}.spm.meeg.averaging.average.userobust.standard = true;
    job{1}.spm.meeg.averaging.average.plv = false;
    job{1}.spm.meeg.averaging.average.prefix = 'm';
    
    spm_jobman('run', job);
end

%% Generate scalp * frequency maps
%==========================================================================
% Setup frequency of interest
%--------------------------------------------------------------------------
mtfs = cellstr(spm_select('FPList', Fmeeg, '^mtf.*\.mat'));

for m = 1:length(mtfs)
    try 
        disp(['File ' num2str(m) '/' num2str(length(mtfs))]);
        clear job
        job{1}.spm.meeg.images.convert2images.D = mtfs(m);
        job{1}.spm.meeg.images.convert2images.mode = 'scalp x frequency'; 
        job{1}.spm.meeg.images.convert2images.channels{1}.all = 'all';
        job{1}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
        job{1}.spm.meeg.images.convert2images.freqwin = [1 30]; 
        job{1}.spm.meeg.images.convert2images.prefix = '';
        spm_jobman('run', job);
    catch 
        disp([mtfs{m} ' did not run']);
    end
end


% Smooth TF Data to allow averaging across individual differences
%--------------------------------------------------------------------------
dirs = cellstr(spm_select('FPList', Fmeeg, 'dir', '^mtf')); 

clear job
for d = 1:length(dirs)
    fnames = cellstr(spm_select('FPList', dirs{d}, '^condition*'));
    for f = 1:length(fnames)
        disp(['Subject ' num2str(d) '/' num2str(length(mtfs)) ', File ' num2str(f)]);
        fname = fnames{f};
        job{1}.spm.spatial.smooth.data = {fname};
        job{1}.spm.spatial.smooth.fwhm = [10 10 0];
        job{1}.spm.spatial.smooth.dtype = 0;
        job{1}.spm.spatial.smooth.im = 0;
        job{1}.spm.spatial.smooth.prefix = 's'; 
        spm_jobman('run', job);
    end
end


%% Second level - effects of experimental condition
%==========================================================================
% Subject specific folders containing smoothed, mean images
%--------------------------------------------------------------------------
dirs = cellstr(spm_select('FPList', Fmeeg, 'dir')); 
filtr = '^s*';

% Loop through folders to define filepaths and regressors
%--------------------------------------------------------------------------
clear job lvl1 K age pal vf split fullset fname scans bgcond fgcond
D       = snd_housekeeping('surround');
Fbase   = D.Fbase;
pheno   = csv2cell([Fdata fs 'MIPDB' fs 'MIPDB_PublicFile.csv'], 'fromfile');

[r agec] = find(~cellfun(@isempty, strfind(pheno, 'Age')));
agecells = pheno(2:end, agec);
for a = 1:length(agecells)
    ages(a) = str2double(agecells{a});
end

[r subc] = find(~cellfun(@isempty, strfind(pheno, 'ID')));
subcells = pheno(2:end, subc);

[r subc] = find(~cellfun(@isempty, strfind(pheno, 'DX_Status')));
diagnosis = pheno(2:end, subc);


count = 0;
for d = 1:length(dirs)
    fnames = cellstr(spm_select('FPList', dirs{d}, filtr));
    for f = 1:length(fnames)
        fname       = fnames{f};
        seppos      = find(fname == fs);
        longname    = fname(seppos(end-1)+1 : seppos(end)-1);
        subname     = longname(find(longname == 'A'):end);
        
        subi        = find(strcmp(subcells, subname));
        thisage     = ages(subi);
        thisdx      = diagnosis{subi};
        
        Floc = find(fname == 'F');  Floc = Floc(end)+1;
        Bloc = find(fname == 'B');  Bloc = Bloc(end)+1;
        
        if thisdx == '0'
            count           = count + 1;
            suball{count}   = subname;
            scans{count}    = fname;
            fgcond(count)   = str2double(fname(Floc)) / 3;
            bgcond(count)   = str2double(fname(Bloc));
            age(count)      = thisage;
        end
    end
end

scans = scans';

% Setup SPM Analysis
%--------------------------------------------------------------------------

clear job
Fspm = [Fanalysis fs 'SPM'];

job{1}.spm.stats.factorial_design.dir = {Fspm};
job{1}.spm.stats.factorial_design.des.mreg.scans    = scans;
job{1}.spm.stats.factorial_design.des.mreg.mcov     = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.des.mreg.incint   = 1;

job{1}.spm.stats.factorial_design.cov(1).c      = fgcond;
job{1}.spm.stats.factorial_design.cov(1).cname  = 'Foreground';
job{1}.spm.stats.factorial_design.cov(1).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(1).iCC    = 1;

job{1}.spm.stats.factorial_design.cov(2).c      = bgcond;
job{1}.spm.stats.factorial_design.cov(2).cname  = 'Background';
job{1}.spm.stats.factorial_design.cov(2).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(2).iCC    = 1;

agecv   = log(age) - mean(log(age));
job{1}.spm.stats.factorial_design.cov(3).c      = agecv;
job{1}.spm.stats.factorial_design.cov(3).cname  = 'Age';
job{1}.spm.stats.factorial_design.cov(3).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(3).iCC    = 1;

fgage = fgcond .* (agecv);
job{1}.spm.stats.factorial_design.cov(4).c      = fgage;
job{1}.spm.stats.factorial_design.cov(4).cname  = 'Age x Foreground';
job{1}.spm.stats.factorial_design.cov(4).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(4).iCC    = 1;

bgage = bgcond .* (agecv);
job{1}.spm.stats.factorial_design.cov(5).c      = bgage;
job{1}.spm.stats.factorial_design.cov(5).cname  = 'Age x Background';
job{1}.spm.stats.factorial_design.cov(5).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(5).iCC    = 1;

fgbg  = fgcond .* bgcond;
job{1}.spm.stats.factorial_design.cov(5).c      = fgbg;
job{1}.spm.stats.factorial_design.cov(5).cname  = 'Foreground x Background';
job{1}.spm.stats.factorial_design.cov(5).iCFI   = 1;
job{1}.spm.stats.factorial_design.cov(5).iCC    = 1;


job{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
job{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
job{1}.spm.stats.factorial_design.masking.im = 1;
job{1}.spm.stats.factorial_design.masking.em = {''};
job{1}.spm.stats.factorial_design.globalc.g_omit = 1;
job{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
job{1}.spm.stats.factorial_design.globalm.glonorm = 1;

job{2}.spm.stats.fmri_est.spmmat = {[Fspm fs 'SPM.mat']};
job{2}.spm.stats.fmri_est.write_residuals = 0;
job{2}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', job);
clear job


% Define Contrasts
%--------------------------------------------------------------------------
job{1}.spm.stats.con.spmmat = {[Fspm fs 'SPM.mat']};

job{1}.spm.stats.con.consess{1}.tcon.name = 'Foreground';
job{1}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0 0 0];
job{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';

job{1}.spm.stats.con.consess{2}.tcon.name = 'Background';
job{1}.spm.stats.con.consess{2}.tcon.weights = [0 0 -1 0 0 0];
job{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

job{1}.spm.stats.con.consess{3}.fcon.name = 'Age';
job{1}.spm.stats.con.consess{3}.fcon.weights = [0 0 0 1 0 0];
job{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';

job{1}.spm.stats.con.consess{4}.tcon.name = 'Age x Foreground';
job{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 0 1 0];
job{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

job{1}.spm.stats.con.consess{5}.tcon.name = 'Age x Background';
job{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 0 1];
job{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

job{1}.spm.stats.con.delete = 1;
spm_jobman('run', job); 

d = job{1}.spm.stats.con.consess;

slist = unique(suball);
save([Fspm fs 'Subjectlist'], 'slist');


%% This is where code goes to die
%==========================================================================

% %% Average first level scalp maps across individuals
% %--------------------------------------------------------------------------
% clear job flist fname
% flist = cellstr(spm_select('FPList', Fmeeg, 'dir', '^BL*'));
% count = 0;
% clear fname
% for f = 1:length(flist)
%     fnames = cellstr(spm_select('FPList', flist{f}, '^s*'));
%     for ff = 1:length(fnames)
%         count = count + 1;
%         fname{count} = fnames{ff};
% 	
%         if count == 1,      cstring = '(i1';
%         else                cstring = [cstring '+i' num2str(count)]; 
%         end
%     end
% end
% 
% fname   = fname';
% cstring = [cstring ')/' num2str(length(flist))];
% 
% job{1}.spm.util.imcalc.input = fname;
% job{1}.spm.util.imcalc.output = 'mean_sEC';
% job{1}.spm.util.imcalc.outdir = {Fmeeg};
% job{1}.spm.util.imcalc.expression = cstring;
% job{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% job{1}.spm.util.imcalc.options.dmtx = 0;
% job{1}.spm.util.imcalc.options.mask = 0;
% job{1}.spm.util.imcalc.options.interp = 1;
% job{1}.spm.util.imcalc.options.dtype = 4;
% 
% spm_jobman('run', job)
% 
% %% Calculate difference between baseline and SSVEP freq
% %--------------------------------------------------------------------------
% clear job flist fname
% flist_BL = cellstr(spm_select('FPList', Fmeeg, 'dir', '^24*'));
% flist_TG = cellstr(spm_select('FPList', Fmeeg, 'dir', '^25*'));
% 
% frq = {num2str(24), num2str(25)};
% 
% for s = 1 %:length(subs)
%     sub = subs{s};
%     BL  = [Fmeeg fs frq{1} 'mtf_ra' sub];
%     TG  = [Fmeeg fs frq{2} 'mtf_ra' sub];
%     if exist(BL, 'dir') == 7, blyes = 1; else disp(['Problem with baseline of ' sub]);   end
%     if exist(TG, 'dir') == 7, tgyes = 1; else disp(['Problem with target frq of ' sub]); end
%     
%     if blyes && tgyes
%         blfiles = cellstr(spm_select('List', BL, '^scondition*'));
%         tgfiles = cellstr(spm_select('List', TG, '^scondition*'));
%        
%         for b = 1:length(blfiles)
%         if strcmp(blfiles{b}, tgfiles{b})
%             disp('Yah')
%             clear job
%             made = mkdir([Fmeeg fs 'ssvep' sub]);
%             job{1}.spm.util.imcalc.input = {[BL fs blfiles{b}]; [TG fs tgfiles{b}]};
%             job{1}.spm.util.imcalc.output = blfiles{b};
%             job{1}.spm.util.imcalc.outdir = {[Fmeeg fs 'ssvep' sub]};
%             job{1}.spm.util.imcalc.expression = 'i2-i1';
%             job{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
%             job{1}.spm.util.imcalc.options.dmtx = 0;
%             job{1}.spm.util.imcalc.options.mask = 0;
%             job{1}.spm.util.imcalc.options.interp = 1;
%             job{1}.spm.util.imcalc.options.dtype = 4;
%             spm_jobman('run', job);
%         end
%         end
%     end
% end

% %% Generate scalp by time maps
% %--------------------------------------------------------------------------
% mtfs = cellstr(spm_select('FPList', Fmeeg, '^mtf.*\.mat'));
% for m = 1:length(mtfs)
%     clear job
%     job{1}.spm.meeg.images.convert2images.D = mtfs(m);
%     job{1}.spm.meeg.images.convert2images.mode = 'scalp x time';
%     job{1}.spm.meeg.images.convert2images.conditions = cell(1,0);
%     job{1}.spm.meeg.images.convert2images.channels{1}.all = 'all';
%     job{1}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
%     job{1}.spm.meeg.images.convert2images.freqwin = [25 25];
%     job{1}.spm.meeg.images.convert2images.prefix = '';
%     spm_jobman('run', job);
% end

