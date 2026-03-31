function [sub_check] = tmfc_PPI(tmfc,ROI_set_number,start_sub)

% ========= Task-Modulated Functional Connectivity (TMFC) toolbox =========
%
% Calculates psychophysiological interactions (PPIs).
% Whitening is applied during deconvolution, consistent with SPM PEB assumptions.
% Mean centering of the psychological regressor (PSY) can be enabled or disabled.
% In the subsequent gPPI model estimation, the raw (not whitened) BOLD signal
% is used for the PHYS regressor to avoid double whitening (see He et al., 2025).
%
%
% FORMAT [sub_check] = tmfc_PPI(tmfc)
% Run a function starting from the first subject in the list.
%
%   tmfc.subjects.path            - Paths to individual SPM.mat files
%   tmfc.subjects.name            - Subject names within the TMFC project
%                                   ('Subject_XXXX' naming will be used by default)
%   tmfc.project_path             - Path where all results will be saved
%   tmfc.defaults.parallel        - 0 or 1 (sequential/parallel computing)
%
%   tmfc.ROI_set                  - List of selected ROIs
%   tmfc.ROI_set.PPI_centering    - Apply mean centering of psychological
%                                   regressor (PSY) prior to deconvolution:
%                                   'with_mean_centering' (default)
%                                   or 'no_mean_centering'
%                                   (see Di, Reynolds & Biswal, 2017; Masharipov et al., 2024)
%
%   tmfc.ROI_set.type             - Type of the ROI set
%   tmfc.ROI_set.set_name         - Name of the ROI set
%   tmfc.ROI_set.ROIs.name        - Name of the selected ROI
%   tmfc.ROI_set.ROIs.path_masked - Paths to the ROI images masked by group
%                                   mean binary mask 
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions        - List of conditions of interest
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.sess   - Session number (as specified in SPM.Sess)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.number - Condition number (as specified in SPM.Sess.U)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.pmod   - Parametric/Time modulator number (see SPM.Sess.U.P)
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.name   - Condition name (as specified in SPM.Sess.U.name(kPmod))
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions.file_name - Condition-specific file names:
%   (['[Sess_' num2str(iSess) ']_[Cond_' num2str(jCond) ']_[' ...
%    regexprep(char(SPM.Sess(iSess).U(jCond).name(1)),' ','_') ']'];)
%
% Session number and condition number must match the original SPM.mat file.
% Consider, for example, a task design with two sessions. Both sessions 
% contain three task regressors for "Cond A", "Cond B" and "Errors". If
% you are only interested in comparing "Cond A" and "Cond B", the following
% structure must be specified (see tmfc_conditions_GUI, nested function:
% [cond_list] = generate_conditions(SPM_path)):
%
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).sess   = 1;   
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).number = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).name = 'Cond_A'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(1).file_name = '[Sess_1]_[Cond_1]_[Cond_A]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).sess   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).name = 'Cond_B'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(2).file_name = '[Sess_1]_[Cond_1]_[Cond_B]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).number = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).pmod   = 1; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).name = 'Cond_A'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(3).file_name = '[Sess_2]_[Cond_1]_[Cond_A]';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).sess   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).pmod   = 1;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).name = 'Cond_B'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(4).file_name = '[Sess_2]_[Cond_2]_[Cond_B]';
%
% If GLMs contain parametric or time modulators, add the following fields:
% e.g. first modulator for fourth condition:
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).sess   = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).number = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).pmod   = 2;
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).name = 'Cond_BxModulator1^1';
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(5).file_name = '[Sess_2]_[Cond_2]_[Cond_BxModulator1^1]'; 
% e.g. second modulator for fourth condition:
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).sess   = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).number = 2; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).pmod = 3; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).name = 'Cond_BxModulator2^1'; 
%   tmfc.ROI_set(ROI_set_number).gPPI.conditions(6).file_name = '[Sess_2]_[Cond_2]_[Cond_BxModulator2^1]'; 
%
% Example of the ROI set (see tmfc_select_ROIs_GUI):
%
%   tmfc.ROI_set(1).set_name = 'two_ROIs';
%   tmfc.ROI_set(1).type = 'binary_images';
%   tmfc.ROI_set(1).ROIs(1).name = 'ROI_1';
%   tmfc.ROI_set(1).ROIs(2).name = 'ROI_2';
%   tmfc.ROI_set(1).ROIs(1).path_masked = 'C:\ROI_set\two_ROIs\ROI_1.nii';
%   tmfc.ROI_set(1).ROIs(2).path_masked = 'C:\ROI_set\two_ROIs\ROI_2.nii';
%
% FORMAT [sub_check] = tmfc_PPI(tmfc,ROI_set_number,start_sub)
% Run the function starting from a specific subject in the path list for
% the selected ROI set.
%
%   tmfc           - As above
%   ROI_set_number - Number of the ROI set in the tmfc structure
%   start_sub      - Subject number in the list to start computations from
%
% =========================================================================
% Copyright (C) 2025 Ruslan Masharipov
% License: GPL-3.0-or-later
% Contact: masharipov@ihb.spb.ru

if nargin == 1
   ROI_set_number = 1;
   start_sub = 1;
elseif nargin == 2
   start_sub = 1;
end

if ~isfield(tmfc.ROI_set(ROI_set_number),'PPI_centering')
    tmfc.ROI_set(ROI_set_number).PPI_centering = 'with_mean_centering';
elseif isempty(tmfc.ROI_set(ROI_set_number).PPI_centering)
    tmfc.ROI_set(ROI_set_number).PPI_centering = 'with_mean_centering';
end

% Check subject names
if ~isfield(tmfc.subjects,'name')
    for iSub = 1:length(tmfc.subjects)
        tmfc.subjects(iSub).name = ['Subject_' num2str(iSub,'%04.f')];
    end
end

% -------------------------------------------------------------------------
% Basic setup
% -------------------------------------------------------------------------
nSub = length(tmfc.subjects);
nROI = length(tmfc.ROI_set(ROI_set_number).ROIs);

cond_list = tmfc.ROI_set(ROI_set_number).gPPI.conditions;
nCond = length(cond_list);

sess = []; sess_num = []; 
for iCond = 1:nCond
    sess(iCond) = cond_list(iCond).sess;
end
sess_num = unique(sess);
maxSess  = max(sess_num);

conds_by_sess = cell(1, maxSess);
for jCond = 1:nCond
    s = cond_list(jCond).sess;
    conds_by_sess{s}(end+1) = jCond;
end

sub_check = zeros(1,nSub);
if start_sub > 1
    sub_check(1:start_sub) = 1;
end

% Initialize SPM
spm('defaults','fmri');

% -------------------------------------------------------------------------
% NOTE ABOUT PARFOR (INTENTIONALLY DISABLED BY DEFAULT)
% -------------------------------------------------------------------------
% spm_PEB (used inside tmfc_PEB_PPI) involves inversion/solves of large
% matrices and is already implicitly multithreaded via high-performance
% linear algebra libraries (e.g., BLAS/LAPACK). In many practical cases,
% using parfor here can be slower due to worker overhead and CPU
% oversubscription. If you want to benchmark parfor speed on your machine,
% you can set the flag below to false.

force_disable_parfor = true;   % <-- set to false to test parfor on your PC

% -------------------------------------------------------------------------
% Prepare CACHE
% -------------------------------------------------------------------------
CACHE = struct();
CACHE.project_path = tmfc.project_path;

CACHE.ROI_set_number = ROI_set_number;
CACHE.ROIset         = tmfc.ROI_set(ROI_set_number);

CACHE.conditions.list       = cond_list;
CACHE.conditions.by_session = conds_by_sess;
CACHE.conditions.nCond      = nCond;

CACHE.sessions.list = sess_num;
CACHE.sessions.max  = maxSess;

CACHE.subjects.sublist = tmfc.subjects;
CACHE.subjects.nSub = nSub; 

% -------------------------------------------------------------------------
% Subject loop 
% -------------------------------------------------------------------------

% Sequential computations
% -------------------------------------------------------------------------
if tmfc.defaults.parallel == 0 || force_disable_parfor

    for iSub = start_sub:nSub
        tmfc_PEB_PPI(CACHE,iSub);
        sub_check(iSub) = 1;
    end

% Parallel computations
% -------------------------------------------------------------------------
else

    % Parallel Loop
    try
        if isempty(gcp('nocreate')), parpool; end
        figure(findobj('Tag','TMFC_GUI'));
    end

    CACHEc = parallel.pool.Constant(CACHE);
    parfor iSub = start_sub:nSub
        tmfc_PEB_PPI(CACHEc.Value, iSub);
        sub_check(iSub) = 1;
    end
end   
