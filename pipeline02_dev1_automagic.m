
function [ score_table ] = pipeline02_dev1_automagic(eeglab_dir, code_dir, data_dir, subj_name_edf)
%function [ score_table, eeg_ICA ] = pipeline01_dev2_happe(eeglab_dir, code_dir, data_dir, subj_name_edf)

%    This function performs a series of preprocessing steps
%    according to the standardize procedure 
%    "AUTOMAGIC" implemented by Pedroni, Langer in 2019

%    https://github.com/methlabUZH/automagic
%    
%    and compute the "scopepochs" after each step of the preprocessing

% INPUTS: 
%   subj_name_edf : file with a single subject data in .edf format
%                    (i.e.: subj_name_edf = 'S003R01.edf'
%   optional (already set in the function):
%       eeglab_dir
%            !!! eeglab versione 2021226 requires this list of PLUGIN: 
%                Biosig3.7.5;
%                Cleanline1.04; 
%                clean_rawdata2.3
%                PrepPipeline0.55.4
%                ICLabel
% 
%
%       code_dir (with scorepochs package)
%       subj_dir (also for saving intermediate steps)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% EXTRA PARAMETERs section to check/set

% OUTPUTs: 
%   score_struct (with score_Xep at each step)

%   andrea.vitale@gmail.com 
%   20210430

%%
    %clear; close all

    % SET DIR
    if ~exist('eeglab_dir') %isempty('eeglab_dir')
        eeglab_dir = 'D:\IIT\EEG\eeglab_20201226'
        cd(eeglab_dir)
        %eeglab
        eeglab('nogui')
    end
    
    if ~exist('code_dir')
        code_dir = 'D:\IIT\_PROJECT\SCORE_epoch\code';
        addpath(genpath(code_dir))
    end

    if ~exist('data_dir')
        data_dir = 'D:\IIT\_PROJECT\SCORE_epoch\data';
        % this folder should contain also the CHANNEL INFO (.txt file)
        cd(data_dir)
    end
    
    if ~exist('subj_name_edf') 
        disp('!!! subj_name is required as INPUT') 
    end
   
    
    %% EXTRA PARAMETERs: - - - - - - - - - - - - - - 
    do_chan_pruning = 1
    do_prep = 1
    
    do_plot_chan = 0
    do_plot_PSD = 0
    
    do_save_prep = 0
    do_save_prep_ICA = 1
    do_save_score = 1
    
    
    %% INPUT: resting state EYES OPEN (R01)
    
    % R01 = resting state eyes open - - - - - - - - - - - 
    %subj_name_edf = 'S003R01.edf';
    %file_name = 'S001R01.edf'; 
    %file_name = 'S002R01.edf'; 
    %file_name = 'S003R01.edf'; % score epoch > 95% already for raw data 
    %file_name = 'S010R01.edf'; % score epoch > 90% already for raw data 

    % R02 = resting state eyes close - - - - - - - - - - - 
    %file_name = 'S001R02.edf';

    chan_file = 'coord_BS_motorEEG.txt';

    
%0) LOAD .edf file - - - - - - -
    eeg_struct = []; 
    eeg_struct = pop_biosig(fullfile(data_dir, subj_name_edf));
    %EEG = pop_biosig('D:\IIT\_PROJECT\SCORE_epoch\data\S001R01.edf');
    % eeg_struct = EEG;

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Nose direction should be set from '-Y' to default +X in EEG.chanlocs
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    eeg_struct=pop_chanedit(eeg_struct, 'load',{fullfile(data_dir, chan_file),'filetype','xyz'},'nosedir','-Y');
    %eeg_struct=pop_chanedit(eeg_struct, 'load',{'D:\\IIT\\_PROJECT\\SCORE_epoch\\data\\coord_BS_motorEEG.txt','filetype','xyz'},'nosedir','-Y');
    %EEG=pop_chanedit(EEG, 'lookup','D:\\IIT\\EEG\\eeglab_20201226\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','load',{'D:\\IIT\\_PROJECT\\SCORE_epoch\\data\\coord_BS_motorEEG.txt','filetype','xyz'},'nosedir','-Y');

    % check PLOT
    if do_plot_chan
        %figure;
        %pop_eegplot( eeg_struct, 1, 1, 1);

        figure; 
        % by NUMBER
        subplot 121
        topoplot([], eeg_struct.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', eeg_struct.chaninfo);
        % by LABEL
        subplot 122
        topoplot([],eeg_struct.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', eeg_struct.chaninfo);
    end
    
    %%
    % SOME CHECKS - - - - - -
    sample_rate = eeg_struct.srate;
    % length in sec of the recording:
    n_sample = eeg_struct.pnts;
    n_sample / sample_rate;  %in sec
    n_chan = eeg_struct.nbchan;

    % number of channel that can be retained for ICA
    %(number of channel)^2 x 20 to 30 data points to perform ICA
    if n_sample > n_chan^2 * 20
        disp([ num2str ' channels can be given as input to ICA'])
    else
        n_chan_max = sqrt(n_sample/20)
        disp([ 'number of channels for ICA should be reduced to ' num2str(n_chan_max)])
    end
    
    
    %CHANNEL PRUNING (to 19 channels):
    chan_toinclude = [ 62, 63, ...
                   48, 52, 56, 53, 49, ...
                   28, 34, 38, 35, 30, ...
                   17, 13, 10, 14, 18, ...
                       3, 4 ]; 
                   % Fpz=64 , Oz=1 not included !!!
                   % Tp9  and Tp10 (mastoids) not included
 
    
% 1) RAW DATA = = = = = = = = = = = = = = = = = = = = = =
    %% remove segment of the data (at the end of the recording) with 0 values
    % !!!! specific for BCI2000 dataset
    i_chan = 1;
    for i_sample = 1:n_sample
        if eeg_struct.data(i_chan,end-i_sample) ~= 0
            zero_last_sample = n_sample-i_sample;
            break
        end
    end
    
    % reject last timepoints
    eeg_struct = eeg_eegrej( eeg_struct, [zero_last_sample n_sample] );
    eeg_raw = eeg_struct;
    
    
    % = = = = = = = = = = = = = = = = = = = = = = == = =
    %% RUN SCOREPOCHS at each step of the preprocessing:
    %INPUT
    %    cfg struct with the following fields
    %           freqRange    - array with the frequency range used to compute the power
    %                          spectrum (see MATLAB pwelch function)
    %           fs           - integer representing sample frequency         
    %           windowL      - integer representing the window length (in seconds)  
    %           smoothFactor - smoothing factor for the power spectrum
    cfg = []; 
    % <<<<<<<<<<<<<<<<<< ENTIRE FREQUENCY RANGE<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %cfg.freqRange = [ 1 : 80 ];
    % <<<<<<<<<<<<<<<<<< ONLY ALPHA BAND <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cfg.freqRange = [ 8 : 13 ]; 
    cfg.fs = eeg_struct.srate;
    cfg.windowL = 5; % in sec <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    cfg.smoothFactor = 0;

       
% ----------------------------------------------------------------
%2) BAD CHANNEL IDENTIFICATION (and REMOVAL) 
%   and (iterative) ROBUST AVERAGE RE-REFERENCING 
%   using the PREP pipeline (Bigdely-Shamlo 2015)
    
    hpf_cutoff = 1;
    %lpf_cutoff = 80;
    lpf_cutoff = eeg_struct.srate/2 -1;
  
    line_noise_freq = 60; %<<<<<<<<<<< TO SET <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  
    chan_vector = [ 1 : n_chan ];
    
    params = [];
    params.name = subj_name_edf(1:end-4);
    params.referenceChannels = chan_vector;
    params.evaluationChannels = chan_vector;
    params.rereferencedChannels = chan_vector; 
    params.lineFrequencies = [ line_noise_freq ];

% For a full listing of the defaults for each step in the pipeline, execute:
%     showPrepDefaults(EEG).

    
    [eeg_prep, params, computationTimes] = prepPipeline(eeg_struct, params)
    
       
%     eeg_prep2 = pop_prepPipeline(eeg_struct, struct('ignoreBoundaryEvents', true, ...
%         'lineNoiseChannels', chan_vector, ...
%         'lineFrequencies', line_noise_freq, 'Fs', sample_rate, 'p', 0.01, ...
%         'fScanBandWidth', 2, 'taperBandWidth', 2, 'taperWindowSize', 4, 'pad', 0, ...
%         'taperWindowStep', 1, 'fPassBand', [0  lpf_cutoff], 'tau', 100, 'maximumIterations', 10, ...
%         'referenceChannels', chan_vector, ...
%         'evaluationChannels', chan_vector, ...
%         'rereferencedChannels', chan_vector, ...
%         'ransacOff', false, 'ransacSampleSize', 50, 'ransacChannelFraction', 0.25, ...
%         'ransacCorrelationThreshold', 0.75, 'ransacUnbrokenTime', 0.4, 'ransacWindowSeconds', 5, ...
%         'srate', sample_rate, 'robustDeviationThreshold', 5, 'correlationWindowSeconds', 1, ...
%         'highFrequencyNoiseThreshold', 5, 'correlationThreshold', 0.4, 'badTimeThreshold', 0.01, ...
%         'maxReferenceIterations', 4, 'referenceType', 'Robust', 'reportingLevel', ...
%         'Verbose', 'interpolationOrder', 'Post-reference', 'meanEstimateType', 'Median', ...
%         'samples', n_sample, 'reportMode', 'normal', 'publishOn', true, ...
%         'consoleFID', 1, 'cleanupReference', false, 'keepFiltered', true, 'removeInterpolatedChannels', true)); 

%         %'sessionFilePath', '.\S001R01Report.pdf', 'summaryFilePath', '.\S001R01Summary.html', ...


    % = = = = = = = = = = = =  = = = = = = = = = = = = = =
    % or thorugh AUTOMAGIC function:
%    params = struct('FilterParams', struct('notch', struct('freq', 50), ...
%                                        'high',  struct('freq', 0.5,...
%                                                        'order', []),...
%                                        'low',   struct('freq', 30,...
%                                                         'order',
%                                                         []))),...
%                   'CRDParams',    struct('ChannelCriterion',   0.85,...
%                                          'LineNoiseCriterion', 4,...
%                                          'BurstCriterion',     5,...
%                                          'WindowCriterion',    0.25, ...
%                                          'Highpass',           [0.25 0.75])
%    [eeg_automagic, plots] = preprocess(eeg_struct, params);
  

% ----------------------------------------------------------------
% ?? in the eeg_prep.etc
    % there's still a field called : stillNoisyChannelNumbers ??
    % what we should do ?? NOT CLEAR IF these noisy CHANNELs should be interpolated 
    
    do_chan_interp2 = 0 % for now NO
    if do_chan_interp2
        bad_chan_prep = eeg_prep.etc.sillNoisyChannelNumbers;
        bad_chan_label = {};
        for ii = 1:length(bad_chan_prep)
            i_chan = bad_chan_prep(ii);
            bad_chan_label{ii} = eeg_prep.chanlocs(i_chan).labels;
        end
    
     
        %if do_interp_badchan
            % INTERPOLATE instead of remove
            % based on the channel location of the initial eeg_struct
            eeg_prep_badchan_interp = pop_interp(eeg_prep, eeg_prep.chanlocs, 'spherical');
        %end
        % ?? do NOT INTERPOLATE before ICA ??
    end

% ----------------------------------------------------------------
%   !!! HIGH- and LOW- pass FILTERS probably NOT necessary
%     eeg_struct = pop_eegfiltnew(eeg_struct, hpf_cutoff, [], [],0,[],0);
%     eeg_hpf = eeg_struct;
%     
%     % - - - 
%     eeg_struct = pop_eegfiltnew(eeg_struct, [], lpf_cutoff, [],0,[],0);
%     eeg_lpf = eeg_struct;
        

%3  NOTCH FILTER seems not necessary
    eeg_prep_notch = pop_eegfiltnew(eeg_prep, 'locutoff',line_noise_freq-2, ...
                              'hicutoff',line_noise_freq+2,'revfilt',1,'plotfreqz',1);

    if do_plot_PSD
        % CHECK the PSD (before and after LINE NOISE removal)                                
        figure; subplot 131; %!!!
        pop_spectopo(eeg_struct, 1, [ ], 'EEG' , 'percent', 50, 'freq', [8 13 20], 'freqrange',[2 lpf_cutoff],'electrodes','off');
        
        subplot 132;
        pop_spectopo(eeg_prep, 1, [ ], 'EEG' , 'percent', 50, 'freq', [8 13 20], 'freqrange',[2 lpf_cutoff],'electrodes','off');
    
        subplot 133;
        pop_spectopo(eeg_prep_notch, 1, [ ], 'EEG' , 'percent', 50, 'freq', [8 13 20], 'freqrange',[2 lpf_cutoff],'electrodes','off');
    end
    
% ----------------------------------------------------------------
%4  RE-REFERENCE to the AVERAGE probably not necesasry !!!!   
    eeg_prep_avgref = pop_reref(eeg_prep_notch, []);
    %eeg_prep_avgref = pop_reref(eeg_prep_badchan_interp, []);
    %eeg_avgref = eeg_struct;
    
% ----------------------------------------------------------------
%5  CHAN PRUNING + IC decomposition (and rejection)    
    
    eeg_prep_red = pop_select(eeg_prep, 'channel', chan_toinclude);
    

    eeg_prep_ICA = pop_runica(eeg_prep_red, 'icatype', 'runica', 'extended',1,'interrupt','on');

    %%(PLOT component topography:)
    % (see: https://github.com/sccn/viewprops)
    %pop_topoplot(eeg_ICA, 0, [1:length(chan_toinclude)] ,'EDF file',[5 5] ,0,'electrodes','on');
    
    eeg_prep_ICA = pop_iclabel(eeg_prep_ICA, 'default');
                     % for component viewing
    %pop_viewprops(eeg_cleanraw_ICA, 0, [1:eeg_cleanraw_ICA.nbchan], [2 50]) % for component properties
    %pop_viewprops( eeg_ICA, 0, [1:length(chan_toinclude)], [2 80], [], 0, eeg_ICA.etc.ic_classification, 'on')
        %spec_opt, erp_opt, scroll_event, classifier_name, fig)

    if do_save_prep_ICA
        cd(fullfile(data_dir))
        pop_saveset(eeg_prep_ICA, 'filename', [subj_name_edf(1:end-4) '_prep_ICA'])
    end
    
    % - - - -  - - - - - - - - - - - - - - - -
    % REMOVE BAD COMPONENT based on ICLabels
    %eeg_nobadICA = pop_icflag(eeg_ICA, [NaN NaN;0.8 1;0.8 1;NaN NaN;0.8 1;0.8 1;0.8 1]);

    % if the % of brain ICA < 0.2 -> then is removed
    eeg_prep_nobadICA = pop_icflag(eeg_prep_ICA, [0 0.2;0.7 1;0.7 1;NaN NaN;0.7 1;0.7 1;0.7 1]);
    %eeg_nobadICA = EEG;

    % or manually: 
    % bad_ICA_idx = []; 
    % eeg_struct = pop_subcomp(eeg_struct, bad_ICA_idx, 0);

end

    
    
    
    % = ==========================================================
    %% COMPUTE SCOREPOCHS at each step:
    % = ==========================================================
    
    prep_step = {
            'eeg_raw';
            'eeg_prep';
            'eeg_prep_notch';
            'eeg_prep_avgref';
            
            'eeg_prep_nobadICA';
                };
             
    % CREATE A SCORE STRUCT for final report:
    % with epoch not sorted !!!
    score_struct = [];
    
    for i_step = 1:length(prep_step)
        eval(['eeg_step_tmp = ' prep_step{i_step} ]);
        
        % reduce to 19 che number of channels 
        if eeg_step_tmp.nbchan > length(chan_toinclude)
            eeg_step_tmp =  pop_select(eeg_step_tmp, 'channel', chan_toinclude);
        end
        
        [idx_best_ep,epoch,score_Xep] = scorEpochs(cfg, eeg_step_tmp.data);
        eval([ 'score_struct.' prep_step{i_step} '= score_Xep' ]);
        disp(mean(score_Xep))
    end
  
    
    if do_save_score
        save_name = [subj_name_edf(1:end-4) '_scorepoch_pipeline02' ];
        save(save_name, 'score_struct')
        
%         % REPORT other METRICS:
%         n_chan; 
%         chan_interpolated = []; 
%         n_badICA =[];
%         Percent_Variance_Kept_of_Post_Waveleted_Data=[];

    end
% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end