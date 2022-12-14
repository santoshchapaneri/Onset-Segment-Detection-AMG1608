% Segment Level features of AMG1608

% What we have so far:
% 1. From wav files, we obtained 70D features per frame with MIRToolbox,
% stored in FeatsX.mat, where X is song number
% 2. We also obtained onsets of each audio file using OnsetDetector_CLR
% python file and the output is stored in X.txt indicating onset times
% 3. Now, we will do block-based segment level features by captuing mean
% and std of each segment, instead of working at frame level

clear; clc;
% Directory of the frame level features
% audioPath = 'D:\MusicPhD\2017\AMG1608_release\amg1608_wav_IDs\';
audioPath = 'D:\Santosh\MusicPhD\Datasets\AMG1608_release\amg1608_wav_IDs\';
namesAudio = dir(fullfile(audioPath,'*.wav') );
namesAudio = {namesAudio(~[namesAudio.isdir]).name};
[namesAudioSorted,idx] = sort_nat(namesAudio);

% onsetPath = 'D:\MusicPhD\2017\AMG1608_release\AMG1608AcousticPosterior\AMG1608Onsets\';
onsetPath = 'D:\Santosh\MusicPhD\Datasets\AMG1608_release\AMG1608OnsetsSegments\AMG1608Onsets\';
namesOnsets = dir(fullfile(onsetPath,'*.txt') );
namesOnsets = {namesOnsets(~[namesOnsets.isdir]).name};
[namesOnsetsSorted,idx1] = sort_nat(namesOnsets);

% featPath = 'D:\MusicPhD\2017\AMG1608_release\AMG1608AcousticPosterior\AMG1608Features\';
featPath = 'D:\Santosh\MusicPhD\Datasets\AMG1608_release\AMG1608AcousticPosterior\AMG1608Features\';
namesFeats = dir(fullfile(featPath,'*.mat') );
namesFeats = {namesFeats(~[namesFeats.isdir]).name};
[namesFeatsSorted,idx2] = sort_nat(namesFeats);

SegmentFeats_cell = cell(numel(namesFeatsSorted),1);

for numSample = 1:numel(namesFeatsSorted)
    fprintf('%d\n',numSample);
    
    % Read audio and find frame indices start and end
    strAudio = sprintf('%s%d.wav',audioPath,numSample);
    [x,sr]=audioread(strAudio);
    w(1) = 1; w(2) = length(x);
    lsz = w(2)-w(1)+1; 
    fr.length.val = 0.05;
    fr.hop.val = 0.5;
    fl = fr.length.val*sr; % fr.length.val = 0.05
    h = fr.hop.val*fl;     % fr.hop.val = 0.5
    n = floor((lsz-fl)/h)+1;   % Number of frames
    st = floor(((1:n)-1)*h)+w(1);
    fp = [st; floor(st+fl-1)];
    fp(:,fp(2,:)>w(2)) = []; % frame indices

    % Import onsets
    strOnsets = sprintf('%s%d.txt',onsetPath,numSample);
    onsd = importdata(strOnsets); % time instants
    onsd = onsd*sr; % convert to sample indices
    startidx = zeros(1,length(onsd)-1);
    endidx = zeros(1,length(onsd)-1);
    for i = 1:length(onsd)-1
        startidx(i) = find(onsd(i)>fp(1,:) & onsd(i)<fp(2,:), 1 );
        tmp = find(onsd(i+1)>fp(1,:) & onsd(i+1)<fp(2,:), 1 );
        if numel(tmp)==0
            endidx(i) = size(fp,2); % edge case at end of frames
        else
            endidx(i) = tmp;
        end
    end

    % Now, compute segment level features
    strFeat = sprintf('%sFeats%d.mat',featPath,numSample);
    feat = load(strFeat); feat = feat.features';
    
    D = size(feat,2); % dimensionality of features
    F = zeros(length(onsd)-1,D*2);
    for i = 1:length(onsd)-1
        M = feat(startidx(i):endidx(i),:);
        if size(M,1) == 1
            M = [M;M]; % hack for just one frame
        end
        F(i,:) = [mean(M) std(M)];
    end
    F = zscore(F);
    
    % Store the segment level features
    SegmentFeats_cell{numSample} = F;
end

save('SegmentFeatures_AMG1608.mat','SegmentFeats_cell');
