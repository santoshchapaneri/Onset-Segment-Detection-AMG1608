% Run OnsetDetector_CLR for all wav files

% Directory of the files
d = 'D:\MusicPhD\2017\Sep2017\OnsetDetectCLR\OnsetDetectorCLR-master\AMG1608Wav\';
% path(path, [pwd, filesep, d]); % add to path
names = dir(fullfile(d,'*.wav') );
names = {names(~[names.isdir]).name};
[namesSorted,idx] = sort_nat(names);

for n = 1:numel(names)
    fprintf('%d\n',n);
    str = sprintf('python onsetDetectorCLR.py %s%d.wav',d,n);
    system(str);
end
