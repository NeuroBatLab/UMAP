clear all;
close all;
%set here desired paths for wav files and prior calculated biosound
%parameters
audiopath='C:\Users\BatLab\Desktop\Juliedata';
biosoundpath='C:\Users\BatLab\Desktop\Juliedata';
%example bin values
bin=25:25:200;
%example overlap value, (should be between .1 and .9)
overlperc=0.3;
%flag gets biosound parameters for each vocalization or snippet
bioflag=1;
%enter here other biosound settings for generating spectrograms, for default see
%prep_audio_umap function
biosoundVals=[];
%flag only gets piezo recordings and piezo biosound parameters
piezoflag=1;
%set dynamic time warping (=1) or time scaling
dtw=0;
if dtw
    dtwstr='dtw';
    IDs{1}={'11648';'14461';'14464';'65696';'71353';
        '14463';'71043';'71047';'71351';'71354'};
else
    dtwstr='tsc';
    %saline
    IDs{1}={'11648';'14461';'14464';'65696';'71353'};
    %deaf
    IDs{2}={'14463';'71043';'71047';'71351';'71354'};
end
%prepares data for all saline and deaf bats
for run=1:length(IDs)
    %run data prep
    [prepDATA]=prep_audio_umap(audiopath,bin,overlperc,IDs{run},dtw,...
        biosoundVals,bioflag,biosoundpath,piezoflag);
    if run==2
        prepDATAd=prepDATA;
    else
        prepDATAs=prepDATA;
    end
end
%example for saving the data
if dtw
    save([pwd '\umap_data\umap_data_saldeaf_spec_' dtwstr '.mat'],...
        'prepDATA','-v7.3')
else    
    save([pwd '\umap_data\Julie_umap_data_saldeaf_spec_' dtwstr '.mat'],...
    'prepDATAd','prepDATAs','bin','overlperc','-v7.3')
end