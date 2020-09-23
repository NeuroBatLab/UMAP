function [prepDATA]=prep_audio_umap(audiopath,bin,overlperc,IDs,dtw,biosoundVals,...
    bioflag,biosoundpath,piezoflag)
%function prepares audio data for UMAP application by either (if dtw=0)
%snipping the spectrogram of a vocalization into defined bin windows with a
%defined overlap; each spectrogram snippet corresponds then to one row of the
%resulting matrix; snippets are time scaled to either the lower or higher
%multiple of the defined bin window; or (if dtw=1) directly calculating the
%Euclidean Distance of each spectrogram to each other through dynamic time
%warping
%INPUTS: 1. path for the audio files (.wav)
%        2. bin window size in ms, can be single value or vector
%        3. overlap of the bin windows (.1-.9)
%        4. cell of IDs to extract audios from individuals, if not declare []
%        5. flag to either time scale snippets (=0) or directly calculate
%        Euclidean Distance with dynamic time warping (=1)
%        6. parameters for calculating the spectrograms, when using default
%        values declare []
%        7. flag to grab corresponding biosound parameters (0 or 1)
%        8. path for the prior generated biosound files (.mat), if not needed declare []
%        9. flag if only use piezo wav files and biosound parameters (0 or
%        1)
%OUTPUTS: structure for each defined bin window
%        1. matrix to input into umap
%        2. the number of rows (=snippets) for each vocalization (if dtw=0)
%        3. the IDs (if defined) for each vocalization
%        4. the mean biosound parameters (if declared) for each snippet (dtw=0) or
%        vocalization (dtw=1)
%
% for this function the SoundAnalysisBats folder needs to bespecified as a path
%(can be found on the Neurobatlab github)

idflag=0;
if ~isempty(IDs)
    idflag=1;
end
pathflag=0;
if strcmp(audiopath,biosoundpath)
    pathflag=1;
end
%default biosound values
DBNOISE=60;
f_high=40e3;
fband=100;
time_res=.001;
if ~isempty(biosoundVals)
    DBNOISE = biosoundVals.dbnoise;
    f_high=biosoundVals.fhigh;
    fband = biosoundVals.band;
    time_res = biosoundVals.timeres;
end
%for dynamic time warping only one bin loop
if dtw
    bin=1;
    overlperc=1;
end
ovl=round(overlperc.*bin);
%check for time scaling that overlap steps match to bin or multiple of bin
if dtw==0
    for ovch=1:length(bin)
        correct=1;
        while correct
            if mod(bin(ovch)/ovl(ovch),1)~=0
                ovl(ovch)=ovl(ovch)-1;
            else
                correct=0;
            end
        end
    end
end
overl=bin-ovl;
allpnf=[];
allIDs=[];
for bn=1:length(bin)
    if bn==1
        %grab IDs for each vocalization if flagged and path dir for all
        if idflag==1
            for i=1:length(IDs)
                if piezoflag
                    %grab only piezo recordings
                    pnf=dir([audiopath '\*' IDs{i} '*Piezo.wav']);
                else
                    if pathflag
                        %grab only mic recordings
                        pnf=dir([audiopath '\*' IDs{i} '*Raw.wav']);
                    else
                        %grab mic recordings
                        pnf=dir([audiopath '\*' IDs{i} '*.wav']);
                    end
                end
                allpnf=[allpnf;pnf];
                ID=IDs{i};
                allIDs(length(allIDs)+1:length(allIDs)+length(pnf))=...
                    ones(1,length(length(allIDs)+1:length(allIDs)+length(pnf))).*str2num(ID);
            end
        else
            if piezoflag
                allpnf=dir([audiopath '\*Piezo.wav']);
            else
                if pathflag
                    allpnf=dir([audiopath '\*Raw.wav']);
                else
                    allpnf=dir([audiopath '\*.wav']);
                end
            end
        end
    end
    ALLbios=[];
    if dtw
        EDs=NaN(length(allpnf),length(allpnf));
        rep=[];
    else
        ALL=[];
        POS=[];
    end
    disp(['Grabbing specs for bin ' num2str(bin(bn))])
    for fln=1:length(allpnf)
        [y,fs]=audioread([audiopath '\' allpnf(fln).name]);
        %calculate spectrogram
        if fln==1 && f_high>fs/2
            f_high=floor(fs/2);
        end
        [to, fo, logB, pg, tError, fError] = spec_only_bats(y,fs,DBNOISE,...
            f_high,fband,'Time_increment',time_res);
        if dtw
            %ED calculation between each vocalization through dynamic time
            %warping of the spectrograms
            [eds,allbios]=dyntimewarp(logB,fln,allpnf,audiopath,DBNOISE,f_high,...
                fband,time_res,bioflag,biosoundpath,allpnf(fln).name,piezoflag,pathflag);
            EDs(fln,:)=[rep cell2mat(eds)];
            if fln+1<=length(allpnf)
                feds=NaN(1,length(EDs(1,2:fln+1)));
                rep=[rep(fln:end) feds];
            end
        else
            %grab spectrogram snippets with prior time scaling
            valdef=round(length(y)/fs*1e3/bin(bn));
            if valdef==0
                valdef=1;
            end
            [MAT,allbios]=timescalespec(logB,valdef,valdef*bin(bn),bin(bn),...
                overl(bn),bioflag,biosoundpath,allpnf(fln).name,piezoflag,pathflag);
            POS(fln)=size(MAT,1);
            ALL=[ALL;MAT];
        end
        if bioflag
            ALLbios=[ALLbios;allbios];
        end
        
    end
    if dtw
        for i=1:size(EDs,2)
            for ii=1:size(EDs,1)
                EDs(ii,i)=EDs(i,ii);
            end
        end
        prepDATA.umap{bn}=EDs;
    else
        prepDATA.umap{bn}=ALL;
        prepDATA.Pos{bn}=POS;
    end
    if bn==1
        prepDATA.IDs=allIDs;
        prepDATA.fhigh=f_high;
    end
    prepDATA.bios{bn}=ALLbios;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MAT,allbios]=timescalespec(logB,valdef,tscale,binv,ovv,bioflag,...
    biosoundpath,fname,piezoflag,pathflag)
%grabs snippets of a time-scaled spectrogram and if flagged the
%corresponding biosound parameters
tw=interp1(1:size(logB,2),logB',linspace(1,size(logB,2),tscale))';
bios=[];
if bioflag
    if piezoflag
        load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
    else
        if pathflag
            load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
        else
            if exist([biosoundpath '\' fname '_Raw_biosound.mat'])
                load([biosoundpath '\' fname '_Raw_biosound.mat'])
            else
                load([biosoundpath '\' fname '_biosound.mat'])
            end
        end
    end
    for bs=1:5
        if bs==1
            vec=BioSoundCall.tAmp;
        elseif bs==2
            vec=BioSoundCall.amp;
        elseif bs==3
            vec=BioSoundCall.SpectralMean;
        elseif bs==4
            if piezoflag
                load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
            else
                if pathflag
                    load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
                else
                    if exist([biosoundpath '\' fname '_Piezo_biosound.mat'])
                        load([biosoundpath '\' fname '_Piezo_biosound.mat'])
                    end
                end
            end
            vec=BioSoundCall.sal;
        else
            if piezoflag
                load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
            else
                if pathflag
                    load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
                else
                    if exist([biosoundpath '\' fname '_Piezo_biosound.mat'])
                        load([biosoundpath '\' fname '_Piezo_biosound.mat'])
                    end
                end
            end
            vec=BioSoundCall.f0;
        end
        if length(vec)>5
            bios(bs,:)=interp1(1:length(vec),vec,linspace(1,length(vec),tscale));
        else
            interin=interp1(1:10,1:10,linspace(1,10,tscale));
            bios(bs,:)=NaN(1,length(interin));
        end
    end
end
MAT=[];
allbios=[];
for i=1:valdef
    grab=tw(:,1+(i-1)*ovv:(i-1)*ovv+binv);
    MAT=[MAT;reshape(grab,1,numel(grab))];
    if bioflag
        if (i-1)*ovv+binv<=size(bios,2)
            allbios=[allbios;nanmean(bios(:,1+(i-1)*ovv:(i-1)*ovv+binv)')];
        else
            allbios=[allbios;NaN(bs,1)];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function [eds,bios]=dyntimewarp(X,fln,allpnf,audiopath,DBNOISE,f_high,fband,time_res,...
    bioflag,biosoundpath,fname,piezoflag,pathflag)
%calculates the ED of the spectrogram of each vocalization to all others
%through dynamic time warping and grabs if flagged the mean biosound
%parameters for each vocalization
eds=cell(1,length(fln:length(allpnf)));
parfor cfln=fln:length(allpnf)
    [y,fs]=audioread([audiopath '\' allpnf(cfln).name]);
    [to, fo, logB, pg, tError, fError] = spec_only_bats(y,fs,DBNOISE,...
        f_high,fband,'Time_increment',time_res);
    eds{1,cfln}=dtw(X,logB,5);
end
bios=[];
if bioflag
    if piezoflag
        load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
    else
        if pathflag
            load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
        else
            if exist([biosoundpath '\' fname '_Raw_biosound.mat'])
                load([biosoundpath '\' fname '_Raw_biosound.mat'])
            else
                load([biosoundpath '\' fname '_biosound.mat'])
            end
        end
    end
    for bs=1:5
        if bs==1
            vec=BioSoundCall.tAmp;
        elseif bs==2
            vec=BioSoundCall.amp;
        elseif bs==3
            vec=BioSoundCall.SpectralMean;
        elseif bs==4
            if piezoflag
                load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
            else
                if pathflag
                    load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
                else
                    if exist([biosoundpath '\' fname '_Piezo_biosound.mat'])
                        load([biosoundpath '\' fname '_Piezo_biosound.mat'])
                    end
                end
            end
            vec=BioSoundCall.sal;
        else
            if piezoflag
                load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
            else
                if pathflag
                    load([biosoundpath '\' fname(1:end-4) '_biosound.mat'])
                else
                    if exist([biosoundpath '\' fname '_Piezo_biosound.mat'])
                        load([biosoundpath '\' fname '_Piezo_biosound.mat'])
                    end
                end
            end
            vec=BioSoundCall.f0;
        end
        bios(bs)=nanmean(vec);
    end
end