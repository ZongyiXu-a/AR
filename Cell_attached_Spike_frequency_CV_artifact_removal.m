clear all; clc

fid=dir('*.mat');%filename
load('dt.mat');

%%-----Adjust if needed----%%
freq=130;
num=1300;
baserange=10;
postime=20;
binwidth=1;
stimlatency=0.003; % remove artifact
pick=[6,11,16,21,26];
xlstitle='Result_APtime_1301300';%filename
%---------------------%%
Recordfilename=[];
Record_basepeakFreq=zeros(length(fid),1);
Record_intrapeakFreq=zeros(length(fid),1);
Record_postpeakFreq=zeros(length(fid),1);
Record_instantFreq=[];
Record_baseCV=zeros(length(fid),1);
Record_intraCV=zeros(length(fid),1);
Record_postCV=zeros(length(fid),1);
Record_binCV=[];
Record_binmeanAP=[];
Record_binmean=zeros(length(fid),30);
Record_binCVCV=zeros(length(fid),30);
Record_CVCV=zeros(length(fid),5);

for k=1:length(fid)
    filename=fid(k).name;
    xlsname=[filename(1:end-4) '.xls'];
    load(filename);
    sti=SMarker;
    baseline=data_act(1:fix(baserange/dtI));
    basemean=mean(baseline);

    %% Remove artifact%%
    artifactidx=zeros(length(sti),fix(stimlatency/dtI));
    artifactidx(:,1)=fix(sti/dtI);
    for a=2:size(artifactidx,2)
        artifactidx(:,a)=artifactidx(:,a-1)+1;
    end
    current_removeA=data_act;
    current_removeA(artifactidx)=basemean;  % Replace artifact with baseline
    time=(1:length(current_removeA))*dtI;   % Time vector

    %% Peak detection on the artifact-removed signal%%
    peaktime=peakfinder(current_removeA, -40, -20, -1, 0) * dtV; % peakfinder on cleaned data
    peaktimereal=diff(peaktime) > 0.0003;
    peakpeak=peaktime(peaktimereal);
    peaktimeplot_idx=fix(peakpeak / dtV);  % Find peak indices

    %%Split the peaks into different categories%%
    basepeakidx=find(peakpeak < baserange);
    basepeak=peakpeak(basepeakidx);
    postpeakidx=find(peakpeak > baserange + 1 / freq * num);
    postpeak=peakpeak(postpeakidx);
    intrapeakidx=find(peakpeak >= baserange & peakpeak < (baserange + 1 / freq * num));
    intrapeak=peakpeak(intrapeakidx);

    %% Calculate instant frequency%%
    instantFreq=[];
    for i=2:length(peakpeak)
        ISI=peakpeak(i) - peakpeak(i - 1);   
        instantfreq=1 / ISI;    
        instantFreq=[instantFreq; instantfreq];
    end
    baseCVidx=instantFreq(1:(length(basepeak) - 1));
    intraCVidx=instantFreq(((length(basepeak) + 1):(length(basepeak) + length(intrapeak) - 1)));
    postCVidx=instantFreq(((length(basepeak) + length(intrapeak) + 1):end));

    %% Bin-based jitter calculation (1 second bin)%%
    bininstantFreq=[];
    for i=1:(baserange + postime)
        binAPidx=find(peakpeak < i & peakpeak > (i - 1));
        if isempty(binAPidx)
            binAPCV=0; 
            binAP=0;
        else
            binAP=peakpeak(binAPidx);
            for m=2:length(binAP)
                binISI=binAP(m) - binAP(m - 1);
                bininstantfreq=1 / binISI;
                bininstantFreq=[bininstantFreq; bininstantfreq];   
            end
            binAPCV=std(bininstantFreq);
        end
        Record_binmeanAP(i, 1)=sum(binAP ~= 0);
        Record_binCV(i, 1)=binAPCV;
        bininstantFreq=[];
    end

    %% Frequency computation for specified bins%%
    binFreq=[];
    for i=1:length(pick)
        APidx=find(peakpeak < (pick(i) + 4) & peakpeak > pick(i));
        if isempty(APidx)
            APCV=0; 
        else
            AP=peakpeak(APidx);
            for m=2:length(AP)
                ISI=AP(m) - AP(m - 1);
                binfreq=1 / ISI;
                binFreq=[binFreq; binfreq];   
            end
            APCV=std(binFreq);
        end
        Record_CV(i, 1)=APCV;
        binFreq=[];
    end

    %% Plotting the results%%
    figure(1), clf;
    set(gcf, 'position', [200, 100, 800, 400]);

    time_V=(1:length(data_actV)) * dtV;  % Time vector for data_actV
    time_stim=(1:length(data_act)) * dtI;  % Time vector for original data

    % Instant Frequency and CV Plot%%
    h1=subplot(311);
    title(filename); axis tight; hold on;
    yyaxis left
    plot(peakpeak(1:end-1), instantFreq, '.', 'markerfacecolor', [60, 176, 106] / 255);
    ylabel('instantFreq (Hz)')
    yyaxis right
    plot(0.5:1:30, Record_binCV);
    ylabel('instantFreq CV/s')

    % Data signal plot with detected peaks%%
    h2=subplot(312);
    axis tight; hold on;
    plot(time_V, current_removeA, 'color', [220, 20, 60] / 255);  % Plot artifact-removed signal
    plot(time_V(peaktimeplot_idx), current_removeA(peaktimeplot_idx), 'o', 'markerfacecolor', [60, 176, 106] / 255);

    % Stimulus plot%%
    h3=subplot(313);
    axis tight; hold on;
    plot(time_stim, data_actV, 'color', [220, 20, 60] / 255);

    pause;

    %% Frequency statistics%%
    basepeaknum=length(basepeak);
    basepeakFreq=basepeaknum / baserange;
    postpeaknum=length(postpeak);
    postpeakFreq=postpeaknum / (postime - (1 / freq * num));
    intrapeaknum=length(intrapeak);
    intrapeakFreq=intrapeaknum / (1 / freq * num);
    baseCV=std(baseCVidx)/mean(baseCVidx);
    intraCV=std(intraCVidx)/mean(intraCVidx);
    postCV=std(intraCVidx)/mean(intraCVidx);

    %% Storing results for output%%
    Recordfilename{k, 1}=filename;
    Record_basepeakFreq(k, 1)=basepeakFreq;
    Record_intrapeakFreq(k, 1)=intrapeakFreq;
    Record_postpeakFreq(k, 1)=postpeakFreq;
    Record_instantFreq(k, 1:length(instantFreq))=instantFreq;
    Record_baseCV(k, 1)=baseCV;
    Record_intraCV(k, 1)=intraCV;
    Record_postCV(k, 1)=postCV;
    Record_binmean(k, :)=Record_binmeanAP;
    Record_binCVCV(k, :)=Record_binCV;
    Record_CVCV(k, :)=Record_CV;
end

%% Write results to Excel%%
xlswrite(xlstitle, {'Name', 'basepeakFreq', 'intrapeakFreq', 'postpeakFreq', 'baseCV', 'intraCV', 'postCV', 'binCV'}, 'para', 'A1');
xlswrite(xlstitle, Recordfilename, 'para', 'A2');
xlswrite(xlstitle, Record_basepeakFreq, 'para', 'B2');
xlswrite(xlstitle, Record_intrapeakFreq, 'para', 'C2');
xlswrite(xlstitle, Record_postpeakFreq, 'para', 'D2');
xlswrite(xlstitle, Record_baseCV, 'para', 'E2');
xlswrite(xlstitle, Record_intraCV, 'para', 'F2');
xlswrite(xlstitle, Record_postCV, 'para', 'G2');
xlswrite(xlstitle, Record_binCVCV, 'para', 'H2');
xlswrite(xlstitle, {'Name', 'binmeanAP'}, 'binmean', 'A1');
xlswrite(xlstitle, Recordfilename, 'binmean', 'A2');
xlswrite(xlstitle, Record_binmean, 'binmean', 'B2');
xlswrite(xlstitle, {'Name', 'binCV4'}, 'binCV4', 'A1');
xlswrite(xlstitle, Recordfilename, 'binCV4', 'A2');
xlswrite(xlstitle, Record_CVCV, 'binCV4', 'B2');
