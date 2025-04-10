% 2024-01-03 modified by DSX,specific for SR alpha fitting
clear all;clc;close all;
load('dt.mat');
fid=dir('*.mat');%filename

%% ------------Adjust if needed----------------
pretime=20;
durtion=30;
stimlatency=0.0011;%remove artifact
SRampthre=5;%SR>5pA
slopdownsample=6;%slop
onsetwindow=2.5;%SR onset, ms
risingwindow=4;%SR peak latency window after stimulus onset: 4–5 ms.;
fitwindow=3.5;%Search range from SR onset (ms);130hz 3.5；67hz 6-8；20hz 8-10
fitonset=1.05;%130 hz 1.1;20 hz 1.05
taumax=5;
c=30;%plot
paraname={'stimlatency','SRampthre','slopdownsample','onsetwindow','risingwindow','fitwindow','fitonset','taumax'};
para=[stimlatency,SRampthre,slopdownsample,onsetwindow,risingwindow,fitwindow,fitonset,taumax];
%%
for k=1:length(fid)
    filename=fid(k).name;
    load(filename);
    datatitle=['DATA_',filename(1:end-4),'_alphaFit_2024'];%filename
    
    sti=SMarker;
    stio=find(sti<=sti(1)+1);
    freq=length(stio);
    num=length(sti);
    plotnum=100;%Number of subplots
    stiwindow=1/freq;%ISI
    
    baseline=data_act(1:fix(pretime/dtI));
    basemean=mean(baseline);
    
    %% remove artifact
    artifactidx=zeros(length(sti),fix(stimlatency/dtI));
    artifactidx(:,1)=fix(sti/dtI);
    for a=2:size(artifactidx,2)
        artifactidx(:,a)=artifactidx(:,a-1)+1;
    end
    current_removeA=data_act;
    current_removeA(artifactidx)=basemean;
    time=(1:length(current_removeA))*dtI;
    
    %% SR_alpha synapse
    for j=1:length(sti)
        chargeduridx=find(time>=sti(j)&time<sti(j)+stiwindow);
        chargedur=current_removeA(chargeduridx);
        chargeoriginal=data_act(chargeduridx);
        time_plot=time(chargeduridx);
        
        tspanI=(time_plot-time_plot(1)).*1000;%ms
        onsetidx=find(chargedur==basemean);
        datasearch=chargedur;
        A=mean(chargeoriginal(1:2));%Baseline (pre-stimulus)
        B=mean(chargeoriginal(length(onsetidx)-1:length(onsetidx)));%Baseline (post-stimulus)
        datasearch(onsetidx)=A;
        datasearch=datasearch-datasearch(1);
        
        sloptimeidx=1:slopdownsample:length(datasearch);
        slopsearch=datasearch(sloptimeidx);
        sloptime=tspanI(sloptimeidx);
        dataslop=zeros(length(slopsearch),1);
        dataslop(2:end)=diff(slopsearch)'./diff(sloptime);
        dataslop(2:end)=diff(dataslop)'./diff(sloptime);
        sloppeakidx=peakfinder(dataslop,5,-100,-1,0);%(x0, amp, base-thresh, extrema, include_endpoints)
        sloppeakidxx=sloptimeidx(sloppeakidx);
        slopC(j,1:length(dataslop))=dataslop;
        sloptimeC(j,1:length(sloptime))=sloptime;
        
        Record_tspanI(j,1:length(tspanI))=tspanI;
        Record_datasearch(j,1:length(datasearch))=datasearch;
        Record_dataoriginal(j,1:length(chargeoriginal))=chargeoriginal;
        
        if isempty(sloppeakidxx)
            Record_peakmin(j)=0;
            Record_peakindex(j)=0;
        else
            firstslop=tspanI(sloppeakidxx(1));
            if firstslop>onsetwindow
                Record_peakmin(j)=0;
                Record_peakindex(j)=0;
            else
                peaksearchidxa=peakfinder(datasearch,SRampthre,5,-1,0);%(x0, amp, base-thresh, extrema, include_endpoints)
                realpeak=find(tspanI(peaksearchidxa)>fitonset+0.1);
                if isempty(realpeak)
                    Record_peakmin(j)=0;
                    Record_peakindex(j)=0;
                elseif tspanI(peaksearchidxa(realpeak(1)))>risingwindow
                    Record_peakmin(j)=0;
                    Record_peakindex(j)=0;            
                else                    
                    Record_peakmin(j)=datasearch(peaksearchidxa(realpeak(1)));
                    Record_peakindex(j)=tspanI(peaksearchidxa(realpeak(1)));
                    
                    FitOnRange=(fitonset:dtI*1000:Record_peakindex(j));% time window of fot onset
                    GOF=zeros(1,length(FitOnRange));
                    Asum=zeros(1,length(FitOnRange));
                   
                    for i=1:length(FitOnRange)
                        FitOn=FitOnRange(i);
                        fitAut_direct=fittype('a*(x-t0)/tau*exp(-(x-t0-tau)/tau)*heaviside(x-t0)',...
                            'coefficients',{'a','t0','tau'});
                        f_direct=fit(tspanI(round(FitOn/dtI/1000):round((FitOn+fitwindow)/dtI/1000))',datasearch(round(FitOn/dtI/1000):round((FitOn+fitwindow)/dtI/1000)),fitAut_direct,'Startpoint',[Record_peakmin(j) firstslop 1.25],'lower',[-2500 fitonset 0.2],'upper',[0 Record_peakindex(j) taumax]);
                        cha=f_direct(tspanI(round(FitOn/dtI/1000):round((FitOn+fitwindow)/dtI/1000)))-datasearch(round(FitOn/dtI/1000):round((FitOn+fitwindow)/dtI/1000));%拟合线与实际值相减
                        GOF(i)=sum(abs(cha))/length(cha);
                        Asum(i)=trapz(tspanI(round(FitOn/dtI/1000):round((FitOn+fitwindow)/dtI/1000)),f_direct(tspanI(round(FitOn/dtI/1000):round((FitOn+fitwindow)/dtI/1000))));%Rui:GOF:pA 绝对值之和除以总点数
                    end
                    temp=[];
                    tempidx=find(FitOnRange<=2);
                    temp(:,1)=GOF(tempidx);
                    temp(:,2)=Asum(tempidx);
                    
%                     [temp,idx]=sortrows(temp,[2,1],'ascend');
                    square=(temp(:,1)-min(temp(:,1)))+(temp(:,2)-min(temp(:,2)));
%                     [GOFmin,GOFindex]=min(temp(:,1));
                    [GOFmin,GOFindex]=min(square);
                    minGOFtime=FitOnRange(GOFindex);
                    
                    figure(100),clf;                    
                    h1=subplot(121);title([num2str(k),' - ',num2str(j)]);hold on;
                    plot(FitOnRange,GOF,'ro');hold on
                    plot(FitOnRange,Asum,'bo');
                    plot([FitOnRange(GOFindex);FitOnRange(GOFindex)],[max(temp(:));min(temp(:))],'k')
                    h2=subplot(122);hold on;
                    plot(GOF,Asum,'ro');hold on
                    plot(GOF(GOFindex),Asum(GOFindex),'k*')
                    
                    fitAut_direct=fittype('a*(x-t0)/tau*exp(-(x-t0-tau)/tau)*heaviside(x-t0)',...
                        'coefficients',{'a','t0','tau'});
                    f_direct=fit(tspanI(round(minGOFtime/dtI/1000):round((minGOFtime+fitwindow)/dtI/1000))',datasearch(round(minGOFtime/dtI/1000):round((minGOFtime+fitwindow)/dtI/1000)),fitAut_direct,'Startpoint',[Record_peakmin(j) minGOFtime 1.25],'lower',[-2500 fitonset 0.2],'upper',[0 risingwindow taumax]);  %fit的时间范围，autapse反应范围，用的fit方法。最终返回三个值：幅度，起始时间，tau
%                     f_direct=fit(tspanI(round(minGOFtime/dtI/1000):end)',datasearch(round(minGOFtime/dtI/1000):end),fitAut_direct,'Startpoint',[Record_peakmin(j) minGOFtime 1.25],'lower',[-700 stimlatency*1000 0.2],'upper',[0 Record_peakindex(j) 4]);  %fit的时间范围，autapse反应范围，用的fit方法。最终返回三个值：幅度，起始时间，tau
                    Fitted=f_direct.a*(tspanI-f_direct.t0)/f_direct.tau.*exp(-(tspanI-f_direct.t0-f_direct.tau)/f_direct.tau).*heaviside(tspanI-f_direct.t0);
                    cha=Fitted'-datasearch;%Subtract the fitted line from the actual values
                   
                    a=f_direct.a;
                    t0=f_direct.t0;
                    tau=f_direct.tau;
                    
                    Record_Fitted(j,1:length(Fitted))=Fitted;
                    Record_cha(j,1:length(cha))=cha;
                    Record_f_directA(j,1)=a;
                    Record_f_directT(j,1)=t0;
                    Record_f_directTau(j,1)=tau;
                    Record_datasearch(j,1:length(datasearch))=datasearch;
                    
                    figure(200),clf
                    h1=subplot(311);title([num2str(k),' - ',num2str(j)]);hold on
                    plot(Record_tspanI(j,:),Record_datasearch(j,:),'color',[60,176,106]/255,'LineWidth',1);hold on
                    h2=subplot(312);title(['Slop']);hold on
                    plot(sloptimeC(j,:),slopC(j,:),'color',[255,106,106]/255,'LineWidth',1);hold on
                    h3=subplot(313);title(['Alpha fitting']);hold on
                    plot(Record_tspanI(j,:),Record_datasearch(j,:),'color',[20,170,220]/255,'LineWidth',1);hold on
                    plot(Record_tspanI(j,:),Record_Fitted(j,:),'g','LineWidth',2);hold on
                    plot(Record_tspanI(j,:),Record_cha(j,:),'r-','LineWidth',1);hold on
                    linkaxes([h3,h2,h1],'x');
                end
            end
        end
        Record_peakmin(abs(Record_peakmin)<=SRampthre)=0;
        SRfackidx=find(Record_peakmin==0);
        Record_peakindex(SRfackidx)=0;
        Record_f_directA(SRfackidx)=0;
        Record_f_directT(SRfackidx)=0;
        Record_f_directTau(SRfackidx)=0;
        Record_Fitted(SRfackidx,:)=0;
        Record_cha(SRfackidx,:)=0;
    end
    
    %% save to matlab
    save([datatitle,'.mat'],'Record_f_directA','Record_f_directT','Record_f_directTau','Record_Fitted','Record_cha','Record_datasearch','Record_tspanI','Record_dataoriginal','Record_peakmin','Record_peakindex','sti','paraname','para')
   

        figure(c),clf
        iset=1;plotnum=100;k=1;
        for i=plotnum*(c-1)+1:plotnum*c
            subplot(10,10,iset);title([num2str(k),' - ',num2str(i)]);hold on;
            iset=iset+1;
            plot(Record_tspanI(i,:),Record_datasearch(i,:),'color',[20,170,220]/255,'LineWidth',1);hold on
            plot(Record_tspanI(i,:),Record_Fitted(i,:),'g','LineWidth',2);hold on
            plot(Record_tspanI(i,:),Record_cha(i,:),'r-','LineWidth',1);hold on
        end

    
    clear Record_f_directA Record_f_directT Record_f_directTau Record_Fitted Record_cha Record_datasearch Record_tspanI Record_dataoriginal Record_peakmin Record_peakindex sti paraname para;
end

