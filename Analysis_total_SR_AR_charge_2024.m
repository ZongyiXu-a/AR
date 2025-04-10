clear all;clc;close all;
 
fid=dir('*.mat');%filename
load('dt.mat');
load('*.mat');%filename
%% ------------Adjust if needed---------------%%
freq=130;
num=3900;
xlstitle='Result_charge_alphaall_2024_1303900';%filename
datatitle='DATA_charge_alphaall_2024_1303900';%filename
choose=[1,100,650,1300,1950,2600,3050,3900];%trail number of stimuli for polt
pretime=20;
durtime=30;
posttime=20;
bin10bin=10;
binPsec=1;
burstidx=(pretime/binPsec+1):1:(pretime+durtime)/binPsec;
stimonset=pretime/bin10bin+1;%Adjust according to pretime and durtime

%% ------------Adjust if needed----------------%%
timewindow=0.02;%Alpha fit window (s)
stiwindow=1/freq;%ISI
stimlatency=0.0011;%remove artifact
bintime=0.001;%The charge area is calculated every 1 ms, save to DATA*.mat

Record_bin10Sum=zeros(length(fid),(pretime+durtime+posttime)/bin10bin);
Record_Charge_current=zeros(length(fid),(pretime+durtime+posttime)/bintime);
Record_binchargesum=zeros(length(fid),(pretime+durtime+posttime)/binPsec);
Record_ChargeBasesum=zeros(length(fid),1);
Record_ChargeSRsum=zeros(length(fid),1);
Record_ChargeINAR=zeros(length(fid),1);
Record_ChargeINsum=zeros(length(fid),1);
Record_ChargePostsum=zeros(length(fid),1);
Record_IntraSR_Ratio=zeros(length(fid),1);
Record_IntraAR_Ratio=zeros(length(fid),1);
Record_IntraARtoSR_Ratio=zeros(length(fid),1);
Record_PostAR20tobase20=zeros(length(fid),1);
Record_PostAR10tobase10=zeros(length(fid),1);
Record_Charge_SR=zeros(length(fid),num);
Record_current_removeA=zeros(length(fid),fix((pretime+durtime+posttime)/dtI)+1);
Record_current_replaceA=zeros(length(fid),fix((pretime+durtime+posttime)/dtI));
Recorde_IntraARtoSR_Ratio_bin1s=zeros(length(fid),num/freq);
Recorde_SR_bin1s=zeros(length(fid),num/freq);
Record_Alphareal=zeros(1,num);
Recordfilename=[];
time_plotC=zeros(length(choose),fix(timewindow/dtI)+1);
chargeoriginalC=zeros(length(choose),fix(timewindow/dtI)+1);
chargedurC=zeros(length(choose),fix(timewindow/dtI)+1);
FitteddurC=zeros(length(choose),fix(timewindow/dtI)+1);
chadurC=zeros(length(choose),fix(timewindow/dtI)+1);
%%
for k=1:length(fid)
    filename=fid(k).name;
    xlsname=[filename(1:end-4) '.xlsx'];
    load(filename);
    sti=SMarker;
    
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
    current_replaceA=data_act;
    current_replaceA(artifactidx)=[];
    time=(1:length(current_removeA))*dtI;
    
    %% based on bintime,total charge analysis
    Y_base=basemean*ones(length(current_removeA),1);
    Charge_current=zeros((pretime+durtime+posttime)/bintime,1);
    for i=1:(pretime+durtime+posttime)/bintime
        if i==1
            Y_base_choose=Y_base(1:fix(i*bintime/dtI));
            current_choose=current_removeA(1:fix(i*bintime/dtI));
            time_choose=time(1:fix(i*bintime/dtI));
            Charge_base=trapz(time_choose,Y_base_choose);
            Charge_c=trapz(time_choose,current_choose);
            Charge_current(i)=abs(Charge_base-Charge_c);
        else
            Y_base_choose=Y_base(fix((i-1)*bintime/dtI):fix(i*bintime/dtI));
            current_choose=current_removeA(fix((i-1)*bintime/dtI):fix(i*bintime/dtI));
            time_choose=time(fix((i-1)*bintime/dtI):fix(i*bintime/dtI));
            Charge_base=trapz(time_choose,Y_base_choose);
            Charge_c=trapz(time_choose,current_choose);
            Charge_current(i)=abs(Charge_base-Charge_c);
        end
    end
    bincharge=reshape(Charge_current,binPsec/bintime,[]);
    binchargemean=mean(bincharge);
    binchargesum=sum(bincharge);
    bin10s=reshape(binchargesum,bin10bin/binPsec,[]);
    bin10mean=mean(bin10s);
    bin10Sum=sum(bin10s);    
    xscale=0.5:binPsec:length(binchargemean);
    
    %% Intra-train basal calculation
    for n=1:length(sti)
        chargeidx=find(time>=sti(n)&time<sti(n)+stiwindow);
        charge=data_act(chargeidx);
        basalcurrent=mean(charge(1:2));
        Record_basal(n)=basalcurrent;
    end    
    Record_basal(Record_basal>0)=0;
    Basal_current=reshape(Record_basal,freq*binPsec,[]);
    Basal_current1s=mean(Basal_current);
    
    figure(1),clf,title(filename);hold on;
    plot(xscale,binchargemean,'color',[20,170,220]/255,'LineWidth',1);
    plot(sti,max(binchargemean)*ones(size(sti)),'b.');
    href_base=refline(0,basemean);
    xlabel('Time(s)');
    ylabel('Charge/1 s bin size (pA*s)');      
 
    figure(2),clf
    subplot(211);title(filename);hold on;
    plot(1:num,Record_peakmin,'color',[20,170,220]/255,'LineWidth',1);hold on;
    plot(1:num,Record_peakindex*100,'color',[255,106,106]/255,'LineWidth',1);hold on;
    plot(1:num,Record_f_directA,'color',[60,176,106]/255,'LineWidth',1);hold on;
    legend('SR peak amp','SR peak time','location','best')
    subplot(212);hold on;
    plot(1:num,Record_basal,'color',[60,176,106]/255,'LineWidth',1);hold on;
    legend('Basal current','location','best')
    ylabel('pA');
    xlabel('Stim num');
    pause(0.1);   
    
    
    %% SR alpha synapse fitting and calculation of charge
    SRrealidx=find(Record_peakmin~=0);
    for i =1:length(SRrealidx)
        chargeduridx=find(time>=sti(SRrealidx(i))&time<sti(SRrealidx(i))+timewindow);
        chargedur=current_removeA(chargeduridx);
        chargeoriginal=data_act(chargeduridx);
        time_plot=time(chargeduridx);
        
        tspanI=(time_plot-time_plot(1)).*1000;%ms
        onsetidx=find(chargedur==basemean);
        datasearch=chargedur;
        datasearch(onsetidx)=mean(chargeoriginal(1:2));
        datasearch=datasearch-datasearch(1);       
               
        A=[Record_peakmin(SRrealidx(i)),Record_f_directA(SRrealidx(i))];
        a=-max(abs(A));
        t0=Record_f_directT(SRrealidx(i));
        tau=Record_f_directTau(SRrealidx(i));
        
        Fitted=a*(tspanI-t0)/tau.*exp(-(tspanI-t0-tau)/tau).*heaviside(tspanI-t0);
        cha=Fitted'-datasearch;
        windownum=fix(timewindow/stiwindow)+1;
        windowstimtimeidx=[];
        for b=1:windownum           
            windowstimtimeidx=[windowstimtimeidx;(SRrealidx(i)+b-1)];
        end
        
        if SRrealidx(i)>=(num-windownum+2)
            windowstimtime=zeros(windownum,1);
            for s=2:length(windowstimtime)
                windowstimtime(s)=stiwindow*(s-1);
            end
        else
            windowstimtime=(sti(windowstimtimeidx)-sti(windowstimtimeidx(1)))*1000;
        end
        
        Record_Alphaneg=[];        
        for m=1:windownum
            negFitteridx=find(tspanI>=windowstimtime(m)&tspanI<=(windowstimtime(m)+stimlatency*1000));
            Alphaneg=trapz(tspanI(negFitteridx),Fitted(negFitteridx));%pA.ms
            Record_Alphaneg=[Record_Alphaneg;Alphaneg];
        end
        
        Alphatotal=trapz(tspanI,Fitted);%pA.ms
        Alphareal=abs(Alphatotal)-abs(sum(Record_Alphaneg));
        Record_Alphareal(SRrealidx(i))=Alphareal;
       
        for o=1:length(choose)
            if SRrealidx(i)==choose(o)
                time_plotC(o,1:length(tspanI))=tspanI;
                chargeoriginalC(o,1:length(chargeoriginal))=chargeoriginal;
                chargedurC(o,1:length(datasearch))=datasearch;
                FitteddurC(o,1:length(Fitted))=Fitted;
                chadurC(o,1:length(cha))=cha;      
                
            end
        end
    end
    
%% choose stim plot
oidx=find(time_plotC(:,size(time_plotC,2)-1)==0);
nego=choose(oidx);

for n=1:length(nego)
    chargeduridx=find(time>=sti(nego(n))&time<(sti(nego(n))+timewindow));
    time_plot=time(chargeduridx);
    tspanI=(time_plot-time_plot(1)).*1000;%ms
    chargeoriginal=data_act(chargeduridx);
    chargedur=current_removeA(chargeduridx);
    onsetidx=find(chargedur<basemean);
    datasearch=chargedur-chargedur(onsetidx(1));
    datasearch(1:onsetidx(1)-1)=0;
    chargeoriginalC(oidx(n),1:length(chargeoriginal))=chargeoriginal;
    chargedurC(oidx(n),1:length(datasearch))=datasearch;
    time_plotC(oidx(n),1:length(time_plot))=tspanI; 
    
end

figure(3),clf
h1=subplot(241);title([num2str(choose(1)),'--',filename]);hold on;
plot(time_plotC(1,1:fix(timewindow/dtI)),chargeoriginalC(1,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h2=subplot(245);hold on;
plot(time_plotC(1,1:fix(timewindow/dtI)),chargedurC(1,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(1,1:fix(timewindow/dtI)),FitteddurC(1,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(1,1:fix(timewindow/dtI)),chadurC(1,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h2,h1],'x');

h3=subplot(242);title([num2str(choose(2)),'--',filename]);hold on;
plot(time_plotC(2,1:fix(timewindow/dtI)),chargeoriginalC(2,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h4=subplot(246);hold on;
plot(time_plotC(2,1:fix(timewindow/dtI)),chargedurC(2,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(2,1:fix(timewindow/dtI)),FitteddurC(2,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(2,1:fix(timewindow/dtI)),chadurC(2,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h4,h3],'x');

h5=subplot(243);title([num2str(choose(3)),'--',filename]);hold on;
plot(time_plotC(3,1:fix(timewindow/dtI)),chargeoriginalC(3,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h6=subplot(247);hold on;
plot(time_plotC(3,1:fix(timewindow/dtI)),chargedurC(3,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(3,1:fix(timewindow/dtI)),FitteddurC(3,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(3,1:fix(timewindow/dtI)),chadurC(3,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h6,h5],'x');

h7=subplot(244);title([num2str(choose(4)),'--',filename]);hold on;
plot(time_plotC(4,1:fix(timewindow/dtI)),chargeoriginalC(4,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h8=subplot(248);hold on;
plot(time_plotC(4,1:fix(timewindow/dtI)),chargedurC(4,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(4,1:fix(timewindow/dtI)),FitteddurC(4,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(4,1:fix(timewindow/dtI)),chadurC(4,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h8,h7],'x');
legend('Trace','Alpha-fitted','Delta','Baseline','location','best')

figure(4),clf
h1=subplot(241);title([num2str(choose(5)),'--',filename]);hold on;
plot(time_plotC(5,1:fix(timewindow/dtI)),chargeoriginalC(5,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h2=subplot(245);hold on;
plot(time_plotC(5,1:fix(timewindow/dtI)),chargedurC(5,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(5,1:fix(timewindow/dtI)),FitteddurC(5,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(5,1:fix(timewindow/dtI)),chadurC(5,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h2,h1],'x');

h3=subplot(242);title([num2str(choose(6)),'--',filename]);hold on;
plot(time_plotC(6,1:fix(timewindow/dtI)),chargeoriginalC(6,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h4=subplot(246);hold on;
plot(time_plotC(6,1:fix(timewindow/dtI)),chargedurC(6,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(6,1:fix(timewindow/dtI)),FitteddurC(6,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(6,1:fix(timewindow/dtI)),chadurC(6,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);

linkaxes([h4,h3],'x');

h5=subplot(243);title([num2str(choose(7)),'--',filename]);hold on;
plot(time_plotC(7,1:fix(timewindow/dtI)),chargeoriginalC(7,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h6=subplot(247);hold on;
plot(time_plotC(7,1:fix(timewindow/dtI)),chargedurC(7,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(7,1:fix(timewindow/dtI)),FitteddurC(7,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(7,1:fix(timewindow/dtI)),chadurC(7,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h6,h5],'x');

h7=subplot(244);title([num2str(choose(8)),'--',filename]);hold on;
plot(time_plotC(8,1:fix(timewindow/dtI)),chargeoriginalC(8,1:fix(timewindow/dtI)),'color',[20,170,220]/255,'LineWidth',1);
h8=subplot(248);hold on;
plot(time_plotC(8,1:fix(timewindow/dtI)),chargedurC(8,1:fix(timewindow/dtI)),'color',[255,106,106]/255,'LineWidth',1);
plot(time_plotC(8,1:fix(timewindow/dtI)),FitteddurC(8,1:fix(timewindow/dtI)),'g','LineWidth',2);
plot(time_plotC(8,1:fix(timewindow/dtI)),chadurC(8,1:fix(timewindow/dtI)),'b','LineWidth',1);
href_base=refline(0,basemean);
linkaxes([h8,h7],'x');
legend('Trace','Alpha-fitted','Delta','Baseline','location','best')
% pause();

%% output
ChargeSRsum=sum(Record_Alphareal)/1000;%Unit: pA*s
Chargedursum=sum(bin10Sum(stimonset:stimonset+durtime/bin10bin-1));
ChargeINAR=Chargedursum-ChargeSRsum;%Unit: pA*s
IntraSR_Ratio=ChargeSRsum/Chargedursum;
IntraAR_Ratio=ChargeINAR/Chargedursum;
IntraARtoSR_Ratio=ChargeINAR/ChargeSRsum;
ChargeBasesum=sum(bin10Sum(1:pretime/bin10bin));
ChargePostsum=sum(bin10Sum(stimonset+durtime/bin10bin:(pretime+durtime+pretime)/bin10bin));
PostAR20tobase20=ChargePostsum/ChargeBasesum;
PostAR10tobase10=bin10Sum(stimonset+durtime/bin10bin)/bin10Sum(pretime/bin10bin);
ChargeSRbinsum=sum(reshape(Record_Alphareal,freq*binPsec,[]))./1000;%Unit: pA*s
IntraARtoSR_Ratio_bin1s=(binchargesum(:,burstidx)-ChargeSRbinsum)./ChargeSRbinsum;
IntraAR_Ratio_bin1s=(binchargesum(:,burstidx)-ChargeSRbinsum)./binchargesum(:,burstidx);

Record_Charge_SR(k,1:length(Record_Alphareal))=Record_Alphareal;
Record_bin10Sum(k,1:length(bin10Sum))=bin10Sum;
Record_Charge_current(k,1:length(Charge_current))=Charge_current;
Record_binchargesum(k,1:length(binchargesum))=binchargesum;
Record_ChargeBasesum(k,1)=bin10Sum(pretime/bin10bin);%pre stim 10s
Record_ChargeSRsum(k,1)=ChargeSRsum;
Record_ChargeINAR(k,1)=ChargeINAR;
Record_ChargeINsum(k,1)=Chargedursum;
Record_ChargePostsum(k,1)=bin10Sum(stimonset+durtime/bin10bin);%post stim 10s
Record_IntraSR_Ratio(k,1)=IntraSR_Ratio;
Record_IntraAR_Ratio(k,1)=IntraAR_Ratio;
Record_IntraARtoSR_Ratio(k,1)=IntraARtoSR_Ratio;
Record_PostAR20tobase20(k,1)=PostAR20tobase20;
Record_PostAR10tobase10(k,1)=PostAR10tobase10;
Record_current_removeA(k,1:length(current_removeA))=current_removeA;
Record_current_replaceA(k,1:length(current_replaceA))=current_replaceA;
Recordfilename{k,1}=filename;
Recorde_IntraARtoSR_Ratio_bin1s(k,1:length(IntraARtoSR_Ratio_bin1s))=IntraARtoSR_Ratio_bin1s;
Recorde_SR_bin1s(k,1:length(ChargeSRbinsum))=ChargeSRbinsum;
Recorde_IntraAR_Ratio_bin1s(k,1:length(IntraAR_Ratio_bin1s))=IntraAR_Ratio_bin1s;
Recorde_Basal_current1s(k,1:length(Basal_current1s))=Basal_current1s;
end
%% mat/excel output
save([datatitle,'.mat'],'Record_current_removeA','Record_current_replaceA','Record_Charge_current','Record_Charge_SR','Recorde_IntraARtoSR_Ratio_bin1s','Recorde_IntraAR_Ratio_bin1s','Recorde_Basal_current1s')

xlswrite(xlstitle,{'Name','ChargeBase10s(pA*s)','ChargeSR(pA*s)','ChargeINAR(pA*s)','ChargeINsum(pA*s)','ChargePost10s(pA*s)','IntraSR_Ratio','IntraAR_Ratio','IntraAR/SR','PostAR20tobase20','PostAR10tobase10'},'para','A1');
xlswrite(xlstitle,Recordfilename,'para','A2');
xlswrite(xlstitle,Record_ChargeBasesum,'para','B2');
xlswrite(xlstitle,Record_ChargeSRsum,'para','C2');
xlswrite(xlstitle,Record_ChargeINAR,'para','D2');
xlswrite(xlstitle,Record_ChargeINsum,'para','E2');
xlswrite(xlstitle,Record_ChargePostsum,'para','F2');
xlswrite(xlstitle,Record_IntraSR_Ratio,'para','G2');
xlswrite(xlstitle,Record_IntraAR_Ratio,'para','H2');
xlswrite(xlstitle,Record_IntraARtoSR_Ratio,'para','I2');
xlswrite(xlstitle,Record_PostAR20tobase20,'para','J2');
xlswrite(xlstitle,Record_PostAR10tobase10,'para','K2');

xlswrite(xlstitle,{'Name','Pre1(pA*s)','Pre2(pA*s)','Dur1(pA*s)','Dur2(pA*s)','Dur3(pA*s)','Post1(pA*s)','Post2(pA*s)'},'bin10s','A1');
xlswrite(xlstitle,Recordfilename,'bin10s','A2');
xlswrite(xlstitle,Record_bin10Sum,'bin10s','B2');

xlswrite(xlstitle,{'Name','bin1s 1-70'},'bin1s','A1');
xlswrite(xlstitle,Recordfilename,'bin1s','A2');
xlswrite(xlstitle,Record_binchargesum,'bin1s','B2');
xlswrite(xlstitle,Recordfilename,'bin1s','A10');
xlswrite(xlstitle,Record_binchargesum(:,(pretime/binPsec+1):(pretime+durtime)/binPsec),'bin1s','B10');

xlswrite(xlstitle,{'Name','ChargeSR_bin1s'},'ChargeSR','A1');
xlswrite(xlstitle,Recordfilename,'ChargeSR','A2');
xlswrite(xlstitle,Recorde_SR_bin1s,'ChargeSR','B2');

xlswrite(xlstitle,{'Name','IntraAR_Ratio_bin1s'},'ARratio','A1');
xlswrite(xlstitle,Recordfilename,'ARratio','A2');
xlswrite(xlstitle,Recorde_IntraAR_Ratio_bin1s,'ARratio','B2');

xlswrite(xlstitle,{'Name','IntraARtoSR_Ratio_bin1s'},'ARtoSR','A1');
xlswrite(xlstitle,Recordfilename,'ARtoSR','A2');
xlswrite(xlstitle,Recorde_IntraARtoSR_Ratio_bin1s,'ARtoSR','B2');

xlswrite(xlstitle,{'Name','Basal_current1s'},'Basal','A1');
xlswrite(xlstitle,Recordfilename,'Basal','A2');
xlswrite(xlstitle,Recorde_Basal_current1s,'Basal','B2');
