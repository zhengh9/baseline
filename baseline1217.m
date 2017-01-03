%get Baseline via BAP, Iwan, New method; progress the test in 海口塔 , X1 main direction
%step1  original 1 directions: NS or EW acceleration
%step2 time-cut acceleration --start
%step3 FFT spectrum with the predominant frequency and the spectrum -fftY
%step4 baseline corrected Acc. in 3 methods:BAP, Iwan, New one,baseline comparing
%step5 corrected Acc with filter-butterworth in different past interval (channel) of frequency
%step6 subplot Acc-Velocity-Displacement
%step7  time-matching via function
% xcorr (Cross-correlation): DispIntegral vs DispObserve to DispIntegral_Reample vs DispObserve_Align
% step8 see the acor,lag
%step9 get the R via function corrcoef (Correlation coefficients ):
% DispIntegral_Reample vs DispObserve_Align
%modified: based on the main direction of pushing, NS，set the two
%direction NS/EW into the same timeDiff with NS. then see the
%corelation coefficient Rns and Rew
% input: original time-cut AccOrig(:,2)-2dir.s, dtAcc=1/256,dtVideo=1/50;
tic
clear  all;
close all
strFolderOrig='./海口塔数据/小震人工波S835-4-XYZ/';
strFolderProc='./海口塔处理结果图片/BLDataS/';
% Acc0(:,1)=load([strFolder '2 s835-4 x _SPECIMEN_123_A123__.TXT']);%point T1  top , unit G
% Acc0(:,2)=load([strFolder '2 s835-4 x _SPECIMEN_124_A124__.TXT']);%point T2  middle ,unit G
% Acc0(:,3)=load([strFolder '2 s835-4 x _SPECIMEN_125_A125__.TXT']);%point T3  bottom , unit G
% Acc=980*Acc0; %unit gal
% save Acc2S835.mat Acc
data= load('./海口塔处理结果图片/BLDataS/Acc2S835.mat');
draw=0%
% strFolder='./海口塔数据/大震人工波L840-7-X/';
% strFolderProc='./海口塔处理结果图片/BLDataL/';
% Acc0(:,1)=load([strFolder '26 L840-7 x _SPECIMEN_123_A123__.TXT']);%point T1  top , unit G
% Acc0(:,2)=load([strFolder '26 L840-7 x _SPECIMEN_124_A124__.TXT']);%point T2  middle ,unit G
% Acc0(:,3)=load([strFolder '26 L840-7 x _SPECIMEN_125_A125__.TXT']);%point T3  bottom , unit G
% Acc=980*Acc0; %unit gal
% save Acc2L840.mat Acc
% data= load('./海口塔处理结果图片/BLDataL/Acc2L840.mat');
fpast= [0.2    25]%[0.41   3.82]  %[0.4103  3.8209]  [0.1874   50]  [0.437   50]       S835T1
strEQ='S835' %'L840'
dir='T2';
% dir='T2';
% dir='T3';
switch dir
    case 'T1'
        AccOrig0=data.Acc(:,1);
%         strFN=[ strFolderOrig 'Displ. 2 s835-4 x 123_A123_.TXT']; %Small Earthqueke
        strFN=[ strFolderOrig '小震人工波S835-4-XYZ-X-T1.TXT'];
%         strFN=[ strFolder 'Displ. 26 L840-7 x 123_A123_.TXT']; %Large Earthqueke
    case 'T2'
        AccOrig0=data.Acc(:,2);
%         strFN=[ strFolderOrig 'Displ. 2 s835-4 x 124_A124_.TXT']; %Small
        strFN=[ strFolderOrig '小震人工波S835-4-XYZ-X-T2.TXT'];
%         strFN=[ strFolder 'Displ. 26 L840-7 x 123_A123_.TXT']; %Large        
     case 'T3'   
        AccOrig0=data.Acc(:,3);
%         strFN=[ strFolderOrig 'Displ. 2 s835-4 x 125_A125_.TXT']; %Small
        strFN=[ strFolderOrig '小震人工波S835-4-XYZ-X-T3.TXT'];
%         strFN=[ strFolder 'Displ. 26 L840-7 x 123_A123_.TXT']; %Large             
end
dtAcc=1/256;
dtVedio=1/50;
NPTS=numel(AccOrig0);

tAOrig0=((0:NPTS-1)*dtAcc)';%time of AccOrig
DispObserve0=load(strFN);%unit mm
DispObserve0=DispObserve0/10;%unit 1mm=1/10 cm
tV0=((0:length(DispObserve0)-1)*dtVedio)';%time of Video
nCutTA=find(tAOrig0==10);
AccOrig=AccOrig0(1:nCutTA);
nCutTV=find(tV0==10.5);
DispObserve=DispObserve0(1:nCutTV);
tA0=tAOrig0(1:nCutTA); %time of AccCut

%step3 FFT spectrum with the predominant frequency and the spectrum -fftY
% [predominantFre,fftY]=FASf(AccOrig,dtAcc);
Accfft=fft(AccOrig);
e=abs(Accfft)*2/NPTS;
f=[0:NPTS-1]/(NPTS*dtAcc);%
per=2;
n=floor(NPTS/per); %the nearest integer less than or equal to that element
fftY=e(1:n);
fre=f(1:n);
[maxY, indexY]=max(fftY);
predominantFre=fre(indexY);
if draw==1
figure
plot(fre,fftY);
hold on
plot([fre(indexY),fre(indexY)],[0,maxY],'--','LineWidth',2,'Color','r');
text(fre(indexY), maxY, [' f = ' num2str(fre(indexY))]);
title('Fourier Amplitude Spectrum');
xlabel('Frequecy(Hz)','FontName','times new roman');
ylabel('Amplitude(Hz^-^1)','FontName','times new roman');
hold off
end
%step4 baseline corrected Acc. and step5 baseline comparing

%BAP
[baselineBAP,AccBAP,preTime]=get_baselin_BAP(AccOrig,dtAcc);
%Iwan
[out] = get_baseline_iwan(AccOrig,dtAcc);
baselineIwan=out.baseline_full;
AccIwan=AccOrig-baselineIwan;
%new way
% [x] = get_baseline(Acc,dtAcc);
% baselineNew=x;
% strfileName=['Acc' dir '.mat'];
strfileName=[strFolderProc 'Acc2' strEQ dir 'BL.mat'];
xzh=load(strfileName);
baselineNew=xzh.x(1:nCutTA);
AccNew=AccOrig-baselineNew;
% Baseline Comparing
if draw==1
figure
title(['Baseline Comparing (' strEQ dir ')']);
xlabel('Time(S)','FontName','times new roman');
ylabel('Acceleration(gal)','FontName','times new roman');
hold on
plot(tA0,baselineBAP,'Color',[0 0.5 0],'LineWidth',1.5);
plot(tA0,baselineIwan,'--','Color',[0 0 1],'LineWidth',1.5);
plot(tA0,baselineNew,'Color',[1 0 0],'LineWidth',1.5);
legend('baselineBAP','baselineIwan','baselineNew',1);
hold off
end
% filter
% fpast=[0.04 25]
% [AccBAPf,D(:,1)]=filterf(AccBAP,fpast,dtAcc,dtVedio);
% [AccIwanf,D(:,2)]=filterf(AccIwan,fpast,dtAcc,dtVedio);
% [AccNewf,D(:,3)]=filterf(AccNew,fpast,dtAcc,dtVedio);
%nyquist frequency   f-采样频率，
% 则ny=f/2,按采样定理信号中的最高成分频率只能为采样频率的一半
ny = 1/(2*dtAcc);
% fpast=[0.04 25]%[0.2004,4.8447]%[0.04 25] NS[0.198 13.25]  [0.2 5.3]  %0.07  0.1 ; 25 50 75
% fpast=[0.035    19.6] %EW [0.1610  2.6000]
% fpast= [2.1964    3.3171]%[0.41   3.82]  %[0.4103  3.8209]  [0.1874   50]  [0.437   50]       S835T1
wn =fpast / ny;
n = 4;%4阶滤波
[b,a]=butter(n,wn);
AccBAPf = filter(b,a,AccBAP);
AccIwanf = filter(b,a,AccIwan);
AccNewf = filter(b,a,AccNew);

%figure4 of Accf and Baseline Comparing,
if draw==1
figure
title(['Acceleration Time History (' num2str(fpast(1)) '-' num2str(fpast(2)) ',' strEQ dir ')']);
xlabel('Time(S)','FontName','times new roman');
ylabel('Acceleration(gal)','FontName','times new roman');
hold on
plot(tA0,AccOrig,'Color','k','LineWidth',0.5)
plot(tA0,AccBAPf,'Color','g','LineWidth',1);
plot(tA0,AccIwanf,'Color','b','LineWidth',0.5);
plot(tA0,AccNewf,'--','Color','r','LineWidth',1.5);
plot(tA0,baselineBAP,'Color',[0 0.5 0],'LineWidth',0.2);
plot(tA0,baselineIwan,'Color',[0 0 0.5],'LineWidth',0.2);
plot(tA0,baselineNew,'--','Color',[0.8 0 0],'LineWidth',0.5);
legend('AccOrig','AccBAPf','AccIwanf','AccNewf','baselineBAP','baselineIwan','baselineNew',1);
hold off
end
%step6 subplot Acc-Velocity-Displacement
%figures of V and D via different Acc baseline corrected relatively  Comparing
[V(:,1),D0(:,1)]=VelocityDisplacement(AccBAPf,dtAcc);
[V(:,2),D0(:,2)]=VelocityDisplacement(AccIwanf,dtAcc);
[V(:,3),D0(:,3)]=VelocityDisplacement(AccNewf,dtAcc);

%step7 - 9
% % strExcelFolder='./data/';
% % fileName=[strExcelFolder '1cut.xls'];
% fileName=['1cut.csv'];
% sheetName='1cut';
% [num,txt,raw]=xlsread(fileName);
% [row clo]=size(num);%获取excel的行列数
% strRow=num2str(row+1)
% DispObserve=xlsread(fileName,sheetName,['C2:C' strRow]);

% DispObserve0=load(strFN);%unit mm
% DispObserve0=DispObserve0/10;%unit 1mm=1/10 cm
% tV0=((0:length(DispObserve0)-1)*dtVedio)';%time of Video
%fix the available time history
% nCutTA=find(tA0==10);
% D(:,1)=D0(1:nCutTA,1);D(:,2)=D0(1:nCutTA,2);D(:,3)=D0(1:nCutTA,3);
% nCutTV=find(tV0==10.5);
% DispObserve=DispObserve0(1:nCutTV);
D=D0;
%figure of Displacement of different Acc baseline corrected , unmatched
tA=((0:length(D(:,1))-1)*dtAcc)';%time of Acc
tV=((0:length(DispObserve)-1)*dtVedio)';%time of Video
figure
title(['Displacement Time History --Unmatched (' num2str(fpast(1)) '-' num2str(fpast(2)) ',' strEQ dir ')']);
xlabel('Time(S)','FontName','times new roman');
ylabel('Displacement(cm)','FontName','times new roman');
hold on
plot(tA,D(:,1),'Color','g','LineWidth',0.5);
plot(tA,D(:,2),'Color','b','LineWidth',0.5);
plot(tA,D(:,3),'Color','k','LineWidth',1.5);
plot(tV,DispObserve,'Color','r','LineWidth',0.5);
legend('DispBAP','DispIwan','DispNew','DispObserve',1);
plot(tA, zeros(1,numel(tA)),'Color',[0.8 0.8 0.8],'LineWidth',0.2);
hold off
%step7  time-matching via function
% xcorr (Cross-correlation): DispIntegral vs DispObserve to DispIntegral_Reample vs DispObserve_Align
% step8 see the acor,lag
%step9 get the R via function corrcoef (Correlation coefficients ):
% DispIntegral_Reample vs DispObserve_Align
[ out1 ] = corrf( D(:,1),DispObserve,fpast,dir,dtAcc,dtVedio)
[ out2 ] = corrf( D(:,2),DispObserve,fpast,dir,dtAcc,dtVedio)
[ out3 ] = corrf( D(:,3),DispObserve,fpast,dir,dtAcc,dtVedio)
%with matching result--by out3

figure
title(['Displacement Time History of 3 Methods with Matching DispNew (' num2str(fpast(1)) '-' num2str(fpast(2)) ',' strEQ dir ')']);
xlabel('Time(S)','FontName','times new roman');
ylabel('Displacement(cm)','FontName','times new roman');
hold on
plot(tA,D(:,1),'Color','g','LineWidth',0.5);
plot(tA,D(:,2),'Color','b','LineWidth',0.5);
plot(tA,D(:,3),'--','Color','k','LineWidth',1.5);
tV=((0:length(out3.DispObserve_AlignCo)-1)*dtVedio)';
if out3.R(2,1)<0 
    Dplot=-out3.DispObserve_AlignCo;
else
    Dplot=out3.DispObserve_AlignCo;
end
plot(tV,Dplot,'Color','r','LineWidth',0.5);
legend('DispBAP','DispIwan','DispNew','DispObserve',1);
text(tV(floor(numel(tV)*0.8)), max(abs(Dplot))/2, [' R=' num2str(out3.R(2,1)) ',' 'timeDiff=' num2str(out3.timeDiff) ],'FontSize',10);
out3.R(2,1)
out3.timeDiff
plot(tA, zeros(1,numel(tA)),'Color',[0.8 0.8 0.8],'LineWidth',0.2);
hold off
toc


%
% figure %DispObserve
% title(['Displacement Time History --Observed NS (' num2str(fpast(1)) '-' num2str(fpast(2)) ',' strEQ dir ')']);
% xlabel('Time(S)','FontName','times new roman');
% ylabel('Displacement(cm)','FontName','times new roman');
% hold on
% % plot(t,D(:,1),'Color','g','LineWidth',0.5);
% % plot(t,D(:,2),'Color','b','LineWidth',0.5);
% % plot(t,D(:,3),'Color','k','LineWidth',1.5);
% tt=((0:length(DispObserve)-1)*dtVedio)';
% plot(tt,DispObserve,'Color','b','LineWidth',1.5);
% legend('DispObserveNS',1);
% plot(tt, zeros(1,length(DispObserve)),'Color',[0.8 0.8 0.8],'LineWidth',0.2);
% hold off
%
% figure
% title(['Original Acceleration Time History (NS)']);
% xlabel('Time(S)','FontName','times new roman');
% ylabel('Acceleration(gal)','FontName','times new roman');
% hold on
% plot( ((0:length(data.AccNS)-1)*dtAcc)' ,data.AccNS,'Color','b','LineWidth',1.5)
% plot( ((70/dtAcc:135/dtAcc)*dtAcc)' ,data.AccNS(70/dtAcc:135/dtAcc),'r','LineWidth',1.0)
% legend('AccOrig','AccOrigCut',1);
% hold off

% VelocityDisplacement(AccBAP,dtAcc);
% VelocityDisplacement(AccIwan,dtAcc);
% VelocityDisplacement(AccNew,dtAcc);
%without filter, get the R
if draw==1
[Vnew,Dnew0]=VelocityDisplacement(AccNew,dtAcc);
Dnew=Dnew0(1:nCutTA);
[ outnew ] = corrf( Dnew,DispObserve,fpast,dir,dtAcc,dtVedio)
figure
title(['Displacement Time History with Matching DispNew (without filter,' strEQ dir ')']);
xlabel('Time(S)','FontName','times new roman');
ylabel('Displacement(cm)','FontName','times new roman');
hold on
plot(tA,Dnew,'--b','LineWidth',1);
% plot(t,D(:,1),'Color','g','LineWidth',0.5);
% plot(t,D(:,2),'Color','b','LineWidth',0.5);
plot(tA,D(:,3),'Color','g','LineWidth',1);
tV=((0:length(outnew.DispObserve_AlignCo)-1)*dtVedio)';
if outnew.R(2,1)<0 
    Dplot=-outnew.DispObserve_AlignCo
else
    Dplot=outnew.DispObserve_AlignCo;
end
plot(tV,Dplot,'Color','r','LineWidth',0.5);
legend('DispNew','DispNewFilter','DispObserve',1);
text(tV(floor(numel(tV)*0.8)), max(abs(Dplot))/2, [' R=' num2str(outnew.R(2,1)) ',' 'timeDiff=' num2str(outnew.timeDiff) ],'FontSize',10);
text(tV(floor(numel(tV)*0.8)), -max(abs(Dplot))/2, [' R_f(' num2str(fpast(1)) '-' num2str(fpast(2)) ')=' num2str(out3.R(2,1)) ',timeDiff_f=' num2str(out3.timeDiff)],'FontSize',10);
outnew.R(2,1)
outnew.timeDiff
plot(tA, zeros(1,numel(tA)),'Color',[0.8 0.8 0.8],'LineWidth',0.2);
hold off
end

%判断速度和位移最后是周期震荡
%思路：1. Vend / Dend tend to be Asin(omiga*x+phi)
%          2. fft(V-last/ D-last) --周期性.
tD=((0:length(Dplot)-1)*dtVedio)';
nt=find(tD==7);
Dend=Dplot(nt:end);
Dendfft=fft(Dend);
nDend=numel(Dendfft)
e=abs(Dendfft)*2/nDend;
fDend=[0:nDend-1]/(nDend*dtAcc);%
per=2;
n=floor(nDend/per); %the nearest integer less than or equal to that element
fftDend=e(1:n);
freDend=fDend(1:n);
[maxDend, indexDend]=max(fftDend);
predominantFre=freDend(indexDend);
% if draw==1
figure
plot(freDend,fftDend);
hold on
plot([freDend(indexDend),freDend(indexDend)],[0,maxDend],'--','LineWidth',2,'Color','r');
text(freDend(indexDend), maxDend, [' f = ' num2str(freDend(indexDend))]);
title('Fourier Amplitude Spectrum of Disp.');
xlabel('Frequecy(Hz)','FontName','times new roman');
ylabel('Amplitude(Hz^-^1)','FontName','times new roman');
hold off
% end



