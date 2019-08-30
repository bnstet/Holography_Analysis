clear variables
close all
clc

filesdir=uigetdir('C:\Users\Gilad\Documents\Experiments\Stim trials');
n=1;
loadmore='Yes';
while strcmp(loadmore,'Yes')
    [name,path]=uigetfile('*.mat','Please choose the first mat file of the analysis',filesdir);
    fname=[path,name];
    newData1 = load('-mat', fname);
    vars = fieldnames(newData1);
    for i = 1:length(vars)
        Data{i,n}=newData1.(vars{i});
    end
    Data{i+1,n}=cell2mat(inputdlg('What is the laser power used in this measurement [mw] ?','Laser power',1));
    
    AmpMean(1:size(newData1.AmpMean,2),n)=newData1.AmpMean;
    AmpSTD(1:size(newData1.AmpSTD,2),n)=newData1.AmpSTD;
    GoodStim(1:size(newData1.GoodStim,2),n)=newData1.GoodStim;
    StimTries(1:size(newData1.StimTries,2),n)=newData1.StimTries; 
    StimDeltaValue(1:size(newData1.StimDeltaValue,1),1:size(newData1.StimDeltaValue,2),n)=newData1.StimDeltaValue; 
    Power(n)=str2num(Data{i+1,n});
    
    loadmore=questdlg('Would you like to load another file ?','Loading files for analysis',...
        'Yes','No','Yes');
    n=n+1;
end
MeanDelta=squeeze(mean(StimDeltaValue,2));

%h1=figure;
%h2=figure;
h3=figure;
h4=figure;
for m=1:size(AmpMean,1)
    StimEff(m,:)=GoodStim(m,:)./StimTries;
    %figure(h1);plot(Power,AmpMean(m,:));hold on
    %figure(h2);plot(Power,AmpSTD(m,:));hold on
    figure(h3);plot(Power,StimEff(m,:));hold on
    figure(h4);plot(Power,MeanDelta(m,:));hold on
end

%figure; plot(Power,mean(StimEff,1),'g','LineWidth',3);
figure; plot(Power,mean(MeanDelta,1),'r','LineWidth',3);
DeltaSEM=std(MeanDelta,0,1)./sqrt(size(MeanDelta,1));
errorshade(Power,mean(MeanDelta,1),DeltaSEM,DeltaSEM,'r','scalar');
xlabel('Average laser power [mw]'); ylabel('Mean (post-pre) signal');grid on;