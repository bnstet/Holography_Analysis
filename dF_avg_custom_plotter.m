% It basically works after the analysis was done and all parameters are
% still in the workspace. you can play with the way dF avg is plotted
figure
PreStimFrames = 30; % the number of frames before stim to include in avg dF plots 
PostStimFrames = 76; % the number of frames after stim to include in avg dF plots 

n=1;
cellschosen=[2 ,3 4, 6, 7, 8];
for plotidx2=cellschosen;
%for plotidx2=1:size(dFvec,1)
    
    %plots2(plotidx2)=subplot(spot_num+1,13,(plotidx2-1)*13+6);
    plots2(plotidx2)=subplot(2,3,n);
    dFcell = [dFvec{plotidx2,:}];
    dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
    dFSEM(plotidx2,1:size(dFvec{1,1},1))=std(dFcell,0,2)/sqrt(size([dFvec{plotidx2,:}],2));
    ActualStims(plotidx2,1:size(dFvec,2)) = cellfun(@(x) dot(dFavg(plotidx2,:),x)./(norm(dFavg(plotidx2,:))*norm(x)),dFvec(plotidx2,:));
    Tavg=-PreStimFrames+Omitpre:PostStimFrames+1-Omitpost; % the time for the average df/f plot
    plot(Tavg,dFavg(plotidx2,:),'r','LineWidth',1);
    hold on
    errorshade(Tavg,dFavg(plotidx2,:),dFSEM(plotidx2,:),dFSEM(plotidx2,:),'r','scalar');
    ylim([mean(dFavg(plotidx2,:))-3*std(dFavg(plotidx2,:)) mean(dFavg(plotidx2,:))+3*std(dFavg(plotidx2,:))])
    xlim([Tavg(1)-1 Tavg(end)+1]);
    
    GoodStim(plotidx2)=sum(StimEff(plotidx2,:));
    AmpMean(plotidx2)=mean(nonzeros(StimAmp(plotidx2,:)));
    AmpSTD(plotidx2)=std(nonzeros(StimAmp(plotidx2,:)));
    
    plot([0,0],[get(plots2(plotidx2),'ylim')],'--b')
%     ylabel({['Cell # ',num2str(plotidx2)],[num2str(GoodStim(plotidx2)),'/',...
%         num2str(size(StimEff,2)),' stimulated'],['Avg Amp ',num2str(AmpMean(plotidx2))]...
%         '\DeltaF/F_0 Avg'},'FontSize',7);
%     ylabel({'\DeltaF/F_0 Avg'},'FontSize',7);
%     xlabel('# Frames from stim onset');
    set(gca,'FontSize',5);
%set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[])

    n=n+1;
end

linkaxes(plots2,'y')
for plotidx2=cellschosen
plot(plots2(plotidx2),[0,0],[get(plots2(plotidx2),'ylim')],'--b')
end
set(plots2(plotidx2),'ylim',[-0.1,0.2]);
set(gcf,'color','w')
