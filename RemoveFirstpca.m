function [Fpca,dFvecpca] = RemoveFirstpca(F,pattern,spot_num,SpotMat,Stack,m)

Fpca=pcares(F,1);% remove first pca component

%% Calculating dF/F signal
    F0=prctile(Fpca,5,1);% We define F0 as the lowest 5% of the fluorescence signal over time.
    dF=(Fpca-F0)./F0;



overlay=zeros(size(pattern,1),size(pattern,2),spot_num);
    h1=figure;
    for idx=1:spot_num
    img(idx)=subplot(spot_num+1,12,(idx-1)*12+4:(idx-1)*12+5);% Includes saving room for sniff and trig time plot
    overlaytemp=(SpotMat{idx}==0).*Stack(:,:,round(size(Stack,3)/2));
    overlaytemp(overlaytemp==0)=max(max(overlaytemp));
    overlay(:,:,idx)=overlaytemp;
    imagesc(overlay(:,:,idx));axis equal;colormap('hot');
    end
    Patidx=cat(3, SpotMat{:});Patidx=sum(Patidx,3);
    totoverlay=(Patidx==0).*(Stack(:,:,round(size(Stack,3)/2)));totoverlay(totoverlay==0)=max(max(totoverlay));
    img(idx+1)=subplot(spot_num+1,12,(idx)*12+4:(idx)*12+5);imagesc(totoverlay);axis equal;
    
    %% plotting dF/F , sniff and triggers time
    for plotidx=1:spot_num
        
        plots(plotidx)=subplot(spot_num+1,12,(plotidx-1)*12+7:(plotidx-1)*12+12);
        
        plot(m.ftrig_shift,dF(:,plotidx),'b');ylim([min(dF(:,plotidx)) max(dF(:,plotidx))]); ylabel('\DeltaF/F_0');hold on;
        for plotline=1:size(m.time(m.shutter_timing>0),2)
            linetime=m.time(m.shutter_timing>0);
            plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx),'ylim')],'--r')
        end
    end
    
    %plotting the sniff measurements
    plots(plotidx+1)=subplot(spot_num+1,12,(plotidx)*12+7:(plotidx)*12+12);
    plot(m.time,m.sniff,'b');ylim([min(m.sniff) max(m.sniff)]); ylabel('Sniff');hold on;
    for plotline=1:size(m.time(m.shutter_timing>0),2)
        linetime=m.time(m.shutter_timing>0);
        plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx+1),'ylim')],'--r')
    end
    linkaxes([plots],'x'); % make sure all the x axis limits are the same
    
    sld = uicontrol('Style', 'slider','Units','normalized','Min',0,'Max',1,'SliderStep',[0.01 0.10],'Value',0,'Position', [0.5254    0.0500    0.3796    0.02],'Callback', {@plotxlim,plots,m}); % generate the slider in the plot

    %% Generating and plotting the dF average of all trials per stimulation spot
    %[dFvec,frameidx] = AvgFluorescenseTrials( m, dF,0 ); % for global f0
    [dFvec,frameidx,Ft0] = AvgFluorescenseTrials( m, F ,1,F0); % for local f0
    
    dFvecmat=cell2mat(dFvec');
    for plotidx2=1:size(dFvec,1)
        
        plots2(plotidx2)=subplot(spot_num+1,12,(plotidx2-1)*12+6);
        for plotidx3 = 1:size(dFvec,2)-1
            plot(-29:57,dFvec{plotidx2,plotidx3},'Color',[0.5 0.5 0.5]);ylabel('\DeltaF/F_0 Avg');
            hold on;
        end
        
        dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
        plot(-29:57,dFavg(plotidx2,:),'r','LineWidth',1); ylim([min(dFvecmat(:,plotidx2)) max(dFvecmat(:,plotidx2))]);   ylabel('\DeltaF/F_0 Avg');
        plot([0,0],[get(plots2(plotidx2),'ylim')],'--b')
    end