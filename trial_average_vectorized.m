frameidx_adj = frameidx(frameidx > preCalcPeriod);
prestimwindow = (-preCalcPeriod:-Omitpre)' + frameidx_adj; prestimlen = size(prestimwindow,1);
poststimwindow = (Omitpost+1:postCalcPeriod)' + frameidx_adj; poststimlen = size(poststimwindow,1);
prestimwindow = reshape(prestimwindow,1,[]);
poststimwindow = reshape(poststimwindow,1,[]);
prestimmean = mean(reshape(F(prestimwindow,:),prestimlen,[],size(F,2)),1);
poststimmean = mean(reshape(F(poststimwindow,:),poststimlen,[],size(F,2)),1);
prestimstd= std(reshape(F(prestimwindow,:),prestimlen,[],size(F,2)),1);
poststimstd = std(reshape(F(poststimwindow,:),poststimlen,[],size(F,2)),1);

newdf = squeeze((poststimmean - prestimmean) ./ prestimmean);
newzscore = squeeze( (poststimmean - prestimmean) ./ sqrt(.5*(prestimstd.^2 + poststimstd.^2)));
newdfmean = squeeze(mean(newdf,1));

figure;
bar(newdfmean);