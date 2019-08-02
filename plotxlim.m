function [] = plotxlim(source,callbackdata,axeshandles,m)


range=m.ftrig_shift(end)-m.ftrig_shift(1);
xlim(axeshandles(1),[range*source.Value range*source.Value+100000]);
end
