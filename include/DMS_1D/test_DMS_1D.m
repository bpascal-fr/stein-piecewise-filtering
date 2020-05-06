S = [((0:.01:2)-3).^2-1, ((2.01:.01:4)-2.5).^2+2, ((4.01:.01:5)-3).^2];

res = dms_1D(S,50,0.0001,'AddNoise',[1 0.5]);

figure(1), clf, colormap gray
plot(res.data,'linewidth',1);   hold on
%plot(res.u,'linewidth',2);
[ax,l1,l2] = plotyy(1:length(res.u),res.u,1:length(res.e),res.e); hold off;
linkaxes(ax,'x');
set(l1,'LineWidth',2);
set(l2,'marker','*','linestyle',':','markersize',5,'Color',[0.9290    0.6940    0.1250]);
axis([-inf +inf -inf +inf])
legend('data z','estimate u','jump set e')