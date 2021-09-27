clear all
Start=100;
Interval=100;
NStep=100000;
column=int8((NStep-Start/Interval)+1);
latticeLength=70;
latticeBreadth=25;
margin=10;     % used for plotting
Frame(column,1)=0;
j=1;
VelField=VideoWriter('forceField','Uncompressed AVI');
VelField.FrameRate= 5;
%VelField.Quality= 100;
open(VelField);
figure(1)
for i=Start:Interval:NStep
%c=a(:,8).*a(:,5);
%d=a(:,8).*a(:,6);
filename = strcat('quiverplot_', int2str(i), '.dat');
a=importdata(filename);
%set(gcf,'Visible','off');
VelPlot=quiver(a(:,1),a(:,2),a(:,11),a(:,12),1.07,'r'); % to see velocity
%aplot=quiver(a(:,1),a(:,2),a(:,5),a(:,6),0,'k'); % to see polarizaiton
%aplot=quiver(a(:,1),a(:,2),c,d,0,'r');          % to see force along Pol
set(VelPlot,'ShowArrowHead','on');
axis([-1*margin latticeLength+margin -1*margin latticeBreadth+margin+5]);
title({'Force Field'; 'k=5, mu1=1, mu2=1, zeta=1, v0=1, h=0.0001' });
xlabel({'Timestep - ', i});
% F(j) = getframe(gcf);
frame=getframe(gcf);
writeVideo(VelField,frame);
j=j+1;
end
%movie(F,1,6); %(Frame, No. of Repitition, fps)
% VelField=VideoWriter('velField','Uncompressed AVI');
% VelField.FrameRate= 5;
% %VelField.Quality= 100;
% open(VelField);
% writeVideo(VelField,F);
close(VelField);
