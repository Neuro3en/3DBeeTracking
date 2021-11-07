%% Batch analysis of 3D tracked Wei Data by Benjamin Paffhausen CC4.0 BY

% put script in base folder were the data2.mat files are and run it, it
% will seach out the data2 files also from subfolders.
% ARROW keys to sync the interactions, the lowest distance should be at timepoint 100 (500 ms)
% then click space
% if there is no curve just click arrow key (LEFT or RIGTH)

%% Import the data 
% all data will be in dataAll struct, with tracks and anglesSelf (percieved oposite angle)
version=2; % if 1 seach for data, with out angles, if 2 seach for data2, with angles
if version ==1
    Folder   = cd;
    FileList = dir(fullfile(Folder, '**', 'data*'));
    cases=length(FileList);
    for i = 1:cases
        load([FileList(i).folder '\data.mat'])
        intermeiiten=cat(3,hornetNew,hornetNewSide(:,:,2));
        dataAll(i).track=intermeiiten;
    end
else % for data2, with angles
    Folder   = cd;
    FileList = dir(fullfile(Folder, '**', 'data2*'));
    cases=length(FileList);
    for i = 1:cases
        load([FileList(i).folder '\data2.mat'])
        intermeiiten=cat(3,hornetNew,hornetNewSide(:,:,2));
        dataAll(i).track=intermeiiten;
        dataAll(i).anglesSelf=cat(3,hornetsAngle,BeeAngle);
        %    dataAll(i).anglesWorld=cat(3,hornetNewO,hornetNewSideO);
    end
end

% some angles are multiples of 360°, this will be fixed here
for badcode=1:10
    for k=1:2
        for n=1:2
            for i=1:cases-1
                for j=1:length(dataAll(i).anglesSelf)
                    if  dataAll(i).anglesSelf(j,k,n) > 180
                        dataAll(i).anglesSelf(j,k,n) = dataAll(i).anglesSelf(j,k,n)-360;
                    end
                    if  dataAll(i).anglesSelf(j,k,n) < -180
                        dataAll(i).anglesSelf(j,k,n) = dataAll(i).anglesSelf(j,k,n)+360;
                    end
                end
            end
        end
    end
end


%% 3D Plot, all of them
figure('units','normalized','outerposition',[0 0 1 1])    % FIG 1
for i = 1:cases
    subplot(ceil(sqrt(cases)),ceil(sqrt(cases)),i)
    plot3(dataAll(i).track(:,1,1),dataAll(i).track(:,1,2),dataAll(i).track(:,1,3),'b.')
    hold on
    plot3(dataAll(i).track(:,2,1),dataAll(i).track(:,2,2),dataAll(i).track(:,2,3),'r.')
end
saveas(gcf,'3D_overview.png')


%% Velocity & Disance
figure                                                  % FIG 2
for i = 1:cases
    %   subplot(9,9,i)
    distanceChangeH=((sqrt((diff(dataAll(i).track(:,1,1)).^2+diff(dataAll(i).track(:,1,2)).^2)+diff(dataAll(i).track(:,1,3)).^2)));
    distanceChangeB=((sqrt((diff(dataAll(i).track(:,2,1))).^2+diff(dataAll(i).track(:,2,2)).^2)+diff(dataAll(i).track(:,2,3)).^2));
    distanceChangeHsmoothed = filloutliers(distanceChangeH,'linear','movmedian',5,'ThresholdFactor',.8);
    distanceChangeBsmoothed = filloutliers(distanceChangeB,'linear','movmedian',5,'ThresholdFactor',.8);
    a=dataAll(i).track(:,1,1)-dataAll(i).track(:,2,1);
    b=dataAll(i).track(:,1,3)-dataAll(i).track(:,2,3);
    c=dataAll(i).track(:,1,2)-dataAll(i).track(:,2,2);
    hold on
    plot(distanceChangeHsmoothed,'b')
    plot(distanceChangeBsmoothed,'r')
    plot(sqrt(a.^2+b.^2+c.^2)*0.02,'g')
    title('velocity over time Hornet(blue) Bee(red)| distance(yellow)')
end


%% sync distance
figure                                               
n=1;
for i = 1:cases
    trackNow=dataAll(i).track(:,:,:);
    buttom=1;
    while buttom<32
        [x, y, buttom]=ginput(1);
        if buttom==29
            trackNow(6:length(trackNow)+5,:,:)=trackNow(:,:,:);
            %     trackNow(6:end,:,:)=trackNow(1:end,:,:);
            trackNow(1:5,:,:)=nan(5,2,3);
            a=trackNow(:,1,1)-trackNow(:,2,1);
            b=trackNow(:,1,3)-trackNow(:,2,3);
            c=trackNow(:,1,2)-trackNow(:,2,2);
            distance=sqrt(a.^2+b.^2+c.^2)*0.02;
            plot(distance)
            hold on
            plot(100,1:.2:10,'k.')
            hold off
        end
        if buttom==30
            distance=distance+10;
            plot(distance)
            hold on
            plot(100,1:.2:10,'k.')
            hold off
        end
        if buttom==28
            trackNow=trackNow(6:end,:,:);
            a=trackNow(:,1,1)-trackNow(:,2,1);
            b=trackNow(:,1,3)-trackNow(:,2,3);
            c=trackNow(:,1,2)-trackNow(:,2,2);
            distance=sqrt(a.^2+b.^2+c.^2)*0.02;
            plot(distance)
            hold on
            plot(100,1:.2:10,'k.')
            hold off
        end
        if buttom==31
            distance=distance-10;
            plot(distance)
            hold on
            plot(100,1:.2:10,'k.')
            hold off
        end
    end
    if buttom ~= 35 % raute; pound; hashtag; #
        syncData(n).track=trackNow;
        n=n+1;
        hold off
    end
    hold off
end
close
syncCases=length(syncData);


%%  Disance SYNCed
figure                                      % FIG 3
hold on
for i = 1:syncCases
    syncData(i).veloH=((sqrt((diff(syncData(i).track(:,1,1)).^2+diff(syncData(i).track(:,1,2)).^2)+diff(syncData(i).track(:,1,3)).^2)));
    syncData(i).veloB=((sqrt((diff(syncData(i).track(:,2,1))).^2+diff(syncData(i).track(:,2,2)).^2)+diff(syncData(i).track(:,2,3)).^2));
    syncData(i).veloH = filloutliers(syncData(i).veloH,'linear','movmedian',5,'ThresholdFactor',.8);
    syncData(i).veloB = filloutliers(syncData(i).veloB,'linear','movmedian',5,'ThresholdFactor',.8);
    a=syncData(i).track(:,1,1)-syncData(i).track(:,2,1);
    b=syncData(i).track(:,1,3)-syncData(i).track(:,2,3);
    c=syncData(i).track(:,1,2)-syncData(i).track(:,2,2);
   % plot(syncData(i).veloH,'b')
   % plot(syncData(i).veloB,'r')
    distance=sqrt(a.^2+b.^2+c.^2)*0.02;
    plot([1:1:length(distance)]/250,distance,'g')
    xlim([0 .8])
    title('distance between honet and bee')
    ylabel('distance [cm]')
    xlabel('time [s]')
end
saveas(gcf,'Distance_syncd.png')
distance2=distance;
dataspeedH=nan(syncCases,1000);
for i=1:syncCases
    dataspeedH(i,1:size(syncData(i).veloH))= syncData(i).veloH;
end

dataspeedB=nan(syncCases,1000);
for i=1:syncCases
    dataspeedB(i,1:size(syncData(i).veloB))= syncData(i).veloB;
end

distance=nan(syncCases,1000);
for i=1:syncCases
    
    a=syncData(i).track(:,1,1)-syncData(i).track(:,2,1);
    b=syncData(i).track(:,1,3)-syncData(i).track(:,2,3);
    c=syncData(i).track(:,1,2)-syncData(i).track(:,2,2);
    
    distance(i,1:length((sqrt(a.^2+b.^2+c.^2)*0.02)))=(sqrt(a.^2+b.^2+c.^2)*0.02);
end

figure('units','normalized','outerposition',[0 0 1 1])  % FIG 4
h1=histogram(5.*dataspeedH(:,1:100)','EdgeColor','none','FaceColor','b');
title('Hornet(blue) and bee(red) speed distribution')
xlim([0 250])
hold on
h2=histogram(5.*dataspeedB(:,1:100)','EdgeColor','none','FaceColor','r');
ylabel('frequency')
xlabel('flight speed [cm/s]')
 h1.Normalization = 'probability';
 h1.BinWidth = 1;
 h2.Normalization = 'probability';
 h2.BinWidth = 1;
saveas(gcf,'speedDistro_till100.png')

figure('units','normalized','outerposition',[0 0 1 1])  % FIG 5
subplot(2,1,1) 
boxplot(dataspeedH.*5)
title('hornet speed distribution over time')
xlim([0 100])
ylim([0 100])
ylabel('flight speed [cm/s]')
xlabel('time [frames]')
% figure
% histogram(5.*dataspeedB',300)
% title('Bees speed distribution')
% xlim([0 250])
% ylabel('frequency')
% xlabel('flight speed [cm/s]')
subplot(2,1,2)                          % FIG 6
boxplot(dataspeedB.*5)
title('Bees speed distribution over time')
xlim([0 100])
ylim([0 100])
ylabel('flight speed [cm/s]')
xlabel('time [frames]')
saveas(gcf,'speedDistro_overTime.png')

figure
%subplot(2,1,2)                         % FIG 7
boxplot(distance)
title('distance distribution over time')
xlim([0 100])
ylim([0 inf])
ylabel('distance [cm]')
xlabel('time [frames]')
saveas(gcf,'distanceOverTime.png')


% acceleration
figure                  % FIG 8
hold on
animal=1; % 1hornet 2bee
alleDiffDiff=zeros(1800,1);
for i=1:syncCases-1
    dataAll(i).DiffDiffH=((diff(dataAll(i).track(1:end-1,animal,1))-diff(dataAll(i).track(2:end,animal,1)))+...
        (diff(dataAll(i).track(1:end-1,animal,2))-diff(dataAll(i).track(2:end,animal,2)))+...
        (diff(dataAll(i).track(1:end-1,animal,3))-diff(dataAll(i).track(2:end,animal,3))));
    plot(dataAll(i).DiffDiffH+30*i,'b')
    alleDiffDiff(1:length(dataAll(i).DiffDiffH))=dataAll(i).DiffDiffH+alleDiffDiff(1:length(dataAll(i).DiffDiffH));
end

figure('units','normalized','outerposition',[0 0 1 1])  % FIG 9
animal=1;
subplot(2,3,1)
for i=1:syncCases-1
    hold on;    plot(dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,2))
end
title('hornet x y [top]')
subplot(2,3,2)
for i=1:syncCases-1
    hold on;    plot(dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,3))
end
title('hornet x z [side]')
set (gca,'Ydir','reverse')
subplot(2,3,3)
for i=1:syncCases-1
    hold on;    plot(dataAll(i).track(:,animal,2),dataAll(i).track(:,animal,3))
end
title('hornet y z [front/imaginary]')
set (gca,'Ydir','reverse')
animal=2;
subplot(2,3,4)
for i=1:syncCases-1
    hold on;    plot(dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,2))
end
title('bee x y [side]')
subplot(2,3,5)
for i=1:syncCases-1
    hold on;    plot(dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,3))
end
title('bee x z [top]')
set (gca,'Ydir','reverse')
subplot(2,3,6)
for i=1:syncCases-1
    hold on;    plot(dataAll(i).track(:,animal,2),dataAll(i).track(:,animal,3))
end
title('bee y z [front/imaginary]')
set (gca,'Ydir','reverse')
saveas(gcf,'trajectoriesAcrossPlanes.png')


figure('units','normalized','outerposition',[0 0 1 1])  % FIG 10
animal=1;
subplot(2,3,1)
allBeesCoor=[0, 0];
for i=1:length(dataAll)
    allBeesCoor=cat(1,allBeesCoor, [dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,2)]);
end
hist3(allBeesCoor,'Nbins',[30 30],'CdataMode','auto','LineStyle','none')
view(2)
title('hornet x y')
pbaspect([1 .5625 1])
xlim([0 1280])
ylim([0 720])
subplot(2,3,2)
allBeesCoor=[0, 0];
for i=1:length(dataAll)
    allBeesCoor=cat(1,allBeesCoor, [dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,3)]);
end
hist3(allBeesCoor,'Nbins',[30 30],'CdataMode','auto','LineStyle','none')
view(2)
title('hornet x z')
pbaspect([1 .5625 1])
xlim([0 1280])
ylim([0 720])
set (gca,'Ydir','reverse')
subplot(2,3,3)
allBeesCoor=[0, 0];
for i=1:length(dataAll)
    allBeesCoor=cat(1,allBeesCoor, [dataAll(i).track(:,animal,2),dataAll(i).track(:,animal,3)]);
end
hist3(allBeesCoor,'Nbins',[30 30],'CdataMode','auto','LineStyle','none')
view(2)
title('hornet y z')
pbaspect([1 1 1])
xlim([0 720])
ylim([0 720])
set (gca,'Ydir','reverse')
animal=2;
subplot(2,3,4)
allBeesCoor=[0, 0];
for i=1:length(dataAll)
    allBeesCoor=cat(1,allBeesCoor, [dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,2)]);
end
hist3(allBeesCoor,'Nbins',[30 30],'CdataMode','auto','LineStyle','none')
view(2)
title('bee x y')
pbaspect([1 .5625 1])
xlim([0 1280])
ylim([0 720])
subplot(2,3,5)
allBeesCoor=[0, 0];
for i=1:length(dataAll)
    allBeesCoor=cat(1,allBeesCoor, [dataAll(i).track(:,animal,1),dataAll(i).track(:,animal,3)]);
end
hist3(allBeesCoor,'Nbins',[30 30],'CdataMode','auto','LineStyle','none')
view(2)
title('bee x z')
xlim([0 1280])
ylim([0 720])
pbaspect([1 .5625 1])
set (gca,'Ydir','reverse')
subplot(2,3,6)
allBeesCoor=[0, 0];
for i=1:length(dataAll)
    allBeesCoor=cat(1,allBeesCoor, [dataAll(i).track(:,animal,2),dataAll(i).track(:,animal,3)]);
end
hist3(allBeesCoor,'Nbins',[30 30],'CdataMode','auto','LineStyle','none')
view(2)
title('bee y z')
pbaspect([1 1 1])
xlim([0 720])
ylim([0 720])
set (gca,'Ydir','reverse')
saveas(gcf,'trajectoriesAcrossPlanes_accumulated.png')


%% angles
%hornet
figure('units','normalized','outerposition',[0 0 1 1])  % FIG 11
subplot(1,2,1)
allBeesAngles=[0, 0];
animal=1;
for i=1:length(dataAll)
    allBeesAngles=cat(1,allBeesAngles, [dataAll(i).anglesSelf(:,animal,1),dataAll(i).anglesSelf(:,animal,2)]);
end
allBeesAngles(allBeesAngles(:,2)==allBeesAngles(:,1))=nan;
for fucksake=1:10
    m=1;
    for i=1:length(allBeesAngles)
        if m==0
            m=1;
        end
        if isnan(allBeesAngles(m,1)) || isnan(allBeesAngles(m,2))
            allBeesAngles(m,:)=[];
            m=m-2;
        end
        m=m+1;
    end
end
hist3(allBeesAngles,'Nbins',[140 140],'CdataMode','auto','EdgeColor','none')
xlim([-180 180])
ylim([-180 180])
title('hornet view the bees')
view(2)
pbaspect([1 1 1])

%bee
subplot(1,2,2)
allBeesAngles=[0, 0];
animal=2;
for i=1:length(dataAll)
    allBeesAngles=cat(1,allBeesAngles, [dataAll(i).anglesSelf(:,animal,1),dataAll(i).anglesSelf(:,animal,2)]);
end
allBeesAngles(allBeesAngles(:,2)==allBeesAngles(:,1))=nan;
for fucksake=1:10
    m=1;
    for i=1:length(allBeesAngles)
        if m==0
            m=1
        end
        if isnan(allBeesAngles(m,1)) || isnan(allBeesAngles(m,2))
            allBeesAngles(m,:)=[];
            m=m-2;
        end
        m=m+1;
    end
end
hist3(allBeesAngles,'Nbins',[70 70],'CdataMode','auto','EdgeColor','none')
xlim([-180 180])
ylim([-180 180])
title('bees view the hornet')
view(2)
pbaspect([1 1 1])
saveas(gcf,'relativeviews.png')




%% test field
% switched off for now
%{
% X Y
for i=1:syncCases-1
    plot(dataAll(i).track(:,1,1),dataAll(i).track(:,1,3),'b.')
    hold on;
    plot(dataAll(i).track(:,2,1),dataAll(i).track(:,2,3),'r.')
end

% x time
for i=1:syncCases-1
    plot(dataAll(i).track(:,1,3)+i*50,'b.')
    hold on;
    plot(dataAll(i).track(:,2,3)+i*50,'r.')
end


% x corsscorr
figure
subplot(1,3,1)
from=1;
too=100;
clear XCF
for i=1:syncCases-1
    if size(dataAll(i).track,1)>too
        hornetX=dataAll(i).track(from:too,1,1);
        beeX=dataAll(i).track(from:too,2,1);
        hornetX(isnan(hornetX))=0;
        beeX(isnan(beeX))=0;
        
        [XCF(:,i),lags,bounds] = crosscorr(hornetX,beeX,'NumLags',50);
        plot(XCF)
        hold on
    end
end
plot(nanmean(XCF'),'*')
%title('X crosscorrelation over time, *mean*')

% Y crosscorr
subplot(1,3,2)
clear XCF
for i=1:syncCases-1
    if size(dataAll(i).track,1)>too
        hornetX=dataAll(i).track(from:too,1,2);
        beeX=dataAll(i).track(from:too,2,2);
        hornetX(isnan(hornetX))=0;
        beeX(isnan(beeX))=0;
        
        [XCF(:,i),lags,bounds] = crosscorr(hornetX,beeX,'NumLags',50);
        plot(XCF)
        hold on
    end
end
plot(nanmean(XCF'),'*')
title('X Y Z crosscorrelation over time, *mean*')

% Z crosscorr
subplot(1,3,3)
clear XCF
for i=1:syncCases-1
    if size(dataAll(i).track,1)>too
        hornetX=dataAll(i).track(from:too,1,3);
        beeX=dataAll(i).track(from:too,2,3);
        hornetX(isnan(hornetX))=0;
        beeX(isnan(beeX))=0;
        
        [XCF(:,i),lags,bounds] = crosscorr(hornetX,beeX,'NumLags',50);
        plot(XCF)
        hold on
    end
end
plot(nanmean(XCF'),'*')


%% speed over time
figure
from=1;
too=100;
for i=1:syncCases-1
    if size(dataAll(i).track,1)>too
        subplot(7,7,i)
        plot(smooth(diff(dataAll(i).track(from:too,1,1)),10))
        hold on
        plot(smooth(diff(dataAll(i).track(from:too,2,1)),10))
    end
end

%}

figure('units','normalized','outerposition',[0 0 1 1])  % FIG 12
from=1;
too=200;
for k=1:2
    for n=1:2
        for i=1:syncCases-1
            if size(dataAll(i).track,1)>too
                if k==1 && n==1
                    m=1;
                elseif k==1 && n==2
                    m=2;
                elseif k==2 && n==1
                    m=3;
                else
                    m=4;
                end
                subplot(2,2,m)
                if m==1
                    title('hornet top self-angle')
                elseif m==2
                    title('hornet side self-angle')
                elseif m==3
                    title('bee top self-angle')
                else 
                    title('bee side self-angle')
                end
                % angels(:,i)=(smooth((dataAll(i).anglesSelf(from:too,1,2)),10));
                plot(((dataAll(i).anglesSelf(from:too,k,n))),'.')
                hold on
                ylim([-180 180])
            end
        end
    end
end
saveas(gcf,'anglesOverTime.png')


figure('units','normalized','outerposition',[0 0 1 1])  % FIG 13
xfrom=50;
xto=100;
yfrom=-100;
yto=100;
subplot(2,2,1)
for i=1:syncCases-1
    if size(dataAll(i).track,1)>99
selfAngleAll1(i,:)=(dataAll(i).anglesSelf(1:100,1,1));
    end
end
boxplot(selfAngleAll1)
title('hornet top selfangle')
xlim([xfrom xto])
ylim([yfrom yto])
subplot(2,2,2)
for i=1:syncCases-1
    if size(dataAll(i).track,1)>99
selfAngleAll2(i,:)=(dataAll(i).anglesSelf(1:100,1,2));
    end
end
boxplot(selfAngleAll2)
title('hornet side selfangle')
xlim([xfrom xto])
ylim([yfrom yto])
subplot(2,2,3)
for i=1:syncCases-1
    if size(dataAll(i).track,1)>99
selfAngleAll3(i,:)=(dataAll(i).anglesSelf(1:100,2,1));
    end
end
boxplot(selfAngleAll3)
title('bee top selfangle')
xlim([xfrom xto])
ylim([yfrom yto])
subplot(2,2,4)
for i=1:syncCases-1
    if size(dataAll(i).track,1)>99
selfAngleAll4(i,:)=(dataAll(i).anglesSelf(1:100,2,2));
    end
end
boxplot(selfAngleAll4)
title('bee side selfangle')
xlim([xfrom xto])
ylim([yfrom yto])
saveas(gcf,'preyAnglesOverTime.png')

% 
% figure
% plot(nanmean((angels')))

save('workspace.mat')




%folder=pwd;
%text(1,1,folder(end-21:end))


