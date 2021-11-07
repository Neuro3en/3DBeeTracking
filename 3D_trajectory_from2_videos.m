% Semi Automatic Hornet Bee Interaction Tracker

%% PRESS F5 AND FOLLOW COMANDS ON SCREEN

figure % tell MatLab what video you want to analyse
plot(1)
title({'first (SIDE) video? left click!';...
    'second (TOP) video? right click!'})
[x, y, buttom] = ginput(1);
videoNr=buttom;
close
[file,path] = uigetfile('*.mp4');
videoText=(file);
video_top=importdata(videoText);
video_top_size=size(video_top,4);
background=squeeze(median(video_top(:,:,:,:),4));
hornet=nan(video_top_size,2);
bee=hornet;
beeH=hornet;
hornetH=hornet;
intervalls=10;
BW_threshold=50;
centroids=nan(100,2,video_top_size);
Area=nan(100,video_top_size);
Orientation=nan(100,video_top_size);
for i=1:video_top_size
    BG_sub_frame=background(:,:,1)-video_top(:,:,1,i);
    BW= BG_sub_frame > BW_threshold;
    s = regionprops(BW,'centroid','Area','Orientation');
    centroids(1:size(cat(1,s.Centroid),1),:,i) = cat(1,s.Centroid);
    Area(1:size(cat(1,s.Centroid),1),i) = cat(1,s.Area);
    Orientation(1:size(cat(1,s.Centroid),1),i) = cat(1,s.Orientation);
end

%% start tracks for hornet and bee
hornetNew=nan(video_top_size,2,2);  % frames, hornet/bee,  X/Y
figure
imshow(video_top(:,:,1,1))
hold on
for i=5:video_top_size
    plot(centroids(:,1,i),centroids(:,2,i),'b.')
end
plot(centroids(:,1,1),centroids(:,2,1),'r*')
plot(centroids(:,1,1),centroids(:,2,1),'r*')
plot(centroids(:,1,1+1),centroids(:,2,1+1),'r*')
plot(centroids(:,1,1+2),centroids(:,2,1+2),'r*')
plot(centroids(:,1,1+3),centroids(:,2,1+3),'r*')
plot(centroids(:,1,1+4),centroids(:,2,1+4),'r*')
title('left click hornet, right click bee. Red points are at start frame')
[x, y]=ginput(2);
plot(x,y,'y*')
pause(.5)
close

[index, distance]=knnsearch(centroids(:,:,1),[x,y]);
for animal = 1:2
    if distance(animal)< 20
        hornetNew(1,animal,:)=centroids(index(animal),:,1);
        hornetNewO(1,animal)=Orientation(index(animal),1);
    else
        hornetNew(1,animal,:)=[x(animal), y(animal)];
    end
end
for i=2:video_top_size
    [index, distance]=knnsearch(centroids(:,:,i),squeeze(hornetNew(i-1,:,:)));
    for animal = 1:2
        if distance(animal)< 20
            hornetNew(i,animal,:)=centroids(index(animal),:,i);
            hornetNewO(i,animal)=Orientation(index(animal),i);
        else
            hornetNew(i,animal,:)=[x(animal), y(animal)];
        end
    end
end

%% attach hornet and bee tracks axis 1
figure('units','normalized','outerposition',[0 0 1 1])
plot(squeeze(centroids(1,1,:)),'k.')
hold on
for i=1:100
    plot(squeeze(centroids(i,1,:)),'k.')
end
plot(hornetNew(:,:,1),'*')
xlim([-10 size(centroids,3)+10])
ylim([-10 max(max(centroids(:,1,:)))+10])
buttom=16;
while buttom ~= 32
    if size(x,1)==1
        if buttom==1
            animal=1;
        elseif buttom==3
            animal=2;
        end
        if buttom < 4
            [index, distance]=knnsearch([squeeze(centroids(:,1,round(x))), zeros(100,1)+x],[y,x]);
            hornetNew(round(x),animal,:)=squeeze(centroids(index,:,round(x)))';
            hornetNewO(round(x),animal)=squeeze(Orientation(index,round(x)))';
            for i=round(x)+1:video_top_size
                [index, distance]=knnsearch(centroids(:,:,i),squeeze(hornetNew(i-1,animal,:))');
                if distance< 50
                    hornetNew(i,animal,:)=centroids(index,:,i);
                    hornetNewO(i,animal)=Orientation(index,i);
                end
            end
        end
        plot(hornetNew(:,animal,1),'*')
        title('find hornet(left click) and bee (right click) start early; quit with spacebar')
    end
    [x, y, buttom] = ginput(1);
end
close

%% attach hornet and bee tracks axis 2
figure('units','normalized','outerposition',[0 0 1 1])
plot(squeeze(centroids(1,2,:)),'k.')
hold on
for i=1:100
    plot(squeeze(centroids(i,2,:)),'k.')
end
plot(hornetNew(:,:,2),'*')
xlim([-10 size(centroids,3)+10])
ylim([-10 max(max(centroids(:,2,:)))+10])
buttom=16;
while buttom ~= 32
    if size(x,1)==1
        if buttom==1
            animal=1;
        elseif buttom==3
            animal=2;
        end
        if buttom < 4
            [index, distance]=knnsearch([squeeze(centroids(:,2,round(x))), zeros(100,1)+x],[y,x]);
            hornetNew(round(x),animal,:)=squeeze(centroids(index,:,round(x)))';
            for i=round(x)+1:video_top_size
                [index, distance]=knnsearch(centroids(:,:,i),squeeze(hornetNew(i-1,animal,:))');
                % for animal = 1:2
                if distance < 50
                    hornetNew(i,animal,:)=centroids(index,:,i);
                    hornetNewO(i,animal)=Orientation(index,i);
                    %    else
                    %        hornetNew(i,animal,:)=[x(animal), y(animal)];
                    % vor zur?ck . ginput....
                end
            end
        end
        plot(hornetNew(:,animal,2),'*')
        title('find hornet(left click) and bee (right click) start early; quit with spacebar')
    end
    [x, y, buttom] = ginput(1);
end
close

%% sort out false positive HORNET X
figure
buttom=2;
animal=1;
while buttom ~= 32
    plot(hornetNew(:,animal,1),'.')
    xlim([-20 20+video_top_size])
    title('HORNET trash rectangle: upper left to lower right, quit with 2x spacebar')
    [x, y, buttom]=ginput(2);
    x(x>video_top_size)=video_top_size;
    x(x<0)=1;
    if buttom(1)==1 && buttom(2)==1
        for i=round(x(1)):round(x(2))
            if hornetNew(i,animal,1)<y(1) && hornetNew(i,animal,1)>y(2)
                hornetNew(i,animal,:)=nan;
                hornetNewO(i,animal)=nan;
            end
        end
    end
end
close

%% sort out false positive HORNET Y
figure
buttom=2;
while buttom ~= 32
    plot(hornetNew(:,animal,2),'.')
    xlim([-20 20+video_top_size])
    title('HORNET trash rectangle: upper left to lower right, quit with 2x spacebar')
    [x, y, buttom]=ginput(2);
    x(x>video_top_size)=video_top_size;
    x(x<0)=1;
    if buttom(1)==1 && buttom(2)==1
        for i=round(x(1)):round(x(2))
            if hornetNew(i,animal,2)<y(1) && hornetNew(i,animal,2)>y(2)
                hornetNew(i,animal,:)=nan;
                hornetNewO(i,animal)=nan;
            end
        end
    end
end
close

%% sort out false positive BEE X
figure
buttom=2;
animal=2;
while buttom ~= 32
    plot(hornetNew(:,animal,1),'.')
    xlim([-20 20+video_top_size])
    title('BEE trash rectangle: upper left to lower right, quit with 2x spacebar')
    [x, y, buttom]=ginput(2);
    x(x>video_top_size)=video_top_size;
    x(x<0)=1;
    
    if buttom(1)==1 && buttom(2)==1
        for i=round(x(1)):round(x(2))
            if hornetNew(i,animal,1)<y(1) && hornetNew(i,animal,1)>y(2)
                hornetNew(i,animal,:)=nan;
                hornetNewO(i,animal)=nan;
            end
        end
    end
end
close

%% sort out false positive BEE Y
figure
buttom=2;
while buttom ~= 32
    plot(hornetNew(:,animal,2),'.')
    xlim([-20 20+video_top_size])
    title('BEE trash rectangle: upper left to lower right, quit with 2x spacebar')
    [x, y, buttom]=ginput(2);
    x(x>video_top_size)=video_top_size;
    x(x<0)=1;
    if buttom(1)==1 && buttom(2)==1
        for i=round(x(1)):round(x(2))
            if hornetNew(i,animal,2)<y(1) && hornetNew(i,animal,2)>y(2)
                hornetNew(i,animal,:)=nan;
                hornetNewO(i,animal)=nan;
            end
        end
    end
end 
close

hornetNewOLD=hornetNew; % to check if the gab filling is fine
for animal=1:2
for i=1:length(hornetNewO)
    if hornetNewO(i,animal)==0
hornetNewO(i,animal)=hornetNewO(i,animal)+rand()*.1;
    end
end
end

% fill gabs
for achsis=1:2
    for animal=1:2
        n=1;
        while n<video_top_size-1
            n=n+1;
            if isnan(hornetNew(n-1,animal,achsis))==0 && isnan(hornetNew(n,animal,achsis))==1
                start=n;
                while isnan(hornetNew(n,animal,achsis))==1 && n<video_top_size-1
                    n=n+1;
                end
                if n+1~=video_top_size
                    if (hornetNew(n,animal,achsis)-hornetNew(start-1,animal,achsis)) == 0
                        hornetNew(start:n-1,animal,achsis)=hornetNew(start-1,animal,achsis);
                        hornetNewO(start:n-1,animal)=     hornetNewO(start-1,animal);
                    else
                    hornetNew(start:n-1,animal,achsis)=[hornetNew(start-1,animal,achsis):...
                        (hornetNew(n,animal,achsis)-hornetNew(start-1,animal,achsis))/((n-1)-start):...
                        hornetNew(n,animal,achsis)];
                    hornetNewO(start:n-1,animal)=[hornetNewO(start-1,animal):...
                        (hornetNewO(n,animal)-hornetNewO(start-1,animal))/((n-1)-start):...
                        hornetNewO(n,animal)];
                    end
                end
            end
        end
    end
end


%%  ONLY THE FIRST / SIDE VIDEO
if videoNr==1
    hornetNewSide=hornetNew;
    hornetNewSideO=hornetNewO;
    centroidsSide=centroids;
else
    
    %% sync that data
    
    figure
    plot(hornetNew(:,1,1))
    hold on
    plot(hornetNewSide(:,1,1))
    title('sync with arrow keys, exit with space bar')
    
    shiftethornetNewSide=hornetNewSide(:,1,1);
    buttom=1;
    while buttom<32
        [x, y, buttom]=ginput(1);
        if buttom==29
            shiftethornetNewSide=cat(1,nan,nan,nan,nan,nan, shiftethornetNewSide);
            plot(hornetNew(:,1,1))
            hold on
            plot(shiftethornetNewSide)
            hold off
        end
        if buttom==30
            shiftethornetNewSide=shiftethornetNewSide+10;
            plot(hornetNew(:,1,1))
            hold on
            plot(shiftethornetNewSide)
            hold off
        end
        if buttom==28
            shiftethornetNewSide=shiftethornetNewSide(5:end);
            plot(hornetNew(:,1,1))
            hold on
            plot(shiftethornetNewSide)
            hold off
        end
        if buttom==31
            shiftethornetNewSide=shiftethornetNewSide-10;
            plot(hornetNew(:,1,1))
            hold on
            plot(shiftethornetNewSide)
            hold off
        end
        title('sync with arrow keys, exit with space bar')
    end
    close
    
    delayTraks=size(shiftethornetNewSide,1)-size(hornetNewSide,1);
    if delayTraks<0
        hornetNewSide=hornetNewSide(abs(delayTraks):end,:,:);
        hornetNewSideO=hornetNewSideO(abs(delayTraks):end,:);
    end
    if delayTraks>0
        hornetNew=hornetNew(abs(delayTraks):end,:,:);
        hornetNewO=hornetNewO(abs(delayTraks):end,:);
    end
    
    endDiffTraks=size(hornetNewSide,1)-size(hornetNew,1);
    if endDiffTraks>0
        hornetNewSide=hornetNewSide(1:-abs(endDiffTraks)+end,:,:);
        hornetNewSideO=hornetNewSideO(1:-abs(endDiffTraks)+end,:);
    end
    if endDiffTraks<0
        hornetNew=hornetNew(1:-abs(endDiffTraks)+end,:,:);
        hornetNewO=hornetNewO(1:-abs(endDiffTraks)+end,:);
    end
    
    hornetNew(:,:,2)=-hornetNew(:,:,2)+720;
    hornetNewSide(:,:,2)=-hornetNewSide(:,:,2)+720;
    
    
    %% 3D tracks Hornet Bule
    
    figure % hornet
    for i=1:length(hornetNewSide)
        plot3([squeeze(hornetNew(i,1,2)) squeeze(hornetNew(i,1,2))-9*sin(deg2rad(+hornetNewO(i,1)))*cos(deg2rad(+hornetNewSideO(i,1)))],...
            [squeeze(hornetNew(i,1,1)) squeeze(hornetNew(i,1,1))-9*cos(deg2rad(hornetNewO(i,1)))*cos(deg2rad(+hornetNewSideO(i,1)))],...
            [squeeze(hornetNewSide(i,1,2)) squeeze(hornetNewSide(i,1,2))-9*sin(deg2rad(+hornetNewSideO(i,1)))],'b')
        hold on
    end
    % bee
    for i=1:length(hornetNewSide)
        plot3([squeeze(hornetNew(i,2,2)) squeeze(hornetNew(i,2,2))+9*sin(deg2rad(+hornetNewO(i,2)))*cos(deg2rad(+hornetNewSideO(i,2)))],...
            [squeeze(hornetNew(i,2,1)) squeeze(hornetNew(i,2,1))+9*cos(deg2rad(+hornetNewO(i,2)))*cos(deg2rad(+hornetNewSideO(i,2)))],...
            [squeeze(hornetNewSide(i,2,2)) squeeze(hornetNewSide(i,2,2))+9*sin(deg2rad(hornetNewSideO(i,2)))],'r')
        hold on
    end
    
    % hive entrance rectangle
    hive=zeros(420,3);
    hive(1:2:end-1,3)=220;
    hive(2:2:end,3)=170;
    hive(:,2)=250;
    hive(1:2:end-1,1)=150:2:568;
    hive(2:2:end,1)=151:2:569;
    plot3(hive(:,1),hive(:,2),hive(:,3),'k')
    hold on
    
    plot3(hornetNew(:,1,2),hornetNew(:,1,1),hornetNewSide(:,1,2),'.')
    plot3(hornetNew(1,1,2),hornetNew(1,1,1),hornetNewSide(1,1,2),'bo','MarkerSize',6)
    plot3(hornetNew(:,2,2),hornetNew(:,2,1),hornetNewSide(:,2,2),'.')
    plot3(hornetNew(find(~isnan(hornetNew(:,2,2)),1),2,2),...
        hornetNew(find(~isnan(hornetNew(:,2,2)),1),2,1),...
        hornetNewSide(find(~isnan(hornetNew(:,2,2)),1),2,2),'ro','MarkerSize',6)
    set(gca, 'YDir','reverse')
    
    xlim([0 720])
    xlabel('X depth from top')
    ylim([0 1200])
    ylabel('Y distance Hive shared')
    zlim([0 720])
    zlabel('Z hight from side')
    
    title('3D - hornet: blue, bee: red')
    % saveas(gcf,[VideoTextChar(1:end-4),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)), '_3D.svg'])
    % vector grafics forpublication
    %saveas(gcf,[VideoTextChar(1:end-4),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)), '_3D.jpg'])
    saveas(gcf,'3D.jpg')
    
    figure
    subplot(2,2,1)
    plot(hornetNewSide(:,1,1),'.')
    xlim([0 video_top_size])
    title('hornet x-time')
    subplot(2,2,2)
    plot(hornetNewSide(:,1,2),'.')
    xlim([0 video_top_size])
    title('hornet y-time')
    subplot(2,2,3)
    plot(hornetNew(:,1,1),'.')
    xlim([0 video_top_size])
    title('hornet x-time 2.cam')
    subplot(2,2,4)
    plot(hornetNew(:,1,2),'.')
    xlim([0 video_top_size])
    title('hornet y-time 2.cam')
    
    figure
    subplot(2,2,1)
    plot(hornetNewSide(:,2,1),'.')
    xlim([0 video_top_size])
    title('bee x-time')
    subplot(2,2,2)
    plot(hornetNewSide(:,2,2),'.')
    xlim([0 video_top_size])
    title('bee y-time')
    subplot(2,2,3)
    plot(hornetNew(:,2,1),'.')
    xlim([0 video_top_size])
    title('bee x-time 2.cam')
    subplot(2,2,4)
    plot(hornetNew(:,2,2),'.')
    xlim([0 video_top_size])
    title('bee y-time 2.cam')
    
    
    %% velocity plot
    distanceChangeH=((sqrt((diff(hornetNewSide(:,1,1)).^2+diff(hornetNewSide(:,1,2)).^2)+diff(hornetNew(:,1,2)).^2)));
    distanceChangeB=((sqrt((diff(hornetNewSide(:,2,1)).^2+diff(hornetNewSide(:,2,2)).^2)+diff(hornetNew(:,2,2)).^2)));
    % distanceChangeHsmoothed = filloutliers(distanceChangeH,'linear','movmedian',5,'ThresholdFactor',.8); % only works in MatLab2019a
    % distanceChangeBsmoothed = filloutliers(distanceChangeB,'linear','movmedian',5,'ThresholdFactor',.8); % only works in MatLab2019a
    distanceChangeHsmoothed = distanceChangeH;
    distanceChangeBsmoothed = distanceChangeB;
    
    figure
    plot(distanceChangeHsmoothed,'.')
    hold on
    plot(distanceChangeBsmoothed,'.')
    title('velocity over time Hornet(blue) Bee(red)| distance(yellow)')
    
    A=hornetNewSide(:,1,1)-hornetNewSide(:,2,1); % distance plot
    B=hornetNewSide(:,1,2)-hornetNewSide(:,2,2);
    C=hornetNew(:,1,2)-hornetNew(:,2,2);
    plot(sqrt(A.^2+B.^2+C.^2)*0.02,'.')   % *0.02 plot factor 
    saveas(gcf,'DYN.jpg')
    
    
    %% percieved angle, hornets View
    figure
    hornetsAngle=nan(length(hornetNew),2); % top angle, side angle
    for i=1:length(hornetNew)
        hornetsAngle(i,1)=rad2deg(atan((hornetNew(i,2,2)-hornetNew(i,1,2))./(    hornetNew(i,2,1)-    hornetNew(i,1,1))))-    hornetNewO(i,1);
        hornetsAngle(i,2)=rad2deg(atan((hornetNewSide(i,2,2)-hornetNewSide(i,1,2))./sqrt((hornetNew(i,2,2)-hornetNew(i,1,2))^2+(hornetNew(i,2,1)-hornetNew(i,1,1))^2)))-hornetNewSideO(i,1);
    end
    for i=1:length(hornetsAngle) -1
        if abs((hornetsAngle(i,1)+hornetNewO(i,1))-(hornetsAngle(i+1,1)+hornetNewO(i+1,1)))>160 && (hornetsAngle(i,1)+hornetNewO(i,1))>0
            hornetsAngle(i+1:end,1)=hornetsAngle(i+1:end,1)+180;
        elseif abs((hornetsAngle(i,1)+hornetNewO(i,1))-(hornetsAngle(i+1,1)+hornetNewO(i+1,1)))>160 && (hornetsAngle(i,1)+hornetNewO(i,1))<0
            hornetsAngle(i+1:end,1)=hornetsAngle(i+1:end,1)-180;
        end
    end
    c=sqrt(A.^2+B.^2+C.^2)*-0.02;
    scatter((hornetsAngle(:,1)),hornetsAngle(:,2),[],c,'filled','MarkerFaceAlpha',.8)
    hold on
    plot((-90:1:90),(zeros(181,1)),'k.')
    plot((zeros(181,1)),(-90:1:90),'k.')
    title('inbound bees from hornets View in °')
    xlabel('color is distance, warm is close, cold is far')
    colormap jet
    set(gca, 'XDir','reverse')
    saveas(gcf,'HornetView.jpg')
    
    
    %% percieved angle, BEE View
    figure
    BeeAngle=nan(length(hornetNew),2); % top angle, side angle
    for i=1:length(hornetNew)
        BeeAngle(i,1)=rad2deg(atan((hornetNew(i,1,2)-hornetNew(i,2,2))./(    hornetNew(i,1,1)-    hornetNew(i,2,1))))-    hornetNewO(i,2);
        BeeAngle(i,2)=rad2deg(atan((hornetNewSide(i,1,2)-hornetNewSide(i,2,2))./sqrt((hornetNew(i,1,2)-hornetNewSide(i,2,2))^2+(hornetNew(i,1,1)-hornetNew(i,2,1))^2)))-hornetNewSideO(i,2);
    end
    for i=1:length(BeeAngle)-1
        if abs((BeeAngle(i,1)+hornetNewO(i,1))-(BeeAngle(i+1,1)+hornetNewO(i+1,1)))>160 && (BeeAngle(i,1)+hornetNewO(i,1))>0
            BeeAngle(i+1:end,1)=BeeAngle(i+1:end,1)+180;
        elseif abs((BeeAngle(i,1)+hornetNewO(i,1))-(BeeAngle(i+1,1)+hornetNewO(i+1,1)))>160 && (BeeAngle(i,1)+hornetNewO(i,1))<0
            BeeAngle(i+1:end,1)=BeeAngle(i+1:end,1)-180;
        end
    end
    A=hornetNewSide(:,1,1)-hornetNewSide(:,2,1); % distance plot
    B=hornetNewSide(:,1,2)-hornetNewSide(:,2,2);
    C=hornetNew(:,1,2)-hornetNew(:,2,2);
    c=sqrt(A.^2+B.^2+C.^2)*-0.02;
    scatter((BeeAngle(:,1)+90),BeeAngle(:,2),[],c,'filled','MarkerFaceAlpha',.8)
    hold on
    plot((-90:1:90),(zeros(181,1)),'k.')
    plot((zeros(181,1)),(-90:1:90),'k.')
    title('hornet from bees View in °')
        xlabel('color is distance, warm is close, cold is far')
    colormap jet
    set(gca, 'XDir','reverse')
    saveas(gcf,'BeeView.jpg')
    
    
    %% save VARIABLES
    % video file name _ month _ day _ hour _ min
    %save([VideoTextChar(1:end-4),'_',num2str(c(2)),'_',num2str(c(3)),'_',num2str(c(4)),'_',num2str(c(5)), '.mat'],'hornetNewSide','hornetNew')
    save('data2.mat','hornetNewSide','hornetNewSideO','hornetNew', 'hornetNewO','BeeAngle','hornetsAngle')
end
