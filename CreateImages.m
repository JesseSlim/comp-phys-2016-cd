%% INIT
clear all;
close all;
clc;
format compact;

%% LOAD
Data = load('output.mat');
Orientation = Data.SaveOrientation;
Location    = Data.SaveLocation;
Radius       = Data.Radius * 2;

%% PLOT IMAGES
minXL  = min(min(Location(:,:,1)));
maxXL = max(max(Location(:,:,1)));
minYL =  min(min(Location(:,:,2)));
maxYL = max(max(Location(:,:,2)));

CenterXL = (maxXL + minXL) / 2;
CenterYL = (maxYL + minYL) / 2;
SizeL = max([maxXL - minXL, maxYL - minYL]) / 2 * 1.2;

parfor (i = 1:1:size(Location,1))
    figure('units','normalized','outerposition',[0 0 1 1], 'Visible','off')
    minX  = min(min(Location(i,:,1)));
    maxX = max(max(Location(i,:,1)));
    minY =  min(min(Location(i,:,2)));
    maxY = max(max(Location(i,:,2)));

    CenterX = (maxX + minX) / 2;
    CenterY = (maxY + minY) / 2;
    Size = max([maxX - minX, maxY - minY]) / 2 * 1.2;
    subplot(1,2,1);
    scatter(Location(i,:,1), Location(i,:,2));
    axis([minXL, maxXL, minYL, maxYL]);
    axis([CenterXL-SizeL, CenterXL+SizeL, CenterYL-SizeL,CenterYL+SizeL]);
    axis square;
    subplot(1,2,2),
    hold on;
    for j = 1:1:size(Location,2)
        rectangle('Position', [Location(i,j,1) Location(i,j,2) Radius(j) Radius(j)], 'Curvature', [1,1]);
        plot([Location(i,j,1) + Radius(j) / 2,Location(i,j,1) + Radius(j) / 2 + 0.5 * cos(Orientation(i,j))],[Location(i,j,2) + Radius(j) / 2 ,Radius(j) / 2 + Location(i,j,2) + 0.5 * sin(Orientation(i,j))], 'k', 'LineWidth', 1.5);
    end


    axis([CenterX-Size, CenterX+Size, CenterY-Size,CenterY+Size]);
    hold off;
    axis square
    fig = gcf();
    saveas(fig, sprintf('images/image%g.jpg', i));
    close fig;
end