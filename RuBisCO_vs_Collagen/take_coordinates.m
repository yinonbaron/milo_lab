clear all 
%%
im = imread('F5.large.jpg');
mega_hum = ((im(:,:,2)>150) & (im(:,:,1)<30));
livestock = (im(:,:,2)>100 & im(:,:,2)<150 & im(:,:,1)<20);

topY = 14;
bottomY = 906;
rightX = 1056;
leftX = 74;

mega_CC = bwconncomp(mega_hum,8);
size_mega_CC = cellfun(@(x) length(x),mega_CC.PixelIdxList);
mega_hum_center = regionprops(mega_CC,'centroid');
mega_hum_box = regionprops(mega_CC,'BoundingBox');
mega_hum_center = mega_hum_center(size_mega_CC>10);
livestock_CC = bwconncomp(livestock,8);
size_livestock_CC = cellfun(@(x) length(x),livestock_CC.PixelIdxList);
livestock_center = regionprops(livestock_CC,'centroid');
livestock_box = regionprops(livestock_CC,'BoundingBox');
livestock_center = livestock_center(size_livestock_CC>10);

livestock_points(1,:) = [74 808];
livestock_points(2,:) = [133 808];
livestock_points(3,:) = [192 813];
livestock_points(4,:) = [212 813];
livestock_points(5,:) = [236 844];
livestock_points(6,:) = [252 844];
livestock_points(7,:) = [271 853];
livestock_points(8,:) = [290 867];
livestock_points(9,:) = [314 865];
livestock_points(10,:) = [349 862];
livestock_points(11,:) = [374 856];
livestock_points(12,:) = [408 843];
livestock_points(13:17,:) = reshape([livestock_center(end-4:end).Centroid],[2 5])';
livestock_points(17,1) = rightX;

mega_hum_points(1,:) = [74 808];
mega_hum_points(2,:) = [133 808];
mega_hum_points(3,:) = [192 813];
mega_hum_points(4,:) = [212 813];
mega_hum_points(5,:) = [236 844];
mega_hum_points(6,:) = [252 844];
mega_hum_points(7,:) = [271 853];
mega_hum_points(8,:) = [290 867];
mega_hum_points(9,:) = [349 867];
mega_hum_points(10,:) = [349 869];
mega_hum_points(11,:) = [374 867];
mega_hum_points(12,:) = [408 863];
mega_hum_points(13:17,:) = reshape([mega_hum_center(end-4:end).Centroid],[2 5])';
mega_hum_points(17,1) = rightX;

mega_hum_corr_points(:,1) = 10.^(-(mega_hum_points(:,1)-rightX)/((rightX-leftX)/5));
livestock_corr_points(:,1) = 10.^(-(livestock_points(:,1)-rightX)/((rightX-leftX)/5));
mega_hum_corr_points(:,2) = 1.6-(mega_hum_points(:,2)-topY)/(bottomY-topY)*1.6;
livestock_corr_points(:,2) = 1.6-(livestock_points(:,2)-topY)/(bottomY-topY)*1.6;
%%
im2 = imread('Human_vs_megafauna.png');
mega_only_points = imfindcircles(im2(:,:,1)>100 & im2(:,:,2)<50,[2 20]);
hum_only_points = imfindcircles(im2(:,:,2)>100 & im2(:,:,1)<100,[5 20]);
%{
imshow(im2),hold on,
plot(mega_only_points(:,1),mega_only_points(:,2),'*g')
plot(hum_only_points(:,1),hum_only_points(:,2),'*g')
%}
topY2 = 124;
bottomY2 = 615;
rightX2 = 1296;
leftX2 = 150;

hum_only_corr_points(:,1) = 10.^(-(hum_only_points(:,1)-rightX2)/((rightX2-leftX2)/5));
mega_only_corr_points(:,1) = 10.^(-(mega_only_points(:,1)-rightX2)/((rightX2-leftX2)/5));
hum_only_corr_points(:,2) = 10.^(12-(hum_only_points(:,2)-topY2)/(bottomY2-topY2)*6);
mega_only_corr_points(:,2) = 10.^(12-(mega_only_points(:,2)-topY2)/(bottomY2-topY2)*3);
[B I] = sort(hum_only_corr_points(:,1));
hum_only_corr_points = hum_only_corr_points(I,:);
hum_only_corr_points = hum_only_corr_points([1:16 18],:);
[B I] = sort(mega_only_corr_points(:,1));
mega_only_corr_points = mega_only_corr_points(I,:);
mega_only_corr_points = mega_only_corr_points([1:5 7:18],:);

% Taken from Hern 1999
cattle_fraction = 663149950000/988565250000;
sheep_fraction = 77871734000/988565250000;
pig_fraction = 95835247000/988565250000;
horses_fraction = 53797867000/988565250000;
buffalo_camel_fraction = 71653831000/988565250000;
sheep_fraction = 77871734000/988565250000;
chicken_fraction = 26082000000/988565250000;
livestock_only_points = flipud(livestock_corr_points(:,2)-mega_hum_corr_points(:,2))*1e12;
cattle_points = livestock_only_points*cattle_fraction;
sheep_points = livestock_only_points*sheep_fraction;
pig_points = livestock_only_points*pig_fraction;
chicken_points = livestock_only_points*chicken_fraction;
other_points = livestock_only_points*(1-(cattle_fraction));
figure();
total_data = [mega_only_corr_points(:,2) hum_only_corr_points(:,2) cattle_points other_points];
area(mega_only_corr_points(:,1),total_data*6/1e12)
xlim([0 100])
set(gca,'Xdir','reverse')
xlabel('Years BP'), ylabel('10^{13} grams')
line([1 100000],[7.5 7.5],'LineWidth',2,'Color','k')

line([300 300],[0 10],'LineStyle','--','LineWidth',2,'Color','k')
text(2000, 7.2,'RuBisCO mass')
text(500, 0.1,'Begining of the industrial revolution','Rotation',90)
legend('Wild Megafauna','Humans','Cattle','Other Livestock')

% The y axis is total mass. to convert to mass of collagen we multiply by
% 0.2*0.3 = 0.06
%{
semilogx(mega_corr_points(:,1),mega_hum_corr_points(:,2),'.g')
hold on
%}
figure()
plot(livestock_corr_points(:,1),livestock_corr_points(:,2)*0.06*100,'or','MarkerFaceColor','r')
set(gca,'Xdir','reverse')
xlabel('Years BP'), ylabel('10^{13} grams')
line([1 100000],[7.5 7.5],'LineWidth',2,'Color','k')
xlim([0 10000])
line([300 300],[0 10],'LineStyle','--','LineWidth',2,'Color','k')
text(50000, 7.2,'RuBisCO mass')
text(500, 0.1,'Begining of the industrial revolution','Rotation',90)
legend('Total collagen')


