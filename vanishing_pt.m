close all;
sw1 = imread('sw_1.jpg');
sw1 = imresize(sw1, 0.1); %make it smaller

numpeaks = 10;
imshow(sw1);
gray = rgb2gray(sw1);
BW = imgaussfilt(gray,'FilterSize',151); %(gray, fspecial('gaussian'));
imshow(BW);
BW = edge(BW, 'Canny');

BW = bwareaopen(BW,75); %Remove all object containing fewer than 30 pixels

[H,theta,rho] = hough(BW);
peaks = houghpeaks(H,numpeaks);
lines = houghlines(BW,theta,rho,peaks,'FillGap',5,'MinLength',7);
figure, imshow(BW), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

%    % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
%    % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end