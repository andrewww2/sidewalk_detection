%Implenting this paper: http://wscg.zcu.cz/wscg2014/Short/N19-full.pdf
%close all;
fName = 'sidewalk_pics/play_2.jpg';
img = imread(fName);
%img = imresize(img, 0.05); %make it smaller
% imwrite(img, fName);
% 2) convert to grayscale
gray = rgb2gray(img);      
% 3) and 4) histogram equalization and median filter
gray_eq = medfilt2(histeq(gray), [9 9]);  
% 5) Canny edge detector
[I_m, I_a] = CannyEdgeDetector(gray);

%Orientation consistency pass filtering (OCPF)
% I_m = I_m(size(I_m,1)/2:end,:);
% I_a = I_a(size(I_a,1)/2:end,:);
thresh_ang = 20; %in degrees
W = find(I_m); %precompute set W -> all nonzero indexes
I_f = zeros(size(I_m)); %output OCPF image
%I_a = abs(I_a);
for p = (1:size(W)) %make sure I_m(x,y) > 0
    ind = W(p);
    n_W = length(W) - 1; %size of W. exclude current element
    setW = I_a(W); 
    setW(p) = []; %exclude current element from W
    diff = I_a(ind) - setW;  
    A_v = sum(diff) / n_W; %angle variation
    
    if (A_v < thresh_ang)
        I_f(ind) = 1;
    else
        I_f(ind) = 0;
    end
end
BW = bwareaopen(I_f,40); %Remove all object containing fewer than 40 pixels
subplot(2,3,1), imshow(img), title('Original');
% subplot(2,3,2), imshow(gray_eq), title('Blurred');
subplot(2,3,2), imshow(I_m), title('Canny Magnitude');
subplot(2,3,3), imshow(I_f), title('After OCPF');
subplot(2,3,4), imshow(BW), title('Remove noise');

% Hough lines
numpeaks = 10;
[H,theta,rho] = hough(BW);
peaks = houghpeaks(H,numpeaks);
lines = houghlines(BW,theta,rho,peaks,'FillGap',5,'MinLength',15);
subplot(2,3,5), imshow(BW), title('Hough'), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   %find equation of line
   m = (xy(4)-xy(3))/(xy(2)-xy(1));
   x = 1:size(img,2);
   b = xy(4) - m*xy(2);
   y = m*x+b;
    plot(x,y, 'r');
    %plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%    % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end