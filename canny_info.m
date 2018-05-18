% my own implementation of Canny edge detector
% based on: http://justin-liang.com/tutorials/canny/

%out_img = final filtered image
%G_mag = magnitude gradient
%G_ang = angle gradient
function [G_mag, G_ang] = canny_info(img)
    sigma = 0.5;
    gauss_img = imgaussfilt(img, sigma, 'filtersize', 5);
    G_y =  fspecial('sobel');
    G_x = G_y';
    x_conv = double(imfilter(gauss_img, G_x))/255;
    y_conv = double(imfilter(gauss_img, G_y))/255;
    G_mag = hypot(x_conv, y_conv);
    G_ang = atan(y_conv./x_conv);
    imshow(G_mag); 
    figure; imshow(G_ang);
end