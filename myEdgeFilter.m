function [Im, Io, Ix, Iy] = myEdgeFilter(img, sigma)
%Your implemention
    ySobel = fspecial('sobel'); %what's the size of sobel filter?
    xSobel = ySobel';
    
    hsize = 2 * ceil(3 * sigma) + 1;
    hGauss = fspecial('gaussian',hsize,sigma);
    img = imfilter(img, hGauss);
%     imshow(imgGauss)
    Ix = double(imfilter(img, xSobel));
    Iy = double(imfilter(img, ySobel));
    Io = atan2d(Iy,Ix); %orientation matrix
    Im = sqrt(Ix.^2 + Iy.^2); %magnitude
    ImCopy = Im; %make a copy
    %imshow(Im)
        
    %nonmaximum suppression - loop over orientation matrix
    for r = 2:(size(img,1) - 1)
        for c = 2:(size(img,2) - 1)
            angle = Io(r,c);
            if (angle < 0)
                angle = angle + 180;
            end
            
            %map to 0 degrees
            if(angle < 22.5) || (angle > 157.5)
                if (Im(r,c) < Im(r, c + 1)) || (Im(r,c) < Im(r, c - 1))
                    ImCopy(r,c) = 0;
                end
            %map to 45 degrees
            elseif(angle >= 22.5) && (angle < 67.5)
                %disp("Hello")
                if (Im(r,c) < Im(r - 1, c - 1)) || (Im(r,c) < Im(r + 1, c + 1))
                    ImCopy(r,c) = 0;
                end
            %map to 90 degrees
            elseif(angle >= 67.5) && (angle < 112.5)
                if (Im(r,c) < Im(r + 1, c)) || (Im(r,c) < Im(r - 1, c))
                    ImCopy(r,c) = 0;
                end
            %map to 135 degrees
            elseif(angle >= 112.5) && (angle < 157.5)
                if (Im(r,c) < Im(r + 1, c - 1)) || (Im(r,c) < Im(r - 1, c + 1))
                    ImCopy(r,c) = 0;
                end
            end
        end
    end
    Im = ImCopy;
end
    
                
        
        
