function [img1] = myImageFilter(img0, h)
    filterH = size(h,1);
    filterW = size(h,2);

    %pad image
    padH = floor(filterH/2);
    padW = floor(filterW/2);
    padImg0 = double(padarray(img0, [padH padW], 'replicate'));
    
    %flip h horizontally and vertically
    h = fliplr(h);
    flipH = double(flipud(h));
    
    %possible vectorization functions: im2col, ind2sub
    img1 = zeros(size(img0,1), size(img0,2));  %init outp to same size as original image
    for r =  (1 + padH):(size(padImg0,1) - padH) %height of image
        for c = (1 + padW):(size(padImg0,2) - padW) %width of image
            conv = flipH.*padImg0((r-padH):(r+padH), (c-padW):(c+padW)); %element-wise multiply
            img1(r - padH,c - padW) = sum(sum(conv));
        end
    end
end