function level = isodata(I)
    %   ISODATA Compute global image threshold using iterative isodata method.
    %   LEVEL = ISODATA(I) computes a global threshold (LEVEL) that can be
    %   used to convert an intensity image to a binary image.
    %   This iterative technique for choosing a threshold was developed by Ridler and Calvard .
    %   The histogram is initially segmented into two parts using a starting threshold value such as 0 = 2B-1, 
    %   half the maximum dynamic range. 
    %   The sample mean (mf,0) of the gray values associated with the foreground pixels and the sample mean (mb,0) 
    %   of the gray values associated with the background pixels are computed. A new threshold value 1 is now computed 
    %   as the average of these two sample means. The process is repeated, based upon the new threshold, 
    %   until the threshold value does not change any more.   
    %
    %   Class Support
    %   -------------
    %   The input image I can be of class uint8, uint16, or double and it
    %   must be nonsparse.  LEVEL is a double scalar.
    %
    %   Example
    %   -------
    %       I = imread('blood1.tif');
    %       level = graythresh(I);
    %       BW = im2bw(I,level);
    %       imshow(BW)
    %
    %   See also IM2BW.
    %
    % Reference :T.W. Ridler, S. Calvard, Picture thresholding using an iterative selection method, 
    %            IEEE Trans. System, Man and Cybernetics, SMC-8 (1978) 630-632.

    % Convert all N-D arrays into a single column.  Convert to uint8 for
    % fastest histogram computation.
    I = I(:);

    % STEP 1: Compute mean intensity of image from histogram, set T=mean(I)
    [counts,N0]=histcounts(I);
    for ii = 2:length(N0)
        N(ii-1) = mean([N0(ii-1) N0(ii)]);
    end

    i=1;
    mu=cumsum(counts);
    T(i)=(sum(N.*counts))/mu(end);
    T(i)=round(T(i));

    % STEP 2: compute Mean above T (MAT) and Mean below T (MBT) using T from
    % step 1
    T_ind = find(abs(N-T(i))==min(abs(N-T(i))));
    mu2=cumsum(counts(1:T_ind));
    MBT=sum(N(1:T_ind).*counts(1:T_ind))/mu2(end);

    mu3=cumsum(counts(T_ind(i):end));
    MAT=sum(N(T_ind:end).*counts(T_ind:end))/mu3(end);
    i=i+1;
    % new T = (MAT+MBT)/2
    T(i)=round((MAT+MBT)/2);

    % STEP 3 to n: repeat step 2 if T(i)~=T(i-1)
    level=[];
    while abs(T(i)-T(i-1))>=1
        T_ind = find(abs(N-T(i))==min(abs(N-T(i))));

        mu2=cumsum(counts(1:T_ind));
        MBT=sum(N(1:T_ind).*counts(1:T_ind))/mu2(end);

        mu3=cumsum(counts(T_ind:end));
        MAT=sum(N(T_ind:end).*counts(T_ind:end))/mu3(end);

        i=i+1;
        T(i)=round((MAT+MBT)/2); 
        level=T(i);
    end   
end