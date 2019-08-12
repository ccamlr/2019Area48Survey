function crop_image(imageFilename)
   % Trim 'whitessace' from the outer boundaries of an image. Intended for
   % removing the often large borders around Matlab-produced figure files.
   % Assumes that the top righthand pixel is the colour of the border that
   % is to be removed.
   %
   % Overwrites the image file with the cropped version.

   A = imread(imageFilename);
   
   bgColour = A(1,1,:);
   
   % trim from top
   topTrim = 1;
   while A(topTrim,:,:) == bgColour & topTrim < size(A, 1)
       topTrim = topTrim + 1;
   end
   bttmTrim = 0;
   while A(end-bttmTrim,:,:) == bgColour & bttmTrim < size(A, 1)
       bttmTrim = bttmTrim + 1;
   end
   leftTrim = 1;
   while A(:,leftTrim,:) == bgColour & leftTrim < size(A, 2)
       leftTrim = leftTrim + 1;
   end
   rightTrim = 0;
   while A(:,end-rightTrim,:) == bgColour & rightTrim < size(A, 2)
       rightTrim = rightTrim + 1;
   end
   
   A = A(topTrim:end-bttmTrim,leftTrim:end-rightTrim,:);
   
   imwrite(A, imageFilename);
   
end


