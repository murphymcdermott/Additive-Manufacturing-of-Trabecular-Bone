function binary_Image = function_Otsu(gray_image)

%% Level 1 Otsu method
% Apply the first Otsu threshold to separate background (background and
% paraffin) from foreground (bone marrow and trabecular bone)
level1 = graythresh(gray_image); % Otsu's method to find the first global threshold
background = imbinarize(gray_image, level1); % Segmenting the image into foreground and background

% Isolate the foreground from the initial image
foreground = gray_image;
foreground(~background) = 0; % Set background pixels to 0
image_Improved = foreground;

%% Level 2 Otsu Method
% Apply the second Otsu threshold to separate trabecular bone (foreground) from the bone marrow (background)
level2 = graythresh(image_Improved(image_Improved > 0)); % Otsu's method on non-background pixels

if level2<0.35 
    level2 = 0.35; % This threshold is used to consider as background the very dark grey gradients found on the top and bottom of the sample 
end

background2 = imbinarize(image_Improved, level2); % Segmenting the image into bone marrow and trabecular bone

% Isolate the trabecular bone from the bone marrow
foreground2 = image_Improved;
foreground2(~background2) = 0; % Set bone marrow pixels to 0
image__further_Improved = foreground2;
binary_Image = imbinarize(image__further_Improved);

% Assess the foreground pixel amount to decide if valid bone is detected
foregroundPixelCount = sum(binary_Image(:));
if foregroundPixelCount < 100 % Example threshold, adjust based on your dataset
    binary_Image = false(size(binary_Image)); % Set all pixels to 0 if below threshold
end

%% %%%%%%%%%%%%%%%%%%%     Visualization code     %%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,2,1)
imshow(gray_image); title("Original Image")
subplot(2,2,2)
imshow(image_Improved) ; title(["Original image with the removed background ","and paraffin","Otsu threshold value:" num2str(level1)])
subplot(2,2,3)
imshow(image__further_Improved) ; title(["Original Image with the background", "paraffin and bone marrow removed","Otsu threshold value:" num2str(level2)]);
subplot(2,2,4)
imshow(binary_Image) ; title("Binary image")

end