function cord=Run_ALL(imagename,nColors,filename, lineTol, dilation_shape,dilation_size, autocontrast)
% input: lineTol is the distance between each adjacent output outline
% points. If the lineTol is big, the total number of points to draw a
% outline is smaller; otherwise, the number of points is bigger but the
% outline is more close to real edge of the shape.
% dilation is the shape, size of dilation you want to perform. 
% Autocontrast can have value 0 or 1. 0 means no auto-adjustment, whereas 1
% means to perform autoadjustment. 
% Suggest: use lineTol=10.
% Suggest: use dilation_shape='disk',dilation_size='5',autocontrast=1;
%%Note AUTOCONTRAST property automatically adjusts contrast of images to
%%optimum level. Acknowledgements: This part is inspired from Auto Contrast 
%%by Divakar Roy. Here is the link. http://www.mathworks.com/matlabcentral/fileexchange/10566-auto-contrast

%% Auto contrast
    if autocontrast==1
        min=0.008;
        max=0.992;
        myImage=imread(imagename);
        [row, col,dim]=size(myImage);
        myImage=double(myImage);
        for i=1:dim
            myArray=sort(reshape(myImage(:,:,i),row*col, 1));
            myMin(i)=myArray(ceil(min*row*col));
            myMax(i)=myArray(ceil(max*row*col));
        end
        if dim==3
            myMin=rgb2ntsc(myMin);
            myMax=rgb2ntsc(myMax);
        end
        myImage=(myImage-myMin(1))/(myMax(1)-myMin(1));
    else
        myImage=imread(imagename);
    end
    
    if nColors<3
        cc=BWseg(myImage,dilation_shape, dilation_size);
    else
        cc=k_edge(myImage, nColors, dilation_shape, dilation_size);
    end
    nrows=cc.ImageSize(1);
    ncols=cc.ImageSize(2);
    tol=input('What is the minimal area of each objects: ');
    cordcut=cut_coord(cc,tol);
    cordcut=minimize_points(cordcut);
    cord=reduce_points(cordcut,lineTol, nrows*ncols); 
    numArray2File(cord,filename);
    figure, imshow(myImage),title('Contour indicated by Green Lines');
    hold on;
    for i=1:length(cord)
        X=cord{i}(:,1);
        Y=cord{i}(:,2);
        plot(X,Y,'Color','green','LineWidth',2);
    end
end

%%Function Segmentation: BWseg--pre-process the image
% input: a black&white image or an image with only two different colors
% output: outline of the image
function outline=BWseg(myImage,dilation_shape, dilation_size)
bwImage=im2bw(myImage);
SE=strel(dilation_shape,dilation_size);
%bwImage=imopen(bwImage,SE);
bwImage=imerode(bwImage,SE);
BWedge=bwperim(bwImage);
% NOTE: you can also add imdiate & imerode combination to smooth the
% objects in order to get better results
outline=bwconncomp(BWedge);
% delete the first connected component, which is the edge of the whole
% picture
outline.PixelIdxList{1}=[];
end

%% Function Segmentation: k_edge--pre-process the image
% input: a color image with more than two different colors, along with
% numbers of different clusters/colors
% output:outline of the image
function outline=k_edge(myImage,nColors, dilation_shape, dilation_size)
% Step 1: convert image from RGB color space to L*a*b color space
cform=makecform('srgb2lab');
lab_myImage=applycform(myImage, cform);
ab=double(lab_myImage(:,:,2:3));
nrows=size(ab, 1);
ncols=size(ab, 2);

% Optional: fill holes
% ab=imfill(ab,'holes');
ab=reshape(ab,nrows*ncols,2);
% Use K-Means method to segment the image into desired clusters
[cluster_idx, ~]=kmeans(ab, nColors, 'distance', 'sqEuclidean',...
                                'Replicates', 5,'emptyaction','drop');
% Step 2: Label Every pixel in the Image using the results from KMEANS
pixel_labels=reshape(cluster_idx,nrows,ncols);

% Step 3: Create image that segment the original image by Color.
segmented_images=cell(1, nColors);
rgb_labels=repmat(pixel_labels, [1 1 3]);
for k=1:nColors
    color=myImage;
    color(rgb_labels~=k)=0;
    segmented_images{k}=color;
end

figure, subplot(nColors, 1, 1),imshow(segmented_images{1}), title('cluster 1');
for m=2:nColors
    subplot(nColors, 1, m), imshow(segmented_images{m}), title(sprintf('cluster %i', m));
end
% Step 4: Choose your desired cluster
cluster_interest=input('Which cluster are you interested in? ');
pixel_reshape=pixel_labels;
pixel_reshape(pixel_labels~=cluster_interest)=0;
SE=strel(dilation_shape,dilation_size);  % It depends on the approximate mean radius of the shapes
pixel_reshape=imerode(pixel_reshape,SE);
pixel_reshape=imdilate(pixel_reshape,SE);
% pixel_reshape=imdilate(pixel_reshape,SE);
% optional steps: further clean up the images and eliminate noises
% NOTE: THESE CLEAN-UP STEPS CAN POTENTIONALLY DISTORT THE IMAGE, BE
% CAREFUL WHEN DOING CLEAN UPS. IF ANY STEP IS NOT NEEDED, COMMENT THEM
% OUT. If you are satisfied by each step, enter 1; otherwise enter 0.
bwFill=imfill(pixel_reshape, 'holes');
figure, imshow(bwFill);
satisfy=input('Are you satisfied by this clean up? Yes 1, No 0: ');
if satisfy
    pixel_reshape=bwFill;
end
bwClean=bwmorph(pixel_reshape, 'clean'); % remove isolated pixels 
%(individual 1s that are surrounded by 0s).
figure, imshow(bwClean);
satisfy=input('Are you satisfied by this clean up? Yes 1, No 0: ');
if satisfy
    pixel_reshape=bwClean;
end
bwMaj=bwmorph(pixel_reshape,'majority'); % Sets a pixel to 1 if five or more 
% pixels in its 3-by-3 neighborhood are 1s; otherwise, it sets the pixel to 0.
figure, imshow(bwMaj);
satisfy=input('Are you satisfied by this clean up? Yes 1, No 0: ');

if satisfy
    pixel_reshape=bwMaj;
end
BWedge=edge(pixel_reshape);
figure, imshow(BWedge), title('Contours of selcted objects');
outline=bwconncomp(BWedge);
end

%% Function Contour-making: cut_coord--create the x,y coordinates and get ready to output into XML file
% input: outline from previous process. A object from bwconncomp() method.
% output: cell array. Each cell represents a connected component. Each row
% of a cell element is a x,y coordinates. The coordinates are in a specific
% order (goes around the object in clockwise or counter-clockwise
% direction. 
function cordcut=cut_coord(cc,tol)
cordcut=cell(cc.NumObjects, 1);
nrows=cc.ImageSize(1);
count=0;
% delete all the objects not qualified
for i=1:cc.NumObjects
    if length(cc.PixelIdxList{i})>=tol
        count=count+1;
        cordcut{count}=cc.PixelIdxList{i};
    end
end
cordcut=cordcut(1:count);
for i=1:count
        list=cordcut{i};
        y=rem(list, nrows); % y represents row number
        x=ceil(list/nrows); % x represents column number
        list=zeros(length(list)+1, 2);
        % calculate all the distances from the rest of points to the first.
        list(1,:)=[x(1) y(1)];
        prev_ind=zeros(length(x),1);
        ind=1;
        prev_ind(ind)=1;
        dist=(x-x(ind)).^2+(y-y(ind)).^2;
        ind=find(dist==min(dist(dist>0)),1);
        list(2,:)=[x(ind) y(ind)];
        prev_ind(ind)=1;
        prev_vector=list(2,:)-list(1,:);
        j=3;
       while(ind~=1&& j<length(list))
         dist=(x-x(ind)).^2+(y-y(ind)).^2;
           if j<= 0.9*length(x) || dist(1)~=min(dist>0)
                candidate= find(dist==min(dist(~prev_ind)));
                cand= candidate(prev_ind(candidate)==0);
                if ~isempty(cand)
                    new_vector=[x(cand)-x(ind) y(cand)-y(ind)];
                    dot_product= sum(repmat(prev_vector,length(cand), 1) .* new_vector, 2);
                    ind=cand(dot_product==max(dot_product));
                    ind=ind(1);
                    prev_ind(ind)=1;
                    list(j,:)=[x(ind) y(ind)];
                    j=j+1;
                else
                    ind=candidate;
                    list(j,:)=[x(ind) y(ind)]; 
                    j=j+1;
                    break;
                end
            else
               j=j+1;
               break;
            end
        end
         list(j,:)=[x(1) y(1)];
         cordcut{i}=list(1:j, :);
end
end

%% Function Contour-Writing: minimize_points-- reduce the number of points needed to form the known contour
% input: the ordered coordinates 
% output: reduced ordered coordinates
% Method: control the angles between each line segments
function newList=minimize_points(cordcut)
    newList=cell(length(cordcut),1);
    for i=1:length(cordcut)
        [nrows, ncols]=size(cordcut{i});
        list=zeros(nrows,ncols);
        list(1,:)=cordcut{i}(1,:);
        number=1;
        j=1;
        while j<=nrows-2
            a=[cordcut{i}(j+1,1)-cordcut{i}(j,1) cordcut{i}(j+1,2)-cordcut{i}(j,2)];
            b=[cordcut{i}(j+2,1)-cordcut{i}(j+1,1) cordcut{i}(j+2,2)-cordcut{i}(j+1,2)];
            cosTheta=dot(a,b)/(norm(a)*norm(b));
            ThetaInDegrees=acos(cosTheta)*180/pi;
            if ThetaInDegrees>=45 && ThetaInDegrees<=135  %% may change often
                list(number+1,:)=cordcut{i}(j+1,:);
                number=number+1;
                j=j+1;
            else
                j=j+1;
            end
        end
        list(number+1,:)=cordcut{i}(j+1,:);
        number=number+1;
        newList{i}=list(1:number,:);
    end
end

%% Function Further reducing number of points
% input: cordcut (x,y coordinates)
% input: tolerance (the minial distance between important points
% output: the reduced coordinates
% Method: control the distances between each line segments
function myList=reduce_points(cordcut, tolerance, area)
List=cell(length(cordcut),1);
for i=1:length(cordcut)
    X=cordcut{i}(:,1);
    Y=cordcut{i}(:,2);
    m=1;
    importance=ones(length(X),1);
    dist=sqrt((X-X(m)).^2+(Y-Y(m)).^2);
    while m<=length(X)
        importance(m)=2;
        importance(dist<tolerance & dist>0 | dist>2*sqrt(area/length(cordcut)))=0;
        m=find(importance==1 & dist>0,1);
        if isempty(m)
            break;
        end
        dist=sqrt((X-X(m)).^2+(Y-Y(m)).^2);
    end
    List{i}=[X(importance==2), Y(importance==2)];
end
myList=cell(length(List),1);
count=0;
for b=1:length(List)
    if length(List{b})>2
       count=count+1;
       myList{count}=List{b};
    end
end
myList=myList(1:count);
end

%% Function Contour output to XML file: NumArray2File(A,fname)
% Write the numeric values in array A to a text file using the numer format
% 10.4f
% fname is a string that names a plain text file. Apends A at the end of the
% original text.
% n is the number of lines written to the file.
function n=numArray2File(A, fname)
fid=fopen(fname,'a');
n=length(A);
if n==0 || ~isnumeric(n)
    return;
else
    % Note: Some Windows text editors, including Microsoft Notepad, require
    % a newline character sequence of '\r'n' instead of '\n'. '\n' is
    % sufficient for Microsoft Word or WordPad.
    str=sprintf('<ShapeCount>%i</ShapeCount>',n);
    fprintf(fid, '%s\r\n', str);
end
for r=1:n
    str=sprintf('<Shape_%i>',r);
    fprintf(fid, '%s\r\n', str);
    [nrows,~]=size(A{r});
    points=sprintf('<PointCount>%i</PointCount>',nrows);
    fprintf(fid,'%s\r\n',points);
    % Here assuming all the coordinates are integers. If they are not,
    % change the format such as '%10.1f' or so.
    for row=1:nrows
        str=sprintf('<X_%i>%d</X_%i>\r\n<Y_%i>%d</Y_%i>',row,A{r}(row,1),row,row,A{r}(row,2),row);
        fprintf(fid, '%s\r\n',str);   
    end
    str=sprintf('</Shape_%i>',r);
    fprintf(fid,'%s\r\n',str);
end
fprintf(fid,'</ImageData>');
fclose(fid);
end