function SalDetectionDemo

tic
clear all;
clc;
addpath('./others/');
%%------------------------set parameters---------------------%%
global imgRoot;
global saldir;
global supdir;
global imnames;
spnumber = 200;
theta = 10;
alpha = 0.99;
imgRoot='./test/';
saldir='./saliencymap/';
supdir='./superpixels/';
mkdir(supdir);
mkdir(saldir);
imnames=dir([imgRoot '*' '.jpg']);

for ii=1:length(imnames)   
    disp(ii);
    imname = [imgRoot imnames(ii).name]; 
    [input_im, w]=removeframe(imname);
    
    [m,n,k] = size(input_im);    
    input_vals=reshape(input_im, m*n, k);
    clear k

    imname=[imname(1:end-4) '.bmp'];
    comm=['SLICSuperpixelSegmentation' ' ' imname ' ' int2str(20) ' ' int2str(spnumber) ' ' supdir];
    system(comm);    
    spname=[supdir imnames(ii).name(1:end-4)  '.dat'];
    superpixels=ReadDAT([m,n],spname);
    spnum=max(superpixels(:));
    adjc = calAdjacentMatrix(superpixels,spnum);
    clear comm spname
    
    % compute the feature (mean color in lab color space) for each superpixels    
    rgb_vals = zeros(spnum,3);
    inds=cell(spnum,1);
    [x, y] = meshgrid(1:1:n, 1:1:m);
    x_vals = zeros(spnum,1);
    y_vals = zeros(spnum,1);
    num_vals = zeros(spnum,1);
    for i=1:spnum
        inds{i}=find(superpixels==i);
        num_vals(i) = length(inds{i});
        rgb_vals(i,:) = mean(input_vals(inds{i},:),1);
        x_vals(i) = sum(x(inds{i}))/num_vals(i);
        y_vals(i) = sum(y(inds{i}))/num_vals(i);
    end  
    seg_vals = colorspace('Lab<-', rgb_vals);
    clear lab_vals input_vals I_gradient Cets superpixels
    
    W_lab = normalize( DistanceZL(seg_vals, seg_vals, 'euclid') ); % euclid
    sal_lab = exp( -theta*W_lab );
    
    % Ranking
    W = sal_lab.*adjc;
    dd = sum(W,2); D = sparse(1:spnum,1:spnum,dd);
    
    P = (D-alpha*W)\eye(spnum); 
    Sal = P*sal_lab;
    Sal = Sal';
    clear W dd D
    
    salSup = Sal.*(ones(spnum,1)*num_vals');
    Isum = sum(salSup,2);
    
    comSal = calfgDistributionSal(salSup,Isum,x_vals,y_vals,m,n,spnum);
    clear Sal salSup Isum
    
    disSup = sqrt((x_vals*ones(1,spnum) - ones(spnum,1)*x_vals').^2 + (y_vals*ones(1,spnum) - ones(spnum,1)*y_vals').^2).*(1-adjc);
    sal_lab = sal_lab.*(1-eye(spnum));
    
    cocomSal = calLocalContrastSal(W_lab,comSal,disSup,sal_lab,spnum,adjc,theta,num_vals,x_vals,y_vals);
    fgSal = normalize( P*cocomSal );
    clear disSup
    
    finSal = normalize( comSal + fgSal );
    clear weight W_lab sal_lab
    
    salall = zeros(m,n);
    for k = 1:spnum
        salall(inds{k}) = finSal(k);
    end
    clear x y x_vals y_vals adjc num_vals spnum inds salbg index bgNum bgSal_lab bgSal bfgSal rgb_vals seg_vals
    
    %% Êä³öÏÔÖøÐÔÍ¼
    mapstage1=zeros(w(1),w(2));
    mapstage1(w(3):w(4),w(5):w(6))=salall;
    mapstage1=uint8(mapstage1*255);
    outname=[saldir imnames(ii).name(1:end-4) '.jpg'];
    imwrite(mapstage1,outname);
    clear mapstage1 outname w salall m n
end

toc