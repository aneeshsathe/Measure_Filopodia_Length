clear
clc

file_path='y:\Naila_filopodia\OLD_STUFF\Filopodia round2\140808 464ATR\';
file_name='140808_XY point 4 464ATR 34min.tif';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               %
%   DO NOT MODIFY BELOW THIS    %
%                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add required paths
addpath(genpath([fileparts(mfilename('fullpath')),'\modules\']))
addpath(genpath([fileparts(mfilename('fullpath')),'\rand_func']))

addpath(genpath(mfilename('fullpath')))

%%
img_write_path=[file_path,'results\',file_name(1:end-4),'_'];
mkdir([file_path,'results\'])
Rimg1=imread([file_path,file_name],'Index',1);%filo
Rimg2=imread([file_path,file_name],'Index',2);%actin

img1=mat2gray(Rimg1);%grey filo
img2=mat2gray(Rimg2);%grey actin

%% reop
imshow(img1+img2,'InitialMagnification',200)
h=imfreehand;
BW_crop=createMask(h);
close all

%%


Mimg1=mat2gray(adapthisteq(Rimg1)+adapthisteq(Rimg2));

Mimg1=noisecomp(Mimg1, 2,6, 3, 6, 0.5);


%%
rat_img=(Rimg1.*2)./Rimg2;

bwrat_img=im2bw(rat_img,graythresh(rat_img));
mImg2=bwrat_img.*img2;
mImg1=bwrat_img.*Mimg1;

MCimg=mat2gray(mImg1+(mImg2.*2)); %add channels
[cell_body, unerod_cell_body]=keep_max_area_obj(an_gray(MCimg));
cell_perim=bwperim(cell_body);
se = strel('diamond', 4);
aft_open=imdilate(cell_perim,se);




%%

INITPSF= fspecial('gaussian',5,5);

[J, P1] = deconvblind(mat2gray(Mimg1),INITPSF,5);
[J2, P2] = deconvblind(J,P1,5);

J2=noisecomp(J2, 2, 6, 2, 6, 0.5);
J2=cell_perim+J2.*~cell_body;
J2=bwareaopen(an_gray(wiener2(J2)),10);

J3=bwareaopen(an_gray(J2+cell_body+(cell_perim+Mimg1.*~cell_body)*2),5);

%%
BW3 = (bwmorph((bwmorph(J3+cell_body,'thin',2)+cell_perim).*double(BW_crop),'skel',Inf).*~cell_body)+cell_perim;
BW4=an_gray(BW3);

%%
J3=bwmorph(imerode(entropyfilt(Mimg1,ones(11)),strel('diamond',4)).*BW_crop,'skel',Inf).*~cell_body+cell_perim;
%
imshow(imoverlay(img1+img2,J3,[1,0,0]))

%%
[edgelist1, ~] = edgelink(J3, 20);
edgelist=cleanedgelist(edgelist1, 10);
edgeim=edgelist2image(edgelist,[size(BW4,1),size(BW4,2)]);

BW4=logical(BW4+logical(edgeim));

%%
bw2 = filledgegaps(BW4, 5);
[edgelist1, labelededgeim] = edgelink(bw2, 20);
edgelist=cleanedgelist(edgelist1, 10);
edgeim=edgelist2image(edgelist,[size(BW4,1),size(BW4,2)]);

BW4=logical(edgeim);
%%
[rj, cj, re, ce] = findendsjunctions(BW4, 0);

pi=find(imdilate(cell_perim,strel('diamond',2))==1);
ji=sub2ind(size(BW4),rj,cj);
LIAi=intersect(pi,ji);
[LIAr,LIAc]=ind2sub(size(BW4),LIAi);

%%
% imshow(img1+img2+BW4), hold on
% 	plot(cj,rj,'y+')
%     plot(LIAc,LIAr,'r+')
% 	plot(ce,re,'g+')
%     hold off
%%%%%%%
%%% FIND COMMON POINTS BETWEEN cell_perim and BW4 to replace perim_R and
%%% perim_C
%%


%%
BWbranchpt=bwmorph(BW4,'branchpoints');
BWbranchendpt=bwmorph(BW4,'endpoints');
%%
BW_sans_branchpt=BW4-BWbranchpt;
BW_sans_branchpt_area_filt = bwmorph(bwmorph(bwareaopen(BW_sans_branchpt, 5,8),'clean'),'spur');

%%
% disp_img=imoverlay(imoverlay(imoverlay(imoverlay(Mimg1,BW3,[1 0 0]),cell_perim,[0 1 0]),BWbranchpt,[0 0 1]),BWbranchendpt,[1 1 0]);
% imshow(disp_img,'InitialMagnification',200)

%% pixels commom to branch end points and cell_perim

BW_perim_and_branch=BWbranchpt & (cell_perim-BWbranchendpt);

[perim_R,perim_C]=find(BW_perim_and_branch==1);
[branch_R,branch_C]=find(BWbranchendpt==1);

%%

% figure,
% imshow(imoverlay(Mimg1, BW4,[1 0 0]))
% hold on
% for pl_count=1:size(perim_R,1)
%     scatter(perim_C(pl_count),perim_R(pl_count),'*y')
%     text(perim_C(pl_count),perim_R(pl_count),num2str(pl_count))
% end
% for pl_count=1:size(branch_C,1)
%     scatter(branch_C(pl_count),branch_R(pl_count),'*c')
%     text(branch_C(pl_count),branch_R(pl_count),num2str(pl_count))
% end
% hold off
%%

disp(size(re,1))
%%
[xls_label{1,1:3}]=deal('No.','ChessDist','EucDist');
tic
for count=1:size(re,1)
    
    BR_C=ce(count);
    BR_R=re(count);
    
    
    parfor count_perim =1:size(LIAc,1)
        path_length(count_perim,1)=an_dist_find_filo( BW4,LIAc(count_perim), LIAr(count_perim),BR_C, BR_R);        
    end
    
    %%
    [path_length, count_perim]=min(path_length);
    %     disp(['path length',num2str(path_length)])
    if path_length<Inf&&path_length>3
        
        D1 = bwdistgeodesic(BW4, LIAc(count_perim), LIAr(count_perim),'chessboard');
        D2=bwdistgeodesic(BW4, ce(count), re(count),'chessboard');
        
        D = D1 + D2;
        D = round(D * 8) / 8;
        
        D(isnan(D)) = inf;
        skeleton_path = imregionalmin(D);
        path_length_out = D(skeleton_path);
        path_length_out = path_length_out(1);
        chess_dist=path_length_out;
        
        P1 = imoverlay(imoverlay([img1+img2,img1+img2], skeleton_path, [1 0 0]),cell_perim,[0,1,0]);
        
        out_P1=insertText(P1,[1,1],['No. ',num2str(count), ' | Dist:', num2str(path_length), ' px. | Method: Cityblock']);
        
        pos   = [ce(count),  re(count);LIAc(count_perim), LIAr(count_perim);...
            ce(count)+512, re(count);LIAc(count_perim)+512, LIAr(count_perim)];
        
        color = {'red', 'red','red', 'red'};
        out_P1 = insertMarker(out_P1, pos, 's', 'color', color, 'size', 3);
        
        %%%%%%%%%%%%%%%%%%%%%
        D1 = bwdistgeodesic(BW4, LIAc(count_perim), LIAr(count_perim),'quasi-euclidean');
        D2 = bwdistgeodesic(BW4, ce(count), re(count), 'quasi-euclidean');
        
        D = D1 + D2;
        D = round(D * 8) / 8;
        
        D(isnan(D)) = inf;
        skeleton_path = imregionalmin(D);
        
        path_length_out = D(skeleton_path);
        path_length_out = path_length_out(1);
        eu_dist=path_length_out;
        %%
        P2 = imoverlay(imoverlay([img1+img2,img1+img2], skeleton_path, [1 0 0]),cell_perim,[0,1,0]);
        out_P2=insertText(P2,[1,1],['No. ',num2str(count), ' | Dist:', num2str(path_length), ' px. | Method: Cityblock']);
        
        pos   = [ce(count),  re(count);LIAc(count_perim), LIAr(count_perim);...
            ce(count)+512, re(count);LIAc(count_perim)+512, LIAr(count_perim)];
        
        color = {'red', 'red','red', 'red'};
        out_P2 = insertMarker(out_P2, pos, 's', 'color', color, 'size', 3);
        
        %% Collect Data
        [pro_xls_label{1,1:3}]=deal(['No ',num2str(count)],chess_dist,eu_dist);
        xls_label=vertcat(xls_label,pro_xls_label);
        
        %% Write Image
        out_P=[out_P1;out_P2];
        if count==1
            imwrite(out_P,[img_write_path,'out_img.tif'])
        else
            imwrite(out_P,[img_write_path,'out_img.tif'],'WriteMode','append')
        end
        disp(['No.',num2str(count),' of ', num2str(size(branch_C,1))])
    end
end
toc
%% write data
g=regionprops(cell_perim,'Area');
[xls_label{1,4:5}]=deal('Perim ',g.Area);
%
xls_file_name=[img_write_path,'out_results.txt'];

dlmcell(xls_file_name,xls_label,',')

%write excel file if OS is windows
if ispc
    xls_file_name2=[img_write_path,  'out_results.xls'];
    xlswrite(xls_file_name2,xls_label,'Sheet1' );
    xlswrite(xls_file_name2,xls_label,'Sheet2' );
end


disp('done writing excel file')
% figure,
% imshow(out_P1)
%http://www.csse.uwa.edu.au/~pk/Research/MatlabFns/#edgelink