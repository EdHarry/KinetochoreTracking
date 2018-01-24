function [MeanMat,StdMat]=MeanInMatrix(Matrix,Colsearch,Colmean)
%MEANINMATRIX Calculates the mean of specified columns.
%
% SYNOPSIS      [MeanMat,StdMat]=MeanInMatrix(Matrix,Colsearch,Colmean);  
%	
% INPUT         Matrix:     matrix to be treated
%               Colsearch:  Column containing non-changing values
%               Colmean:    Column to ba averaged
%
% OUTPUT        MeanMat:    Reduced matrix with means
%               StdMat:     Reduced matrix with standard deviations
%
% Example:      A=[1 1 1;1 1 2;1 2 1;1 2 4]
%               B=meaninmatrix(A,1,[2,3]) and C=meaninmatrix(A,2,3)
%               Results: B=[1 1.5 2] and C=[1 1 1.5;1 2 2.5]
%
% CB, 01-06-02 (M.24)

nanout=sum(isnan(Matrix(:,Colmean)),2);
Matrix=Matrix(find(nanout(:)==0),:);
Colfind=[0;find([Matrix(2:end,Colsearch)-Matrix(1:end-1,Colsearch);1])];


for i=1:(size(Colfind,1)-1)
    MeanMat(i,:)=Matrix(Colfind(i)+1,:);
    MeanMat(i,Colmean)=mean(Matrix(Colfind(i)+1:Colfind(i+1),Colmean),1);
    StdMat(i,1:size(Colmean,2))=std(Matrix(Colfind(i)+1:Colfind(i+1),Colmean),1);
end