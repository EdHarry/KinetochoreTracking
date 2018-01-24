function spy3d(matIn)

%The spy function, but in 3d!!

if ndims(matIn) ~= 3
    error('Input matrix must be 3d!!')
end
    

[M,N,P] = size(matIn);

figure
hold on

for p = 1:P
    
    [row,col] = find(matIn(:,:,p));
    
    plot3(row,col,p * ones(1,length(row)),'.','color',[ (P-p)/P 0 p/P ]  );
    
    
end

axis vis3d