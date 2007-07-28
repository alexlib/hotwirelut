function [mat2 mat2dist] = myNearestNeighbors(mat1,k)
% 
% Program to find the k - nearest neighbors (kNN) within a set of points. Given a data matrix with 
% each row corresponding to a vector, the program outputs the indexes of the k-NNs for every vector 
% in the matrix. It also outputs the corresponding Euclidean distances.
% 
% Usage:
% [neighbors distances] = kNearestNeighbors(dataMatrix,k);
% dataMatrix = (N x D) - N vectors each having a dimensionality D
% k - Number of nearest neighbors desired
% 
% Example:
% a = [ 1 2 3 4 5; 1 2 2 4 6]';
% [neighbors distances] = myNearestNeighbors(a,2);
% 
% Output:
% neighbours =      
%      2     3
%      3     1
%      2     1
%      3     5
%      4     3
%      
% distances = 
%     1.4142    2.2361
%     1.0000    1.4142
%     1.0000    2.2361
%     2.2361    2.2361
%     2.2361    4.4721


% mat2 = zeros(size(mat1,1),k);
% mat2dist = mat2;

numVectors = size(mat1,1);
for i=numVectors %1:numVectors,
    dist = sum((repmat(mat1(i,:),numVectors,1)-mat1).^2,2);
    [sortval sortpos] = sort(dist,'ascend');
    mat2(1,:) = sortpos(2:2+k-1);
    mat2dist(1,:) = sqrt(sortval(2:2+k-1));
end
    
    