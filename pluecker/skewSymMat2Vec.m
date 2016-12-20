function [ vec ] = skewSymMat2Vec( ssMatNoisy )
%SKEWSYMMAT2VEC Returns a 3-vector from a 3x3 skew-symmetric matrix.

	if(size(ssMatNoisy,1) ~= 3 || size(ssMatNoisy,2) ~= 3)
		error('skewSymMat2Vec:input', 'The input has to be a 3x3 matrix.');
	end
	
	% in the case the skew-symmetric matrix is noisy, compute the nearest
	% truly skew-symmetric matrix in the sense of the Frobenius norm
	% reference: G. H. Golub: Matrix Computations
	ssMat = (ssMatNoisy - ssMatNoisy') / 2;
	
	vec = [ ...
		ssMat(3,2); ...
		ssMat(1,3); ...
		ssMat(2,1)  ...
	];
end

