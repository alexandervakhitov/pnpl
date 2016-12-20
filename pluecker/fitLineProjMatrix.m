function [ P, M ] = fitLineProjMatrix( lines_3D, lines_2D, w )
%FITLINEPROJMATRIX Finds 3x6 projection matrix for 3D lines in Plucker
%coordinates and their 2D images.
%
% algorithm based on: http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/EPSRC_SSAZ/node5.html
%

%% Input checks
	if (size(lines_3D,2) ~= size(lines_2D,2))
		error('fitLineProjMatrix:linecount', 'number of image lines and world lines is not equal');
	end
	
	if (size(lines_3D,2) ~= length(w))
		error('fitLineProjMatrix:weightcount', 'number of lines and weights is not equal');
	end
	
%% Least-squares fitting

	% assemble the system of linear equations Mp = 0;  p = P(:)
	M = getMeasurementMatrix(lines_3D, lines_2D);
	
	% diagonal matrix of weights
	W = diag([w; w]);
	
	% the solution is the right singular vector with least singular value
	[~, ~, V] = svd(W * M, 'econ');
	proj_mat_lines = V(:,end);

	% form a 3x6 matrix from the eigenvector
	P = reshape(proj_mat_lines, 3, 6);
	
	return;
end

