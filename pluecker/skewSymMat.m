function [ ssMat ] = skewSymMat( vec )
%SKEWSYMMAT Returns a 3x3 skew-symmetric matrix for a 3-vector.

	if(length(vec) ~= 3)
		error('skewSymMat:input', 'The input has to be a 3-vector.');
	end

	ssMat = [ ... 
		   0    -vec(3)  vec(2); ...
		 vec(3)    0    -vec(1); ...
		-vec(2)  vec(1)    0   ...
	];

end

