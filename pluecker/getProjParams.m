function [ x, y, z, angleX, angleY, angleZ ] = getProjParams( P )
%GETPROJPARAMS Extracts position and orientation of the camera from a line projection matrix.

	% input checks

	if(size(P,1) ~= 3)
		error('getProjParams:input', 'The input matrix must have 3 rows.');
	end
	
	if(size(P,2) ~= 6)
		error('getProjParams:input', 'The input matrix must have 6 columns.');
	end


	%%
	
	% scale the projection matrix so that det(P1) == det(R) == 1
	P1 = P(:, 1:3);
	s = nthroot(1/det(P1), 3);
	Pscale = s * P;
	P2scale = Pscale(:, 4:6);
	
	% method inspired by decomposition of the Essential matrix = R.[t]x by
	% Hartley: MVGinCV, p.258
	[U, Sigma, V] = svd(P2scale);

	Z = [ ...
		0  1  0; ...
	 -1  0  0; ...
		0  0  0  ...
	];
	W = [ ...
		0 -1  0; ...
		1  0  0; ...
		0  0  1  ...
	];
	avgSGValue = (Sigma(1,1) + Sigma(2,2)) / 2;

	% 2 possible solutions (explained in Y. Ma, S. Soatto, J. Kosecka, S. S. Sastry: An Invitation to 3D Vision, ch. 5, p. 83)
	da = det(U * W  * V'); % either +1 or -1
	db = det(U * W' * V'); % either +1 or -1
	
	Ra = U * W  *  diag([1 1 da])  * V';
	Rb = U * W' *  diag([1 1 db])  * V';

	tXa = avgSGValue * V * Z  * V';
	tXb = avgSGValue * V * Z' * V';

	ta = skewSymMat2Vec(tXa);
	tb = skewSymMat2Vec(tXb);
	[angleXa, angleYa, angleZa] = rotMatrix2EulerAngles(Ra);
	[angleXb, angleYb, angleZb] = rotMatrix2EulerAngles(Rb);

	x = num2cell([-ta(1) -tb(1)]);
	y = num2cell([-ta(2) -tb(2)]);
	z = num2cell([-ta(3) -tb(3)]);
	
	angleX = num2cell([angleXa angleXb]);
	angleY = num2cell([angleYa angleYb]);
	angleZ = num2cell([angleZa angleZb]);

	return;
end

