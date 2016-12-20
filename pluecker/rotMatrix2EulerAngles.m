function [ angleX, angleY, angleZ ] = rotMatrix2EulerAngles ( varargin )

	% Code based on Gregory G. Slabaugh: Computing Euler angles from a
	% rotation matrix, https://truesculpt.googlecode.com/hg-history/38000e9dfece971460473d5788c235fbbe82f31b/Doc/rotation_matrix_to_euler.pdf
	% Permutation of rotation matrix elements based on http://en.wikipedia.org/wiki/Euler_angles
	
	% process parameters
	switch nargin
		case 1
			rotOrder = 'zyx';
		case 2
			rotOrder = varargin{2};
		otherwise
			error('rotMatrix2EulerAngles:nargin', 'rotMatrix2EulerAngles expects 1 - 2 input parameters, %d given.', nargin);
	end
	
	R = varargin{1};
	
	if (~checkRotOrder(rotOrder))
		error('rotMatrix2EulerAngles:rotOrder', 'The 2nd argument has to contain exactly 3 characters, ''x'', ''y'', and ''z''. The ordering is arbitrary.');
	end
	
	if(~strcmp(rotOrder, 'zyx'))
		error('rotMatrix2EulerAngles:rotOrder', 'Any other rotOrder than ''zyx'' is not allowed.');
	end
	
	
	% for getRotationMatrix(rotX, rotY, rotZ, 'zyx')
	if (R(1,3) ~= 1 && R(1,3) ~= -1)
		angleY1 = -asin(R(1,3));
		angleY2 = pi - angleY1;
		angleZ1 = atan2(R(1,2) / cos(angleY1), R(1,1) / cos(angleY1));
		angleZ2 = atan2(R(1,2) / cos(angleY2), R(1,1) / cos(angleY2));
		angleX1 = atan2(R(2,3) / cos(angleY1), R(3,3) / cos(angleY1));
		angleX2 = atan2(R(2,3) / cos(angleY2), R(3,3) / cos(angleY2));
		
		% convert angles into the interval [-pi; pi)
		angleX = mod((angleX1 + pi), 2*pi) - pi;
		angleY = mod((angleY1 + pi), 2*pi) - pi;
		angleZ = mod((angleZ1 + pi), 2*pi) - pi;
				
	else
		angleZ = 0.0; % anything
		if (R(1,3) == -1)
			angleY =  pi/2;
			angleX =  angleZ + atan2(R(3,2), R(2,2));
		else
			angleY = -pi/2;
			angleX = -angleZ + atan2(-R(3,2), -R(3,1));
		end
		
		% convert angles into the interval [-pi; pi)
		angleX = mod((angleX + pi), 2*pi) - pi;
		angleY = mod((angleY + pi), 2*pi) - pi;
		angleZ = mod((angleZ + pi), 2*pi) - pi;
	end
	
end


function [ valid ] = checkRotOrder( rotOrder )
	if (length(rotOrder) ~= 3)
		valid = false;
		return;
	end
	
	xPos = strfind(rotOrder, 'x');
	yPos = strfind(rotOrder, 'y');
	zPos = strfind(rotOrder, 'z');
	
	if (length(xPos) ~= 1 || length(yPos) ~= 1 || length(zPos) ~= 1)
		valid = false;
		return;
	end
	
	valid = true;
	return;
end

