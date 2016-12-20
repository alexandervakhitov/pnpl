function [ R ] = getRotationMatrix( varargin )
%GETROTATIONMATRIX Returns a 3x3 rotation matrix based on input rotation angles.

	% camera orientation (rotation/Euler angles -> rotation matrix)

	% process parameters
	switch nargin
		case 3
			rotOrder = 'zyx';
		case 4
			rotOrder = varargin{4};
		otherwise
			error('getRotationMatrix expects 3 - 4 input parameters, %d given.', nargin);
	end
	
	rotX_rad = varargin{1};
	rotY_rad = varargin{2};
	rotZ_rad = varargin{3};

	if (~checkRotOrder(rotOrder))
		error('getRotationmatrix: the 4th argument has to contain exactly 3 characters, ''x'', ''y'', and ''z''. The ordering is arbitrary.');
	end


	% build individual axis rotation matrices
	Rx = [ ...
		1        0             0;       ...
		0  cos(rotX_rad) sin(rotX_rad); ...
		0 -sin(rotX_rad) cos(rotX_rad)  ...
	];

	Ry = [                             ...
		cos(rotY_rad) 0 -sin(rotY_rad); ...
				0        1         0;       ...
		sin(rotY_rad) 0  cos(rotY_rad)  ...
	];

	Rz = [                             ...
		cos(rotZ_rad) sin(rotZ_rad) 0;  ...
	 -sin(rotZ_rad) cos(rotZ_rad) 0;  ...
				 0              0        1  ...
	];

	
	% Evaluate order of rotations and build the final R matrix
	R = eye(3);
	
	for i = 1:3
		rotChar = rotOrder(i);
		switch (rotChar)
			case 'x'
				Ri = Rx;
			case 'y'
				Ri = Ry;
			case 'z'
				Ri = Rz;
		end
		R = Ri * R;
	end
	
	return;

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
