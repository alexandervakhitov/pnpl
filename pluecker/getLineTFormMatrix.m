function [ LT ] = getLineTFormMatrix( shiftX, shiftY, shiftZ, rotX_rad, rotY_rad, rotZ_rad )
%GETLINETFORMMATRIX  Returns a 6x6 transformation matrix for 
% transforming Plucker lines based on input rotation angles and shift
% distances.

	R = getRotationMatrix(rotX_rad, rotY_rad, rotZ_rad);
	t = getTranslationVector(shiftX, shiftY, shiftZ);
	
	tX = skewSymMat(t);

	LT = [ ...
		   R      R*tX; ...
		zeros(3)    R   ...
	];

end

