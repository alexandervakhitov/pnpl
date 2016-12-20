function [ t ] = getTranslationVector( x, y, z )
%GETTRANSLATIONVECTOR Returns a 3x1 translation vector based on input parameters.

	% camera position (position of world origin with respect to camera = translaction vector)

	t = [-x; -y; -z];

end

