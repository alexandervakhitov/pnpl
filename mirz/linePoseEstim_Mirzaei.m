function [ R, t, varargout ] = linePoseEstim_Mirzaei( line3DEndPts, line2DEndPts, zPM )
%LINEPOSEESTIM_MIRZAEI Camera pose estimation from line correspondences based on
%the paper Mirzaei, Roumeliotis: Globally Optimal Pose Estimation from Line
%Correspondences, ICRA 2011 and on the VPtoolbox from Mirzaei.
%   line3DEndPts - 4x(2N) matrix of 3D line endpoints [x; y; z; w]
%   line2DEndPts - 3x(2N) matrix of 2D line endpoints [x; y; w]
%   focalL - focal length in pixels
%   principalX - X coordinates of the principal point of the image
%   principalY - Y coordinates of the principal point of the image

	% input checks
	if (rem(size(line3DEndPts,2), 2))
		error('Number of 3D line endpoints has to be an even number.');
	end
	
	if (size(line3DEndPts,1) ~= 4)
		error('3D line endpoints have to be homogeneous coordinates - 4-tuples [x; y; z; w].');
	end
	
	if (size(line3DEndPts,2) ~= size(line2DEndPts,2))
		error('Number of 3D and 2D line endpoints has to be equal.');
	end;
	
	if (size(line2DEndPts,1) ~= 3)
		error('2D line endpoints have to be homogeneous coordinates - 3-tuples [x; y; w].');
	end

	NLINES = size(line3DEndPts,2)/2;
	
	if (NLINES < 3)
		error('At least 3 lines has to be supplied.');
	end
	
	%% Normalize 3D line endpoints
	for i = 1:4
		line3DEndPts(i,:) = line3DEndPts(i,:) ./ line3DEndPts(4,:);
	end

	%% Create 3D lines (as directions)
	lines3Ddir = line3DEndPts(1:3, 2:2:end) - line3DEndPts(1:3, 1:2:end);
	lines3Ddir = normc(lines3Ddir);
	
	%% Change the camera coordinate system to fit Matlab img. coord. system
	% Change the direction of Y-axis to conform to the Matlab image coord. system,
	% and after that also change the direction of the Z-axis (going now out of the
	% camera into the scene) in order to keep the system right-handed. This is
	% done by multiplying the 2D coords by diag([-1; -1; 1]) ... I don't fully
	% understand why.
	line2DEndPts = diag([-1; -1; 1]) * line2DEndPts;
	
	%% Create 2D lines
	
	clear 'lines2D';
	lines2D(NLINES).point1 = [];
	lines2D(NLINES).point2 = [];
	for lineNr=1:NLINES
		idx = 2*(lineNr-1)+1 : 2*(lineNr-1)+2;
		lines2D(lineNr).point1 = line2DEndPts(1:2,idx(1))';
		lines2D(lineNr).point2 = line2DEndPts(1:2,idx(2))';
	end

	% compute line moments (and remove radial distortion)
	fc = [1;1];
	cc = [0;0];
	alpha_c = 0;
	lines2D = normalize_lines(lines2D, fc, cc, alpha_c);

	%% Estimate VP (i.e. camera orientation) using lines
	
	if (NLINES == 3)
		opts = EstimateVanishingPointsDefaultOpts('minimal');
% 		opts.kSolutionInClassTol = 1; % default .99
% 		opts.kNormalizationTol = 1e-32; % default 1e-16;
		
		tStart_R = tic;
		[bestQuat, bestRes] = EstimateVanishingPoints([lines2D.nmoment], lines3Ddir, 'minimal', opts);
		time_R = toc(tStart_R);
	else
		opts = EstimateVanishingPointsDefaultOpts('relaxed');
		opts.kSolutionInClassTol = .99; % default .99
		opts.kNormalizationTol = 1e-16; % default 1e-16;
		opts.kMulMatrixCondToAbort = 1e14; % default 1e14; do not set lower otherwsise no solutions will be returned
		
		tStart_R = tic;
		[bestQuat, bestRes] = EstimateVanishingPoints([lines2D.nmoment], lines3Ddir, 'relaxed', opts);
		time_R = toc(tStart_R);
	end
	
	% convert quaternions to rotation matrices
	if (size(bestQuat, 2) > 0)
		R_Mirzaei = SpinCalc('QtoDCM', bestQuat', eps, false);
	else
		R_Mirzaei = NaN(3, 3);
	end
	
	err = bestRes;
	
	%% Estimate camera position
	% according to Eq. (20) of the paper "Mirzaei, Roumeliotis:
	% Globally Optimal Pose Estimation from Line Correspondences, ICRA 2011"
	%
	% Notation regarding Eq. (20):
	% C{w}i'.. a point on an image line ............................ line2Dpoints
	% G{l}i .. direction vector of the 3D line in the global frame . lines3Ddir
	% G{m}i .. moment of the 3D line in the global frame ........... lines3Dmoment
	% G{p}C .. camera position which is being estimated ............ t
	
% 	% direction vectors of the 3D lines
% 	lines3Ddir = line3DEndPts(:,2:2:end) - line3DEndPts(:,1:2:end);
% 	lines3Ddir = normc(lines3Ddir);
	
	% moments of the 3D lines
	lines3Dmoment = cross(line3DEndPts(1:3, 1:2:end), lines3Ddir);
	
	% points on the image lines
	line2DEndPtsA = [[lines2D.pnt1_n]; ones(1,NLINES)];
	line2DEndPtsB = [[lines2D.pnt2_n]; ones(1,NLINES)];
	thetas = -atan((line2DEndPtsB(1,:) - line2DEndPtsA(1,:)) ./ (line2DEndPtsB(2,:) - line2DEndPtsA(2,:)));
	rhos = -line2DEndPtsA(1,:) .* cos(thetas) - line2DEndPtsA(2,:) .* sin(thetas);
	line2Dpoints = [-rhos .* cos(thetas);   -rhos .* sin(thetas);   ones(1,NLINES)];

	% compute a position-solution for each rotation-matrix-solution
	NSOLUTIONS = size(R_Mirzaei,3);
	t_rough_Mirzaei = zeros(3, 1, NSOLUTIONS);
	valid = false(NSOLUTIONS, 1);
	time_t = 0;
	
	for i = 1 : NSOLUTIONS
		
		tStart_t = tic;
		
		% construct left and right side of equation (20)
		left = zeros(NLINES, 3);
		for j = 1:NLINES
			left(j,:) =  line2Dpoints(:,j)' * R_Mirzaei(:,:,i)' * skewSymMat(lines3Ddir(:,j));
		end

		right = zeros(NLINES, 1);
		for j = 1:NLINES
			right(j) =  -line2Dpoints(:,j)' * R_Mirzaei(:,:,i)' * lines3Dmoment(:,j);
		end

		camPos = left \ right; %!V!|
		t_rough_Mirzaei(:,:,i) = getTranslationVector(camPos(1), camPos(2), camPos(3));
		
		time_t = time_t + toc(tStart_t);
	end

	time = time_R + time_t;
	
	% Convert from coord. system of Mirzaei to my coord. system (reversed Z axis in the camera frame)
	R = zeros(3, 3, NSOLUTIONS);
	for i = 1:NSOLUTIONS
		R(:,:,i) = getRotationMatrix(0, 0, pi) * R_Mirzaei(:,:,i)';
	end
	
	t = t_rough_Mirzaei;
	
	
	% If multiple solutions exist, choose the best (in front of the camera and least reprojection error)
	if (NSOLUTIONS > 1)
		reprojErrors = zeros(NSOLUTIONS, 1);
		for i = 1:NSOLUTIONS
			T = R(:,:,i) * [eye(3) t(:,:,i)];
			line3DImgEndPts = T * line3DEndPts;
			
			% choose a solution which is physically plausible, i.e. scene in front of camera
			center3D = mean(line3DEndPts, 2);
			
			% this test needs to be -dot(.) for VGG and dot(.) for other experiments
			test = sign(zPM) * dot(R(3,:,i), (center3D(1:3) - t(:,:,i))/norm(center3D(1:3) - t(:,:,i)));
			if (test > 0)
				valid(i) = true;
			end
			
			% reprojection errors
			lines        = cross(line3DImgEndPts(:, 1:2:end), line3DImgEndPts(:, 2:2:end));
			lineSegments = diag([-1; -1; 1]) * line2DEndPts;
			reprojErrors(i) = reprojErr_lines(lines, lineSegments);
		end
		
		% keep only solutions in front of the camera
		if (~any(valid))
			valid = true(NSOLUTIONS, 1);
		end

		R         = R(:,:,valid);
		t         = t(:,:,valid);
		
		err       = err(valid);
		reprojErrors    = reprojErrors(valid);		
		
		[~, bestSolInd] = min(reprojErrors);
		
		R = R(:,:,bestSolInd);
		t = t(:,:,bestSolInd);
		
		err = err(bestSolInd);
	end
	

	if(nargout > 2)
		varargout{1} = err;
	end
	
	if(nargout > 3)
		varargout{2} = time;
	end
	
end

