function [ reproj_err, varargout ] = reprojErr_lines( lines_2D, line_2D_end_pts )
%REPROJERR_LINES Computes reprojection error of lines / line segments.
%   This function computes a reprojection error of line segments to lines,
%   based on: C.J.Taylor and D.J.Kriegman: "Structure and Motion from Line
%   Segments in Multiple Images".
%
%		The return value is a scalar (number) - a sum of individual
%		reprojection errors between the pairs lines2D(i) - line2DEndPts(j,j+1).

	if(size(lines_2D, 1) ~= 3)
		error('reprojErr_lines:argin', 'The 1st input argument has to be of size 3xN.');
	end

	if((size(line_2D_end_pts, 1) ~= 3) || rem(size(line_2D_end_pts, 2), 2))
		error('reprojErr_lines:argin', 'The 2nd input argument has to be of size 3xN, where N is even.');
	end
	
	if((2 * size(lines_2D, 2)) ~= size(line_2D_end_pts, 2))
		error('reprojErr_lines:argin', 'The 2nd argument has to contain twice as much elements as the 1st argument (endpoints of line segments).');
	end
	
	% normalize line2DEndPts to [x; y; 1]
	for i = 1:3
		line_2D_end_pts(i,:) = line_2D_end_pts(i,:) ./ line_2D_end_pts(3,:);
	end

	NLINES = size(lines_2D, 2);
	
	% compute the reprojection error
	reproj_errs = zeros(NLINES, 1);
	
	for line_nr = 1:NLINES
		
		idx = 2*(line_nr-1)+1 : 2*(line_nr-1)+2;
		
		m = lines_2D(:,line_nr);
		A = line_2D_end_pts(:,idx)';
		
		line_len = sqrt((A(1,1) - A(2,1))^2 + (A(1,2) - A(2,2))^2);
		
		B = (line_len / (3 * (m(1)^2 + m(2)^2))) * [1.0  0.5; 0.5  1.0];
		
		reproj_errs(line_nr) = m' * (A' * B * A) * m;
	end
	
	reproj_err = sum(reproj_errs);
	
	if(nargout > 1)
		varargout{1} = reproj_errs;
	end

end
