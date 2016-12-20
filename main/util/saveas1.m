function saveas( h, name, format )
%SAVEAS Save Figure or Simulink block diagram in desired output format
%   SAVEAS(H,'FILENAME')
%   Will save the Figure or Simulink block diagram with handle H to file 
%   called FILENAME. 
%   The format of the file is determined from the extension of FILENAME.
%
%   SAVEAS(H,'FILENAME','FORMAT')
%   Will save the Figure or Simulink block diagram  with handle H to file 
%   called FILENAME in the format specified by FORMAT. FORMAT can be the 
%   same values as extensions of FILENAME. 
%   The FILENAME extension does not have to be the same as FORMAT.  
%   The specified FORMAT overrides FILENAME extension.
%
%   Valid options for FORMAT are:
%
%   'fig'  - save figure to a single binary FIG-file.  Reload using OPEN. 
%   'm'    - save figure to binary FIG-file, and produce callable
%            MATLAB code file for reload.
%   'mfig' - same as M.
%   'mmat' - save figure to callable MATLAB code file as series of creation
%            commands with param-value pair arguments.  Large data is saved
%            to MAT-file.  
%            Note: MMAT Does not support some newer graphics features. Use
%                  this format only when code inspection is the primary goal.
%                  FIG-files support all features, and load more quickly. 
%
%   Additional FORMAT options include devices allowed by PRINT.
%
%   NOTE: not all format options are allowed for Simulink block diagrams.
%   See the online help for more information.
%
%   Examples:
%
%   Write current figure to MATLAB fig file
%
%       saveas(gcf, 'output', 'fig')
%
%   Write current figure to windows bitmap file
%
%       saveas(gcf, 'output', 'bmp')
%
%   Write block diagram named 'demo' to an Encapsulated Postscript file
%
%       saveas(get_param('demo', 'Handle'), 'output', 'eps')
%
%   In the following list, SAVE_SYSTEM is available for Simulink users only. 
%   See also LOAD, SAVE, OPEN, PRINT, SAVE_SYSTEM

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.16.4.12 $ 

%Input validation
if nargin < 2 
    error(message('MATLAB:saveas:tooFewInputs'))
end

if ~ishandle(h)
    error('MATLAB:saveas:invalidHandle','%s',...
        getString(message('MATLAB:saveas:invalidHandle')))
end

if ishghandle(h)
    while ~isempty(h) &&  ~isfigure(h)  
        h = get(h,'Parent');
    end
    if isempty(h)
        error('MATLAB:saveas:invalidHandle','%s',...
           getString(message('MATLAB:saveas:invalidFigureHandle')))
    end
else
    if ~strcmp( 'block_diagram', get_param( h, 'type' ) ) ...
            && ~strcmp( 'SubSystem', get_param( h, 'blocktype' ) )
        error('MATLAB:saveas:invalidHandle','%s',...
            getString(message('MATLAB:saveas:invalidSimulinkHandle')))
    end
end


if ~ischar(name) || isempty(name)
    error('MATLAB:saveas:invalidFilename','%s',...
    getString(message('MATLAB:saveas:invalidFilename')));
end

% since this is a callback from the figure menu, it is possible that
% uimenufcn will give a bogus filename (one with a * in it) if we get
% any *'s in the filename, error out gracefully.
% might want to generalize this to the same set of chars windows prevents
% in filenames - could be crippling on unix though.
if any(name == '*')
    error('MATLAB:saveas:invalidFilename','%s',...
         getString(message('MATLAB:saveas:invalidFilenameAsterisk', name)))
end

% Make sure we can write given file. Note fileparts returns leaf directory in
% name return argument if given filename is only a path, i.e. 'c:\temp'.
[fp,fn,fe]=fileparts(name);

% NOTE: this does not allow to say:
% saveas(gcf, 'foo.', '-fig')
% maybe that's OK...
% william - 11/98
if ~isempty(fe) && strcmp(fe, '.')
    error('MATLAB:saveas:invalidFilename','%s',...
        getString(message('MATLAB:saveas:invalidFilenameExtension',name)));
end
if isempty(fe)
   fe = '.fig';
end
if isempty(fn) 
    error('MATLAB:saveas:invalidFilename','%s',...
        getString(message('MATLAB:saveas:invalidFilenameFile',name)));
end
if ~isempty(fp) && ~exist(fp, 'dir')
    error('MATLAB:saveas:invalidFilename','%s',...
        getString(message('MATLAB:saveas:invalidFilenameBadPath' ,name)))
end

if nargin == 2
    format = fe(2:end); %Do not want the '.'
end

%For some formats we have helper SAVEAS... functions
% make sure format is defined to prevent recursion into this file
if ~isempty(format) && any( exist( ['saveas' format]) == [2 3 5 6] ) %#ok
    feval( ['saveas' format], h, name )
    
else    
    %If FORMAT is specified, look first to see if it is an extension
    %we know is that of a PRINT output format. If not, look to see if
    %it is a PRINT supported device format specifier.
    [ops,dev,ext] = printtables(printjob(h)); %#ok
    i = strmatch( format, ext ); %#ok
    
    if length(i) >= 1
        %Handle special cases, more then one device, i.e. format='ps'
        i = i(1);
        
    elseif isempty(i)
        i = strmatch( format, dev, 'exact'  );
        if isempty(i)
            i = strmatch( format, dev  ); %#ok
            if ~isempty(i)
                %Need to handle cases of multiple devices, i.e. format='ljet'
                i = i(1);
            end
        end
    end
    
    %FORMAT is a PRINT support ext or driver
    if isempty(i)
        error(message('MATLAB:saveas:badFormat', format))
    else
        print( h, name, ['-d' dev{i} '-r0'] )
        return
    end
    
end

% LocalWords:  mfig mmat online blocktype uimenufcn william ps ljet
