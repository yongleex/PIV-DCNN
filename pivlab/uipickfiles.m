function out = uipickfiles(varargin)

% Version: 1.3, 5 February 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})



%uipickfiles: GUI program to select file(s) and/or directories.
%
% Syntax:
%   files = uipickfiles('PropertyName',PropertyValue,...)
%
% The current directory can be changed by operating in the file navigator:
% double-clicking on a directory in the list to move further down the tree,
% using the popup menu to move up the tree, typing a path in the box to
% move to any directory or right-clicking on the path box to revisit a
% previously-listed directory.
%
% Files can be added to the list by double-clicking or selecting files
% (non-contiguous selections are possible with the control key) and
% pressing the Add button.  Files in the list can be removed or re-ordered.
% When finished, a press of the Done button will return the full paths to
% the selected files in a cell array, structure array or character array.
% If the Cancel button is pressed then zero is returned.
%
% The following optional property/value pairs can be specified as arguments
% to control the indicated behavior:
%
%   Property    Value
%   ----------  ----------------------------------------------------------
%   FilterSpec  String to specify starting directory and/or file filter.
%               Ex:  'C:\bin' will start up in that directory.  '*.txt'
%               will list only files ending in '.txt'.  'c:\bin\*.txt' will
%               do both.  Default is to start up in the current directory
%               and list all files.  Can be changed with the GUI.
%
%   Type        Two-column cell array where the first column contains file
%               suffixes and the second column contains descriptions.
%
%   REFilter    String containing a regular expression used to filter the
%               file list.  Ex: '\.m$|\.mat$' will list files ending in
%               '.m' and '.mat'.  Default is empty string.  Can be used
%               with FilterSpec and both filters are applied.  Can be
%               changed with the GUI.
%
%   Prompt      String containing a prompt appearing in the title bar of
%               the figure.  Default is 'Select files'.
%
%   NumFiles    Scalar or vector specifying number of files that must be
%               selected. A scalar specifies an exact value; a two-element
%               vector can be used to specify a range, [min max].  The
%               function will not return unless the specified number of
%               files have been chosen.  Default is [] which accepts any
%               number of files.
%
%   Output      String specifying the data type of the output: 'cell',
%               'struct' or 'char'.  Specifying 'cell' produces a cell
%               array of strings, the strings containing the full paths of
%               the chosen files.  'Struct' returns a structure array like
%               the result of the dir function except that the 'name' field
%               contains a full path instead of just the file name.  'Char'
%               returns a character array of the full paths.  This is most
%               useful when you have just one file and want it in a string
%               instead of a cell array containing just one string.  The
%               default is 'cell'.
%
% All properties and values are case-insensitive and need only be
% unambiguous.  For example,
%
%   files = uipickfiles('num',1,'out','ch')
%
% is valid usage.

% Version: 1.3, 5 February 2010
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Define properties and set default values.
prop.filterspec = '*';
prop.type = {};
prop.refilter = '';
prop.prompt = 'Select files';
prop.numfiles = [];
prop.output = 'cell';

% Process inputs and set prop fields.
properties = fieldnames(prop);
arg_index = 1;
while arg_index <= nargin
	arg = varargin{arg_index};
	if ischar(arg)
		prop_index = find(strncmpi(arg,properties,length(arg)));
		if length(prop_index) == 1
			prop.(properties{prop_index}) = varargin{arg_index + 1};
		else
			error('Property ''%s'' does not exist or is ambiguous.',arg)
		end
		arg_index = arg_index + 2;
	elseif isstruct(arg)
		arg_fn = fieldnames(arg);
		for i = 1:length(arg_fn)
			prop_index = find(strncmpi(arg_fn{i},properties,...
				length(arg_fn{i})));
			if length(prop_index) == 1
				prop.(properties{prop_index}) = arg.(arg_fn{i});
			else
				error('Property ''%s'' does not exist or is ambiguous.',...
					arg_fn{i})
			end
		end
		arg_index = arg_index + 1;
	else
		error(['Properties must be specified by property/value pairs',...
			' or structures.'])
	end
end

% Validate FilterSpec property.
if isempty(prop.filterspec)
	prop.filterspec = '*';
end
if ~ischar(prop.filterspec)
	error('FilterSpec property must contain a string.')
end

% Validate Type property.
if isempty(prop.type)
elseif iscellstr(prop.type) && isvector(prop.type)
	prop.type = repmat(prop.type(:),1,2);
elseif iscellstr(prop.type) && size(prop.type,2) == 2
else
	error(['Type property must be empty or a cellstr vector or ',...
		'a 2-column cellstr matrix.'])
end

% Validate REFilter property.
if ~ischar(prop.refilter)
	error('REFilter property must contain a string.')
end

% Validate Prompt property.
if ~ischar(prop.prompt)
	error('Prompt property must contain a string.')
end

% Validate NumFiles property.
if numel(prop.numfiles) > 2 || any(prop.numfiles < 0)
	error('NumFiles must be empty, a scalar or two-element vector.')
end
prop.numfiles = unique(prop.numfiles);
if isequal(prop.numfiles,1)
	numstr = 'Select exactly 1 file.';
elseif length(prop.numfiles) == 1
	numstr = sprintf('Select exactly %d files.',prop.numfiles);
else
	numstr = sprintf('Select %d to %d files.',prop.numfiles);
end

% Validate Output property.
legal_outputs = {'cell','struct','char'};
out_idx = find(strncmpi(prop.output,legal_outputs,length(prop.output)));
if length(out_idx) == 1
	prop.output = legal_outputs{out_idx};
else
	error(['Value of ''Output'' property, ''%s'', is illegal or '...
		'ambiguous.'],prop.output)
end


% Initialize file lists.
[current_dir,f,e] = fileparts(prop.filterspec);
filter = [f,e];
if isempty(current_dir)
	current_dir = pwd;
end
if isempty(filter)
	filter = '*';
end
re_filter = prop.refilter;
full_filter = fullfile(current_dir,filter);
path_cell = path2cell(current_dir);
fdir = filtered_dir(full_filter,re_filter);
filenames = {fdir.name}';
filenames = annotate_file_names(filenames,fdir);

% Initialize some data.
file_picks = {};
full_file_picks = {};
%dir_picks = dir(' ');  % Create empty directory structure. --> adapt for
%r2015a
dir_picks = repmat(dir(char(127)),0,1); % Create empty directory structure.
show_full_path = false;
nodupes = true;
history = {current_dir};

% Create figure.
gray = get(0,'DefaultUIControlBackgroundColor');
fig = figure('Position',[0 0 740 580],...
	'Color',gray,...
	'WindowStyle','modal',...
	'Resize','off',...
	'NumberTitle','off',...
	'Name',prop.prompt,...
	'IntegerHandle','off',...
	'CloseRequestFcn',@cancel,...
	'CreateFcn',{@movegui,'center'});

% Create uicontrols.
uicontrol('Style','frame',...
	'Position',[255 290 110 40])
uicontrol('Style','frame',...
	'Position',[275 135 110 100])

navlist = uicontrol('Style','listbox',...
	'Position',[10 10 250 320],...
	'String',filenames,...
	'Value',[],...
	'BackgroundColor','w',...
	'Callback',@clicknav,...
	'Max',2,'tooltipstring','Current folder contents');
pickslist = uicontrol('Style','listbox',...
	'Position',[380 10 350 240],...
	'String',{},...
	'BackgroundColor','w',...
	'Callback',@clickpicks,...
	'Max',2,'tooltipstring','Selected images');

%{
openbut = uicontrol('Style','pushbutton',...
	'Position',[270 300 80 20],...
	'String','Open',...
	'Enable','off',...
	'Callback',@open);
%}

arrow = [2 2 2 2 2 2 2 2 1 2 2 2;...
         2 2 2 2 2 2 2 2 2 0 2 2;...
	     2 2 2 2 2 2 2 2 2 2 0 2;...
		 0 0 0 0 0 0 0 0 0 0 0 0;...
		 2 2 2 2 2 2 2 2 2 2 0 2;...
		 2 2 2 2 2 2 2 2 2 0 2 2;...
		 2 2 2 2 2 2 2 2 1 2 2 2];
arrow(arrow == 2) = NaN;
arrow_im = NaN*ones(16,76);
arrow_im(6:12,45:56) = arrow/2;
im = repmat(arrow_im,[1 1 3]);
addbut = uicontrol('Style','pushbutton',...
	'Position',[270 300 80 20],...
	'String','Add    ',...
	'Enable','off',...
	'CData',im,...
	'Callback',@add,'tooltipstring','Add images to selection','fontweight','bold','fontsize',10);

removebut = uicontrol('Style','pushbutton',...
	'Position',[290 205 80 20],...
	'String','Remove',...
	'Enable','off',...
	'Callback',@remove,'tooltipstring','Remove images from selection');
moveupbut = uicontrol('Style','pushbutton',...
	'Position',[290 175 80 20],...
	'String','Move Up',...
	'Enable','off',...
	'Callback',@moveup,'tooltipstring','Rearrange images');
movedownbut = uicontrol('Style','pushbutton',...
	'Position',[290 145 80 20],...
	'String','Move Down',...
	'Enable','off',...
	'Callback',@movedown,'tooltipstring','Rearrange images');

uicontrol('Position',[10 384 250 16],...
	'Style','text',...
	'String','Current Directory',...
	'HorizontalAlignment','center')
dir_popup = uicontrol('Style','popupmenu',...
	'Position',[10 335 250 20],...
	'BackgroundColor','w',...
	'String',path_cell(end:-1:1),...
	'Value',1,...
	'Callback',@dirpopup);
hist_cm = uicontextmenu;
pathbox = uicontrol('Position',[10 360 250 24],...
	'Style','edit',...
	'BackgroundColor','w',...
	'String',current_dir,...
	'HorizontalAlignment','left',...
	'Callback',@change_path,...
	'UIContextMenu',hist_cm,'tooltipstring','Current directory');
hist_menus = [];
hist_cb = @history_cb;
hist_menus = make_history_cm(hist_cb,hist_cm,hist_menus,history);

uicontrol('Position',[10 425+30 80 16],...
	'Style','text',...
	'String','File Filter',...
	'HorizontalAlignment','left',...
    'Visible', 'off')
showallfiles = uicontrol('Position',[10 405+40 200 20],...
	'Style','checkbox',...
	'String','Show All Files in Current Directory',...
	'Value',0,...
	'HorizontalAlignment','left',...
	'Callback',@togglefilter,'tooltipstring','Show also non-image files');
PIVlab_select = uicontrol('Position',[10 405+20 200 20],...
	'Style','radiobutton',...
	'String','Sequencing style 1-2, 2-3, 3-4, ...',...
	'Value',0,...
	'HorizontalAlignment','left',...
	'Callback',@PIVlab_select_Callback,'tooltipstring','Image sequencing style');
PIVlab_select2 = uicontrol('Position',[10 405-0 200 20],...
	'Style','radiobutton',...
	'String','Sequencing style 1-2, 3-4, 5-6, ...',...
	'Value',1,...
	'HorizontalAlignment','left',...
	'Callback',@PIVlab_select_Callback2,'tooltipstring','Image sequencing style');
filter_ed = uicontrol('Position',[10 405+30 80 20],...
	'Style','edit',...
	'BackgroundColor','w',...
	'String',filter,...
	'HorizontalAlignment','left',...
	'Callback',@setfilspec,...
     'Visible', 'off');
uicontrol('Position',[10 460+30 160 16],...
	'Style','text',...
	'String','Reg. Exp. Filter',...
	'HorizontalAlignment','left',...
    'Visible', 'on')
 
 refilter_ed = uicontrol('Position',[10 440+30 160 20],...
	'Style','edit',...
	'BackgroundColor','w',...
	'String',re_filter,...
	'HorizontalAlignment','left',...
	'Callback',@setrefilter,...
    'Visible', 'on','tooltipstring','Example:  \.m$|\.mat$  will list files ending in ''.m'' and ''.mat''.');
%z.b. insert  46 (| 47 ( for my files...

%new folder selection
dirbrowse = uicontrol('Style','pushbutton',...
	'Position',[263 362 30 20],...
	'String','...',...
	'Enable','on',...
	'Callback',@selectdire,'tooltipstring','Select directory','fontweight','bold');


if ~isempty(prop.type)
	type_value = 1;
	set(filter_ed,'String',prop.type{type_value,1},...
		'Enable','inactive')
	setfilspec()
	type_popup = uicontrol('Position',[10 405 250 20],...
		'Style','popupmenu',...
		'String',prop.type(:,2),...
		'Value',type_value,...
		'Callback',@filter_type_callback);
end

	function filter_type_callback(varargin)
		type_value = get(type_popup,'Value');
		set(filter_ed,'String',prop.type{type_value,1})
		setfilspec()
	end

viewfullpath = uicontrol('Style','checkbox',...
	'Position',[620 270 230 20],...
	'String','Show full paths',...
	'Value',show_full_path,...
	'HorizontalAlignment','left',...
	'Callback',@showfullpath,'tooltipstring','Show the full file path');
remove_dupes = uicontrol('Style','checkbox',...
	'Position',[380 270 230 20],...
	'String','Remove duplicates (as per full path)',...
	'Value',nodupes,...
	'HorizontalAlignment','left',...
	'Callback',@removedupes,'tooltipstring','Prevent images to appear twice in the list');
uicontrol('Position',[480 255 150 15],...
	'Style','text',...
	'String','Selected Files',...
	'HorizontalAlignment','center','FontWeight','bold')
uicontrol('Position',[280 80 80 30],'String','Import',...
	'Callback',@done,'tooltipstring','Finish selection and return to PIVlab','fontweight','bold','fontsize',10);
uicontrol('Position',[280 30 80 30],'String','Cancel',...
	'Callback',@cancel,'tooltipstring','Cancel image selection');

if ~isempty(prop.numfiles)
	uicontrol('Position',[620 255 100 16],...
		'Style','text',...
		'String',numstr,...
		'ForegroundColor','r',...
		'HorizontalAlignment','center')
end
    

ah1 = axes('units', 'pixels', 'Position',[380 295 350 270], 'xtick', [],'ytick',[]);
imshow(zeros(27,35)+0.75)%,'parent',ah1);

uicontrol('Position',[380 566 150 14],...
	'Style','text',...
	'String','Image preview',...
	'HorizontalAlignment','left','FontWeight','bold','tooltipstring','Image preview')


hgui=getappdata(0,'hgui');
setappdata(hgui, 'sequencer', 1); %default style

%set(fig,'HandleVisibility','off')

uiwait(fig)

% Compute desired output.
switch prop.output
	case 'cell'
		out = full_file_picks;
	case 'struct'
		out = dir_picks(:);
	case 'char'
		out = char(full_file_picks);
	case 'cancel'
		out = 0;
end


% -------------------- Callback functions --------------------

	function add(varargin)
		values = get(navlist,'Value');
        for i = 1:length(values)
            if fdir(values(i)).isdir==0
                dir_pick = fdir(values(i));
                pick = dir_pick.name;
                pick_full = fullfile(current_dir,pick);
                dir_pick.name = pick_full;
                if ~nodupes || ~any(strcmp(full_file_picks,pick_full))
                    file_picks{end + 1} = pick;
                    full_file_picks{end + 1} = pick_full;
                    dir_picks(end + 1) = dir_pick;
                end
            end
		end
		if show_full_path
			set(pickslist,'String',full_file_picks,'Value',[]);
		else
			set(pickslist,'String',file_picks,'Value',[]);
		end
		set([removebut,moveupbut,movedownbut],'Enable','off');
    end

    function selectdire(varargin)
        folder_name = uigetdir(get(pathbox,'String'),'Select folder');
        if folder_name ~= 0
            set(pathbox,'String',folder_name);
            change_path();
        end
    end

	function remove(varargin)
		values = get(pickslist,'Value');
		file_picks(values) = [];
		full_file_picks(values) = [];
		dir_picks(values) = [];
		top = get(pickslist,'ListboxTop');
		num_above_top = sum(values < top);
		top = top - num_above_top;
		num_picks = length(file_picks);
		new_value = min(min(values) - num_above_top,num_picks);
		if num_picks == 0
			new_value = [];
			set([removebut,moveupbut,movedownbut],'Enable','off')
		end
		if show_full_path
			set(pickslist,'String',full_file_picks,'Value',new_value,...
				'ListboxTop',top)
		else
			set(pickslist,'String',file_picks,'Value',new_value,...
				'ListboxTop',top)
		end
	end

	function open(varargin)
		values = get(navlist,'Value');
		if fdir(values).isdir
			if strcmp(fdir(values).name,'.')
				return
			elseif strcmp(fdir(values).name,'..')
				set(dir_popup,'Value',min(2,length(path_cell)))
				dirpopup();
				return
			end
			current_dir = fullfile(current_dir,fdir(values).name);
			history{end+1} = current_dir;
			history = unique(history);
			hist_menus = make_history_cm(hist_cb,hist_cm,hist_menus,...
				history);
			full_filter = fullfile(current_dir,filter);
			path_cell = path2cell(current_dir);
			fdir = filtered_dir(full_filter,re_filter);
			filenames = {fdir.name}';
			filenames = annotate_file_names(filenames,fdir);
			set(dir_popup,'String',path_cell(end:-1:1),'Value',1)
			set(pathbox,'String',current_dir)
			set(navlist,'ListboxTop',1,'Value',[],'String',filenames)
			set(addbut,'Enable','off')
			%set(openbut,'Enable','off')
        end
    end

    function clicknav(varargin)
        %hier preview einbauen...
        value = get(navlist,'Value');
        nval = length(value);
        dbl_click_fcn = @add;
        switch nval
            case 0
                set(addbut,'Enable','off')
            case 1
                set(addbut,'Enable','on');
                if fdir(value).isdir
                    %set(openbut,'Enable','on')
                    set(addbut,'Enable','off')
                    dbl_click_fcn = @open;
                else
%                    set(openbut,'Enable','off')
                    try
                    dir_pick = fdir(value(1));
                    pick = dir_pick.name;
                    pick_full = fullfile(current_dir,pick);
                    %dir_pick.name = pick_full;
                    
                    imshow (imread(pick_full))%,'parent',ah1);
                    catch 
                    end
                end
            otherwise %more than 1 file selected
                %if fdir(value(1)).isdir==1 || fdir(value(2)).isdir==1
                %    set(addbut,'Enable','off')
                %else
                    set(addbut,'Enable','on')
                %end

        end
        if strcmp(get(fig,'SelectionType'),'open')
            dbl_click_fcn();
        end

        
        
    end

	function clickpicks(varargin)
        value = get(pickslist,'Value');
		if isempty(value)
			set([removebut,moveupbut,movedownbut],'Enable','off')
		else
			set(removebut,'Enable','on')
			if min(value) == 1
				set(moveupbut,'Enable','off')
			else
                set(moveupbut,'Enable','on')
            end
            if max(value) == length(file_picks)
                set(movedownbut,'Enable','off')
            else
                set(movedownbut,'Enable','on')
            end
try
            pick = file_picks(value(1));
            pick_full = fullfile(current_dir,pick{1,1});
            figure(fig)
            
            imshow (imread(pick_full))%,'parent',ah1);

catch

end

		end
		if strcmp(get(fig,'SelectionType'),'open')
			remove();
		end
	end

	function dirpopup(varargin)
		value = get(dir_popup,'Value');
		len = length(path_cell);
		path_cell = path_cell(1:end-value+1);
		if ispc && value == len
			current_dir = '';
			full_filter = filter;
			drives = getdrives;
			num_drives = length(drives);
			temp = tempname;
			mkdir(temp)
			dir_temp = dir(temp);
			rmdir(temp)
			fdir = repmat(dir_temp(1),num_drives,1);
			[fdir.name] = deal(drives{:});
		else
			current_dir = cell2path(path_cell);
			history{end+1} = current_dir;
			history = unique(history);
			hist_menus = make_history_cm(hist_cb,hist_cm,hist_menus,...
				history);
			full_filter = fullfile(current_dir,filter);
			fdir = filtered_dir(full_filter,re_filter);
		end
		filenames = {fdir.name}';
		filenames = annotate_file_names(filenames,fdir);
		set(dir_popup,'String',path_cell(end:-1:1),'Value',1)
		set(pathbox,'String',current_dir)
		set(navlist,'String',filenames,'Value',[])
		set(addbut,'Enable','off')
	end

	function change_path(varargin)
		proposed_path = get(pathbox,'String');
		% Process any directories named '..'.
		proposed_path_cell = path2cell(proposed_path);
		ddots = strcmp(proposed_path_cell,'..');
		ddots(find(ddots) - 1) = true;
		proposed_path_cell(ddots) = [];
		proposed_path = cell2path(proposed_path_cell);
		% Check for existance of directory.
		if ~exist(proposed_path,'dir')
			uiwait(errordlg(['Directory "',proposed_path,...
				'" does not exist.'],'','modal'))
			return
		end
		current_dir = proposed_path;
		history{end+1} = current_dir;
		history = unique(history);
		hist_menus = make_history_cm(hist_cb,hist_cm,hist_menus,history);
		full_filter = fullfile(current_dir,filter);
		path_cell = path2cell(current_dir);
		fdir = filtered_dir(full_filter,re_filter);
		filenames = {fdir.name}';
		filenames = annotate_file_names(filenames,fdir);
		set(dir_popup,'String',path_cell(end:-1:1),'Value',1)
		set(pathbox,'String',current_dir)
		set(navlist,'String',filenames,'Value',[])
		set(addbut,'Enable','off')
		%set(openbut,'Enable','off')
	end

	function showfullpath(varargin)
		show_full_path = get(viewfullpath,'Value');
		if show_full_path
			set(pickslist,'String',full_file_picks)
		else
			set(pickslist,'String',file_picks)
		end
	end

	function removedupes(varargin)
		nodupes = get(remove_dupes,'Value');
		if nodupes
			num_picks = length(full_file_picks);
			[unused,rev_order] = unique(full_file_picks(end:-1:1));
			order = sort(num_picks + 1 - rev_order);
			full_file_picks = full_file_picks(order);
			file_picks = file_picks(order);
			dir_picks = dir_picks(order);
			if show_full_path
				set(pickslist,'String',full_file_picks,'Value',[])
			else
				set(pickslist,'String',file_picks,'Value',[])
			end
			set([removebut,moveupbut,movedownbut],'Enable','off')
		end
	end

	function moveup(varargin)
		value = get(pickslist,'Value');
		set(removebut,'Enable','on')
		n = length(file_picks);
		omega = 1:n;
		index = zeros(1,n);
		index(value - 1) = omega(value);
		index(setdiff(omega,value - 1)) = omega(setdiff(omega,value));
		file_picks = file_picks(index);
		full_file_picks = full_file_picks(index);
		dir_picks = dir_picks(index);
		value = value - 1;
		if show_full_path
			set(pickslist,'String',full_file_picks,'Value',value)
		else
			set(pickslist,'String',file_picks,'Value',value)
		end
		if min(value) == 1
			set(moveupbut,'Enable','off')
		end
		set(movedownbut,'Enable','on')
	end

	function movedown(varargin)
		value = get(pickslist,'Value');
		set(removebut,'Enable','on')
		n = length(file_picks);
		omega = 1:n;
		index = zeros(1,n);
		index(value + 1) = omega(value);
		index(setdiff(omega,value + 1)) = omega(setdiff(omega,value));
		file_picks = file_picks(index);
		full_file_picks = full_file_picks(index);
		dir_picks = dir_picks(index);
		value = value + 1;
		if show_full_path
			set(pickslist,'String',full_file_picks,'Value',value)
		else
			set(pickslist,'String',file_picks,'Value',value)
		end
		if max(value) == n
			set(movedownbut,'Enable','off')
		end
		set(moveupbut,'Enable','on')
	end

	function togglefilter(varargin)
		value = get(showallfiles,'Value');
		if value
			filter = '*';
			re_filter = '';
			set([filter_ed,refilter_ed],'Enable','off')
		else
			filter = get(filter_ed,'String');
			re_filter = get(refilter_ed,'String');
			set([filter_ed,refilter_ed],'Enable','on')
		end
		full_filter = fullfile(current_dir,filter);
		fdir = filtered_dir(full_filter,re_filter);
		filenames = {fdir.name}';
		filenames = annotate_file_names(filenames,fdir);
		set(navlist,'String',filenames,'Value',[])
		set(addbut,'Enable','off')
    end

    function PIVlab_select_Callback(varargin)
        hgui=getappdata(0,'hgui');
        if get(PIVlab_select,'Value')==1
            setappdata(hgui, 'sequencer', 0);
            set(PIVlab_select2, 'Value', 0)
        else
            setappdata(hgui, 'sequencer', 0);
            set(PIVlab_select, 'Value', 1)
            set(PIVlab_select2, 'Value', 0)
        end
        
    end
    function PIVlab_select_Callback2(varargin)
        hgui=getappdata(0,'hgui');
        if get(PIVlab_select2,'Value')==1
            setappdata(hgui, 'sequencer', 1);
            set(PIVlab_select, 'Value', 0)
        else
            setappdata(hgui, 'sequencer', 1);
            set(PIVlab_select2, 'Value', 1)
            set(PIVlab_select, 'Value', 0)
        end
        
    end

    function setfilspec(varargin)
		filter = get(filter_ed,'String');
		if isempty(filter)
			filter = '*';
			set(filter_ed,'String',filter)
		end
		% Process file spec if a subdirectory was included.
		[p,f,e] = fileparts(filter);
		if ~isempty(p)
			newpath = fullfile(current_dir,p,'');
			set(pathbox,'String',newpath)
			filter = [f,e];
			if isempty(filter)
				filter = '*';
			end
			set(filter_ed,'String',filter)
			change_path();
		end
		full_filter = fullfile(current_dir,filter);
		fdir = filtered_dir(full_filter,re_filter);
		filenames = {fdir.name}';
		filenames = annotate_file_names(filenames,fdir);
		set(navlist,'String',filenames,'Value',[])
		set(addbut,'Enable','off')
	end

	function setrefilter(varargin)
		re_filter = get(refilter_ed,'String');
		fdir = filtered_dir(full_filter,re_filter);
		filenames = {fdir.name}';
		filenames = annotate_file_names(filenames,fdir);
		set(navlist,'String',filenames,'Value',[])
		set(addbut,'Enable','off')
	end

	function done(varargin)
		% Optional shortcut: click on a file and press 'Done'.
 		if isempty(full_file_picks) && strcmp(get(addbut,'Enable'),'on')
 			add();
            drawnow;
            pause(0.75);
 		end
		numfiles = length(full_file_picks);
		if ~isempty(prop.numfiles)
			if numfiles < prop.numfiles(1)
				msg = {'Too few files selected.',numstr};
				uiwait(errordlg(msg,'','modal'))
				return
			elseif numfiles > prop.numfiles(end)
				msg = {'Too many files selected.',numstr};
				uiwait(errordlg(msg,'','modal'))
				return
			end
		end
		delete(fig)
	end

	function cancel(varargin)
		prop.output = 'cancel';
		delete(fig)
	end

	function history_cb(varargin)
		current_dir = history{varargin{3}};
		full_filter = fullfile(current_dir,filter);
		path_cell = path2cell(current_dir);
		fdir = filtered_dir(full_filter,re_filter);
		filenames = {fdir.name}';
		filenames = annotate_file_names(filenames,fdir);
		set(dir_popup,'String',path_cell(end:-1:1),'Value',1)
		set(pathbox,'String',current_dir)
		set(navlist,'ListboxTop',1,'Value',[],'String',filenames)
		set(addbut,'Enable','off')
		%set(openbut,'Enable','off')
	end
end


% -------------------- Subfunctions --------------------

function c = path2cell(p)
% Turns a path string into a cell array of path elements.
c = strread(p,'%s','delimiter','\\/');
if ispc
	c = [{'My Computer'};c];
else
	c = [{filesep};c(2:end)];
end
end


function p = cell2path(c)
% Turns a cell array of path elements into a path string.
if ispc
	p = fullfile(c{2:end},'');
else
	p = fullfile(c{:},'');
end
end


function d = filtered_dir(full_filter,re_filter)
% Like dir, but applies filters and sorting.
p = fileparts(full_filter);
if isempty(p) && full_filter(1) == '/'
	p = '/';
end
if exist(full_filter,'dir')
	dfiles = dir(' ');
else
	dfiles = dir(full_filter);
end
if ~isempty(dfiles)
	dfiles([dfiles.isdir]) = [];
end
ddir = dir(p);
ddir = ddir([ddir.isdir]);
% Additional regular expression filter.
if nargin > 1 && ~isempty(re_filter)
	if ispc || ismac
		no_match = cellfun('isempty',regexpi({dfiles.name},re_filter));
	else
		no_match = cellfun('isempty',regexp({dfiles.name},re_filter));
	end
	dfiles(no_match) = [];
end
% Set navigator style:
%	1 => mix file and directory names
%	2 => means list all files before all directories
%	3 => means list all directories before all files
%	4 => same as 2 except put . and .. directories first
if isunix
	style = 3;
else
	style = 4;
end
switch style
	case 1
		d = [dfiles;ddir];
		[unused,index] = sort({d.name});
		d = d(index);
	case 2
		[unused,index1] = sort({dfiles.name});
		[unused,index2] = sort({ddir.name});
		d = [dfiles(index1);ddir(index2)];
	case 3
		[unused,index1] = sort({dfiles.name});
		[unused,index2] = sort({ddir.name});
		d = [ddir(index2);dfiles(index1)];
	case 4
		[unused,index1] = sort({dfiles.name});
		dot1 = find(strcmp({ddir.name},'.'));
		dot2 = find(strcmp({ddir.name},'..'));
		ddot1 = ddir(dot1);
		ddot2 = ddir(dot2);
		ddir([dot1,dot2]) = [];
		[unused,index2] = sort({ddir.name});
		d = [ddot1;ddot2;dfiles(index1);ddir(index2)];
end
end


function drives = getdrives
% Returns a cell array of drive names on Windows.
letters = char('A':'Z');
num_letters = length(letters);
drives = cell(1,num_letters);
for i = 1:num_letters
	if exist([letters(i),':\'],'dir');
		drives{i} = [letters(i),':'];
	end
end
drives(cellfun('isempty',drives)) = [];
end


function filenames = annotate_file_names(filenames,dir_listing)
% Adds a trailing filesep character to directory names.
fs = filesep;
for i = 1:length(filenames)
	if dir_listing(i).isdir
		filenames{i} = [filenames{i},fs];
% 		filenames{i} = ['<html><font color="6495ed">',filenames{i},fs,'</font></html>'];
% 		filenames{i} = ['<html>',filenames{i},'<b style="color:bf0000">',fs,'</b></html>'];
% 		filenames{i} = ['<html>',filenames{i},...
% 			'<span style="font-family:Symbol">',fs,'</span></html>'];
% 		filenames{i} = ['<html>',filenames{i},'<font face="Times">',char(300),'</font></html>'];
% 		filenames{i} = ['<html><hr>',filenames{i},fs,'</hr></html>'];
% 		filenames{i} = ['<html>a<img src="icon.png">b',filenames{i},fs,'</html>'];
	end
end
end


function hist_menus = make_history_cm(cb,hist_cm,hist_menus,history)
% Make context menu for history.
if ~isempty(hist_menus)
	delete(hist_menus)
end
num_hist = length(history);
hist_menus = zeros(1,num_hist);
for i = 1:num_hist
	hist_menus(i) = uimenu(hist_cm,'Label',history{i},...
		'Callback',{cb,i});
end
end
