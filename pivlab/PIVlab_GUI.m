%PIVlab - Digital Particle Image Velocimetry Tool for MATLAB
%{
developed by Dr. William Thielicke and Prof. Dr. Eize J. Stamhuis
programmed with MATLAB Version 7.10 (R2010a)
March 09, 2010

third party content:
3-point gaussian sub-pixel estimator by Uri Shavit, Roi Gurka, Alex Liberzon
inpaint_nans by John D'Errico
uipickfiles by Douglas Schwarz
smoothn by Damien Garcia
dctn, idctn by Damien Garcia
Ellipse by D.G. Long
NaN Suite by Jan Gläscher
Exportfig by Ben Hinkle
mmstream2 by Duane Hanselman 
%}

%{
ToDo in next releases:
Load masks from another session?
plot correlation result in a Popup window?
subtract a flow field: Make sure that it has the same size, and give
information on the size, e.g. shoulb be 8x12 but is 10x12.?
colored vectors (I think there was a good submission in FEX)?
transparency, so that derivatives are displayed over the original image? http://www.mathworks.es/company/newsletters/articles/image-overlay-using-transparency.html ?
replace DECEV by q-criterion. http://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&ved=0CFQQFjAD&url=http%3A%2F%2Fwww.terpconnect.umd.edu%2F~lbravo%2Fdocs%2FQ-criterion.pdf&ei=1wrqUp3VNaOw4QSFzoCABQ&usg=AFQjCNGr6Ge9_H_e1in1niHiGuGBAfOFqw&sig2=Ig5eA0RNVRF4HOdnFExsAg&bvm=bv.60444564,d.bGE&cad=rja?

transparency: Convert the particle image and the derivatives image to RGB.
dann einfach beide mit hold... übereinander plotten. und dann transparenz
des oberen verringern. Sollte gehen, oder?
A=uint8(rand(100,100)*255);
map=colormap(jet(255))
RGB = ind2rgb(A,map)
image(RGB)
figure;imagesc(A);colormap('jet');
%}
function varargout = PIVlab_GUI(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PIVlab_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @PIVlab_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function PIVlab_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
%The opening function contains code that is executed just before the GUI is made visible to the user. You can access all the components for the GUI in the opening function, because all objects in the GUI are created before the opening function is called. You can add code to the opening function to perform tasks that need to be done before the user has access to the GUI -- for example,
% It is executed just before the GUI is made visible to the user, but after
% all the components have been created, i.e., after the components' CreateFcn callbacks, if any, have been run.
handles.output = hObject;
guidata(hObject, handles);
handles=guihandles(hObject);
setappdata(0,'hgui',gcf);
%clc
try
    if verLessThan('matlab', '7.10.0') == 0
        disp('-> Matlab version check ok.')
    else
        disp('WARNING: Your Matlab version might be too old for running PIVlab.')
    end
catch
    disp('MATLAB version could not be checked automatically. You need at least version 7.10.0 (R2010a) to run PIVlab.')
end
try
    result=license('checkout','Image_Toolbox');
    if result == 1
        try
            J = adapthisteq(rand(8,8));
            disp('-> Image Processing Toolbox found.')
        catch
            disp('ERROR: Image Processing Toolbox not found! PIVlab won''t work like this.')
            result=0;
        end
    else
        disp('ERROR: Image Processing Toolbox not found! PIVlab won''t work like this.')
    end
    
    
    ctr=0;
    pivFiles = {'PIVlab_GUI.fig' 'dctn.m' 'idctn.m' 'inpaint_nans.m' 'piv_DCC.m' 'piv_FFTmulti.m' 'PIVlab_preproc.m' 'PIVlablogo.jpg' 'smoothn.m' 'uipickfiles.m' 'PIVlab_settings_default.mat' 'hsbmap.mat' 'parula.mat' 'ellipse.m' 'nanmax.m' 'nanmin.m' 'nanstd.m' 'nanmean.m' 'exportfig.m' 'fastLICFunction.m' 'icons.mat' 'mmstream2.m' 'PIVlab_citing.fig' 'PIVlab_citing.m'};
    for i=1:size(pivFiles,2)
        if exist(pivFiles{1,i})~=2
            disp(['ERROR: A required file was not found: ' pivFiles{1,i}]);
            beep;
        else
            ctr=ctr+1;
        end
    end
    if ctr==size(pivFiles,2)
        disp('-> Required additional files found.')
    end
catch
    result=0;
    disp('Toolboxes could not be checked automatically. You need the Image Processing Toolbox')
end
set (hObject, 'position', [329.25 148.5 605 470]); %standard size in points
drawnow;
movegui(hObject,'center')
%Variable initialization
put('PIVver', '1.41');
put ('toggler',0);
put('caluv',1);
put('calxy',1);
put('subtr_u', 0);
put('subtr_v', 0);
put('displaywhat',1);%vectors
put('imgproctoolbox',result);
set(gcf, 'Name',['PIVlab ' retr('PIVver') ' by W. Thielicke and E.J. Stamhuis'])
if result==1
    msg={'Welcome to PIVlab!';...
        ['version: ' retr('PIVver')];...
        '';...
        'Start by selecting';...
        '"File" -> "New session"';...
        'from the menu. Load';...
        'your PIV images by';...
        'clicking "Load images" on the';...
        'right hand side.';...
        '';...
        'Then, work your way through';...
        'the menu from left to right.';...
        };
else
    msg={'!!! WARNING !!! The Image Processing toolbox was not detected on your system. PIVlab will most likely not work like this. You definetively need this toolbox to run PIVlab!'};
end
set(handles.text6,'String', msg);
%read current and last directory.....:
try
    lastdir=importdata('last.nf','\t');
    put('homedir',lastdir{1});
    put('pathname',lastdir{2});
catch
    try
        lastdir{1}=pwd;
        lastdir{2}=pwd;
        dlmwrite('last.nf', lastdir{1}, 'delimiter', '', 'precision', 6, 'newline', 'pc')
        dlmwrite('last.nf', lastdir{2}, '-append', 'delimiter', '', 'precision', 6, 'newline', 'pc')
        put('pathname',lastdir{2});
    catch
    end
end
try
    %XP Wu modification:
    psdfile=which('PIVlab_settings_default.mat');
    dindex=strfind(psdfile,filesep);
    read_settings ('PIVlab_settings_default.mat',psdfile(1:(dindex(end)-1)));
catch
    disp('Could not load default settings.')
end

load icons.mat
set(handles.zoomon, 'cdata',zoompic);
set(handles.panon, 'cdata',panpic);

set(findobj('tag','uipanel28'),'visible','off'); %disable the "skip" functionality. I think noone has ever used it...

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>

function varargout = PIVlab_GUI_OutputFcn(hObject, eventdata, handles) %#ok<*INUSL>
varargout{1} = handles.output;
displogo(1)

function displogo(zoom)
logoimg=imread('PIVlablogo.jpg');
if zoom==1
    h=image(logoimg+255, 'parent', gca);
    axis image;
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    set(gca, 'xlim', [1 size(logoimg,2)]);
    set(gca, 'ylim', [1 size(logoimg,1)]);
    set(gca, 'ydir', 'reverse');
    set(gca, 'xcolor', [0.94 0.94 0.94], 'ycolor', [0.94 0.94 0.94]) ;
    for i=255:-10:0
        RGB2=logoimg+i;
        try
            set (h, 'cdata', RGB2);
        catch %#ok<*CTCH>
            disp('.')
        end
        drawnow expose;
    end
end
%get(gca,'position')
image(logoimg, 'parent', gca);
set(gca, 'xcolor', [0.94 0.94 0.94], 'ycolor', [0.94 0.94 0.94]) ;

axis image;
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca, 'xlim', [1 size(logoimg,2)]);
set(gca, 'ylim', [1 size(logoimg,1)]);

set(gca, 'ydir', 'reverse');
text (440,450,['version: ' retr('PIVver')], 'fontsize', 8,'fontangle','italic');
imgproctoolbox=retr('imgproctoolbox');
put('imgproctoolbox',[]);
if imgproctoolbox==0
    text (90,200,'Image processing toolbox not found!', 'fontsize', 16, 'color', [1 0 0], 'backgroundcolor', [0 0 0]);
end

function switchui (who)
handles=guihandles(getappdata(0,'hgui')); %#ok<*NASGU>

if get(handles.zoomon,'Value')==1
    set(handles.zoomon,'Value',0);
    zoomon_Callback(handles.zoomon)
end
if get(handles.panon,'Value')==1
    set(handles.panon,'Value',0);
    panon_Callback(handles.panon)
end

turnoff=findobj('-regexp','Tag','multip');
set(turnoff, 'visible', 'off');
turnon=findobj('-regexp','Tag',who);
set(turnon, 'visible', 'on');
drawnow;

function put(name, what)
hgui=getappdata(0,'hgui');
setappdata(hgui, name, what);

function var = retr(name)
hgui=getappdata(0,'hgui');
var=getappdata(hgui, name);

function handles=gethand
hgui=getappdata(0,'hgui');
handles=guihandles(hgui);

function sliderdisp
handles=gethand;
toggler=retr('toggler');
selected=2*floor(get(handles.fileselector, 'value'))-(1-toggler);
filepath=retr('filepath');
%if the images are not found on the current path, then let user choose new path
%not found: assign new path to all following elements.
%check next file. not found -> assign new path to all following.
%and so on...
%checking if all files exist takes 0.5 s each time... need for optimization
%e.g. do this only one time at the start.
if isempty(filepath) == 0 && exist(filepath{selected},'file') ~=2
    for i=1:size(filepath,1)
        while exist(filepath{i,1},'file') ~=2
            errordlg(['The image ' sprintf('\n') filepath{i,1} sprintf('\n') '(and probably some more...) could not be found.' sprintf('\n') 'Please select the path where the images are located.'],'File not found!','on')
            uiwait
            new_dir = uigetdir(pwd,'Please specify the path to all the images');
            if new_dir==0
                break
            else
                for j=i:size(filepath,1) %apply new path to all following imgs.
                    if ispc==1
                        zeichen=findstr('\',filepath{j,1});
                    else
                        zeichen=findstr('/',filepath{j,1});
                    end
                    currentobject=filepath{j,1};
                    currentpath=currentobject(1:(zeichen(1,size(zeichen,2))));
                    currentfile=currentobject(zeichen(1,size(zeichen,2))+1:end);
                    if ispc==1
                        filepath{j,1}=[new_dir '\' currentfile];
                    else
                        filepath{j,1}=[new_dir '/' currentfile];
                    end
                end
            end
            put('filepath',filepath);
        end
        if new_dir==0
            break
        end
    end
end

currentframe=2*floor(get(handles.fileselector, 'value'))-1;
%display derivatives if available and desired...
displaywhat=retr('displaywhat');
delete(findobj('tag', 'derivhint'));
if size(filepath,1)>0
    
    
    if get(handles.zoomon,'Value')==1
        set(handles.zoomon,'Value',0);
        zoomon_Callback(handles.zoomon)
    end
    if get(handles.panon,'Value')==1
        set(handles.panon,'Value',0);
        panon_Callback(handles.panon)
    end
    xzoomlimit=retr('xzoomlimit');
    yzoomlimit=retr('yzoomlimit');
    
    derived=retr('derived');
    if isempty(derived)==0   %derivatives were calculated
        %derived=retr('derived');
        %1=vectors only
        if displaywhat==1 %vectors only
            currentimage=imread(filepath{selected});
            image(currentimage, 'parent',gca, 'cdatamapping', 'scaled');
            colormap('gray');
            
            vectorcolor=[str2double(get(handles.validr,'string')) str2double(get(handles.validg,'string')) str2double(get(handles.validb,'string'))];
            %vectorcolor='g';
            %end
        else %displaywhat>1
            if size(derived,2)>=(currentframe+1)/2 && numel(derived{displaywhat-1,(currentframe+1)/2})>0 %derived parameters requested and existant
                currentimage=derived{displaywhat-1,(currentframe+1)/2};
                %is currentimage 3d? That would cause problems.-....
                
                %pcolor(resultslist{1,(currentframe+1)/2},resultslist{2,(currentframe+1)/2},currentimage);shading interp;
                image(rescale_maps(currentimage), 'parent',gca, 'cdatamapping', 'scaled');
                if displaywhat ~=10 %10 is LIC
                    
                    avail_maps=get(handles.colormap_choice,'string');
                    selected_index=get(handles.colormap_choice,'value');
                    if selected_index == 4 %HochschuleBremen map
                        load hsbmap.mat;
                        colormap(hsb);
                    elseif selected_index== 1 %rainbow
                        %load rainbow.mat;
                        %colormap (rainbow);
                        load parula.mat;
                        colormap (parula);
                    else
                        colormap(avail_maps{selected_index});
                    end
                else
                    colormap('gray');
                end
                if get(handles.autoscaler,'value')==1
                    minscale=min(min(currentimage));
                    maxscale=max(max(currentimage));
                    set (handles.mapscale_min, 'string', num2str(minscale))
                    set (handles.mapscale_max, 'string', num2str(maxscale))
                else
                    minscale=str2double(get(handles.mapscale_min, 'string'));
                    maxscale=str2double(get(handles.mapscale_max, 'string'));
                end
                caxis([minscale maxscale])
                vectorcolor=[str2double(get(handles.validdr,'string')) str2double(get(handles.validdg,'string')) str2double(get(handles.validdb,'string'))];
                %vectorcolor='k';
                if get(handles.displ_colorbar,'value')==1
                    name=get(handles.derivchoice,'string');
                    posichoice = get(handles.colorbarpos,'String');
                    coloobj=colorbar (posichoice{get(handles.colorbarpos,'Value')},'FontWeight','bold','Fontsize',12);
                    if strcmp(posichoice{get(handles.colorbarpos,'Value')},'East')==1 | strcmp(posichoice{get(handles.colorbarpos,'Value')},'West')==1
                        set(coloobj,'YTickLabel',num2str(get(coloobj,'YTick')','%5.5g'))
                        ylabel(coloobj,name{retr('displaywhat')},'fontsize',9,'fontweight','normal');
                    end
                    if strcmp(posichoice{get(handles.colorbarpos,'Value')},'North')==1 | strcmp(posichoice{get(handles.colorbarpos,'Value')},'South')==1
                        set(coloobj,'XTickLabel',num2str(get(coloobj,'XTick')','%5.5g'))
                        xlabel(coloobj,name{retr('displaywhat')},'fontsize',9,'fontweight','normal');
                    end
                end
            else %no deriv available
                currentimage=imread(filepath{selected});
                image(currentimage, 'parent',gca, 'cdatamapping', 'scaled');
                colormap('gray');
                vectorcolor=[str2double(get(handles.validr,'string')) str2double(get(handles.validg,'string')) str2double(get(handles.validb,'string'))];
                %vectorcolor='g';
                text(10,10,'This parameter needs to be calculated for this frame first. Go to Plot -> Derive Parameters and click "Apply to current frame".','color','r','fontsize',9, 'BackgroundColor', 'k', 'tag', 'derivhint')
            end
        end
    else %not in derivatives panel
        try
            currentimage=imread(filepath{selected});
        catch
            disp(['Error: ' filepath{selected} ' --> Image could not be found!']);
            resultslist=retr('resultslist');
            maximgx=max(max(resultslist{1,1}))+min(min(resultslist{1,1}));
            maximgy=max(max(resultslist{2,1}))+min(min(resultslist{2,1}));
            currentimage=zeros(maximgy,maximgx);
        end
        image(currentimage, 'parent',gca, 'cdatamapping', 'scaled');
        colormap('gray');
        vectorcolor=[str2double(get(handles.validr,'string')) str2double(get(handles.validg,'string')) str2double(get(handles.validb,'string'))];
        %vectorcolor='g';
    end
    axis image;
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    filename=retr('filename');
    
    ismean=retr('ismean');
    if size(ismean,1)>=(currentframe+1)/2
        if ismean((currentframe+1)/2,1) ==1
            currentwasmean=1;
        else
            currentwasmean=0;
        end
    else
        currentwasmean=0;
    end
    
    if currentwasmean==1
        set (handles.filenameshow,'BackgroundColor',[0.65 0.65 1]);
    else
        set (handles.filenameshow,'BackgroundColor',[0.9412 0.9412 0.9412]);
    end
        
    set (handles.filenameshow, 'string', ['Frame (' int2str(floor(get(handles.fileselector, 'value'))) '/' int2str(size(filepath,1)/2) '):' sprintf('\n') filename{selected}]);
    set (handles.filenameshow, 'tooltipstring', filepath{selected});
    if strmatch(get(handles.multip01, 'visible'), 'on');
        set(handles.imsize, 'string', ['Image size: ' int2str(size(currentimage,2)) '*' int2str(size(currentimage,1)) 'px' ])
    end
    maskiererx=retr('maskiererx');
    if size(maskiererx,2)>=currentframe
        ximask=maskiererx{1,currentframe};
        if size(ximask,1)>1
            if displaywhat == 1 %%when vectors only are display: transparent mask
                dispMASK(0.333)
            else %otherwise: 100% opaque mask.
                dispMASK(1)
            end
        end
    end
    roirect=retr('roirect');
    if size(roirect,2)>1
        dispROI
    end
    resultslist=retr('resultslist');
    delete(findobj('tag', 'smoothhint'));
    
    
    %manualmarkers
    if get(handles.displmarker,'value')==1
        manmarkersX=retr('manmarkersX');
        manmarkersY=retr('manmarkersY');
        delete(findobj('tag','manualmarker'));
        if numel(manmarkersX)>0
            hold on
            plot(manmarkersX,manmarkersY, 'o','MarkerEdgeColor','k','MarkerFaceColor',[.2 .2 1], 'MarkerSize',9, 'tag', 'manualmarker');
            plot(manmarkersX,manmarkersY, '*','MarkerEdgeColor','w', 'tag', 'manualmarker');
            hold off
        end
    end
    
    
    if size(resultslist,2)>=(currentframe+1)/2 && numel(resultslist{1,(currentframe+1)/2})>0
        x=resultslist{1,(currentframe+1)/2};
        y=resultslist{2,(currentframe+1)/2};
        if size(resultslist,1)>6 %filtered exists
            if size(resultslist,1)>10 && numel(resultslist{10,(currentframe+1)/2}) > 0 %smoothed exists
                u=resultslist{10,(currentframe+1)/2};
                v=resultslist{11,(currentframe+1)/2};
                typevector=resultslist{9,(currentframe+1)/2};
                text(3,size(currentimage,1)-4, 'Smoothed dataset','tag', 'smoothhint', 'backgroundcolor', 'k', 'color', 'y','fontsize',6);
                if numel(typevector)==0 %happens if user smoothes sth without NaN and without validation
                    typevector=resultslist{5,(currentframe+1)/2};
                end
            else
                u=resultslist{7,(currentframe+1)/2};
                if size(u,1)>1
                    v=resultslist{8,(currentframe+1)/2};
                    typevector=resultslist{9,(currentframe+1)/2};
                else %filter was applied for other frames but not for this one
                    u=resultslist{3,(currentframe+1)/2};
                    v=resultslist{4,(currentframe+1)/2};
                    typevector=resultslist{5,(currentframe+1)/2};
                end
            end
        else
            u=resultslist{3,(currentframe+1)/2};
            v=resultslist{4,(currentframe+1)/2};
            typevector=resultslist{5,(currentframe+1)/2};
        end
        if get(handles.highp_vectors, 'value')==1 & strmatch(get(handles.multip08, 'visible'), 'on') %#ok<AND2>
            strength=54-round(get(handles.highpass_strength, 'value'));
            h = fspecial('gaussian',strength,strength) ;
            h2= fspecial('gaussian',3,3);
            ubg=imfilter(u,h,'replicate');
            vbg=imfilter(v,h,'replicate');
            ufilt=u-ubg;
            vfilt=v-vbg;
            u=imfilter(ufilt,h2,'replicate');
            v=imfilter(vfilt,h2,'replicate');
        end
        autoscale_vec=get(handles.autoscale_vec, 'Value');
        vecskip=str2double(get(handles.nthvect,'String'));
        if autoscale_vec == 1
            autoscale=1;
            %from quiver autoscale function:
            if min(size(x))==1, n=sqrt(numel(x)); m=n; else [m,n]=size(x); end
            delx = diff([min(x(:)) max(x(:))])/n;
            dely = diff([min(y(:)) max(y(:))])/m;
            del = delx.^2 + dely.^2;
            if del>0
                len = sqrt((u.^2 + v.^2)/del);
                maxlen = max(len(:));
            else
                maxlen = 0;
            end
            if maxlen>0
                autoscale = autoscale/ maxlen * vecskip;
            else
                autoscale = autoscale; %#ok<*ASGSL>
            end
            vecscale=autoscale;
        else %autoscale off
            vecscale=str2num(get(handles.vectorscale,'string')); %#ok<*ST2NM>
        end
        hold on;
        
        vectorcolorintp=[str2double(get(handles.interpr,'string')) str2double(get(handles.interpg,'string')) str2double(get(handles.interpb,'string'))];
        if vecskip==1
            q=quiver(x(typevector==1),y(typevector==1),...
                (u(typevector==1)-(retr('subtr_u')/retr('caluv')))*vecscale,...
                (v(typevector==1)-(retr('subtr_v')/retr('caluv')))*vecscale,...
                'Color', vectorcolor,'autoscale', 'off','linewidth',str2double(get(handles.vecwidth,'string')));
            q2=quiver(x(typevector==2),y(typevector==2),...
                (u(typevector==2)-(retr('subtr_u')/retr('caluv')))*vecscale,...
                (v(typevector==2)-(retr('subtr_v')/retr('caluv')))*vecscale,...
                'Color', vectorcolorintp,'autoscale', 'off','linewidth',str2double(get(handles.vecwidth,'string')));
            scatter(x(typevector==0),y(typevector==0),'rx') %masked
        else
            typevector_reduced=typevector(1:vecskip:end,1:vecskip:end);
            x_reduced=x(1:vecskip:end,1:vecskip:end);
            y_reduced=y(1:vecskip:end,1:vecskip:end);
            u_reduced=u(1:vecskip:end,1:vecskip:end);
            v_reduced=v(1:vecskip:end,1:vecskip:end);
            q=quiver(x_reduced(typevector_reduced==1),y_reduced(typevector_reduced==1),...
                (u_reduced(typevector_reduced==1)-(retr('subtr_u')/retr('caluv')))*vecscale,...
                (v_reduced(typevector_reduced==1)-(retr('subtr_v')/retr('caluv')))*vecscale,...
                'Color', vectorcolor,'autoscale', 'off','linewidth',str2double(get(handles.vecwidth,'string')));
            q2=quiver(x_reduced(typevector_reduced==2),y_reduced(typevector_reduced==2),...
                (u_reduced(typevector_reduced==2)-(retr('subtr_u')/retr('caluv')))*vecscale,...
                (v_reduced(typevector_reduced==2)-(retr('subtr_v')/retr('caluv')))*vecscale,...
                'Color', vectorcolorintp,'autoscale', 'off','linewidth',str2double(get(handles.vecwidth,'string')));
            scatter(x_reduced(typevector_reduced==0),y_reduced(typevector_reduced==0),'rx') %masked
        end
        hold off;
        %streamlines:
        streamlinesX=retr('streamlinesX');
        streamlinesY=retr('streamlinesY');
        delete(findobj('tag','streamline'));
        if numel(streamlinesX)>0
            ustream=u*retr('caluv')-retr('subtr_u');
            vstream=v*retr('caluv')-retr('subtr_v');
            ustream(typevector==0)=nan;
            vstream(typevector==0)=nan;
            h=streamline(mmstream2(x,y,ustream,vstream,streamlinesX,streamlinesY,'on'));
            set (h,'tag','streamline');
            contents = get(handles.streamlcolor,'String');
            set(h,'LineWidth',get(handles.streamlwidth,'value'),'Color', contents{get(handles.streamlcolor,'Value')});
        end
        
        if verLessThan('matlab','8.4')
            set(q, 'ButtonDownFcn', @veclick, 'hittestarea', 'on');
            set(q2, 'ButtonDownFcn', @veclick, 'hittestarea', 'on');
           else
            % >R2014a
            set(q, 'ButtonDownFcn', @veclick, 'PickableParts', 'visible');
            set(q2, 'ButtonDownFcn', @veclick, 'PickableParts', 'visible'); 
        end
        
        if strmatch(get(handles.multip14, 'visible'), 'on'); %statistics panel visible
            update_Stats (x,y,u,v);
        end
        if strmatch(get(handles.multip06, 'visible'), 'on'); %validation panel visible
            manualdeletion=retr('manualdeletion');
            frame=floor(get(handles.fileselector, 'value'));
            framemanualdeletion=[];
            if numel(manualdeletion)>0
                if size(manualdeletion,2)>=frame
                    if isempty(manualdeletion{1,frame}) ==0
                        framemanualdeletion=manualdeletion{frame};
                    end
                end
            end
            if isempty(framemanualdeletion)==0
                
                
                hold on;
                for i=1:size(framemanualdeletion,1)
                    scatter (x(framemanualdeletion(i,1),framemanualdeletion(i,2)),y(framemanualdeletion(i,1),framemanualdeletion(i,2)), 'rx', 'tag','manualdot')
                end
                hold off;
            end
        end
        
        %{
        figure;
        [Vx2,Vy2] = pppiv(u,v);
        quiver(Vx2,Vy2)
        %}
    end
    
    if isempty(xzoomlimit)==0
        set(gca,'xlim',xzoomlimit)
        set(gca,'ylim',yzoomlimit)
    end
    if verLessThan('matlab','8.4')
        %do nothing
    else
        % >R2014a
        set(gca,'YlimMode','manual');set(gca,'XlimMode','manual') %in r2014b, vectors are not clipped when set to auto... (?!?)
    end
    drawnow;
end

function update_Stats(x,y,u,v)
handles=gethand;
caluv=retr('caluv');
calxy=retr('calxy');
x=reshape(x,size(x,1)*size(x,2),1);
y=reshape(y,size(y,1)*size(y,2),1);
u=reshape(u,size(u,1)*size(u,2),1);
v=reshape(v,size(v,1)*size(v,2),1);
if retr('caluv')==1 && retr('calxy')==1
    set (handles.meanu,'string', [num2str(nanmean(u*caluv)) ' ± ' num2str(nanstd(u*caluv)) ' [px/frame]'])
    set (handles.meanv,'string', [num2str(nanmean(v*caluv)) ' ± ' num2str(nanstd(v*caluv)) ' [px/frame]'])
else
    set (handles.meanu,'string', [num2str(nanmean(u*caluv)) ' ± ' num2str(nanstd(u*caluv)) ' [m/s]'])
    set (handles.meanv,'string', [num2str(nanmean(v*caluv)) ' ± ' num2str(nanstd(v*caluv)) ' [m/s]'])
end

function veclick(src,eventdata)
%only active if vectors are displayed.
handles=gethand;
currentframe=2*floor(get(handles.fileselector, 'value'))-1;
resultslist=retr('resultslist');
x=resultslist{1,(currentframe+1)/2};
y=resultslist{2,(currentframe+1)/2};
pos=get(gca,'CurrentPoint');
xposition=round(pos(1,1));
yposition=round(pos(1,2));
findx=abs(x/xposition-1);
[trash, imagex]=find(findx==min(min(findx)));
findy=abs(y/yposition-1);
[imagey, trash]=find(findy==min(min(findy)));
info(1,1)=imagey(1,1);
info(1,2)=imagex(1,1);
%LOAD INTERPOLATED RESULT IF EXISTENT
if size(resultslist,1)>6 %filtered exists
    u=resultslist{7,(currentframe+1)/2};
    typevector=resultslist{5,(currentframe+1)/2};
    if numel(u)>0
        v=resultslist{8,(currentframe+1)/2};
    else
        u=resultslist{3,(currentframe+1)/2};
        v=resultslist{4,(currentframe+1)/2};
    end
else
    u=resultslist{3,(currentframe+1)/2};
    v=resultslist{4,(currentframe+1)/2};
    typevector=resultslist{5,(currentframe+1)/2};
end
if typevector(info(1,1),info(1,2)) ~=0
    delete(findobj('tag', 'infopoint'));
    %here, the calibration matters...
    
    if retr('caluv')==1 && retr('calxy')==1
        set(handles.u_cp, 'String', ['u:' num2str(round((u(info(1,1),info(1,2))*retr('caluv')-retr('subtr_u'))*10000)/10000) ' [px/fr]']);
        set(handles.v_cp, 'String', ['v:' num2str(round((v(info(1,1),info(1,2))*retr('caluv')-retr('subtr_v'))*10000)/10000) ' [px/fr]']);
        set(handles.x_cp, 'String', ['x:' num2str(round((x(info(1,1),info(1,2))*retr('calxy'))*1000)/1000) ' [px]']);
        set(handles.y_cp, 'String', ['y:' num2str(round((y(info(1,1),info(1,2))*retr('calxy'))*1000)/1000) ' [px]']);
    else
        set(handles.u_cp, 'String', ['u:' num2str(round((u(info(1,1),info(1,2))*retr('caluv')-retr('subtr_u'))*10000)/10000) ' [m/s]']);
        set(handles.v_cp, 'String', ['v:' num2str(round((v(info(1,1),info(1,2))*retr('caluv')-retr('subtr_v'))*10000)/10000) ' [m/s]']);
        set(handles.x_cp, 'String', ['x:' num2str(round((x(info(1,1),info(1,2))*retr('calxy'))*1000)/1000) ' [m]']);
        set(handles.y_cp, 'String', ['y:' num2str(round((y(info(1,1),info(1,2))*retr('calxy'))*1000)/1000) ' [m]']);
    end
    derived=retr('derived');
    displaywhat=retr('displaywhat');
    if displaywhat>1
        if size (derived,2) >= (currentframe+1)/2
            if numel(derived{displaywhat-1,(currentframe+1)/2})>0
                map=derived{displaywhat-1,(currentframe+1)/2};
                name=get(handles.derivchoice,'string');
                try
                    set(handles.scalar_cp, 'String', [name{displaywhat} ': ' num2str(round(map(info(1,1),info(1,2))*10000)/10000)]);
                catch
                    plot_derivs_Callback
                    name=get(handles.derivchoice,'string');
                    set(handles.scalar_cp, 'String', [name{displaywhat} ': ' num2str(round(map(info(1,1),info(1,2))*10000)/10000)]);
                end
            else
                set(handles.scalar_cp, 'String','N/A');
            end
        else
            set(handles.scalar_cp, 'String','N/A');
        end
    else
        set(handles.scalar_cp, 'String','N/A');
    end
    
    hold on;
    plot(x(info(1,1),info(1,2)),y(info(1,1),info(1,2)), 'yo', 'tag', 'infopoint','linewidth', 1, 'markersize', 10);
    hold off;
end

function toolsavailable(inpt);
%0: disable all tools
%1: re-enable tools that were previously also enabled
hgui=getappdata(0,'hgui');
handles=gethand;
if inpt==0
    if get(handles.zoomon,'Value')==1
        set(handles.zoomon,'Value',0);
        zoomon_Callback(handles.zoomon)
    end
    if get(handles.panon,'Value')==1
        set(handles.panon,'Value',0);
        panon_Callback(handles.panon)
    end
end

elementsOfCrime=findobj(hgui, 'type', 'uicontrol');
elementsOfCrime2=findobj(hgui, 'type', 'uimenu');
statuscell=get (elementsOfCrime, 'enable');
wasdisabled=zeros(size(statuscell),'uint8');

if inpt==0
    set(elementsOfCrime, 'enable', 'off');
    for i=1:size(statuscell,1)
        if strmatch(statuscell{i,1}, 'off') ==1
            wasdisabled(i)=1;
        end
    end
    put('wasdisabled', wasdisabled);
    set(elementsOfCrime2, 'enable', 'off');
else
    wasdisabled=retr('wasdisabled');
    set(elementsOfCrime, 'enable', 'on');
    set(elementsOfCrime(wasdisabled==1), 'enable', 'off');
    set(elementsOfCrime2, 'enable', 'on');
end
set(handles.progress, 'enable', 'on');
set(handles.overall, 'enable', 'on');
set(handles.totaltime, 'enable', 'on');
set(handles.messagetext, 'enable', 'on');

function overlappercent
handles=gethand;
perc=100-str2double(get(handles.step,'string'))/str2double(get(handles.intarea,'string'))*100;
set (handles.steppercentage, 'string', ['= ' int2str(perc) '%']);

function figure1_ResizeFcn(hObject, eventdata, handles)
Figure_Size = get(hObject, 'Position');
%{
minimalwidth=801;
minimalheight=610;
panelwidth=188;
panelheight=484;
toolheight=120;
%}

%in points
minimalwidth=605;
minimalheight=470;
panelwidth=141;
panelheight=363;
toolheight=110;

if  Figure_Size(4)<minimalheight
    try
        set(hObject,'position', [Figure_Size(1) Figure_Size(2) Figure_Size(3) minimalheight+3]);
    catch
    end
end

handles=guihandles(hObject);
try
    %set (findobj('-regexp','Tag','multip'), 'position', [Figure_Size(3)-panelwidth Figure_Size(4)-panelheight-5 panelwidth panelheight]);
    set (findobj('-regexp','Tag','multip'), 'position', [Figure_Size(3)-panelwidth Figure_Size(4)-panelheight-3.75 panelwidth panelheight]);
    set (handles.tools, 'position', [Figure_Size(3)-panelwidth Figure_Size(4)-panelheight-toolheight panelwidth toolheight]);
    
    %set (gca, 'position', [5 5 Figure_Size(3)-198 Figure_Size(4)-10]);
    set (gca, 'position', [3.75 3.75 Figure_Size(3)-148.5 Figure_Size(4)-7.5]);
catch
end

function loadimgsbutton_Callback(hObject, eventdata, handles)
if ispc==1
    pathname=[retr('pathname') '\'];
else
    pathname=[retr('pathname') '/'];
end

handles=gethand;


displogo(0)
if ispc==1
    path=uipickfiles ('FilterSpec', pathname, 'REFilter', '\.bmp$|\.jpg$|\.tif$', 'numfiles', [2 inf], 'output', 'struct', 'prompt', 'Select images. Images from one set should have identical dimensions to avoid problems.');
else
    path=uipickfiles ('FilterSpec', pathname, 'numfiles', [2 inf], 'output', 'struct', 'prompt', 'Select images. Images from one set should have identical dimensions to avoid problems.');
end
if isequal(path,0) ==0
    
    if get(handles.zoomon,'Value')==1
        set(handles.zoomon,'Value',0);
        zoomon_Callback(handles.zoomon)
    end
    if get(handles.panon,'Value')==1
        set(handles.panon,'Value',0);
        panon_Callback(handles.zoomon)
    end
    put('xzoomlimit',[]);
    put('yzoomlimit',[]);
    
    sequencer=retr('sequencer');
    if sequencer==1
        for i=1:size(path,1)
            if path(i).isdir == 0 %remove directories from selection
                if exist('filepath','var')==0 %first loop
                    filepath{1,1}=path(i).name;
                else
                    filepath{size(filepath,1)+1,1}=path(i).name;
                end
            end
        end
    else %sequencer=0
        for i=1:size(path,1)
            if path(i).isdir == 0 %remove directories from selection
                if exist('filepath','var')==0 %first loop
                    filepath{1,1}=path(i).name;
                else
                    filepath{size(filepath,1)+1,1}=path(i).name;
                    filepath{size(filepath,1)+1,1}=path(i).name;
                end
            end
        end
    end
    if size(filepath,1) > 1
        if mod(size(filepath,1),2)==1
            cutoff=size(filepath,1);
            filepath(cutoff)=[];
        end
        filename=cell(1);
        for i=1:size(filepath,1)
            if ispc==1
                zeichen=strfind(filepath{i,1},'\');
            else
                zeichen=strfind(filepath{i,1},'/');
            end
            currentpath=filepath{i,1};
            if mod(i,2) == 1
                filename{i,1}=['A: ' currentpath(zeichen(1,size(zeichen,2))+1:end)];
            else
                filename{i,1}=['B: ' currentpath(zeichen(1,size(zeichen,2))+1:end)];
            end
        end
        %extract path:
        pathname=currentpath(1:zeichen(1,size(zeichen,2))-1);
        put('pathname',pathname); %last path
        put ('filename',filename); %only for displaying
        put ('filepath',filepath); %full path and filename for analyses
        sliderrange
        set (handles.filenamebox, 'string', filename);
        put ('resultslist', []); %clears old results
        put ('derived',[]);
        put('displaywhat',1);%vectors
        put('ismean',[]);
        put('framemanualdeletion',[]);
        put('manualdeletion',[]);
        put('streamlinesX',[]);
        put('streamlinesY',[]);
        set(handles.fileselector, 'value',1);
        %Clear all things
        clear_vel_limit_Callback %clear velocity limits
        clear_roi_Callback
        %clear_mask_Callback:
        delete(findobj(gca,'tag', 'maskplot'));
        put ('maskiererx',{});
        put ('maskierery',{});
        set(handles.mask_hint, 'String', 'Mask inactive', 'backgroundcolor', [0.9411764705882353 0.9411764705882353 0.9411764705882353]);
        
        %reset zoom
        set(handles.panon,'Value',0);
        set(handles.zoomon,'Value',0);
        put('xzoomlimit', []);
        put('yzoomlimit', []);
        
        
        %clear_cali_Callback do not clear calibration anymore.....
        %Problems...?
        sliderdisp %displays raw image when slider moves
        zoom reset
        set(handles.skipper, 'enable', 'on');
        set(handles.applyskipper, 'enable', 'on');
    else
        errordlg('Please select at least two images ( = 1 pair of images)','Error','on')
    end
end

function sliderrange
filepath=retr('filepath');
handles=gethand;
if size(filepath,1)>2
    sliderstepcount=size(filepath,1)/2;
    set(handles.fileselector, 'enable', 'on');
    set (handles.fileselector,'value',1, 'min', 1,'max',sliderstepcount,'sliderstep', [1/(sliderstepcount-1) 1/(sliderstepcount-1)*10]);
else
    sliderstepcount=1;
    set(handles.fileselector, 'enable', 'off');
    set (handles.fileselector,'value',1, 'min', 1,'max',2,'sliderstep', [0.5 0.5]);
end

function fileselector_Callback(hObject, eventdata, handles)

filepath=retr('filepath');
if size(filepath,1) > 1
    sliderdisp
end

function togglepair_Callback(hObject, eventdata, handles)
toggler=get(gco, 'value');
put ('toggler',toggler);
filepath=retr('filepath');
if size(filepath,1) > 1
    sliderdisp
    handles=gethand;
    if strmatch(get(handles.multip03, 'visible'), 'on')
        preview_preprocess_Callback
    end
end

%{
set(findobj('-not','type','uimenu','-not','type', 'image'),'units','points');
get(findobj('-regexp','tag','tools'),'position') %459 3.75 141 90
get(findobj('-regexp','tag','multip01'),'position') % 459 90 141 363
get(gcf,'position')
get(gcf,'units')
%}

function loadimgs_Callback(hObject, eventdata, handles)
switchui('multip01')
delete(findobj('tag','hinting'))
test1=get(gca,'xlim');
test2=get(gca,'ylim');
if test1(2)==603 && test2(2)==580 %%only display hint when logo is shown
text(590,30,'import your image pairs by clicking on ''Load images'' \rightarrow','horizontalalignment','right','verticalalignment','middle','fontsize',14,'tag','hinting')
end

function img_mask_Callback(hObject, eventdata, handles)
switchui('multip02')

function pre_proc_Callback(hObject, eventdata, handles)
switchui('multip03')

function piv_sett_Callback(hObject, eventdata, handles)
switchui('multip04')
pause(0.01) %otherwise display isn't updated... ?!?
drawnow;drawnow;
dispinterrog
handles=gethand;
if get(handles.dcc,'value')==1
    countparticles
    set(handles.recommendation,'visible','on');
else
    set(handles.recommendation,'visible','off');
end
overlappercent

function do_analys_Callback(hObject, eventdata, handles)
handles=gethand;
set(handles.progress, 'String','Frame progress: N/A');
set(handles.overall, 'String','Total progress: N/A');
set(handles.totaltime, 'String','Time left: N/A');
set(handles.messagetext, 'String','');
switchui('multip05')

function vector_val_Callback(hObject, eventdata, handles)
switchui('multip06')

function cal_actual_Callback(hObject, eventdata, handles)
switchui('multip07')
pointscali=retr('pointscali');
if numel(pointscali>0)
    xposition=pointscali(:,1);
    yposition=pointscali(:,2);
    caliimg=retr('caliimg');
    if numel(caliimg)>0
        image(caliimg, 'parent',gca, 'cdatamapping', 'scaled');
        colormap('gray');
        axis image;
        set(gca,'ytick',[])
        set(gca,'xtick',[])
    else
        sliderdisp
    end
    hold on;
    plot (xposition,yposition,'ro-', 'markersize', 10,'LineWidth',3 , 'tag', 'caliline');
    plot (xposition,yposition,'y+:', 'tag', 'caliline');
    hold off;
    for j=1:2
        text(xposition(j)+10,yposition(j)+10, ['x:' num2str(round(xposition(j)*10)/10) sprintf('\n') 'y:' num2str(round(yposition(j)*10)/10) ],'color','y','fontsize',7, 'BackgroundColor', 'k', 'tag', 'caliline')
    end
    text(mean(xposition),mean(yposition), ['s = ' num2str(round((sqrt((xposition(1)-xposition(2))^2+(yposition(1)-yposition(2))^2))*100)/100) ' px'],'color','k','fontsize',7, 'BackgroundColor', 'r', 'tag', 'caliline','horizontalalignment','center')
end

function plot_derivs_Callback(hObject, eventdata, handles)
handles=gethand;
switchui('multip08');
if retr('caluv')==1 && retr('calxy')==1
    set(handles.derivchoice,'String',{'Vectors [px/frame]';'Vorticity [1/frame]';'Velocity magnitude [px/frame]';'u component [px/frame]';'v component [px/frame]';'Divergence [1/frame]';'Vortex locator [1]';'Simple shear rate [1/frame]';'Simple strain rate [1/frame]';'Line integral convolution (LIC) [1]'});
    set(handles.text35,'String','u [px/frame]:')
    set(handles.text36,'String','v [px/frame]:')
else
    set(handles.derivchoice,'String',{'Vectors [m/s]';'Vorticity [1/s]';'Velocity magnitude [m/s]';'u component [m/s]';'v component [m/s]';'Divergence [1/s]';'Vortex locator [1]';'Simple shear rate [1/s]';'Simple strain rate [1/s]';'Line integral convolution (LIC) [1]'});
    set(handles.text35,'String','u [m/s]:')
    set(handles.text36,'String','v [m/s]:')
end
derivchoice_Callback(handles.derivchoice)


function modif_plot_Callback(hObject, eventdata, handles)
switchui('multip09');

function ascii_chart_Callback(hObject, eventdata, handles)
switchui('multip10')

function matlab_file_Callback(hObject, eventdata, handles)
switchui('multip11')

function poly_extract_Callback(hObject, eventdata, handles)
handles=gethand;
switchui('multip12')
if retr('caluv')==1 && retr('calxy')==1
    set(handles.extraction_choice,'string', {'Vorticity [1/frame]';'Velocity magnitude [px/frame]';'u component [px/frame]';'v component [px/frame]';'Divergence [1/frame]';'Vortex locator [1]';'Shear rate [1/frame]';'Strain rate [1/frame]';'Tangent velocity [px/frame]'});
else
    set(handles.extraction_choice,'string', {'Vorticity [1/s]';'Velocity magnitude [m/s]';'u component [m/s]';'v component [m/s]';'Divergence [1/s]';'Vortex locator [1]';'Shear rate [1/s]';'Strain rate [1/s]';'Tangent velocity [m/s]'});
end

function dist_angle_Callback(hObject, eventdata, handles)
switchui('multip13')

function statistics_Callback(hObject, eventdata, handles)
switchui('multip14')
filepath=retr('filepath');
if size(filepath,1) > 1
    sliderdisp
end

function part_img_sett_Callback(hObject, eventdata, handles)
switchui('multip15')

function save_movie_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
if size(resultslist,2)>=1
    startframe=0;
    endframe=0;
    for i=1:size(resultslist,2)
        if numel(resultslist{1,i})>0 && startframe==0
            startframe=i;
        end
        if numel(resultslist{1,i})>0
            endframe=i;
        end
    end
    set(handles.firstframe, 'String',int2str(startframe));
    set(handles.lastframe, 'String',int2str(endframe));
    if strmatch(get(handles.multip08, 'visible'), 'on');
        put('p8wasvisible',1)
    else
        put('p8wasvisible',0)
    end
    switchui('multip16');
else
    msgbox('No analyses yet...')
end

function area_extract_Callback(hObject, eventdata, handles)
handles=gethand;
switchui('multip17');
if retr('caluv')==1 && retr('calxy')==1
    set(handles.area_para_select,'string', {'Vorticity [1/frame]';'Velocity magnitude [px/frame]';'u component [px/frame]';'v component [px/frame]';'Divergence [1/frame]';'Vortex locator [1]';'Shear rate [1/frame]';'Strain rate [1/frame]'});
else
    set(handles.area_para_select,'string', {'Vorticity [1/s]';'Velocity magnitude [m/s]';'u component [m/s]';'v component [m/s]';'Divergence [1/s]';'Vortex locator [1]';'Shear rate [1/s]';'Strain rate [1/s]'});
end

function scatterplotter_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    if size(resultslist,1)>6 %filtered exists
        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
            u=resultslist{10,currentframe};
            v=resultslist{11,currentframe};
        else
            u=resultslist{7,currentframe};
            if size(u,1)>1
                v=resultslist{8,currentframe};
            else
                %filter was applied to some other frame than this
                %load unfiltered results
                u=resultslist{3,currentframe};
                v=resultslist{4,currentframe};
            end
        end
    else
        u=resultslist{3,currentframe};
        v=resultslist{4,currentframe};
    end
    caluv=retr('caluv');
    u=reshape(u,size(u,1)*size(u,2),1);
    v=reshape(v,size(v,1)*size(v,2),1);
    h=figure;
    screensize=get( 0, 'ScreenSize' );
    %rect = [screensize(3)/2-300, screensize(4)/2-250, 600, 500];
    rect = [screensize(3)/4-300, screensize(4)/2-250, 600, 500];
    set(h,'position', rect);
    set(h,'numbertitle','off','menubar','none','toolbar','figure','dockcontrols','off','name',['Scatter plot u & v, frame ' num2str(currentframe)],'tag', 'derivplotwindow');
    h2=scatter(u*caluv-retr('subtr_u'),v*caluv-retr('subtr_v'),'r.');
    set (gca, 'xgrid', 'on', 'ygrid', 'on', 'TickDir', 'in')
    if retr('caluv')==1 && retr('calxy')==1
        xlabel('u [px/frame]');
        ylabel('v [px/frame]');
    else
        xlabel('u [m/s]');
        ylabel('v [m/s]');
    end
end

function autocrop (file,fmt)
A=imread(file);
B=rgb2gray(A);
for i=1:ceil(size(B,1)/2)
    val(i)=mean(B(i,:));
end
startcropy=max([find(val==255) 1]);
for i=size(B,1):-1:ceil(size(B,1)/2)
    val2(i)=mean(B(i,:));
end

endcropy=min(find(val2==255));
clear val val2
for i=1:ceil(size(B,2)/2)
    val(i)=mean(B(:,i));
end
startcropx=max([find(val==255) 1]);
for i=size(B,2):-1:ceil(size(B,2)/2)
    val2(i)=mean(B(:,i));
end
endcropx=min(find(val2==255));
A=A(startcropy:endcropy,startcropx:endcropx,:);

if fmt==1 %jpg
    imwrite(A,file,'quality', 100);
else
    imwrite(A,file);
end


function file_save (currentframe,FileName,PathName,type)
handles=gethand;
resultslist=retr('resultslist');
derived=retr('derived');
filename=retr('filename');
caluv=retr('caluv');
calxy=retr('calxy');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    x=resultslist{1,currentframe};
    y=resultslist{2,currentframe};
    if size(resultslist,1)>6 %filtered exists
        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
            u=resultslist{10,currentframe};
            v=resultslist{11,currentframe};
            typevector=resultslist{9,currentframe};
            if numel(typevector)==0%happens if user smoothes sth without NaN and without validation
                typevector=resultslist{5,currentframe};
            end
        else
            u=resultslist{7,currentframe};
            if size(u,1)>1
                v=resultslist{8,currentframe};
                typevector=resultslist{9,currentframe};
            else
                %filter was applied to some other frame than this
                %load unfiltered results
                u=resultslist{3,currentframe};
                v=resultslist{4,currentframe};
                typevector=resultslist{5,currentframe};
            end
        end
    else
        u=resultslist{3,currentframe};
        v=resultslist{4,currentframe};
        typevector=resultslist{5,currentframe};
    end
end
u(typevector==0)=NaN;
v(typevector==0)=NaN;
subtract_u=retr('subtr_u');
subtract_v=retr('subtr_v');

if type==1 %ascii file
    delimiter=get(handles.delimiter, 'value');
    if delimiter==1
        delimiter=',';
    elseif delimiter==2
        delimiter='\t';
    elseif delimiter==3
        delimiter=' ';
    end
    if get(handles.addfileinfo, 'value')==1
        header1=['PIVlab by W.Th. & E.J.S., ASCII chart output - ' date];
        header2=['FRAME: ' int2str(currentframe) ', filenames: ' filename{currentframe*2-1} ' & ' filename{currentframe*2} ', conversion factor xy (px -> m): ' num2str(calxy) ', conversion factor uv (px/frame -> m/s): ' num2str(caluv)];
    else
        header1=[];
        header2=[];
    end
    if get(handles.add_header, 'value')==1
        if retr('calxy')==1 && retr('caluv')==1
            if get(handles.export_vort, 'Value') == 1
                header3=['x [px]' delimiter 'y [px]' delimiter 'u [px/frame]' delimiter 'v [px/frame]' delimiter 'vorticity [1/frame]'];%delimiter 'magnitude[m/s]' delimiter 'divergence[1]' delimiter 'vorticity[1/s]' delimiter 'dcev[1]']
            else
                header3=['x [px]' delimiter 'y [px]' delimiter 'u [px/frame]' delimiter 'v [px/frame]'];%delimiter 'magnitude[m/s]' delimiter 'divergence[1]' delimiter 'vorticity[1/s]' delimiter 'dcev[1]']
            end
        else
            if get(handles.export_vort, 'Value') == 1
                header3=['x [m]' delimiter 'y [m]' delimiter 'u [m/s]' delimiter 'v [m/s]' delimiter 'vorticity [1/s]'];%delimiter 'magnitude[m/s]' delimiter 'divergence[1]' delimiter 'vorticity[1/s]' delimiter 'dcev[1]']
            else
                header3=['x [m]' delimiter 'y [m]' delimiter 'u [m/s]' delimiter 'v [m/s]'];%delimiter 'magnitude[m/s]' delimiter 'divergence[1]' delimiter 'vorticity[1/s]' delimiter 'dcev[1]']
            end
        end
    else
        header3=[];
    end
    if isempty(header1)==0
        fid = fopen(fullfile(PathName,FileName), 'w');
        fprintf(fid, [header1 '\r\n']);
        fclose(fid);
    end
    if isempty(header2)==0
        fid = fopen(fullfile(PathName,FileName), 'a');
        fprintf(fid, [header2 '\r\n']);
        fclose(fid);
    end
    if isempty(header3)==0
        fid = fopen(fullfile(PathName,FileName), 'a');
        fprintf(fid, [header3 '\r\n']);
        fclose(fid);
    end
    if get(handles.export_vort, 'Value') == 1
        derivative_calc(currentframe,2,1);
        derived=retr('derived');
        vort=derived{1,currentframe};
        wholeLOT=[reshape(x*calxy,size(x,1)*size(x,2),1) reshape(y*calxy,size(y,1)*size(y,2),1) reshape(u*caluv-subtract_u,size(u,1)*size(u,2),1) reshape(v*caluv-subtract_v,size(v,1)*size(v,2),1) reshape(vort,size(vort,1)*size(vort,2),1)];
    else
        wholeLOT=[reshape(x*calxy,size(x,1)*size(x,2),1) reshape(y*calxy,size(y,1)*size(y,2),1) reshape(u*caluv-subtract_u,size(u,1)*size(u,2),1) reshape(v*caluv-subtract_v,size(v,1)*size(v,2),1)];
    end
    dlmwrite(fullfile(PathName,FileName), wholeLOT, '-append', 'delimiter', delimiter, 'precision', 10, 'newline', 'pc');
end %type==1
if type==2 %matlab file
    u=u*caluv-subtract_u;
    v=v*caluv-subtract_v;
    x=x*calxy;
    y=y*calxy;
    if retr('calxy')==1 && retr('caluv')==1
        info='units are px respectively px/frame';
    else
        info='units are m respectively m/s';
    end
    if get(handles.export_vort2, 'Value') == 1
        derivative_calc(currentframe,2,1);
        derived=retr('derived');
        vort=derived{1,currentframe};
        
        save(fullfile(PathName,FileName), 'u', 'v', 'x', 'y', 'typevector', 'calxy', 'caluv', 'vort' ,'info');
    else
        save(fullfile(PathName,FileName), 'u', 'v', 'x', 'y', 'typevector', 'calxy', 'caluv', 'info');
    end
end
if type==3
    u=u*caluv-subtract_u;
    v=v*caluv-subtract_v;
    x=x*calxy;
    y=y*calxy;


nr_of_elements=numel(x);
fid = fopen(fullfile(PathName,FileName), 'w'); 
if retr('calxy')==1 && retr('caluv')==1
    info='[px/frame]';
else
    info='[m/s]';
end
%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, ['VTK from PIVlab ' info '\n']);
fprintf(fid, 'BINARY\n\n');
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ' num2str(size(x,1)) ' ' num2str(size(x,2)) ' ' num2str(size(x,3)) '\n']);
fprintf(fid, ['POINTS ' num2str(nr_of_elements) ' float\n']);
fclose(fid);

%append binary x,y,z data
fid = fopen(fullfile(PathName,FileName), 'a'); 
fwrite(fid, [reshape(x,1,nr_of_elements);  reshape(y,1,nr_of_elements); reshape(y,1,nr_of_elements)*0],'float','b');

%append another ASCII sub header
fprintf(fid, ['\nPOINT_DATA ' num2str(nr_of_elements) '\n']);
fprintf(fid, 'VECTORS velocity_vectors float\n');

%append binary u,v,w data
fwrite(fid, [reshape(u,1,nr_of_elements);  reshape(v,1,nr_of_elements); reshape(v,1,nr_of_elements)*0],'float','b');

fclose(fid);

end
% --------------------------------------------------------------------
function roi_select_Callback(hObject, eventdata, handles)
filepath=retr('filepath');
handles=gethand;
if size(filepath,1) > 1
    delete(findobj('tag','warning'));
    toolsavailable(0);
    toggler=retr('toggler');
    selected=2*floor(get(handles.fileselector, 'value'))-(1-toggler);
    filepath=retr('filepath');
    roirect = round(getrect(gca));
    if roirect(1,3)~=0 && roirect(1,4)~=0
        imagesize(1)=size(imread(filepath{selected}),1);
        imagesize(2)=size(imread(filepath{selected}),2);
        if roirect(1)<1
            roirect(1)=1;
        end
        if roirect(2)<1
            roirect(2)=1;
        end
        if roirect(3)>imagesize(2)-roirect(1)
            roirect(3)=imagesize(2)-roirect(1);
        end
        if roirect(4)>imagesize(1)-roirect(2)
            roirect(4)=imagesize(1)-roirect(2);
        end
        put ('roirect',roirect);
        dispROI
        
        set(handles.roi_hint, 'String', 'ROI active' , 'backgroundcolor', [0.5 1 0.5]);
    else
        text(50,50,'Invalid selection: Click and hold left mouse button to create a rectangle.','color','r','fontsize',8, 'BackgroundColor', 'k','tag','warning');
    end
    toolsavailable(1);
end

% --- Executes on button press in clear_roi.
function clear_roi_Callback(hObject, eventdata, handles)
handles=gethand;
delete(findobj(gca,'tag', 'roiplot'));
delete(findobj(gca,'tag', 'roitext'));
delete(findobj('tag','warning'));
put ('roirect',[]);
set(handles.roi_hint, 'String', 'ROI inactive', 'backgroundcolor', [0.9411764705882353 0.9411764705882353 0.9411764705882353]);
set(handles.ROI_Man_x,'String','');
set(handles.ROI_Man_y,'String','');
set(handles.ROI_Man_w,'String','');
set(handles.ROI_Man_h,'String','');

function dispROI
handles=gethand;
roirect=retr('roirect');
x=[roirect(1)  roirect(1)+roirect(3) roirect(1)+roirect(3)  roirect(1)            roirect(1) ];
y=[roirect(2)  roirect(2)            roirect(2)+roirect(4)  roirect(2)+roirect(4) roirect(2) ];
delete(findobj(gca,'tag', 'roiplot'));
delete(findobj(gca,'tag', 'roitext'));
rectangle('Position',roirect,'LineWidth',1,'LineStyle','-','edgecolor','b','tag','roiplot')
rectangle('Position',roirect,'LineWidth',1,'LineStyle',':','edgecolor','y','tag','roiplot')
set(handles.ROI_Man_x,'String',int2str(roirect(1)));
set(handles.ROI_Man_y,'String',int2str(roirect(2)));
set(handles.ROI_Man_w,'String',int2str(roirect(3)));
set(handles.ROI_Man_h,'String',int2str(roirect(4)));

function dispMASK(opaqueness)
if opaqueness == 1
    maskcolor = [0.3 0.1 0.1];
else
    maskcolor = [1 0 0];
end
handles=gethand;
currentframe=2*floor(get(handles.fileselector, 'value'))-1;
maskiererx=retr('maskiererx');
maskierery=retr('maskierery');
delete(findobj(gca,'tag', 'maskplot'));
hold on;
for j=1:size(maskiererx,1)
    if isempty(maskiererx{j,currentframe})==0
        ximask=maskiererx{j,currentframe};
        yimask=maskierery{j,currentframe};
        if verLessThan('matlab','8.4')
            h=fill(ximask,yimask,'r','facecolor', [0.3 0.1 0.1],'linestyle','none','tag','maskplot');
        else
            % >R2014a
            h=fill(ximask,yimask,'r','facecolor', maskcolor,'linestyle','none','tag','maskplot','Facealpha',opaqueness);
        end
        %h=area(ximask,yimask,'facecolor', [0.3 0.1 0.1],'linestyle', 'none','tag','maskplot');
    else
        break;
    end
end
hold off;

function draw_mask_Callback(hObject, eventdata, handles)
filepath=retr('filepath');
handles=gethand;
if size(filepath,1) > 1
    toolsavailable(0);
    currentframe=2*floor(get(handles.fileselector, 'value'))-1;
    filepath=retr('filepath');
    amount=size(filepath,1);
    %currentframe and currentframe+1 =is a pair with identical mask.
    %maskiererx&y contains masks. 3rd dimension is frame nr.
    maskiererx=retr('maskiererx');
    maskierery=retr('maskierery');
    [mask,ximask,yimask]=roipoly;
    insertion=1;
    for j=size(maskiererx,1):-1:1
        try
            if isempty(maskiererx{j,currentframe})==0
                insertion=j+1;
                break
            end
        catch
            maskiererx{1,currentframe}=[];
            maskierery{1,currentframe}=[];
            insertion=1;
        end
    end
    maskiererx{insertion,currentframe}=ximask;
    maskiererx{insertion,currentframe+1}=ximask;
    maskierery{insertion,currentframe}=yimask;
    maskierery{insertion,currentframe+1}=yimask;
    put('maskiererx' ,maskiererx);
    put('maskierery' ,maskierery);
    dispMASK(0.333)
    set(handles.mask_hint, 'String', 'Mask active', 'backgroundcolor', [0.5 1 0.5]);
    toolsavailable(1);
end



% --- Executes on button press in external_mask.
function external_mask_Callback(hObject, eventdata, handles)
filepath=retr('filepath');
handles=gethand;
if size(filepath,1) > 1
    h=helpdlg(['Mask MAT files have to be in the following format' sprintf('\n') 'MAT file must contain two variables: xmask and ymask' sprintf('\n') 'Each column contains a number of polygon coordinates' sprintf('\n') 'Column 1 is for frame 1, column 2 for frame 2 etc.' sprintf('\n') 'You could try this code in Matlab to understand how polygons are created:' sprintf('\n') '[junk,xmask,ymask] = roipoly (rand(100,100))'],'External masks');
    uiwait(h)
    [FileName,PathName] = uigetfile('*.mat','Select the mask MAT-file');
    if isequal(FileName,0) | isequal(PathName,0)
    else
        
        %oad('maskemaske.mat','maskiererx','maskierery','xmask','ymask')
        %dann werden nur diese vars geladen. wenns eine davon nicht gibt wird sie nicht geladen.
        %im anschluss abfragen ob maskiererx existiert oder xmask.
        %danach entscheiden wie masken geladen werden.
        load (fullfile(PathName,FileName));
        for i=1:size(xmask,2)
            %1   2  3  4
            %12 34 56 78
            maskiererx{1,i*2-1}=xmask(:,i);
            maskiererx{1,i*2}=xmask(:,i);
            maskierery{1,i*2-1}=ymask(:,i);
            maskierery{1,i*2}=ymask(:,i);
        end
        put('maskiererx' ,maskiererx);
        put('maskierery' ,maskierery);
        set(handles.mask_hint, 'String', 'Mask active', 'backgroundcolor', [0.5 1 0.5]);
        dispMASK(0.333)
    end
end


% --- Executes on button press in clear_mask.
function clear_mask_Callback(hObject, eventdata, handles)
button = questdlg('Do you want to remove all masks?','Delete?','Yes','Cancel','Cancel');
if strmatch(button,'Yes')==1
    handles=gethand;
    delete(findobj(gca,'tag', 'maskplot'));
    put ('maskiererx',{});
    put ('maskierery',{});
    set(handles.mask_hint, 'String', 'Mask inactive', 'backgroundcolor', [0.9411764705882353 0.9411764705882353 0.9411764705882353]);
end
% --- Executes on button press in clear_current_mask.
function clear_current_mask_Callback(hObject, eventdata, handles)
filepath=retr('filepath');
handles=gethand;
if size(filepath,1) > 1
    delete(findobj(gca,'tag', 'maskplot'));
    currentframe=2*floor(get(handles.fileselector, 'value'))-1;
    maskiererx=retr('maskiererx');
    maskierery=retr('maskierery');
    for i=1:size(maskiererx,1)
        maskiererx{i,currentframe}=[];
        maskiererx{i,currentframe+1}=[];
        maskierery{i,currentframe}=[];
        maskierery{i,currentframe+1}=[];
    end
    try
        emptycells=cellfun('isempty',maskiererx);
    catch
        disp('Problems with old Matlab version... Please update Matlab or unexpected things might happen...')
    end
    if mean(double(emptycells))==1 %not very sophisticated way to determine if all cells are empty
        set(handles.mask_hint, 'String', 'Mask inactive', 'backgroundcolor', [0.9411764705882353 0.9411764705882353 0.9411764705882353]);
    end
    put('maskiererx' ,maskiererx);
    put('maskierery' ,maskierery);
end

% --- Executes on button press in mask_applytoall.


% --- Executes on button press in maskToSelected.
function maskToSelected_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
if size(filepath,1) > 1
    str = strrep(get(handles.maskapplyselect,'string'),'-',':');
    endinside=findstr(str, 'end');
    if isempty(endinside)==0
        str = strrep(str,'end',num2str(size(filepath,1)/2));
        
    end
    selectionok=1;
    strnum=str2num(str);
    if isempty(strnum)==1 || isempty(findstr(str,'.'))==0 || isempty(findstr(str,';'))==0
        msgbox(['Error in frame selection syntax. Please use the following syntax (examples):' sprintf('\n') '1:3' sprintf('\n') '1,3,7,9' sprintf('\n') '1:3,7,8,9,11:13' ],'Error','error','modal')
        selectionok=0;
    end
    amount=max(strnum);
    if amount*2>size(filepath,1)
        msgbox(['Selected frames out of range.'],'Error','error','modal')
        selectionok=0;
    end
    
    
    mini=min(strnum);
    
    %checken ob nicht größer als geladene frame anzahl.
    
    if selectionok==1
        currentframe=2*floor(get(handles.fileselector, 'value'))-1;
        %amount=size(filepath,1);
        
        maskiererx=retr('maskiererx');
        maskierery=retr('maskierery');
        
        for i=1:size(strnum,2)
            for j=1:size(maskiererx,1)
                %keyboard
                %in frame 1=maskiererx 1und2
                maskiererx{j,strnum(i)*2-1}=maskiererx{j,currentframe};
                maskiererx{j,strnum(i)*2}=maskiererx{j,currentframe+1};
                maskierery{j,strnum(i)*2-1}=maskierery{j,currentframe};
                maskierery{j,strnum(i)*2}=maskierery{j,currentframe+1};
                
            end
        end
        put('maskiererx' ,maskiererx);
        put('maskierery' ,maskierery);
    end
end




function preview_preprocess_Callback(hObject, eventdata, handles)
filepath=retr('filepath');
if size(filepath,1) >1
    handles=gethand;
    toggler=retr('toggler');
    filepath=retr('filepath');
    selected=2*floor(get(handles.fileselector, 'value'))-(1-toggler);
    img=imread(filepath{selected});
    clahe=get(handles.clahe_enable,'value');
    highp=get(handles.enable_highpass,'value');
    %clip=get(handles.enable_clip,'value');
    intenscap=get(handles.enable_intenscap, 'value');
    clahesize=str2double(get(handles.clahe_size, 'string'));
    highpsize=str2double(get(handles.highp_size, 'string'));
    wienerwurst=get(handles.wienerwurst, 'value');
    wienerwurstsize=str2double(get(handles.wienerwurstsize, 'string'));
    %clipthresh=str2double(get(handles.clip_thresh, 'string'));
    roirect=retr('roirect');
    if size (roirect,2)<4
        roirect=[1,1,size(img,2)-1,size(img,1)-1];
    end
    out = PIVlab_preproc (img,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize);
    image(out, 'parent',gca, 'cdatamapping', 'scaled');
    colormap('gray');
    axis image;
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    roirect=retr('roirect');
    if size(roirect,2)>1
        dispROI
    end
    currentframe=2*floor(get(handles.fileselector, 'value'))-1;
    maskiererx=retr('maskiererx');
    if size(maskiererx,2)>=currentframe
        ximask=maskiererx{currentframe};
        if size(ximask,1)>1
            dispMASK(0.333)
        end
    end
end

function dispinterrog
handles=gethand;
selected=2*floor(get(handles.fileselector, 'value'))-1;
filepath=retr('filepath');
if numel(filepath)>1
    size_img(1)=size(imread(filepath{selected}),2)/2;
    size_img(2)=size(imread(filepath{selected}),1)/2;
    step=str2double(get(handles.step,'string'));
    delete(findobj(gca,'Type','hggroup')); %=vectors and scatter markers
    delete(findobj(gca,'tag','intareadispl'));
    centre(1)= size(imread(filepath{selected}),1)/2; %y
    centre(2)= size(imread(filepath{selected}),2)/2; %x
    
    intarea1=str2double(get(handles.intarea,'string'))/2;
    x1=[centre(2)-intarea1 centre(2)+intarea1 centre(2)+intarea1 centre(2)-intarea1 centre(2)-intarea1];
    y1=[centre(1)-intarea1 centre(1)-intarea1 centre(1)+intarea1 centre(1)+intarea1 centre(1)-intarea1];
    hold on;
    plot(x1,y1,'c-', 'linewidth', 1, 'linestyle', ':','tag','intareadispl');
    if get(handles.fftmulti,'value')==1
        text(x1(1),y1(1), ['pass 1'],'color','c','fontsize',8,'tag','intareadispl','HorizontalAlignment','right','verticalalignment','bottom')
        if get(handles.checkbox26,'value')==1
            intarea2=str2double(get(handles.edit50,'string'))/2;
            x2=[centre(2)-intarea2 centre(2)+intarea2 centre(2)+intarea2 centre(2)-intarea2 centre(2)-intarea2];
            y2=[centre(1)-intarea2 centre(1)-intarea2 centre(1)+intarea2 centre(1)+intarea2 centre(1)-intarea2];
            plot(x2,y2,'y-', 'linewidth', 1, 'linestyle', ':','tag','intareadispl');
            text(x2(2),y2(1), ['pass 2'],'color','y','fontsize',8,'tag','intareadispl','HorizontalAlignment','left','verticalalignment','bottom')
        end
        if get(handles.checkbox27,'value')==1
            intarea3=str2double(get(handles.edit51,'string'))/2;
            x3=[centre(2)-intarea3 centre(2)+intarea3 centre(2)+intarea3 centre(2)-intarea3 centre(2)-intarea3];
            y3=[centre(1)-intarea3 centre(1)-intarea3 centre(1)+intarea3 centre(1)+intarea3 centre(1)-intarea3];
            plot(x3,y3,'g-', 'linewidth', 1, 'linestyle', ':','tag','intareadispl');
            text(x3(2),y3(3), ['pass 3'],'color','g','fontsize',8,'tag','intareadispl','HorizontalAlignment','left','verticalalignment','top')
        end
        if get(handles.checkbox28,'value')==1
            intarea4=str2double(get(handles.edit52,'string'))/2;
            x4=[centre(2)-intarea4 centre(2)+intarea4 centre(2)+intarea4 centre(2)-intarea4 centre(2)-intarea4];
            y4=[centre(1)-intarea4 centre(1)-intarea4 centre(1)+intarea4 centre(1)+intarea4 centre(1)-intarea4];
            plot(x4,y4,'r-', 'linewidth', 1, 'linestyle', ':','tag','intareadispl');
            text(x4(1),y4(3), ['pass 4'],'color','r','fontsize',8,'tag','intareadispl','HorizontalAlignment','right','verticalalignment','top')
        end
    end
    hold off;
    text(x1(1),y1(1)-20, ['Interrogation area(s) example'],'color','y','fontsize',9,'tag','intareadispl')
end

function countparticles
%{
handles=gethand;
selected=2*floor(get(handles.fileselector, 'value'))-1;
filepath=retr('filepath');
if numel(filepath)>1
    A=imread(filepath{selected});
    A=A(round(0.25*size(A,1)):round(0.75*size(A,1)),round(0.25*size(A,2)):round(0.75*size(A,2)));
    clahe=get(handles.clahe_enable,'value');
    highp=get(handles.enable_highpass,'value');
    %clip=get(handles.enable_clip,'value');
    intenscap=get(handles.enable_intenscap, 'value');
    clahesize=str2double(get(handles.clahe_size, 'string'))*2; % faster...
    highpsize=str2double(get(handles.highp_size, 'string'));
        wienerwurst=get(handles.wienerwurst, 'value');
        wienerwurstsize=str2double(get(handles.wienerwurstsize, 'string'));
    %clipthresh=str2double(get(handles.clip_thresh, 'string'));
    roirect=retr('roirect');
    A = PIVlab_preproc (A,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize);
    A(A<=80)=0;
    A(A>80)=255;
    [spots,numA]=bwlabeln(A,8);
    XA=numA/(size(A,1)*size(A,2));
    YA=8/XA;
    Y1A=16/XA;
    recommendedMIN=round(sqrt(YA)); % 8 peaks are in Z*Z area
    recommendedMAX=round(sqrt(Y1A));
    set (handles.recommendation, 'String', ['Minimal int area size: ' int2str(recommendedMIN) 'px to ' int2str(recommendedMAX) 'px']);
end
%}

function intarea_Callback(hObject, eventdata, handles)
overlappercent
dispinterrog

function step_Callback(hObject, eventdata, handles)
overlappercent

function AnalyzeAll_Callback(hObject, eventdata, handles)
ok=checksettings;
if ok==1
    handles=gethand;
    filepath=retr('filepath');
    filename=retr('filename');
    resultslist=cell(0); %clear old results
    toolsavailable(0);
    set (handles.cancelbutt, 'enable', 'on');
    ismean=retr('ismean');
    maskiererx=retr('maskiererx');
    maskierery=retr('maskierery');
    for i=size(ismean,1):-1:1 %remove averaged results
        if ismean(i,1)==1
            filepath(i*2,:)=[];
            filename(i*2,:)=[];
            
            filepath(i*2-1,:)=[];
            filename(i*2-1,:)=[];
            if size(maskiererx,2)>=i*2
                maskiererx(:,i*2)=[];
                maskierery(:,i*2)=[];
                maskiererx(:,i*2-1)=[];
                maskierery(:,i*2-1)=[];
            end
        end
    end
    put('filepath',filepath);
    put('filename',filename);
    put('ismean',[]);
    sliderrange
    
  
  %for ensemble correltaion experiments
  %{
    setappdata(0,'cormap1',[]);
    setappdata(0,'cormap2',[]);
    setappdata(0,'cormap3',[]);
    setappdata(0,'cormap4',[]);
    %}

    
    for i=1:2:size(filepath,1)
        if i==1
            tic
        end
        cancel=retr('cancel');
        if isempty(cancel)==1 || cancel ~=1
            image1=imread(filepath{i});
            image2=imread(filepath{i+1});
            if size(image1,3)>1
                image1=uint8(mean(image1,3));
                image2=uint8(mean(image2,3));
                disp('Warning: To optimize speed, your images should be grayscale, 8 bit!')
            end
            set(handles.progress, 'string' , ['Frame progress: 0%']);drawnow; %#ok<*NBRAK>
            clahe=get(handles.clahe_enable,'value');
            highp=get(handles.enable_highpass,'value');
            %clip=get(handles.enable_clip,'value');
            intenscap=get(handles.enable_intenscap, 'value');
            clahesize=str2double(get(handles.clahe_size, 'string'));
            highpsize=str2double(get(handles.highp_size, 'string'));
            wienerwurst=get(handles.wienerwurst, 'value');
            wienerwurstsize=str2double(get(handles.wienerwurstsize, 'string'));            
            %clipthresh=str2double(get(handles.clip_thresh, 'string'));
            roirect=retr('roirect');
            image1 = PIVlab_preproc (image1,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize);
            image2 = PIVlab_preproc (image2,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize);
            maskiererx=retr('maskiererx');
            maskierery=retr('maskierery');
            ximask={};
            yimask={};
            if size(maskiererx,2)>=i
                for j=1:size(maskiererx,1);
                    if isempty(maskiererx{j,i})==0
                        ximask{j,1}=maskiererx{j,i}; %#ok<*AGROW>
                        yimask{j,1}=maskierery{j,i};
                    else
                        break
                    end
                end
                if size(ximask,1)>0
                    mask=[ximask yimask];
                else
                    mask=[];
                end
            else
                mask=[];
            end
            interrogationarea=str2double(get(handles.intarea, 'string'));
            step=str2double(get(handles.step, 'string'));
            subpixfinder=get(handles.subpix,'value');
            if get(handles.dcc,'Value')==1
                [x y u v typevector] = piv_DCC (image1,image2,interrogationarea, step, subpixfinder, mask, roirect);
            elseif get(handles.fftmulti,'Value')==1
                passes=1;
                if get(handles.checkbox26,'value')==1
                    passes=2;
                end
                if get(handles.checkbox27,'value')==1
                    passes=3;
                end
                if get(handles.checkbox28,'value')==1
                    passes=4;
                end
                int2=str2num(get(handles.edit50,'string'));
                int3=str2num(get(handles.edit51,'string'));
                int4=str2num(get(handles.edit52,'string'));
                contents = get(handles.popupmenu16,'string');
                imdeform=contents{get(handles.popupmenu16,'Value')};

                [x y u v typevector] = piv_FFTmulti (image1,image2,interrogationarea, step, subpixfinder, mask, roirect,passes,int2,int3,int4,imdeform);
                %u=real(u)
                %v=real(v)
            end
            resultslist{1,(i+1)/2}=x;
            resultslist{2,(i+1)/2}=y;
            resultslist{3,(i+1)/2}=u;
            resultslist{4,(i+1)/2}=v;
            resultslist{5,(i+1)/2}=typevector;
            resultslist{6,(i+1)/2}=[];
            put('resultslist',resultslist);
            set(handles.fileselector, 'value', (i+1)/2);
            set(handles.progress, 'string' , ['Frame progress: 100%'])
            set(handles.overall, 'string' , ['Total progress: ' int2str((i+1)/2/(size(filepath,1)/2)*100) '%'])
            put('subtr_u', 0);
            put('subtr_v', 0);
            sliderdisp
            xpos=size(image1,2)/2-40;
            text(xpos,50, ['Analyzing... ' int2str((i+1)/2/(size(filepath,1)/2)*100) '%' ],'color', 'r','FontName','FixedWidth','fontweight', 'bold', 'fontsize', 20, 'tag', 'annoyingthing')
            zeit=toc;
            done=(i+1)/2;
            tocome=(size(filepath,1)/2)-done;
            zeit=zeit/done*tocome;
            hrs=zeit/60^2;
            mins=(hrs-floor(hrs))*60;
            secs=(mins-floor(mins))*60;
            hrs=floor(hrs);
            mins=floor(mins);
            secs=floor(secs);
            set(handles.totaltime,'string', ['Time left: ' sprintf('%2.2d', hrs) 'h ' sprintf('%2.2d', mins) 'm ' sprintf('%2.2d', secs) 's']);
        end %cancel==0
    end
    delete(findobj('tag', 'annoyingthing'));
    set(handles.overall, 'string' , ['Total progress: ' int2str(100) '%'])
    set(handles.totaltime, 'String','Time left: N/A');
    put('cancel',0);
    %{
    try
    if exist(fullfile(pwd,'PIVlab_temp.mat'),'file') ==2
    delete(fullfile(pwd,'PIVlab_temp.mat'))
    end
    savesessionfuntion (pwd,'PIVlab_temp.mat')
    set(handles.messagetext, 'String',['-> temporary results saved as' sprintf('\n') '   ''PIVlab_temp.mat''.']);
    catch
    set(handles.messagetext, 'String','-> could not save temporary results.');
    end
    %}
end
toolsavailable(1);

function AnalyzeSingle_Callback(hObject, eventdata, handles)
handles=gethand;
ok=checksettings;
if ok==1
    resultslist=retr('resultslist');
    set(handles.progress, 'string' , ['Frame progress: 0%']);drawnow;
    handles=gethand;
    filepath=retr('filepath');
    selected=2*floor(get(handles.fileselector, 'value'))-1;
    ismean=retr('ismean');
    if size(ismean,1)>=(selected+1)/2
        if ismean((selected+1)/2,1) ==1
            currentwasmean=1;
        else
            currentwasmean=0;
        end
    else
        currentwasmean=0;
    end
    if currentwasmean==0
        tic;
        image1=imread(filepath{selected});
        image2=imread(filepath{selected+1});
        if size(image1,3)>1
            image1=uint8(mean(image1,3));
            image2=uint8(mean(image2,3));
            disp('Warning: To optimize speed, your images should be grayscale, 8 bit!')
        end
        clahe=get(handles.clahe_enable,'value');
        highp=get(handles.enable_highpass,'value');
        %clip=get(handles.enable_clip,'value');
        intenscap=get(handles.enable_intenscap, 'value');
        clahesize=str2double(get(handles.clahe_size, 'string'));
        highpsize=str2double(get(handles.highp_size, 'string'));
        wienerwurst=get(handles.wienerwurst, 'value');
        wienerwurstsize=str2double(get(handles.wienerwurstsize, 'string'));
        %clipthresh=str2double(get(handles.clip_thresh, 'string'));
        roirect=retr('roirect');
        image1 = PIVlab_preproc (image1,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize);
        image2 = PIVlab_preproc (image2,roirect,clahe, clahesize,highp,highpsize,intenscap,wienerwurst,wienerwurstsize);
        maskiererx=retr('maskiererx');
        maskierery=retr('maskierery');
        ximask={};
        yimask={};
        if size(maskiererx,2)>=selected
            for i=1:size(maskiererx,1);
                if isempty(maskiererx{i,selected})==0
                    ximask{i,1}=maskiererx{i,selected};
                    yimask{i,1}=maskierery{i,selected};
                else
                    break
                end
            end
            if size(ximask,1)>0
                mask=[ximask yimask];
            else
                mask=[];
            end
        else
            mask=[];
        end
        interrogationarea=str2double(get(handles.intarea, 'string'));
        step=str2double(get(handles.step, 'string'));
        subpixfinder=get(handles.subpix,'value');
        if get(handles.dcc,'Value')==1
            [x y u v typevector] = piv_DCC (image1,image2,interrogationarea, step, subpixfinder, mask, roirect);
        elseif get(handles.fftmulti,'Value')==1
            passes=1;
            if get(handles.checkbox26,'value')==1
                passes=2;
            end
            if get(handles.checkbox27,'value')==1
                passes=3;
            end
            if get(handles.checkbox28,'value')==1
                passes=4;
            end
            int2=str2num(get(handles.edit50,'string'));
            int3=str2num(get(handles.edit51,'string'));
            int4=str2num(get(handles.edit52,'string'));
            contents = get(handles.popupmenu16,'string');
            imdeform=contents{get(handles.popupmenu16,'Value')};
            [x y u v typevector] = piv_FFTmulti (image1,image2,interrogationarea, step, subpixfinder, mask, roirect,passes,int2,int3,int4,imdeform);
            %u=real(u)
            %v=real(v)
        end
        resultslist{1,(selected+1)/2}=x;
        resultslist{2,(selected+1)/2}=y;
        resultslist{3,(selected+1)/2}=u;
        resultslist{4,(selected+1)/2}=v;
        resultslist{5,(selected+1)/2}=typevector;
        resultslist{6,(selected+1)/2}=[];
        %clear previous interpolation results
        resultslist{7, (selected+1)/2} = [];
        resultslist{8, (selected+1)/2} = [];
        resultslist{9, (selected+1)/2} = [];
        resultslist{10, (selected+1)/2} = [];
        resultslist{11, (selected+1)/2} = [];
        put('derived', [])
        put('resultslist',resultslist);
        set(handles.progress, 'string' , ['Frame progress: 100%'])
        set(handles.overall, 'string' , ['Total progress: 100%'])
        time1frame=toc;
        set(handles.totaltime, 'String',['Analysis time: ' num2str(round(time1frame*100)/100) ' s']);
        set(handles.messagetext, 'String','');
        put('subtr_u', 0);
        put('subtr_v', 0);
        sliderdisp
    end
    
end

function ok=checksettings
handles=gethand;
mess={};
filepath=retr('filepath');
if size(filepath,1) <2
    mess{size(mess,2)+1}='No images were loaded';
end
if get(handles.clahe_enable, 'value')==1
    if isnan(str2double(get(handles.clahe_size, 'string')))==1
        mess{size(mess,2)+1}='CLAHE window size contains NaN';
    end
end
if get(handles.enable_highpass, 'value')==1
    if isnan(str2double(get(handles.highp_size, 'string')))==1
        mess{size(mess,2)+1}='Highpass filter size contains NaN';
    end
end
if get(handles.wienerwurst, 'value')==1
    if isnan(str2double(get(handles.wienerwurstsize, 'string')))==1
        mess{size(mess,2)+1}='Wiener2 filter size contains NaN';
    end
end
%if get(handles.enable_clip, 'value')==1
%    if isnan(str2double(get(handles.clip_thresh, 'string')))==1
%        mess{size(mess,2)+1}='Clipping threshold contains NaN';
%    end
%end
if isnan(str2double(get(handles.intarea, 'string')))==1
    mess{size(mess,2)+1}='Interrogation area size contains NaN';
end
if isnan(str2double(get(handles.step, 'string')))==1
    mess{size(mess,2)+1}='Step size contains NaN';
end
if size(mess,2)>0 %error somewhere
    msgbox(['Errors found:' mess],'Errors detected.','warn','modal')
    ok=0;
else
    ok=1;
end

function cancelbutt_Callback(hObject, eventdata, handles)
put('cancel',1);
drawnow;
toolsavailable(1);

function load_settings_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*.mat','Load PIVlab settings','PIVlab_settings.mat');
if isequal(FileName,0)==0
    read_settings (FileName,PathName)
end

function read_settings (FileName,PathName)
handles=gethand;
load(fullfile(PathName,FileName));
set(handles.clahe_enable,'value',clahe_enable);
set(handles.clahe_size,'string',clahe_size);
set(handles.enable_highpass,'value',enable_highpass);
set(handles.highp_size,'string',highp_size);
set(handles.wienerwurst,'value',wienerwurst);
set(handles.wienerwurstsize,'string',wienerwurstsize);
%set(handles.enable_clip,'value',enable_clip);
%set(handles.clip_thresh,'string',clip_thresh);
set(handles.enable_intenscap,'value',enable_intenscap);
set(handles.intarea,'string',intarea);
set(handles.step,'string',stepsize);
set(handles.subpix,'value',subpix);  %popup
set(handles.stdev_check,'value',stdev_check);
set(handles.stdev_thresh,'string',stdev_thresh);
set(handles.loc_median,'value',loc_median);
set(handles.loc_med_thresh,'string',loc_med_thresh);
set(handles.epsilon,'string',epsilon);
set(handles.interpol_missing,'value',interpol_missing);
set(handles.vectorscale,'string',vectorscale);
set(handles.colormap_choice,'value',colormap_choice); %popup
set(handles.addfileinfo,'value',addfileinfo);
set(handles.add_header,'value',add_header);
set(handles.delimiter,'value',delimiter);%popup
set(handles.img_not_mask,'value',img_not_mask);
set(handles.autoscale_vec,'value',autoscale_vec);

set(handles.popupmenu16, 'value',imginterpol);
set(handles.dcc, 'value',dccmark);
set(handles.fftmulti, 'value',fftmark);
if fftmark==1
    set (handles.uipanel42,'visible','on')
else
    set (handles.uipanel42,'visible','off')
end
set(handles.checkbox26, 'value',pass2);
set(handles.checkbox27, 'value',pass3);
set(handles.checkbox28, 'value',pass4);
set(handles.edit50, 'string',pass2val);
set(handles.edit51, 'string',pass3val);
set(handles.edit52, 'string',pass4val);
set(handles.text126, 'string',step2);
set(handles.text127, 'string',step3);
set(handles.text128, 'string',step4);
set(handles.holdstream, 'value',holdstream);
set(handles.streamlamount, 'string',streamlamount);
set(handles.streamlcolor, 'value',streamlcolor);
set(handles.streamlwidth, 'value',streamlcolor);

set(handles.realdist, 'string',realdist);
set(handles.time_inp, 'string',time_inp);

set(handles.nthvect, 'string',nthvect);
set(handles.validr,'string',validr);
set(handles.validg,'string',validg);
set(handles.validb,'string',validb);
set(handles.validdr,'string',validdr);
set(handles.validdg,'string',validdg);
set(handles.validdb,'string',validdb);
set(handles.interpr,'string',interpr);
set(handles.interpg,'string',interpg);
set(handles.interpb,'string',interpb);

if caluv~=1 || calxy ~=1
    set(handles.calidisp, 'string', ['1 px = ' num2str(round(calxy*100000)/100000) ' m' sprintf('\n') '1 px/frame = ' num2str(round(caluv*100000)/100000) ' m/s'],  'backgroundcolor', [0.5 1 0.5]);
end

put('calxy',calxy);
put('caluv',caluv);

function curr_settings_Callback(hObject, eventdata, handles)
handles=gethand;
clahe_enable=get(handles.clahe_enable,'value');
clahe_size=get(handles.clahe_size,'string');
enable_highpass=get(handles.enable_highpass,'value');
highp_size=get(handles.highp_size,'string');
wienerwurst=get(handles.wienerwurst,'value');
wienerwurstsize=get(handles.wienerwurstsize,'string');

%enable_clip=get(handles.enable_clip,'value');
%clip_thresh=get(handles.clip_thresh,'string');
enable_intenscap=get(handles.enable_intenscap,'value');
intarea=get(handles.intarea,'string');
stepsize=get(handles.step,'string');
subpix=get(handles.subpix,'value');  %popup
stdev_check=get(handles.stdev_check,'value');
stdev_thresh=get(handles.stdev_thresh,'string');
loc_median=get(handles.loc_median,'value');
loc_med_thresh=get(handles.loc_med_thresh,'string');
epsilon=get(handles.epsilon,'string');
interpol_missing=get(handles.interpol_missing,'value');
vectorscale=get(handles.vectorscale,'string');
colormap_choice=get(handles.colormap_choice,'value'); %popup
addfileinfo=get(handles.addfileinfo,'value');
add_header=get(handles.add_header,'value');
delimiter=get(handles.delimiter,'value');%popup
img_not_mask=get(handles.img_not_mask,'value');
autoscale_vec=get(handles.autoscale_vec,'value');

imginterpol=get(handles.popupmenu16, 'value');
dccmark=get(handles.dcc, 'value');
fftmark=get(handles.fftmulti, 'value');
pass2=get(handles.checkbox26, 'value');
pass3=get(handles.checkbox27, 'value');
pass4=get(handles.checkbox28, 'value');
pass2val=get(handles.edit50, 'string');
pass3val=get(handles.edit51, 'string');
pass4val=get(handles.edit52, 'string');
step2=get(handles.text126, 'string');
step3=get(handles.text127, 'string');
step4=get(handles.text128, 'string');
holdstream=get(handles.holdstream, 'value');
streamlamount=get(handles.streamlamount, 'string');
streamlcolor=get(handles.streamlcolor, 'value');
streamlcolor=get(handles.streamlwidth, 'value');
realdist=get(handles.realdist, 'string');
time_inp=get(handles.time_inp, 'string');

nthvect=get(handles.nthvect, 'string');
validr=get(handles.validr,'string');
validg=get(handles.validg,'string');
validb=get(handles.validb,'string');
validdr=get(handles.validdr,'string');
validdg=get(handles.validdg,'string');
validdb=get(handles.validdb,'string');
interpr=get(handles.interpr,'string');
interpg=get(handles.interpg,'string');
interpb=get(handles.interpb,'string');

calxy=retr('calxy');
caluv=retr('caluv');

if ispc==1
    [FileName,PathName] = uiputfile('*.mat','Save current settings as...',['PIVlab_set_' getenv('USERNAME') '.mat']);
else
    try
        [FileName,PathName] = uiputfile('*.mat','Save current settings as...',['PIVlab_set_' getenv('USER') '.mat']);
    catch
        [FileName,PathName] = uiputfile('*.mat','Save current settings as...','PIVlab_set.mat');
    end
end
clear handles hObject eventdata
if isequal(FileName,0)==0
    save('-v6', fullfile(PathName,FileName))
end

function vel_limit_Callback(hObject, eventdata, handles)
toolsavailable(0)
%if analys existing
resultslist=retr('resultslist');
handles=gethand;
currentframe=2*floor(get(handles.fileselector, 'value'))-1;
if size(resultslist,2)>=(currentframe+1)/2 %data for current frame exists
    x=resultslist{1,(currentframe+1)/2};
    if size(x,1)>1
        if get(handles.meanofall,'value')==1 %calculating mean doesn't mae sense...
            index=1;
            foundfirst=0;
            for i = 1:size(resultslist,2)
                x=resultslist{1,i};
                if isempty(x)==0 && foundfirst==0
                    firstsizex=size(x,1);
                    secondsizex=size(x,2);
                    foundfirst==1;
                end
                if size(x,1)>1 && size(x,1)==firstsizex && size(x,2) == secondsizex
                    u(:,:,index)=resultslist{3,i};
                    v(:,:,index)=resultslist{4,i};
                    index=index+1;
                end
            end
        else
            y=resultslist{2,(currentframe+1)/2};
            u=resultslist{3,(currentframe+1)/2};
            v=resultslist{4,(currentframe+1)/2};
            typevector=resultslist{5,(currentframe+1)/2};
        end
        velrect=retr('velrect');
        caluv=retr('caluv');
        if numel(velrect>0)
            %user already selected window before...
            %"filter u+v" and display scatterplot
            %problem: if user selects limits and then wants to refine vel
            %limits, all data is filterd out...
            umin=velrect(1);
            umax=velrect(3)+umin;
            vmin=velrect(2);
            vmax=velrect(4)+vmin;
            %check if all results are nan...
            u_backup=u;
            v_backup=v;
            u(u*caluv<umin)=NaN;
            u(u*caluv>umax)=NaN;
            v(u*caluv<umin)=NaN;
            v(u*caluv>umax)=NaN;
            v(v*caluv<vmin)=NaN;
            v(v*caluv>vmax)=NaN;
            u(v*caluv<vmin)=NaN;
            u(v*caluv>vmax)=NaN;
            if mean(mean(mean((isnan(u)))))>0.9 || mean(mean(mean((isnan(v)))))>0.9
                disp('User calibrated after selecting velocity limits. Discarding limits.')
                u=u_backup;
                v=v_backup;
            end
        end
        
        %problem: wenn nur ein frame analysiert, dann gibts probleme wenn display all frames in scatterplot an.
        datau=reshape(u*caluv,1,size(u,1)*size(u,2)*size(u,3));
        datav=reshape(v*caluv,1,size(v,1)*size(v,2)*size(v,3));
        
        if size(datau,2)>20000 %more than 20000 value pairs are too slow in scatterplot.
            pos=unique(ceil(rand(21000,1)*(size(datau,2)-1))); %select random entries...
            scatter(datau(pos),datav(pos), 'b.'); %.. and plot them
        else
            scatter(datau,datav, 'b.');
        end
        
        %skipper=ceil(size(datau,2)/8000);
        %scatter(datau(:,1:skipper:end),datav(:,1:skipper:end), 'b.');
        
        oldsize=get(gca,'outerposition');
        newsize=[0 0 oldsize(3)*0.87 oldsize(4)*0.87];
        set(gca,'outerposition', newsize)
        %%{
        if retr('caluv')==1 && retr('calxy')==1
            xlabel(gca, 'u velocity [px/frame]', 'fontsize', 12)
            ylabel(gca, 'v velocity [px/frame]', 'fontsize', 12)
        else
            xlabel(gca, 'u velocity [m/s]', 'fontsize', 12)
            ylabel(gca, 'v velocity [m/s]', 'fontsize', 12)
        end
        
        grid on
        %axis equal;
        set (gca, 'tickdir', 'in');
        rangeu=nanmax(nanmax(nanmax(u*caluv)))-nanmin(nanmin(nanmin(u*caluv)));
        rangev=nanmax(nanmax(nanmax(v*caluv)))-nanmin(nanmin(nanmin(v*caluv)));
        %set(gca,'xlim',[nanmin(nanmin(nanmin(u*caluv)))-rangeu*0.15 nanmax(nanmax(nanmax(u*caluv)))+rangeu*0.15])
        %set(gca,'ylim',[nanmin(nanmin(nanmin(v*caluv)))-rangev*0.15 nanmax(nanmax(nanmax(v*caluv)))+rangev*0.15])
        %=range of data +- 15%
        %%}
        velrect = getrect(gca);
        if velrect(1,3)~=0 && velrect(1,4)~=0
            put('velrect', velrect);
            set (handles.vel_limit_active, 'String', 'Limit active', 'backgroundcolor', [0.5 1 0.5]);
            umin=velrect(1);
            umax=velrect(3)+umin;
            vmin=velrect(2);
            vmax=velrect(4)+vmin;
            if retr('caluv')==1 && retr('calxy')==1
                set (handles.limittext, 'String', ['valid u: ' num2str(round(umin*100)/100) ' to ' num2str(round(umax*100)/100) ' [px/frame]' sprintf('\n') 'valid v: ' num2str(round(vmin*100)/100) ' to ' num2str(round(vmax*100)/100) ' [px/frame]']);
            else
                set (handles.limittext, 'String', ['valid u: ' num2str(round(umin*100)/100) ' to ' num2str(round(umax*100)/100) ' [m/s]' sprintf('\n') 'valid v: ' num2str(round(vmin*100)/100) ' to ' num2str(round(vmax*100)/100) ' [m/s]']);
            end
            sliderdisp
            delete(findobj(gca,'Type','text','color','r'));
            text(50,50,'Result will be shown after applying vector validation','color','r','fontsize',10, 'fontweight','bold', 'BackgroundColor', 'k')
            set (handles.vel_limit, 'String', 'Refine velocity limits');
        else
            sliderdisp
            text(50,50,'Invalid selection: Click and hold left mouse button to create a rectangle.','color','r','fontsize',8, 'BackgroundColor', 'k')
        end
        
    end
end
toolsavailable(1)
figure1_ResizeFcn(gcf)

function apply_filter_current_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
put('derived', []); %clear derived parameters if user modifies source data
filtervectors(currentframe)
%put('manualdeletion',[]); %only valid one time, why...? Could work without this line.
sliderdisp;

function apply_filter_all_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
toolsavailable(0)
%put('manualdeletion',[]); %not available for filtering several images
put('derived', []); %clear derived parameters if user modifies source data
for i=1:floor(size(filepath,1)/2)+1
    filtervectors(i)
    set (handles.apply_filter_all, 'string', ['Please wait... (' int2str((i-1)/size(filepath,1)*200) '%)']);
    drawnow;
end
set (handles.apply_filter_all, 'string', 'Apply to all frames');
toolsavailable(1)
sliderdisp;

function restore_all_Callback(hObject, eventdata, handles)
%clears resultslist at 7,8,9
resultslist=retr('resultslist');

if size(resultslist,1) > 6
    resultslist(7:9,:)={[]};
    if size(resultslist,1) > 9
        resultslist(10:11,:)={[]};
    end
    put('resultslist', resultslist);
    sliderdisp
end
put('manualdeletion',[]);

% --- Executes on button press in clear_vel_limit.
function clear_vel_limit_Callback(hObject, eventdata, handles)
put('velrect', []);
handles=gethand;
set (handles.vel_limit_active, 'String', 'Limit inactive', 'backgroundcolor', [0.9411764705882353 0.9411764705882353 0.9411764705882353]);
set (handles.limittext, 'String', '');
set (handles.vel_limit, 'String', 'Select velocity limits');

function filtervectors(frame)
%executes filters one after another, writes results to resultslist 7,8,9
handles=gethand;
resultslist=retr('resultslist');
resultslist{10,frame}=[]; %remove smoothed results when user modifies original data
resultslist{11,frame}=[];
if size(resultslist,2)>=frame
    caluv=retr('caluv');
    u=resultslist{3,frame};
    v=resultslist{4,frame};
    typevector_original=resultslist{5,frame};
    typevector=typevector_original;
    if numel(u>0)
        %velocity limits
        velrect=retr('velrect');
        
        if numel(velrect>0) %velocity limits were activated
            umin=velrect(1);
            umax=velrect(3)+umin;
            vmin=velrect(2);
            vmax=velrect(4)+vmin;
            u(u*caluv<umin)=NaN;
            u(u*caluv>umax)=NaN;
            v(u*caluv<umin)=NaN;
            v(u*caluv>umax)=NaN;
            v(v*caluv<vmin)=NaN;
            v(v*caluv>vmax)=NaN;
            u(v*caluv<vmin)=NaN;
            u(v*caluv>vmax)=NaN;
        end
        %manual point deletion
        manualdeletion=retr('manualdeletion');
        
        if numel(manualdeletion)>0
            if size(manualdeletion,2)>=frame
                if isempty(manualdeletion{1,frame}) ==0
                    framemanualdeletion=manualdeletion{frame};
                    for i=1:size(framemanualdeletion,1)
                        u(framemanualdeletion(i,1),framemanualdeletion(i,2))=NaN;
                        v(framemanualdeletion(i,1),framemanualdeletion(i,2))=NaN;
                        
                        %                u(manualdeletion(i,1),manualdeletion(i,2))=NaN;
                        %                v(manualdeletion(i,1),manualdeletion(i,2))=NaN;
                    end
                end
            end
        end
        %% stddev check
        if get(handles.stdev_check, 'value')==1
            stdthresh=str2double(get(handles.stdev_thresh, 'String'));
            meanu=nanmean(nanmean(u));
            meanv=nanmean(nanmean(v));
            std2u=nanstd(reshape(u,size(u,1)*size(u,2),1));
            std2v=nanstd(reshape(v,size(v,1)*size(v,2),1));
            minvalu=meanu-stdthresh*std2u;
            maxvalu=meanu+stdthresh*std2u;
            minvalv=meanv-stdthresh*std2v;
            maxvalv=meanv+stdthresh*std2v;
            u(u<minvalu)=NaN;
            u(u>maxvalu)=NaN;
            v(v<minvalv)=NaN;
            v(v>maxvalv)=NaN;
        end
        %local median check
        %Westerweel & Scarano (2005): Universal Outlier detection for PIV data
        if get(handles.loc_median, 'value')==1
            epsilon=str2double(get(handles.epsilon,'string'));
            thresh=str2double(get(handles.loc_med_thresh,'string'));
            [J,I]=size(u);
            medianres=zeros(J,I);
            normfluct=zeros(J,I,2);
            b=1;
            eps=0.1;
            for c=1:2
                if c==1; velcomp=u;else;velcomp=v;end %#ok<*NOSEM>
                for i=1+b:I-b
                    for j=1+b:J-b
                        neigh=velcomp(j-b:j+b,i-b:i+b);
                        neighcol=neigh(:);
                        neighcol2=[neighcol(1:(2*b+1)*b+b);neighcol((2*b+1)*b+b+2:end)];
                        med=median(neighcol2);
                        fluct=velcomp(j,i)-med;
                        res=neighcol2-med;
                        medianres=median(abs(res));
                        normfluct(j,i,c)=abs(fluct/(medianres+epsilon));
                    end
                end
            end
            info1=(sqrt(normfluct(:,:,1).^2+normfluct(:,:,2).^2)>thresh);
            u(info1==1)=NaN;
            v(info1==1)=NaN;
        end
        %0=mask
        %1=normal
        %2=manually filtered
        u(isnan(v))=NaN;
        v(isnan(u))=NaN;
        typevector(isnan(u))=2;
        typevector(isnan(v))=2;
        typevector(typevector_original==0)=0; %restores typevector for mask
        %interpolation using inpaint_NaNs
        if get(handles.interpol_missing, 'value')==1
            u=inpaint_nans(u,4);
            v=inpaint_nans(v,4);
        end
        resultslist{7, frame} = u;
        resultslist{8, frame} = v;
        resultslist{9, frame} = typevector;
        put('resultslist', resultslist);
    else
    end
end
%sliderdisp

function rejectsingle_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
frame=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=frame %2nd dimesnion = frame
    x=resultslist{1,frame};
    y=resultslist{2,frame};
    u=resultslist{3,frame};
    v=resultslist{4,frame};
    typevector_original=resultslist{5,frame};
    typevector=typevector_original;
    manualdeletion=retr('manualdeletion');
    framemanualdeletion=[];
    if numel(manualdeletion)>0
        if size(manualdeletion,2)>=frame
            if isempty(manualdeletion{1,frame}) ==0
                framemanualdeletion=manualdeletion{frame};
            end
        end
    end
    
    if numel(u>0)
        delete(findobj(gca,'tag','manualdot'));
        text(50,10,'Right mouse button exits manual validation mode.','color','g','fontsize',8, 'BackgroundColor', 'k', 'tag', 'hint')
        toolsavailable(0);
        button = 1;
        while button == 1
            [xposition,yposition,button] = ginput(1);
            if button~=1
                break
            end
            if numel (xposition>0) %will be 0 if user presses enter
                xposition=round(xposition);
                yposition=round(yposition);
                %manualdeletion=zeros(size(xposition,1),2);
                findx=abs(x/xposition-1);
                [trash, imagex]=find(findx==min(min(findx)));
                findy=abs(y/yposition-1);
                [imagey, trash]=find(findy==min(min(findy)));
                idx=size(framemanualdeletion,1);
                %manualdeletion(idx+1,1)=imagey(1,1);
                %manualdeletion(idx+1,2)=imagex(1,1);
                
                framemanualdeletion(idx+1,1)=imagey(1,1);
                framemanualdeletion(idx+1,2)=imagex(1,1);
                
                hold on;
                plot (x(framemanualdeletion(idx+1,1),framemanualdeletion(idx+1,2)),y(framemanualdeletion(idx+1,1),framemanualdeletion(idx+1,2)), 'yo', 'markerfacecolor', 'r', 'markersize', 10,'tag','manualdot')
                hold off;
            end
        end
        manualdeletion{frame}=framemanualdeletion;
        put('manualdeletion',manualdeletion);
        
        delete(findobj(gca,'Type','text','color','r'));
        delete(findobj(gca,'tag','hint'));
        text(50,50,'Result will be shown after applying vector validation','color','r','fontsize',10, 'fontweight','bold', 'BackgroundColor', 'k')
    end
end
toolsavailable(1);

function draw_line_Callback(hObject, eventdata, handles)
filepath=retr('filepath');
caliimg=retr('caliimg');
if numel(caliimg)==0 && size(filepath,1) >1
    sliderdisp
end
if size(filepath,1) >1 || numel(caliimg)>0
    handles=gethand;
    toolsavailable(0)
    delete(findobj('tag', 'caliline'))
    for i=1:2
        [xposition(i),yposition(i)] = ginput(1);
        if numel(caliimg)==0
            sliderdisp
        end
        hold on;
        plot (xposition,yposition,'ro-', 'markersize', 10,'LineWidth',3, 'tag', 'caliline');
        plot (xposition,yposition,'y+:', 'tag', 'caliline');
        hold off;
        for j=1:i
            text(xposition(j)+10,yposition(j)+10, ['x:' num2str(round(xposition(j)*10)/10) sprintf('\n') 'y:' num2str(round(yposition(j)*10)/10) ],'color','y','fontsize',7, 'BackgroundColor', 'k', 'tag', 'caliline')
        end
        
        put('pointscali',[xposition' yposition']);
    end
    text(mean(xposition),mean(yposition), ['s = ' num2str(round((sqrt((xposition(1)-xposition(2))^2+(yposition(1)-yposition(2))^2))*100)/100) ' px'],'color','k','fontsize',7, 'BackgroundColor', 'r', 'tag', 'caliline','horizontalalignment','center')
    toolsavailable(1)
end

function calccali
put('derived',[]) %calibration makes previously derived params incorrect
handles=gethand;
pointscali=retr('pointscali');
if numel(pointscali)>0
    xposition=pointscali(:,1);
    yposition=pointscali(:,2);
    dist=sqrt((xposition(1)-xposition(2))^2 + (yposition(1)-yposition(2))^2);
    realdist=str2double(get(handles.realdist, 'String'));
    time=str2double(get(handles.time_inp, 'String'));
    calxy=(realdist/1000)/dist; %m/px %realdist=realdistance in m; dist=distance in px
    caluv=calxy/(time/1000);
    put('caluv',caluv);
    put('calxy',calxy);
    set(handles.calidisp, 'string', ['1 px = ' num2str(round(calxy*100000)/100000) ' m' sprintf('\n') '1 px/frame = ' num2str(round(caluv*100000)/100000) ' m/s'],  'backgroundcolor', [0.5 1 0.5]);
end
sliderdisp

function clear_cali_Callback(hObject, eventdata, handles)
handles=gethand;
put('pointscali',[]);
put('caluv',1);
put('calxy',1);
put('caliimg', []);
filepath=retr('filepath');
set(handles.calidisp, 'string', ['inactive'], 'backgroundcolor', [0.9411764705882353 0.9411764705882353 0.9411764705882353]);
delete(findobj('tag', 'caliline'));
set(handles.realdist, 'String','1');
set(handles.time_inp, 'String','1');
if size(filepath,1) >1
    sliderdisp
else
    displogo(0)
end


function load_ext_img_Callback(hObject, eventdata, handles)
cali_folder=retr('cali_folder');
if isempty (cali_folder)==1
    if ispc==1
        cali_folder=[retr('pathname') '\'];
    else
        cali_folder=[retr('pathname') '/'];
    end
end
try
    [filename, pathname, filterindex] = uigetfile({'*.bmp;*.tif;*.jpg;','Image Files (*.bmp,*.tif,*.jpg)'; '*.tif','tif'; '*.jpg','jpg'; '*.bmp','bmp'; },'Select calibration image',cali_folder);
catch
    [filename, pathname, filterindex] = uigetfile({'*.bmp;*.tif;*.jpg;','Image Files (*.bmp,*.tif,*.jpg)'; '*.tif','tif'; '*.jpg','jpg'; '*.bmp','bmp'; },'Select calibration image'); %unix/mac system may cause problems, can't be checked due to lack of unix/mac systems...
end
if isequal(filename,0)==0
    
    %caliimg=adapthisteq(imread(fullfile(pathname, filename)));
    caliimg=imread(fullfile(pathname, filename));
    if size(caliimg,3)>1 == 0
        caliimg=adapthisteq(caliimg);
    else
        try
        caliimg=adapthisteq(rgb2gray(caliimg));
        catch
        end
    end
    image(caliimg, 'parent',gca, 'cdatamapping', 'scaled');
    colormap('gray');
    axis image;
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    put('caliimg', caliimg);
    put('cali_folder', pathname);
end

function write_workspace_Callback(hObject, eventdata, handles)

resultslist=retr('resultslist');
if isempty(resultslist)==0
    derived=retr('derived');
    calxy=retr('calxy');
    caluv=retr('caluv');
    nrframes=size(resultslist,2);
    if size(resultslist,1)< 11
        resultslist{11,nrframes}=[]; %make sure resultslist has cells for all params
    end
    if isempty(derived)==0
        if size(derived,1)< 9 || size(derived,2) < nrframes
            derived{9,nrframes}=[]; %make sure derived has cells for all params
        end
    else
        derived=cell(9,nrframes);
    end
    
    if calxy==1 && caluv==1
        units='[px] respectively [px/frame]';
    else
        units='[m] respectively [m/s]';
    end
    %ohne alles: 6 hoch
    %mit filtern: 11 hoch
    %mit smoothed, 11 hoch und inhalt...
    u_original=cell(size(resultslist,2),1);
    v_original=u_original;
    x=u_original;
    y=u_original;
    typevector_original=u_original;
    u_filtered=u_original;
    v_filtered=v_original;
    typevector_filtered=u_original;
    u_smoothed=u_original;
    v_smoothed=u_original;
    vorticity=cell(size(derived,2),1);
    velocity_magnitude=vorticity;
    u_component=vorticity;
    v_component=vorticity;
    divergence=vorticity;
    vortex_locator=vorticity;
    shear_rate=vorticity;
    strain_rate=vorticity;
    LIC=vorticity;
    
    for i=1:nrframes
        x{i,1}=resultslist{1,i}*calxy;
        y{i,1}=resultslist{2,i}*calxy;
        u_original{i,1}=resultslist{3,i}*caluv;
        v_original{i,1}=resultslist{4,i}*caluv;
        typevector_original{i,1}=resultslist{5,i};
        u_filtered{i,1}=resultslist{7,i}*caluv;
        v_filtered{i,1}=resultslist{8,i}*caluv;
        typevector_filtered{i,1}=resultslist{9,i};
        u_smoothed{i,1}=resultslist{10,i}*caluv;
        v_smoothed{i,1}=resultslist{11,i}*caluv;
        
        vorticity{i,1}=derived{1,i};
        velocity_magnitude{i,1}=derived{2,i};
        u_component{i,1}=derived{3,i};
        v_component{i,1}=derived{4,i};
        divergence{i,1}=derived{5,i};
        vortex_locator{i,1}=derived{6,i};
        shear_rate{i,1}=derived{7,i};
        strain_rate{i,1}=derived{8,i};
        LIC{i,1}=derived{9,i};
    end
    
    assignin('base','x',x);
    assignin('base','y',y);
    assignin('base','u_original',u_original);
    assignin('base','v_original',v_original);
    assignin('base','typevector_original',typevector_original);
    assignin('base','u_filtered',u_filtered);
    assignin('base','v_filtered',v_filtered);
    assignin('base','typevector_filtered',typevector_filtered);
    assignin('base','u_smoothed',u_smoothed);
    assignin('base','v_smoothed',v_smoothed);
    
    assignin('base','vorticity',vorticity);
    
    assignin('base','velocity_magnitude',velocity_magnitude);
    assignin('base','u_component',u_component);
    assignin('base','v_component',v_component);
    assignin('base','divergence',divergence);
    assignin('base','vortex_locator',vortex_locator);
    assignin('base','shear_rate',shear_rate);
    assignin('base','strain_rate',strain_rate);
    assignin('base','LIC',LIC);
    
    assignin('base','calxy',calxy);
    assignin('base','caluv',caluv);
    assignin('base','units',units);
    
    
    clc
    disp('EXPLANATIONS:')
    disp(' ')
    disp('The first dimension of the variables is the frame number.')
    disp('The variables contain all data that was calculated in the PIVlab GUI.')
    disp('If some data was not calculated, the corresponding cell is empty.')
    disp('Typevector is 0 for masked vector, 1 for regular vector, 2 for filtered vector')
    disp(' ')
    disp('u_original and v_original are the unmodified velocities from the cross-correlation.')
    disp('u_filtered and v_filtered is the above incl. your data validation selection.')
    disp('u_smoothed and v_smoothed is the above incl. your smoothing selection.')
end

function mean_u_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0 %analysis exists
    if size(resultslist,1)>6 && numel(resultslist{7,currentframe})>0 %filtered exists
        u=resultslist{7,currentframe};
    else
        u=resultslist{3,currentframe};
    end
    caluv=retr('caluv');
    set(handles.subtr_u, 'string', num2str(nanmean(nanmean(u*caluv))));
else
    set(handles.subtr_u, 'string', '0');
end

function mean_v_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0 %analysis exists
    if size(resultslist,1)>6 && numel(resultslist{7,currentframe})>0 %filtered exists
        v=resultslist{8,currentframe};
    else
        v=resultslist{4,currentframe};
    end
    caluv=retr('caluv');
    set(handles.subtr_v, 'string', num2str(nanmean(nanmean(v*caluv))));
else
    set(handles.subtr_v, 'string', '0');
end

function derivative_calc (frame,deriv,update)
handles=gethand;
resultslist=retr('resultslist');
if size(resultslist,2)>=frame && numel(resultslist{1,frame})>0 %analysis exists
    filenames=retr('filenames');
    filepath=retr('filepath');
    derived=retr('derived');
    caluv=retr('caluv');
    calxy=retr('calxy');
    currentimage=imread(filepath{2*frame-1});
    x=resultslist{1,frame};
    y=resultslist{2,frame};
    %subtrayct mean u
    subtr_u=str2double(get(handles.subtr_u, 'string'));
    if isnan(subtr_u)
        subtr_u=0;set(handles.subtr_u, 'string', '0');
    end
    subtr_v=str2double(get(handles.subtr_v, 'string'));
    if isnan(subtr_v)
        subtr_v=0;set(handles.subtr_v, 'string', '0');
    end
    if size(resultslist,1)>6 && numel(resultslist{7,frame})>0 %filtered exists
        u=resultslist{7,frame};
        v=resultslist{8,frame};
        typevector=resultslist{9,frame};
    else
        u=resultslist{3,frame};
        v=resultslist{4,frame};
        typevector=resultslist{5,frame};
    end
    if get(handles.interpol_missing,'value')==1
        if any(any(isnan(u)))==1 || any(any(isnan(v)))==1
            if isempty(strfind(get(handles.apply_deriv_all,'string'), 'Please'))==1 && isempty(strfind(get(handles.ascii_all,'string'), 'Please'))==1 && isempty(strfind(get(handles.save_mat_all,'string'), 'Please'))==1%not in batch
                drawnow;
                msgbox('Your dataset contains NaNs. A vector interpolation will be performed automatically to interpolate missing vectors.', 'modal')
                uiwait
            end
            typevector_original=typevector;
            u(isnan(v))=NaN;
            v(isnan(u))=NaN;
            typevector(isnan(u))=2;
            typevector(typevector_original==0)=0;
            u=inpaint_nans(u,4);
            v=inpaint_nans(v,4);
            resultslist{7, frame} = u;
            resultslist{8, frame} = v;
            resultslist{9, frame} = typevector;
        end
    else
        if isempty(strfind(get(handles.apply_deriv_all,'string'), 'Please'))==1 && isempty(strfind(get(handles.ascii_all,'string'), 'Please'))==1 && isempty(strfind(get(handles.save_mat_all,'string'), 'Please'))==1%not in batch
            drawnow;
            msgbox('Your dataset contains NaNs. Derived parameters will have a lot of missing data. Redo the vector validation with the option to interpolate missing data turned on.', 'modal')
            uiwait
        end
    end
    if get(handles.smooth, 'Value') == 1
        smoothfactor=floor(get(handles.smoothstr, 'Value'));
        try
            
            u = smoothn(u,smoothfactor/10); %not supported in prehistoric Matlab versions like the one I have to use :'-(
            smoothfactor %#ok<*NOPRT>
            v = smoothn(v,smoothfactor/10); %not supported in prehistoric Matlab versions like the one I have to use :'-(
            clc
            
            disp ('Using smoothn.m from Damien Garcia for data smoothing.')
            
        catch
            h=fspecial('gaussian',smoothfactor+2,(smoothfactor+2)/7);
            u=imfilter(u,h,'replicate');
            v=imfilter(v,h,'replicate');
            clc
            disp ('Using Gaussian kernel for data smoothing (your Matlab version is pretty old btw...).')
        end
        resultslist{10,frame}=u; %smoothed u
        resultslist{11,frame}=v; %smoothed v
    else
        %careful if more things are added, [] replaced by {[]}
        resultslist{10,frame}=[]; %remove smoothed u
        resultslist{11,frame}=[]; %remove smoothed v
    end
    if deriv==1 %vectors only
        %do nothing
    end
    if deriv==2 %vorticity
        [curlz,cav]= curl(x*calxy,y*calxy,u*caluv,v*caluv);
        derived{1,frame}=curlz;
    end
    if deriv==3 %magnitude
        %andersrum, (u*caluv)-subtr_u
        derived{2,frame}=sqrt((u*caluv-subtr_u).^2+(v*caluv-subtr_v).^2);
    end
    if deriv==4
        derived{3,frame}=u*caluv-subtr_u;
    end
    if deriv==5
        derived{4,frame}=v*caluv-subtr_v;
    end
    if deriv==6
        derived{5,frame}=divergence(x*calxy,y*calxy,u*caluv,v*caluv);
    end
    if deriv==7
        derived{6,frame}=dcev(x*calxy,y*calxy,u*caluv,v*caluv);
    end
    if deriv==8
        derived{7,frame}=shear(x*calxy,y*calxy,u*caluv,v*caluv);
    end
    if deriv==9
        derived{8,frame}=strain(x*calxy,y*calxy,u*caluv,v*caluv);
    end
    if deriv==10
%{
        A=rescale_maps(LIC(v*caluv-subtr_v,u*caluv-subtr_u,frame));
        [curlz,cav]= curl(x*calxy,y*calxy,u*caluv,v*caluv);
        B= rescale_maps(curlz);
        
        C=B-min(min(B));
        C=C/max(max(C));
        RGB_B = ind2rgb(uint8(C*255),colormap('jet'));
        RGB_A = ind2rgb(uint8(A*255),colormap('gray'));
%}
 
        
        %EDITED for williams visualization
        %Original: 
        derived{9,frame}=LIC(v*caluv-subtr_v,u*caluv-subtr_u,frame);
        
        
        
    end
    put('subtr_u', subtr_u);
    put('subtr_v', subtr_v);
    put('resultslist', resultslist);
    put ('derived',derived);
    if update==1
        put('displaywhat', deriv);
    end
end
function out=LIC(vx,vy,frame)
handles=gethand;
LICreso=round(get (handles.licres, 'value')*10)/10;
resultslist=retr('resultslist');
x=resultslist{1,frame};
y=resultslist{2,frame};
text(mean(x(1,:)/1.5),mean(y(:,1)), ['Please wait. LIC in progress...' sprintf('\n') 'If this message stays here for > 20s,' sprintf('\n') 'check MATLABs command window.' sprintf('\n') 'The function might need to be compiled first.'],'tag', 'waitplease', 'backgroundcolor', 'k', 'color', 'r','fontsize',10);
drawnow;
iterations=2;
set(gca,'units','pixels');
axessize=get(gca,'position');
set(gca,'units','points');
axessize=axessize(3:4);
%was ist größer, x oder y. dann entsprechend die x oder y größe der axes nehemn
xextend=size(vx,2);
yextend=size(vx,1);
if yextend<xextend
    scalefactor=axessize(1)/xextend;
else
    scalefactor=axessize(2)/yextend;
end

vx=imresize(vx,scalefactor*LICreso,'bicubic');
vy=imresize(vy,scalefactor*LICreso,'bicubic');


%{
this function is from:
Matlab VFV Toolbox 1.0
By Nima Bigdely Shamlo (email: bigdelys-vfv@yahoo.com)
Computational Science Research Center
San Diego State University
Under GNU GENERAL PUBLIC LICENSE
%}

[width,height] = size(vx);
LIClength = round(max([width,height]) / 30);

kernel = ones(2 * LIClength);
LICImage = zeros(width, height);
intensity = ones(width, height); % array containing vector intensity

% Making white noise
rand('state',0) % reset random generator to original state
noiseImage=rand(width,height);

% Making LIC Image
try
for m = 1:iterations
    [LICImage, intensity,normvx,normvy] = fastLICFunction(vx,vy,noiseImage,kernel); % External Fast LIC implemennted in C language
    LICImage = imadjust(LICImage); % Adjust the value range
    noiseImage = LICImage;
end
out=LICImage;
delete(findobj('tag', 'waitplease'));
catch
    h=errordlg(['Could not run the LIC tool.' sprintf('\n') 'Probably the tool is not compiled correctly.' sprintf('\n')  'Please execute the following command in Matlab:' sprintf('\n') sprintf('\n') '     mex fastLICFunction.c     ' sprintf('\n') sprintf('\n') 'Then try again.'],'Error','on');

    uiwait(h);
    out=zeros(size(vx));
end

function apply_deriv_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
deriv=get(handles.derivchoice, 'value');
derivative_calc (currentframe,deriv,1)
sliderdisp

function out=dcev(x,y,u,v);
dUdX=conv2(u,[ 0, 0, 0;-1, 0, 1; 0, 0, 0],'valid')./...
    conv2(x,[ 0, 0, 0;-1, 0, 1; 0, 0, 0],'valid');
dVdX=conv2(v,[ 0, 0, 0;-1, 0, 1; 0, 0, 0],'valid')./...
    conv2(x,[ 0, 0, 0;-1, 0, 1; 0, 0, 0],'valid');
dUdY=conv2(u,[ 0,-1, 0; 0, 0, 0; 0, 1, 0],'valid')./...
    conv2(y,[ 0,-1, 0; 0, 0, 0; 0, 1, 0],'valid');
dVdY=conv2(v,[ 0,-1, 0; 0, 0, 0; 0, 1, 0],'valid')./...
    conv2(y,[ 0,-1, 0; 0, 0, 0; 0, 1, 0],'valid');
res=(dUdX+dVdY)/2+sqrt(0.25*(dUdX+dVdY).^2+dUdY.*dVdX);
d=zeros(size(x));
d(2:end-1,2:end-1)=imag(res);
out=((d/(max(max(d))-(min(min(d)))))+abs(min(min(d))))*255;%normalize


function out=strain(x,y,u,v)
hx = x(1,:);
hy = y(:,1);
[px junk] = gradient(u, hx, hy);
[junk qy] = gradient(v, hx, hy); %#ok<*ASGLU>
out = px-qy;

function out=shear(x,y,u,v)
hx = x(1,:);
hy = y(:,1);
[junk py] = gradient(u, hx, hy);
[qx junk] = gradient(v, hx, hy);
out= qx+py;

function out=rescale_maps(in)
%input has same dimensions as x,y,u,v,
%output has size of the piv image
handles=gethand;
filepath=retr('filepath');
currentframe=floor(get(handles.fileselector, 'value'));
currentimage=imread(filepath{2*currentframe-1});
resultslist=retr('resultslist');
x=resultslist{1,currentframe};
y=resultslist{2,currentframe};
out=zeros(size(currentimage));
if size(out,3)>1
    out(:,:,2:end)=[];
end
out(:,:)=mean(mean(in));
step=x(1,2)-x(1,1)+1;
minx=(min(min(x))-step/2);
maxx=(max(max(x))+step/2);
miny=(min(min(y))-step/2);
maxy=(max(max(y))+step/2);
width=maxx-minx;
height=maxy-miny;
if size(in,3)>1 %why would this actually happen...?
    in(:,:,2:end)=[];
end
dispvar = imresize(in,[height width],'bilinear');
if miny<1
    miny=1;
end
if minx<1
    minx=1;
end;
try
    out(floor(miny):floor(maxy-1),floor(minx):floor(maxx-1))=dispvar;
catch
    disp('temp workaround')
    A=out(floor(miny):floor(maxy-1),floor(minx):floor(maxx-1));
    out(floor(miny):floor(maxy-1),floor(minx):floor(maxx-1))=dispvar(1:size(A,1),1:size(A,2));
end
maskiererx=retr('maskiererx');
if numel(maskiererx)>0
    if get(handles.img_not_mask, 'value')==1 && numel(maskiererx{currentframe*2-1})>0
        maskierery=retr('maskierery');
        ximask=maskiererx{currentframe*2-1};
        yimask=maskierery{currentframe*2-1};
        BW=poly2mask(ximask,yimask,size(out,1),size(out,2));
        max_img=double(max(max(currentimage)));
        max_map=max(max(out));
        currentimage=double(currentimage)/max_img*max_map;
        out(BW==1)=currentimage(BW==1);
    end
end

function out=rescale_maps_nan(in)
%input has same dimensions as x,y,u,v,
%output has size of the piv image
handles=gethand;
filepath=retr('filepath');
currentframe=floor(get(handles.fileselector, 'value'));
currentimage=imread(filepath{2*currentframe-1});
resultslist=retr('resultslist');
x=resultslist{1,currentframe};
y=resultslist{2,currentframe};
out=zeros(size(currentimage));
if size(out,3)>1
    out(:,:,2:end)=[];
end
out(:,:)=nan;
step=x(1,2)-x(1,1)+1;
minx=(min(min(x))-step/2);
maxx=(max(max(x))+step/2);
miny=(min(min(y))-step/2);
maxy=(max(max(y))+step/2);
width=maxx-minx;
height=maxy-miny;
if size(in,3)>1 %why would this actually happen...?
    in(:,:,2:end)=[];
end
dispvar = imresize(in,[height width],'bilinear');
if miny<1
    miny=1;
end
if minx<1
    minx=1;
end;
try
    out(floor(miny):floor(maxy-1),floor(minx):floor(maxx-1))=dispvar;
catch
    disp('temp workaround')
    A=out(floor(miny):floor(maxy-1),floor(minx):floor(maxx-1));
    out(floor(miny):floor(maxy-1),floor(minx):floor(maxx-1))=dispvar(1:size(A,1),1:size(A,2));
end

maskiererx=retr('maskiererx');
if numel(maskiererx)>0
    try
        if numel(maskiererx{currentframe*2-1})>0
            maskierery=retr('maskierery');
            ximask=maskiererx{currentframe*2-1};
            yimask=maskierery{currentframe*2-1};
            BW=poly2mask(ximask,yimask,size(out,1),size(out,2));
            out(BW==1)=nan;
        end
    catch
    end
end


function apply_cali_Callback(hObject, eventdata, handles)
calccali

function apply_deriv_all_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
toolsavailable(0)
for i=1:floor(size(filepath,1)/2)+1
    deriv=get(handles.derivchoice, 'value');
    derivative_calc(i,deriv,1)
    set (handles.apply_deriv_all, 'string', ['Please wait... (' int2str((i-1)/size(filepath,1)*200) '%)']);
    drawnow;
end
set (handles.apply_deriv_all, 'string', 'Apply to all frames');
toolsavailable(1)
sliderdisp

function autoscaler_Callback(hObject, eventdata, handles)
handles=gethand;
if get(handles.autoscaler, 'value')==1
    set (handles.mapscale_min, 'enable', 'off')
    set (handles.mapscale_max, 'enable', 'off')
else
    set (handles.mapscale_min, 'enable', 'on')
    set (handles.mapscale_max, 'enable', 'on')
end

function figure1_CloseRequestFcn(hObject, eventdata, handles)
button = questdlg('Do you want to quit PIVlab?','Quit?','Yes','Cancel','Cancel');
toolsavailable(1)
if strmatch(button,'Yes')==1
    try
        homedir=retr('homedir');
        pathname=retr('pathname');
        dlmwrite([homedir '/last.nf'], homedir, 'delimiter', '', 'precision', 6, 'newline', 'pc')
        dlmwrite([homedir '/last.nf'], pathname, '-append', 'delimiter', '', 'precision', 6, 'newline', 'pc')
    catch
    end
    delete(hObject);
end

function vectorscale_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    sliderdisp
end

function ascii_current_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    [FileName,PathName] = uiputfile('*.txt','Save vector data as...','PIVlab.txt'); %framenummer in dateiname
    if isequal(FileName,0) | isequal(PathName,0) %#ok<*OR2>
    else
        file_save(currentframe,FileName,PathName,1);
    end
end

function ascii_all_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
resultslist=retr('resultslist');
[FileName,PathName] = uiputfile('*.txt','Save vector data as...','PIVlab.txt'); %framenummer in dateiname
if isequal(FileName,0) | isequal(PathName,0)
else
    toolsavailable(0)
    for i=1:floor(size(filepath,1)/2)
        %if analysis exists
        if size(resultslist,2)>=i && numel(resultslist{1,i})>0
            [Dir Name Ext] = fileparts(FileName);
            FileName_nr=[Name sprintf('_%.4d', i) Ext];
            file_save(i,FileName_nr,PathName,1)
            set (handles.ascii_all, 'string', ['Please wait... (' int2str((i-1)/size(filepath,1)*200) '%)']);
            drawnow;
        end
    end
    toolsavailable(1)
    set (handles.ascii_all, 'string', 'Export all frames');
end

function save_mat_current_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    [FileName,PathName] = uiputfile('*.mat','Save MATLAB file as...','PIVlab.mat'); %framenummer in dateiname
    if isequal(FileName,0) | isequal(PathName,0)
    else
        file_save(currentframe,FileName,PathName,2);
    end
end

function save_mat_all_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
resultslist=retr('resultslist');
[FileName,PathName] = uiputfile('*.mat','Save MATLAB file as...','PIVlab.mat'); %framenummer in dateiname
if isequal(FileName,0) | isequal(PathName,0)
else
    toolsavailable(0)
    for i=1:floor(size(filepath,1)/2)
        %if analysis exists
        if size(resultslist,2)>=i && numel(resultslist{1,i})>0
            [Dir Name Ext] = fileparts(FileName);
            FileName_nr=[Name sprintf('_%.4d', i) Ext];
            file_save(i,FileName_nr,PathName,2)
            set (handles.save_mat_all, 'string', ['Please wait... (' int2str((i-1)/size(filepath,1)*200) '%)']);
            drawnow;
        end
    end
    toolsavailable(1)
    set (handles.save_mat_all, 'string', 'Save all frames');
end

% --- Executes on button press in set_points.
function set_points_Callback(hObject, eventdata, handles)
sliderdisp
hold on;
toolsavailable(0)
delete(findobj('tag', 'measure'));
n=0;
for i=1:2
    [xi,yi,but] = ginput(1);
    n=n+1;
    xposition(n)=xi;
    yposition(n)=yi;
    plot(xposition(n),yposition(n), 'r*','Color', [0.55,0.75,0.9], 'tag', 'measure');
    line(xposition,yposition,'LineWidth',3, 'Color', [0.05,0,0], 'tag', 'measure');
    line(xposition,yposition,'LineWidth',1, 'Color', [0.05,0.75,0.05], 'tag', 'measure');
end
line([xposition(1,1) xposition(1,2)],[yposition(1,1) yposition(1,1)], 'LineWidth',3, 'Color', [0.05,0.0,0.0], 'tag', 'measure');
line([xposition(1,1) xposition(1,2)],[yposition(1,1) yposition(1,1)], 'LineWidth',1, 'Color', [0.95,0.05,0.01], 'tag', 'measure');
line([xposition(1,2) xposition(1,2)], yposition,'LineWidth',3, 'Color',[0.05,0.0,0], 'tag', 'measure');
line([xposition(1,2) xposition(1,2)], yposition,'LineWidth',1, 'Color',[0.35,0.35,1], 'tag', 'measure');
hold off;
toolsavailable(1)
deltax=abs(xposition(1,1)-xposition(1,2));
deltay=abs(yposition(1,1)-yposition(1,2));
length=sqrt(deltax^2+deltay^2);
alpha=(180/pi) *(acos(deltax/length));
beta=(180/pi) *(asin(deltax/length));
handles=gethand;
calxy=retr('calxy');
if retr('caluv')==1 && retr('calxy')==1
    set (handles.deltax, 'String', [num2str(deltax*calxy) ' [px]']);
    set (handles.deltay, 'String', [num2str(deltay*calxy) ' [px]']);
    set (handles.length, 'String', [num2str(length*calxy) ' [px]']);
    
else
    set (handles.deltax, 'String', [num2str(deltax*calxy) ' [m]']);
    set (handles.deltay, 'String', [num2str(deltay*calxy) ' [m]']);
    set (handles.length, 'String', [num2str(length*calxy) ' [m]']);
end
set (handles.alpha, 'String', num2str(alpha));
set (handles.beta, 'String', num2str(beta));
% --- Executes on button press in draw_stuff.
function draw_stuff_Callback(hObject, eventdata, handles)
sliderdisp;
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    toolsavailable(0);
    xposition=[];
    yposition=[];
    n = 0;
    but = 1;
    hold on;
    if get(handles.draw_what,'value')==1 %polyline
        while but == 1
            [xi,yi,but] = ginput(1);
            if but~=1
                break
            end
            delete(findobj('tag', 'extractpoint'))
            plot(xi,yi,'r+','tag','extractpoint')
            n = n+1;
            xposition(n)=xi;
            yposition(n)=yi;
            delete(findobj('tag', 'extractline'))
            delete(findobj('tag','areaint'));
            line(xposition,yposition,'LineWidth',3, 'Color', [0,0,0.95],'tag','extractline');
            line(xposition,yposition,'LineWidth',1, 'Color', [0.95,0.5,0.01],'tag','extractline');
        end
    elseif get(handles.draw_what,'value')==2 %circle
        for i=1:2
            [xi,yi,but] = ginput(1);
            if i==1;delete(findobj('tag', 'extractpoint'));end
            n=n+1;
            xposition_raw(n)=xi;
            yposition_raw(n)=yi;
            plot(xposition_raw(n),yposition_raw(n), 'r+', 'MarkerSize',10,'tag','extractpoint');
        end
        deltax=abs(xposition_raw(1,1)-xposition_raw(1,2));
        deltay=abs(yposition_raw(1,1)-yposition_raw(1,2));
        radius=sqrt(deltax^2+deltay^2);
        valtable=linspace(0,2*pi,361);
        for i=1:size(valtable,2)
            xposition(1,i)=sin(valtable(1,i))*radius+xposition_raw(1,1);
            yposition(1,i)=cos(valtable(1,i))*radius+yposition_raw(1,1);
        end
        delete(findobj('tag', 'extractline'))
        line(xposition,yposition,'LineWidth',3, 'Color', [0,0,0.95],'tag','extractline');
        line(xposition,yposition,'LineWidth',1, 'Color', [0.95,0.5,0.01],'tag','extractline');
        text(xposition(1,1),yposition(1,1),'\leftarrow start/end','FontSize',7, 'Rotation', 90, 'BackgroundColor',[1 1 1],'tag','extractline')
        text(xposition(1,1),yposition(1,1)+8,'\rightarrow','FontSize',7, 'BackgroundColor',[1 1 1],'tag','extractline')
        text(xposition(1,1),yposition(1,1)-8-radius*2,'\leftarrow','FontSize',7, 'BackgroundColor',[1 1 1],'tag','extractline')
        text(xposition(1,1)-radius-8,yposition(1,1)-radius,'\leftarrow','FontSize',7, 'BackgroundColor',[1 1 1], 'Rotation', 90,'tag','extractline')
        text(xposition(1,1)+radius+8,yposition(1,1)-radius,'\rightarrow','FontSize',7, 'BackgroundColor',[1 1 1], 'Rotation', 90,'tag','extractline')
    elseif get(handles.draw_what,'value')==3 %circle series
        set(handles.extraction_choice,'Value',9);
        for i=1:2
            [xi,yi,but] = ginput(1);
            if i==1;delete(findobj('tag', 'extractpoint'));end
            n=n+1;
            xposition_raw(n)=xi;
            yposition_raw(n)=yi;
            plot(xposition_raw(n),yposition_raw(n), 'r+', 'MarkerSize',10,'tag','extractpoint');
        end
        deltax=abs(xposition_raw(1,1)-xposition_raw(1,2));
        deltay=abs(yposition_raw(1,1)-yposition_raw(1,2));
        radius=sqrt(deltax^2+deltay^2);
        valtable=linspace(0,2*pi,361);
        for m=1:30
            for i=1:size(valtable,2)
                xposition(m,i)=sin(valtable(1,i))*(radius-((30-m)/30)*radius)+xposition_raw(1,1);
                yposition(m,i)=cos(valtable(1,i))*(radius-((30-m)/30)*radius)+yposition_raw(1,1);
            end
        end
        delete(findobj('tag', 'extractline'))
        for m=1:30
            line(xposition(m,:),yposition(m,:),'LineWidth',1.5, 'Color', [0.95,0.5,0.01],'tag','extractline');
        end
        text(xposition(30,1),yposition(30,1),'\leftarrow start/end','FontSize',7, 'Rotation', 90, 'BackgroundColor',[1 1 1],'tag','extractline')
        text(xposition(30,1),yposition(30,1)+8,'\rightarrow','FontSize',7, 'BackgroundColor',[1 1 1],'tag','extractline')
        text(xposition(30,1),yposition(30,1)-8-radius*2,'\leftarrow','FontSize',7, 'BackgroundColor',[1 1 1],'tag','extractline')
        text(xposition(30,1)-radius-8,yposition(30,1)-radius,'\leftarrow','FontSize',7, 'BackgroundColor',[1 1 1], 'Rotation', 90,'tag','extractline')
        text(xposition(30,1)+radius+8,yposition(30,1)-radius,'\rightarrow','FontSize',7, 'BackgroundColor',[1 1 1], 'Rotation', 90,'tag','extractline')
    end
    hold off;
    put('xposition',xposition)
    put('yposition',yposition)
    toolsavailable(1);
end

% --- Executes on button press in save_data.
function save_data_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));

if get(handles.extractLineAll, 'value')==0
    startfr=currentframe;
    endfr=currentframe;
else
    startfr=1;
    endfr=size(resultslist,2);
end
selected=0;
for i=startfr:endfr
    set(handles.fileselector, 'value',i)
    %sliderdisp
    currentframe=floor(get(handles.fileselector, 'value'));
    if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
        delete(findobj('tag', 'derivplotwindow'));
        plot_data_Callback %make sure that data was calculated
        %close figure...
        delete(findobj('tag', 'derivplotwindow'));
        extractwhat=get(handles.extraction_choice,'Value');
        current=get(handles.extraction_choice,'string');
        current=current{extractwhat};
        if selected==0
            imgsavepath=retr('imgsavepath');
            if isempty(imgsavepath)
                imgsavepath=retr('pathname');
            end
            %find '\', replace with 'per'
            part1= current(1:strfind(current,'/')-1) ;
            part2= current(strfind(current,'/')+1:end);
            if isempty(part1)==1
                currentED=current;
            else
                currentED=[part1 ' per ' part2];
            end
            [FileName,PathName] = uiputfile('*.txt','Save extracted data as...',fullfile(imgsavepath,['PIVlab_Extr_' currentED '.txt'])); %framenummer in dateiname
            selected=1;
        end
        if isequal(FileName,0) | isequal(PathName,0)
            %exit for
            break;
        else
            put('imgsavepath',PathName);
            pointpos=strfind(FileName, '.');
            pointpos=pointpos(end);
            FileName_final=[FileName(1:pointpos-1) '_' num2str(currentframe) '.' FileName(pointpos+1:end)];
            c=retr('c');
            distance=retr('distance');
            if size(c,2)>1
                if retr('calxy')==1 && retr('caluv')==1
                    header=['circle nr.,' 'Distance on line [px],' current];
                else
                    header=['circle nr.,' 'Distance on line [m],' current];
                end
                wholeLOT=[];
                for z=1:30
                    wholeLOT=[wholeLOT;[linspace(z,z,size(c,2))' distance(z,:)' c(z,:)']]; %anders.... untereinander
                end
            else
                if retr('calxy')==1 && retr('caluv')==1
                    header=['Distance on line [px],' current];
                else
                    header=['Distance on line [m],' current];
                end
                wholeLOT=[distance c];
            end
            fid = fopen(fullfile(PathName,FileName_final), 'w');
            fprintf(fid, [header '\r\n']);
            fclose(fid);
            dlmwrite(fullfile(PathName,FileName_final), wholeLOT, '-append', 'delimiter', ',', 'precision', 10, 'newline', 'pc');
        end
    end
end
sliderdisp

% --- Executes on button press in plot_data.
function plot_data_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    x=resultslist{1,currentframe};
    y=resultslist{2,currentframe};
    xposition=retr('xposition'); %not conflicting...?
    yposition=retr('yposition'); %not conflicting...?
    if numel(xposition)>1
        for i=1:size(xposition,2)-1
            %length of one segment:
            laenge(i)=sqrt((xposition(1,i+1)-xposition(1,i))^2+(yposition(1,i+1)-yposition(1,i))^2);
        end
        length=sum(laenge);
        percentagex=xposition/max(max(x));
        xaufderivative=percentagex*size(x,2);
        percentagey=yposition/max(max(y));
        yaufderivative=percentagey*size(y,1);
        nrpoints=str2num(get(handles.nrpoints,'string'));
        if get(handles.draw_what,'value')==3 %circle series
            set(handles.extraction_choice,'Value',9); %set to tangent
        end
        extractwhat=get(handles.extraction_choice,'Value');
        switch extractwhat
            case {1,2,3,4,5,6,7,8}
                derivative_calc(currentframe,extractwhat+1,0);
                derived=retr('derived');
                maptoget=derived{extractwhat,currentframe};
                
                maptoget=rescale_maps_nan(maptoget);
                [cx, cy, c] = improfile(maptoget,xposition,yposition,round(nrpoints),'bicubic');
                
                distance=linspace(0,length,size(c,1))';
                
            case 9 %tangent
                if size(xposition,1)<=1 %user did not choose circle series
                    if size(resultslist,1)>6 %filtered exists
                        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
                            u=resultslist{10,currentframe};
                            v=resultslist{11,currentframe};
                            typevector=resultslist{9,currentframe};
                            if numel(typevector)==0%happens if user smoothes sth without NaN and without validation
                                typevector=resultslist{5,currentframe};
                            end
                        else
                            u=resultslist{7,currentframe};
                            if size(u,1)>1
                                v=resultslist{8,currentframe};
                                typevector=resultslist{9,currentframe};
                            else
                                u=resultslist{3,currentframe};
                                v=resultslist{4,currentframe};
                                typevector=resultslist{5,currentframe};
                            end
                        end
                    else
                        u=resultslist{3,currentframe};
                        v=resultslist{4,currentframe};
                        typevector=resultslist{5,currentframe};
                    end
                    caluv=retr('caluv');
                    u=u*caluv-retr('subtr_u');
                    v=v*caluv-retr('subtr_v');
                    
                    u=rescale_maps_nan(u);
                    v=rescale_maps_nan(v);
                    
                    [cx, cy, cu] = improfile(u,xposition,yposition,round(nrpoints),'bicubic');
                    cv = improfile(v,xposition,yposition,round(nrpoints),'bicubic');
                    cx=cx';
                    cy=cy';
                    deltax=zeros(1,size(cx,2)-1);
                    deltay=zeros(1,size(cx,2)-1);
                    laenge=zeros(1,size(cx,2)-1);
                    alpha=zeros(1,size(cx,2)-1);
                    sinalpha=zeros(1,size(cx,2)-1);
                    cosalpha=zeros(1,size(cx,2)-1);
                    for i=2:size(cx,2)
                        deltax(1,i)=cx(1,i)-cx(1,i-1);
                        deltay(1,i)=cy(1,i)-cy(1,i-1);
                        laenge(1,i)=sqrt(deltax(1,i)*deltax(1,i)+deltay(1,i)*deltay(1,i));
                        alpha(1,i)=(acos(deltax(1,i)/laenge(1,i)));
                        if deltay(1,i) < 0
                            sinalpha(1,i)=sin(alpha(1,i));
                        else
                            sinalpha(1,i)=sin(alpha(1,i))*-1;
                        end
                        cosalpha(1,i)=cos(alpha(1,i));
                    end
                    sinalpha(1,1)=sinalpha(1,2);
                    cosalpha(1,1)=cosalpha(1,2);
                    cu=cu.*cosalpha';
                    cv=cv.*sinalpha';
                    c=cu-cv;
                    cx=cx';
                    cy=cy';
                    distance=linspace(0,length,size(cu,1))';
                else %user chose circle series
                    for m=1:30
                        for i=1:size(xposition,2)-1
                            %length of one segment:
                            laenge(m,i)=sqrt((xposition(m,i+1)-xposition(m,i))^2+(yposition(m,i+1)-yposition(m,i))^2);
                        end
                        length(m)=sum(laenge(m,:));
                    end
                    if size(resultslist,1)>6 %filtered exists
                        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
                            u=resultslist{10,currentframe};
                            v=resultslist{11,currentframe};
                            typevector=resultslist{9,currentframe};
                            if numel(typevector)==0%happens if user smoothes sth without NaN and without validation
                                typevector=resultslist{5,currentframe};
                            end
                        else
                            u=resultslist{7,currentframe};
                            if size(u,1)>1
                                v=resultslist{8,currentframe};
                                typevector=resultslist{9,currentframe};
                            else
                                u=resultslist{3,currentframe};
                                v=resultslist{4,currentframe};
                                typevector=resultslist{5,currentframe};
                            end
                        end
                    else
                        u=resultslist{3,currentframe};
                        v=resultslist{4,currentframe};
                        typevector=resultslist{5,currentframe};
                    end
                    caluv=retr('caluv');
                    u=u*caluv-retr('subtr_u');
                    v=v*caluv-retr('subtr_v');
                    u=rescale_maps_nan(u);
                    v=rescale_maps_nan(v);
                    min_y=floor(min(min(yposition)))-1;
                    max_y=ceil(max(max(yposition)))+1;
                    min_x=floor(min(min(xposition)))-1;
                    max_x=ceil(max(max(xposition)))+1;
                    if min_y<1
                        min_y=1;
                    end
                    if max_y>size(u,1)
                        max_y=size(u,1);
                    end
                    if min_x<1
                        min_x=1;
                    end
                    if max_x>size(u,2)
                        max_x=size(u,2);
                    end
                    
                    uc=u(min_y:max_y,min_x:max_x);
                    vc=v(min_y:max_y,min_x:max_x);
                    for m=1:30
                        [cx(m,:),cy(m,:),cu(m,:)] = improfile(uc,xposition(m,:)-min_x,yposition(m,:)-min_y,100,'nearest');
                        cv(m,:) =improfile(vc,xposition(m,:)-min_x,yposition(m,:)-min_y,100,'nearest');
                    end
                    deltax=zeros(1,size(cx,2)-1);
                    deltay=zeros(1,size(cx,2)-1);
                    laenge=zeros(1,size(cx,2)-1);
                    alpha=zeros(1,size(cx,2)-1);
                    sinalpha=zeros(1,size(cx,2)-1);
                    cosalpha=zeros(1,size(cx,2)-1);
                    for m=1:30
                        for i=2:size(cx,2)
                            deltax(m,i)=cx(m,i)-cx(m,i-1);
                            deltay(m,i)=cy(m,i)-cy(m,i-1);
                            laenge(m,i)=sqrt(deltax(m,i)*deltax(m,i)+deltay(m,i)*deltay(m,i));
                            alpha(m,i)=(acos(deltax(m,i)/laenge(m,i)));
                            if deltay(m,i) < 0
                                sinalpha(m,i)=sin(alpha(m,i));
                            else
                                sinalpha(m,i)=sin(alpha(m,i))*-1;
                            end
                            cosalpha(m,i)=cos(alpha(m,i));
                        end
                        sinalpha(m,1)=sinalpha(m,2); %ersten winkel füllen
                        cosalpha(m,1)=cosalpha(m,2);
                    end
                    cu=cu.*cosalpha;
                    cv=cv.*sinalpha;
                    c=cu-cv;
                    for m=1:30
                        distance(m,:)=linspace(0,length(m),size(cu,2))'; %in pixeln...
                    end
                end
                
        end
        if size(xposition,1)<=1 %user did not choose circle series
            h=figure;
            screensize=get( 0, 'ScreenSize' );
            rect = [screensize(3)/4-300, screensize(4)/2-250, 600, 500];
            set(h,'position', rect);
            current=get(handles.extraction_choice,'string');
            current=current{extractwhat};
            set(h,'numbertitle','off','menubar','none','toolbar','figure','dockcontrols','off','name',[current ', frame ' num2str(currentframe)],'tag', 'derivplotwindow');
            calxy=retr('calxy');
            
            %removing nans for integral!
            distance2=distance(isnan(c)==0);
            c2=c(isnan(c)==0);
            
            integral=trapz(distance2*calxy,c2);
            h2=plot(distance*calxy,c);
            
            %get units
            if retr('caluv')==1 && retr('calxy')==1
                distunit='px^2';
            else
                distunit='m^2';
            end
            
            unitpar=get(handles.extraction_choice,'string');
            unitpar=unitpar{get(handles.extraction_choice,'value')};
            unitpar=unitpar(strfind(unitpar,'[')+1:end-1);
            
            %text(0+0.05*max(distance*calxy),min(c)+0.05*max(c),['Integral = ' num2str(integral) ' [' unitpar '*' distunit ']'], 'BackgroundColor', 'w','fontsize',7)
            set (gca, 'xgrid', 'on', 'ygrid', 'on', 'TickDir', 'in')
            
            
            if retr('caluv')==1 && retr('calxy')==1
                distunit_2=' [px]';
            else
                distunit_2=' [m]';
            end
            
            currentstripped=current(1:strfind(current,'[')-1);
            
            %modified units...
            xlabel(['Distance on line' distunit_2 sprintf('\n') 'Integral of ' currentstripped ' = ' num2str(integral) ' [' unitpar '*' distunit_2 ']']);
            
            ylabel(current);
            put('distance',distance*calxy);
            put('c',c);
            h_extractionplot=retr('h_extractionplot');
            h_extractionplot(size(h_extractionplot,1)+1,1)=h;
            put ('h_extractionplot', h_extractionplot);
        else %user chose circle series
            calxy=retr('calxy');
            for m=1:30
                integral(m)=trapz(distance(m,:)*calxy,c(m,:));
            end
            %highlight circle with highest circ
            delete(findobj('tag', 'extractline'))
            for m=1:30
                line(xposition(m,:),yposition(m,:),'LineWidth',1.5, 'Color', [0.95,0.5,0.01],'tag','extractline');
            end
            [r,col]=find(max(abs(integral))==abs(integral)); %find absolute max of integral
            line(xposition(col,:),yposition(col,:),'LineWidth',2.5, 'Color', [0.2,0.5,0.7],'tag','extractline');
            h=figure;
            screensize=get( 0, 'ScreenSize' );
            rect = [screensize(3)/4-300, screensize(4)/2-250, 600, 500];
            set(h,'position', rect);
            current=get(handles.extraction_choice,'string');
            current=current{extractwhat};
            set(h,'numbertitle','off','menubar','none','toolbar','figure','dockcontrols','off','name',[current ', frame ' num2str(currentframe)],'tag', 'derivplotwindow');
            hold on;
            for m=1:30
                h2=plot(distance(m,:)*calxy,c(m,:), 'color',[m/30, rand(1)/4+0.5, 1-m/30]);
            end
            hold off;
            set (gca, 'xgrid', 'on', 'ygrid', 'on', 'TickDir', 'in')
            if retr('caluv')==1 && retr('calxy')==1
                xlabel('Distance on line [px]');
            else
                xlabel('Distance on line [m]');
            end
            ylabel(current);
            h3=figure;
            screensize=get( 0, 'ScreenSize' );
            rect = [screensize(3)/4-300, screensize(4)/2-250, 600, 500];
            set(h3,'position', rect);
            current=get(handles.extraction_choice,'string');
            current=current{extractwhat};
            set(h3,'numbertitle','off','menubar','none','toolbar','figure','dockcontrols','off','name',[current ', frame ' num2str(currentframe)],'tag', 'derivplotwindow');
            calxy=retr('calxy');
            %user can click on point, circle will be displayed in main window
            plot (1:30, integral);
            hold on;
            scattergroup1=scatter(1:30, integral, 80, 'ko');
            hold off;
            
            if verLessThan('matlab','8.4')
                set(scattergroup1, 'ButtonDownFcn', @hitcircle, 'hittestarea', 'off');
            else
                % >R2014a
                set(scattergroup1, 'ButtonDownFcn', @hitcircle, 'pickableparts', 'visible');
            end
            
            
            
            title('Click the points of the graph to highlight it''s corresponding circle.')
            set (gca, 'xgrid', 'on', 'ygrid', 'on', 'TickDir', 'in')
            xlabel('circle series nr. (circle with max. circulation highlighted in blue)');
            if retr('caluv')==1 && retr('calxy')==1
                ylabel('tangent velocity loop integral (circulation) [px^2/frame]');
            else
                ylabel('tangent velocity loop integral (circulation) [m^2/s]');
            end
            put('distance',distance*calxy);
            put('c',c);
            put('h3plot', h3);
            put('integral', integral);
            h_extractionplot=retr('h_extractionplot');
            h_extractionplot(size(h_extractionplot,1)+1,1)=h;
            put ('h_extractionplot', h_extractionplot);
            h_extractionplot2=retr('h_extractionplot2');
            h_extractionplot2(size(h_extractionplot2,1)+1,1)=h3;
            put ('h_extractionplot2', h_extractionplot2);
        end
    end
end

function hitcircle(src,eventdata)
posreal=get(gca,'CurrentPoint');
delete(findobj('tag','circstring'));
pos=round(posreal(1,1));
xposition=retr('xposition');
yposition=retr('yposition');
integral=retr('integral');
hgui=getappdata(0,'hgui');
h3plot=retr('h3plot');
figure(hgui);
delete(findobj('tag', 'extractline'))
for m=1:30
    line(xposition(m,:),yposition(m,:),'LineWidth',1.5, 'Color', [0.95,0.5,0.01],'tag','extractline');
end
line(xposition(pos,:),yposition(pos,:),'LineWidth',2.5, 'Color',[0.2,0.5,0.7],'tag','extractline');
figure(h3plot);
marksize=linspace(80,80,30)';
marksize(pos)=150;
set(gco, 'SizeData', marksize);
if retr('caluv')==1 && retr('calxy')==1
    units='px^2/frame';
else
    units='m^2/s';
end
text(posreal(1,1)+0.75,posreal(2,2),['\leftarrow ' num2str(integral(pos)) ' ' units],'tag','circstring','BackgroundColor', [1 1 1], 'margin', 0.01, 'fontsize', 7, 'HitTest', 'off')

% --- Executes on button press in clear_plot.
function clear_plot_Callback(hObject, eventdata, handles)
h_extractionplot=retr('h_extractionplot');
h_extractionplot2=retr('h_extractionplot2');
for i=1:size(h_extractionplot,1)
    try
        close (h_extractionplot(i));
    catch
    end
    try
        close (h_extractionplot2(i));
    catch
    end
end
put ('h_extractionplot', []);
put ('h_extractionplot2', []);
delete(findobj('tag', 'extractpoint'));
delete(findobj('tag', 'extractline'));
delete(findobj('tag', 'circstring'));

% --- Executes on button press in histdraw.
function histdraw_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    x=resultslist{1,currentframe};
    y=resultslist{2,currentframe};
    if size(resultslist,1)>6 %filtered exists
        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
            u=resultslist{10,currentframe};
            v=resultslist{11,currentframe};
        else
            u=resultslist{7,currentframe};
            if size(u,1)>1
                v=resultslist{8,currentframe};
            else
                u=resultslist{3,currentframe};
                v=resultslist{4,currentframe};
            end
        end
    else
        u=resultslist{3,currentframe};
        v=resultslist{4,currentframe};
    end
    caluv=retr('caluv');
    calxy=retr('calxy');
    x=reshape(x,size(x,1)*size(x,2),1);
    y=reshape(y,size(y,1)*size(y,2),1);
    u=reshape(u,size(u,1)*size(u,2),1);
    v=reshape(v,size(v,1)*size(v,2),1);
    choice_plot=get(handles.hist_select,'value');
    current=get(handles.hist_select,'string');
    current=current{choice_plot};
    h=figure;
    screensize=get( 0, 'ScreenSize' );
    rect = [screensize(3)/4-300, screensize(4)/2-250, 600, 500];
    set(h,'position', rect);
    set(h,'numbertitle','off','menubar','none','toolbar','figure','dockcontrols','off','name',['Histogram ' current ', frame ' num2str(currentframe)],'tag', 'derivplotwindow');
    nrofbins=str2double(get(handles.nrofbins, 'string'));
    if choice_plot==1
        [n xout]=hist(u*caluv-retr('subtr_u'),nrofbins);
        xdescript='velocity (u)';
    elseif choice_plot==2
        [n xout]=hist(v*caluv-retr('subtr_v'),nrofbins);
        xdescript='velocity (v)';
    elseif choice_plot==3
        [n xout]=hist(sqrt((u*caluv-retr('subtr_u')).^2+(v*caluv-retr('subtr_v')).^2),nrofbins);
        xdescript='velocity magnitude';
    end
    h2=bar(xout,n);
    set (gca, 'xgrid', 'on', 'ygrid', 'on', 'TickDir', 'in')
    if retr('caluv')==1 && retr('calxy')==1
        xlabel([xdescript ' [px/frame]']);
    else
        xlabel([xdescript ' [m/s]']);
    end
    ylabel('frequency');
end

% --- Executes on button press in generate_it.
function generate_it_Callback(hObject, eventdata, handles)
handles=gethand;
flow_sim=get(handles.flow_sim,'value');
switch flow_sim
    case 1 %rankine
        v0 = str2double(get(handles.rank_displ,'string')); %max velocity
        vortexplayground=[str2double(get(handles.img_sizex,'string')),str2double(get(handles.img_sizey,'string'))]; %width, height)
        center1=[str2double(get(handles.rankx1,'string')),str2double(get(handles.ranky1,'string'))]; %x,y
        center2=[str2double(get(handles.rankx2,'string')),str2double(get(handles.ranky2,'string'))]; %x,y
        [x,y]=meshgrid(-center1(1):vortexplayground(1)-center1(1)-1,-center1(2):vortexplayground(2)-center1(2)-1);
        [o,r] = cart2pol(x,y);
        uo=zeros(size(x));
        R0 = str2double(get(handles.rank_core,'string')); %radius %35
        uoin = (r <= R0);
        uout = (r > R0);
        uo = uoin+uout;
        uo(uoin) =  v0*r(uoin)/R0;
        uo(uout) =  v0*R0./r(uout);
        uo(isnan(uo))=0;
        u = -uo.*sin(o);
        v = uo.*cos(o);
        if get(handles.singledoublerankine,'value')==2
            [x,y]=meshgrid(-center2(1):vortexplayground(1)-center2(1)-1,-center2(2):vortexplayground(2)-center2(2)-1);
            [o,r] = cart2pol(x,y);
            uo=zeros(size(x));
            R0 = str2double(get(handles.rank_core,'string')); %radius %35
            uoin = (r <= R0);
            uout = (r > R0);
            uo = uoin+uout;
            uo(uoin) =  v0*r(uoin)/R0;
            uo(uout) =  v0*R0./r(uout);
            uo(isnan(uo))=0;
            u2 = -uo.*sin(o);
            v2 = uo.*cos(o);
            u=u-u2;
            v=v-v2;
        end
    case 2 %oseen
        v0 = str2double(get(handles.oseen_displ,'string'))*3; %max velocity
        vortexplayground=[str2double(get(handles.img_sizex,'string')),str2double(get(handles.img_sizey,'string'))]; %width, height)
        center1=[str2double(get(handles.oseenx1,'string')),str2double(get(handles.oseeny1,'string'))]; %x,y
        center2=[str2double(get(handles.oseenx2,'string')),str2double(get(handles.oseeny2,'string'))]; %x,y
        [x,y]=meshgrid(-center1(1):vortexplayground(1)-center1(1)-1,-center1(2):vortexplayground(2)-center1(2)-1);
        [o,r] = cart2pol(x,y);
        uo=zeros(size(x));
        zaeh=1;
        t=str2double(get(handles.oseen_time,'string'));
        r=r/100;
        
        %uo wird im zwentrum NaN!!
        uo=(v0./(2*pi*r)).*(1-exp(-r.^2/(4*zaeh*t)));
        uo(isnan(uo))=0;
        u = -uo.*sin(o);
        v = uo.*cos(o);
        if get(handles.singledoubleoseen,'value')==2
            [x,y]=meshgrid(-center2(1):vortexplayground(1)-center2(1)-1,-center2(2):vortexplayground(2)-center2(2)-1);
            [o,r] = cart2pol(x,y);
            r=r/100;
            uo=(v0./(2*pi*r)).*(1-exp(-r.^2/(4*zaeh*t)));
            uo(isnan(uo))=0;
            u2 = -uo.*sin(o);
            v2 = uo.*cos(o);
            u=u-u2;
            v=v-v2;
         end
    case 3 %linear
        u=zeros(str2double(get(handles.img_sizey,'string')),str2double(get(handles.img_sizex,'string')));
        v(1:str2double(get(handles.img_sizey,'string')),1:str2double(get(handles.img_sizex,'string')))=str2double(get(handles.shiftdisplacement,'string'));
    case 4 % rotation
        [v,u] = meshgrid(-(str2double(get(handles.img_sizex,'string')))/2:1:(str2double(get(handles.img_sizex,'string')))/2-1,-(str2double(get(handles.img_sizey,'string')))/2:1:(str2double(get(handles.img_sizey,'string')))/2-1);
        
        u=u/max(max(u));
        v=-v/max(max(v));
        u=u*str2double(get(handles.rotationdislacement,'string'));
        v=v*str2double(get(handles.rotationdislacement,'string'));
        [x,y]=meshgrid(1:1:str2double(get(handles.img_sizex,'string'))+1);
    case 5 %membrane
        [x,y]=meshgrid(linspace(-3,3,str2double(get(handles.img_sizex,'string'))),linspace(-3,3,str2double(get(handles.img_sizey,'string'))));
        u = peaks(x,y)/3;
        v = peaks(y,x)/3;
end
%% Create Particle Image
set(handles.status_creation,'string','Calculating particles...');drawnow;
i=[];
j=[];
sizey=str2double(get(handles.img_sizey,'string'));
sizex=str2double(get(handles.img_sizex,'string'));
noise=str2double(get(handles.part_noise,'string'));
A=zeros(sizey,sizex);
B=A;
partAm=str2double(get(handles.part_am,'string'));
Z=str2double(get(handles.sheetthick,'string')); %0.25 sheet thickness
dt=str2double(get(handles.part_size,'string')); %particle diameter
ddt=str2double(get(handles.part_var,'string')); %particle diameter variation

z0_pre=randn(partAm,1); %normal distributed sheet intensity
randn('state', sum(100*clock)); %#ok<*RAND>
z1_pre=randn(partAm,1); %normal distributed sheet intensity

z0=z0_pre*(str2double(get(handles.part_z,'string'))/200+0.5)+z1_pre*(1-((str2double(get(handles.part_z,'string'))/200+0.5)));
z1=z1_pre*(str2double(get(handles.part_z,'string'))/200+0.5)+z0_pre*(1-((str2double(get(handles.part_z,'string'))/200+0.5)));

%z0=abs(randn(partAm,1)); %normal distributed sheet intensity
I0=255*exp(-(Z^2./(0.125*z0.^2))); %particle intensity
I0(I0>255)=255;
I0(I0<0)=0;

I1=255*exp(-(Z^2./(0.125*z1.^2))); %particle intensity
I1(I1>255)=255;
I1(I1<0)=0;

randn('state', sum(100*clock));
d=randn(partAm,1)/2; %particle diameter distribution
d=dt+d*ddt;
d(d<0)=0;
rand('state', sum(100*clock));
x0=rand(partAm,1)*sizex;
y0=rand(partAm,1)*sizey;
rd = -8.0 ./ d.^2;
offsety=v;
offsetx=u;

xlimit1=floor(x0-d/2); %x min particle extent image1
xlimit2=ceil(x0+d/2); %x max particle extent image1
ylimit1=floor(y0-d/2); %y min particle extent image1
ylimit2=ceil(y0+d/2); %y max particle extent image1
xlimit2(xlimit2>sizex)=sizex;
xlimit1(xlimit1<1)=1;
ylimit2(ylimit2>sizey)=sizey;
ylimit1(ylimit1<1)=1;

%calculate particle extents for image2 (shifted image)
x0integer=round(x0);
x0integer(x0integer>sizex)=sizex;
x0integer(x0integer<1)=1;
y0integer=round(y0);
y0integer(y0integer>sizey)=sizey;
y0integer(y0integer<1)=1;

xlimit3=zeros(partAm,1);
xlimit4=xlimit3;
ylimit3=xlimit3;
ylimit4=xlimit3;
for n=1:partAm
    xlimit3(n,1)=floor(x0(n)-d(n)/2-offsetx((y0integer(n)),(x0integer(n)))); %x min particle extent image2
    xlimit4(n,1)=ceil(x0(n)+d(n)/2-offsetx((y0integer(n)),(x0integer(n)))); %x max particle extent image2
    ylimit3(n,1)=floor(y0(n)-d(n)/2-offsety((y0integer(n)),(x0integer(n)))); %y min particle extent image2
    ylimit4(n,1)=ceil(y0(n)+d(n)/2-offsety((y0integer(n)),(x0integer(n)))); %y max particle extent image2
end
xlimit3(xlimit3<1)=1;
xlimit4(xlimit4>sizex)=sizex;
ylimit3(ylimit3<1)=1;
ylimit4(ylimit4>sizey)=sizey;

set(handles.status_creation,'string','Placing particles...');drawnow;
for n=1:partAm
    r = rd(n);
    for j=xlimit1(n):xlimit2(n)
        rj = (j-x0(n))^2;
        for i=ylimit1(n):ylimit2(n)
            A(i,j)=A(i,j)+I0(n)*exp((rj+(i-y0(n))^2)*r);
        end
    end
    for j=xlimit3(n):xlimit4(n)
        for i=ylimit3(n):ylimit4(n)
            B(i,j)=B(i,j)+I1(n)*exp((-(j-x0(n)+offsetx(i,j))^2-(i-y0(n)+offsety(i,j))^2)*-rd(n)); %place particle with gaussian intensity profile
        end
    end
end

A(A>255)=255;
B(B>255)=255;

gen_image_1=imnoise(uint8(A),'gaussian',0,noise);
gen_image_2=imnoise(uint8(B),'gaussian',0,noise);

set(handles.status_creation,'string','...done')
figure;imshow(gen_image_1,'initialmagnification', 100);
figure;imshow(gen_image_2,'initialmagnification', 100);
put('gen_image_1',gen_image_1);
put('gen_image_2',gen_image_2);

% --- Executes on button press in save_imgs.
function save_imgs_Callback(hObject, eventdata, handles)
gen_image_1=retr('gen_image_1');
gen_image_2=retr('gen_image_2');
if isempty(gen_image_1)==0
    [FileName,PathName] = uiputfile('*.tif','Save generated images as...',['PIVlab_gen.tif']);
    if isequal(FileName,0) | isequal(PathName,0)
    else
        [Dir Name Ext] = fileparts(FileName);
        FileName_1=[Name '_01' Ext];
        FileName_2=[Name '_02' Ext];
        if exist(fullfile(PathName,FileName_1),'file') >0 || exist(fullfile(PathName,FileName_2),'file') >0
            butt = questdlg(['Warning: File ' FileName_1 ' already exists.'],'File exists','Overwrite','Cancel','Overwrite')
            if strmatch(butt, 'Overwrite') == 1
                write_it=1;
            else
                write_it=0;
            end
        else
            write_it=1;
        end
        if write_it==1
            imwrite(gen_image_1,fullfile(PathName,FileName_1),'Compression','none')
            imwrite(gen_image_2,fullfile(PathName,FileName_2),'Compression','none')
        end
    end
end

% --- Executes on button press in dummy.
function dummy_Callback(hObject, eventdata, handles)
sliderdisp

function applyskipper_Callback(hObject, eventdata, handles)
filename=retr('filename');
filepath=retr('filepath');
handles=gethand;
skipper=str2num(get(handles.skipper, 'string'))+2;
filepathnew(1,1)=filepath(1,1);
filepathnew(2,1)=filepath(2,1);
filenamenew(1,1)=filename(1,1);
filenamenew(2,1)=filename(2,1);
countr=3;
for i=skipper+1:skipper:size(filepath,1)-2
    filepathnew(countr,1)=filepath(i,1);
    filepathnew(countr+1,1)=filepath(i+1,1);
    filenamenew(countr,1)=filename(i,1);
    filenamenew(countr+1,1)=filename(i+1,1);
    countr=countr+2;
end
%user alters source -> results have to be removed.
put('maskiererx' ,[]);
put('maskierery' ,[]);
put ('derived',[]);
put ('resultslist',[]);
filename=filenamenew;
filepath=filepathnew;
put ('filename',filename); %only for displaying
put ('filepath',filepath); %full path and filename for analyses
if size(filepath,1)>2
    sliderstepcount=size(filepath,1)/2;
    set(handles.fileselector, 'enable', 'on');
    set (handles.fileselector,'value',1, 'min', 1,'max',sliderstepcount,'sliderstep', [1/(sliderstepcount-1) 1/(sliderstepcount-1)*10]);
else
    sliderstepcount=1;
    set(handles.fileselector, 'enable', 'off');
    set (handles.fileselector,'value',1, 'min', 1,'max',2,'sliderstep', [0.5 0.5]);
end
set (handles.filenamebox, 'string', filename);
set(handles.skipper, 'enable', 'off');
set(handles.applyskipper, 'enable', 'off');
sliderdisp

function pivlabhelp_Callback(hObject, eventdata, handles)
try
web('http://pivlab.blogspot.de/p/blog-page_19.html','-browser')
catch
    %why does 'web' not work in v 7.1.0.246 ...?
    disp('Ooops, MATLAB couldn''t open the website.')
    disp('You''ll have to open the website manually:')
    disp('http://pivlab.blogspot.de/p/blog-page_19.html')
end

% --------------------------------------------------------------------
function aboutpiv_Callback(hObject, eventdata, handles)
string={...
    'PIVlab - Time-Resolved Digital Particle Image Velocimetry Tool for MATLAB';...
    ['version: ' retr('PIVver')];...
    '';...
    'developed by Dr. William Thielicke and Prof. Dr. Eize J. Stamhuis';...
    'published under the BSD and CC-BY licence.';...
    '';...
    '';...
    'programmed with MATLAB Version 7.10 (R2010a)';...
    'first released March 09, 2010';...
    '';...
    'http://PIVlab.blogspot.com';...
    'contact: PIVlab@gmx.com';...
    };
helpdlg(string,'About')

% --- Executes on button press in clear_everything.
function clear_everything_Callback(hObject, eventdata, handles)
put ('resultslist', []); %clears old results
put ('derived', []);
handles=gethand;
set(handles.progress, 'String','Frame progress: N/A');
set(handles.overall, 'String','Total progress: N/A');
set(handles.totaltime, 'String','Time left: N/A');
set(handles.messagetext, 'String','');
sliderdisp

% --- Executes on button press in autoscale_vec.
function autoscale_vec_Callback(hObject, eventdata, handles)
handles=gethand;
if get(handles.autoscale_vec, 'value')==1
    set(handles.vectorscale,'enable', 'off');
else
    set(handles.vectorscale,'enable', 'on');
end

% --------------------------------------------------------------------
function save_session_Callback(hObject, eventdata, handles)
sessionpath=retr('sessionpath');
if isempty(sessionpath)
    sessionpath=retr('pathname');
end
[FileName,PathName] = uiputfile('*.mat','Save current session as...',fullfile(sessionpath,'PIVlab_session.mat'));
if isequal(FileName,0) | isequal(PathName,0)
else
    put('sessionpath',PathName );
    savesessionfuntion (PathName,FileName)
end

function savesessionfuntion (PathName,FileName)
hgui=getappdata(0,'hgui');
handles=gethand;
app=getappdata(hgui);
text(50,50,'Please wait, saving session...','color','y','fontsize',15, 'BackgroundColor', 'k','tag','savehint')
drawnow;
%Newer versions of Matlab do really funny things when the following vars are not empty...:
app.GUIDEOptions =[];
app.GUIOnScreen  =[];
app.Listeners  =[];
app.SavedVisible  =[];
app.ScribePloteditEnable  =[];
app.UsedByGUIData_m  =[];
app.lastValidTag =[];
iptPointerManager=[];
app.ZoomObject=[]; %Matlab crashes if this is not empty. Weird...
app.ZoomFigureState=[];
app.ZoomOnState=[];
app.PanFigureState=[];
app.uitools_FigureToolManager=[];

try
    iptPointerManager(gcf, 'disable')
catch
end
clear hgui iptPointerManager GUIDEOptions GUIOnScreen Listeners SavedVisible ScribePloteditEnable UsedByGUIData_m ZoomObject

save('-v6', fullfile(PathName,FileName), '-struct', 'app')
clear app %hgui iptPointerManager
clear hgui iptPointerManager GUIDEOptions GUIOnScreen Listeners SavedVisible ScribePloteditEnable UsedByGUIData_m
clahe_enable=get(handles.clahe_enable,'value');
clahe_size=get(handles.clahe_size,'string');
enable_highpass=get(handles.enable_highpass,'value');
highp_size=get(handles.highp_size,'string');
wienerwurst=get(handles.wienerwurst,'value');
wienerwurstsize=get(handles.wienerwurstsize,'string');
%enable_clip=get(handles.enable_clip,'value');
%clip_thresh=get(handles.clip_thresh,'string');
enable_intenscap=get(handles.enable_intenscap,'value');
intarea=get(handles.intarea,'string');
stepsize=get(handles.step,'string');
subpix=get(handles.subpix,'value');  %popup
stdev_check=get(handles.stdev_check,'value');
stdev_thresh=get(handles.stdev_thresh,'string');
loc_median=get(handles.loc_median,'value');
loc_med_thresh=get(handles.loc_med_thresh,'string');
epsilon=get(handles.epsilon,'string');
interpol_missing=get(handles.interpol_missing,'value');
vectorscale=get(handles.vectorscale,'string');
colormap_choice=get(handles.colormap_choice,'value'); %popup
addfileinfo=get(handles.addfileinfo,'value');
add_header=get(handles.add_header,'value');
delimiter=get(handles.delimiter,'value');%popup
img_not_mask=get(handles.img_not_mask,'value');
autoscale_vec=get(handles.autoscale_vec,'value');
calxy=retr('calxy');
caluv=retr('caluv');
pointscali=retr('pointscali');
realdist_string=get(handles.realdist, 'String');
time_inp_string=get(handles.time_inp, 'String');

imginterpol=get(handles.popupmenu16, 'value');
dccmark=get(handles.dcc, 'value');
fftmark=get(handles.fftmulti, 'value');
pass2=get(handles.checkbox26, 'value');
pass3=get(handles.checkbox27, 'value');
pass4=get(handles.checkbox28, 'value');
pass2val=get(handles.edit50, 'string');
pass3val=get(handles.edit51, 'string');
pass4val=get(handles.edit52, 'string');
step2=get(handles.text126, 'string');
step3=get(handles.text127, 'string');
step4=get(handles.text128, 'string');
holdstream=get(handles.holdstream, 'value');
streamlamount=get(handles.streamlamount, 'string');
streamlcolor=get(handles.streamlcolor, 'value');
save('-v6', fullfile(PathName,FileName), '-append');
delete(findobj('tag','savehint'));
drawnow;

% --------------------------------------------------------------------
function load_session_Callback(hObject, eventdata, handles)
sessionpath=retr('sessionpath');
if isempty(sessionpath)
    sessionpath=retr('pathname');
end
[FileName,PathName, filterindex] = uigetfile({'*.mat','MATLAB Files (*.mat)'; '*.mat','mat'},'Load PIVlab session',fullfile(sessionpath, 'PIVlab_session.mat'));
if isequal(FileName,0) | isequal(PathName,0)
else
    clear iptPointerManager
    put('sessionpath',PathName );
    put('derived',[]);
    put('resultslist',[]);
    put('maskiererx',[]);
    put('maskierery',[]);
    put('roirect',[]);
    put('velrect',[]);
    put('filename',[]);
    put('filepath',[]);
    hgui=getappdata(0,'hgui');
    warning off all
    vars=load(fullfile(PathName,FileName),'yposition', 'FileName', 'PathName', 'add_header', 'addfileinfo', 'autoscale_vec', 'caliimg', 'caluv', 'calxy', 'cancel', 'clahe_enable', 'clahe_size', 'colormap_choice', 'delimiter', 'derived', 'displaywhat', 'distance', 'enable_highpass', 'enable_intenscap', 'epsilon', 'filename', 'filepath', 'highp_size', 'homedir', 'img_not_mask', 'intarea', 'interpol_missing', 'loc_med_thresh', 'loc_median', 'manualdeletion', 'maskiererx', 'maskierery', 'pathname', 'pointscali', 'resultslist', 'roirect', 'sequencer', 'sessionpath', 'stdev_check', 'stdev_thresh', 'stepsize', 'subpix', 'subtr_u', 'subtr_v', 'toggler', 'vectorscale', 'velrect', 'wasdisabled', 'xposition','realdist_string','time_inp_string','streamlinesX','streamlinesY','manmarkersX','manmarkersY','imginterpol','dccmark','fftmark','pass2','pass3','pass4','pass2val','pass3val','pass4val','step2','step3','step4','holdstream','streamlamount','streamlcolor','ismean','wienerwurst','wienerwurstsize');
    names=fieldnames(vars);
    for i=1:size(names,1)
        setappdata(hgui,names{i},vars.(names{i}))
    end
    sliderrange
    handles=gethand;
    
    set(handles.clahe_enable,'value',retr('clahe_enable'));
    set(handles.clahe_size,'string',retr('clahe_size'));
    set(handles.enable_highpass,'value',retr('enable_highpass'));
    set(handles.highp_size,'string',retr('highp_size'));
    
 set(handles.wienerwurst,'value',retr('wienerwurst'));
 set(handles.wienerwurstsize,'string',retr('wienerwurstsize'));
    
    %set(handles.enable_clip,'value',retr('enable_clip'));
    %set(handles.clip_thresh,'string',retr('clip_thresh'));
    set(handles.enable_intenscap,'value',retr('enable_intenscap'));
    set(handles.intarea,'string',retr('intarea'));
    set(handles.step,'string',retr('stepsize'));
    set(handles.subpix,'value',retr('subpix'));  %popup
    set(handles.stdev_check,'value',retr('stdev_check'));
    set(handles.stdev_thresh,'string',retr('stdev_thresh'));
    set(handles.loc_median,'value',retr('loc_median'));
    set(handles.loc_med_thresh,'string',retr('loc_med_thresh'));
    set(handles.epsilon,'string',retr('epsilon'));
    set(handles.interpol_missing,'value',retr('interpol_missing'));
    set(handles.vectorscale,'string',retr('vectorscale'));
    set(handles.colormap_choice,'value',retr('colormap_choice')); %popup
    set(handles.addfileinfo,'value',retr('addfileinfo'));
    set(handles.add_header,'value',retr('add_header'));
    set(handles.delimiter,'value',retr('delimiter'));%popup
    set(handles.img_not_mask,'value',retr('img_not_mask'));
    set(handles.autoscale_vec,'value',retr('autoscale_vec'));
    
    set(handles.popupmenu16, 'value',vars.imginterpol);
    set(handles.dcc, 'value',vars.dccmark);
    set(handles.fftmulti, 'value',vars.fftmark);
    if vars.fftmark==1
        set (handles.uipanel42,'visible','on')
    else
        set (handles.uipanel42,'visible','off')
    end
    set(handles.checkbox26, 'value',vars.pass2);
    set(handles.checkbox27, 'value',vars.pass3);
    set(handles.checkbox28, 'value',vars.pass4);
    set(handles.edit50, 'string',vars.pass2val);
    set(handles.edit51, 'string',vars.pass3val);
    set(handles.edit52, 'string',vars.pass4val);
    set(handles.text126, 'string',vars.step2);
    set(handles.text127, 'string',vars.step3);
    set(handles.text128, 'string',vars.step4);
    set(handles.holdstream, 'value',vars.holdstream);
    set(handles.streamlamount, 'string',vars.streamlamount);
    set(handles.streamlcolor, 'value',vars.streamlcolor);
    set(handles.streamlwidth, 'value',vars.streamlcolor);
    
    try
        set(handles.realdist, 'String',vars.realdist_string);
        set(handles.time_inp, 'String',vars.time_inp_string);
        if isempty(vars.pointscali)==0
            caluv=retr('caluv');
            set(handles.calidisp, 'string', ['1 px = ' num2str(round(calxy*100000)/100000) ' m' sprintf('\n') '1 px/frame = ' num2str(round(caluv*100000)/100000) ' m/s'],  'backgroundcolor', [0.5 1 0.5]);
        end
    catch
        disp('.')
    end
    
    %reset zoom
    set(handles.panon,'Value',0);
    set(handles.zoomon,'Value',0);
    put('xzoomlimit', []);
    put('yzoomlimit', []);
    
    sliderdisp
    zoom reset
end

% --- Executes on button press in save_only_one.
function save_only_one_Callback(hObject, eventdata, handles)
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
filepath=retr('filepath');
set(handles.firstframe,'string',int2str(currentframe));
set(handles.lastframe,'string',int2str(currentframe));
put('only_single_frame',1);
saveavi_Callback

% --- Executes on button press in saveavi.
function saveavi_Callback(hObject, eventdata, handles)
only_single_frame=retr('only_single_frame');
handles=gethand;
filepath=retr('filepath');
if only_single_frame==1
    startframe=str2num(get(handles.firstframe,'string'));
    endframe=str2num(get(handles.lastframe,'string'));
    %formattype=2; %single frame --> only jpg
    %set (handles.jpgfilesave, 'value',1)
    %set (handles.avifilesave, 'value',0)
    drawnow;
else
    set(handles.fileselector, 'value',1)
    startframe=str2num(get(handles.firstframe,'string'));
    if startframe <1
        startframe=1;
    elseif startframe>size(filepath,1)/2
        startframe=size(filepath,1)/2;
    end
    set(handles.firstframe,'string',int2str(startframe));
    endframe=str2num(get(handles.lastframe,'string'));
    if endframe <startframe
        endframe=startframe;
    elseif endframe>size(filepath,1)/2
        endframe=size(filepath,1)/2;
    end
    set(handles.lastframe,'string',int2str(endframe));
end
if get (handles.avifilesave, 'value')==1
    formattype=1;
end
if get (handles.jpgfilesave, 'value')==1
    formattype=2;
end
if get (handles.bmpfilesave, 'value')==1
    formattype=3;
end
if get (handles.epsfilesave, 'value')==1
    formattype=4;
end
if get (handles.pdffilesave, 'value')==1
    formattype=5;
end

put('only_single_frame',0);

p8wasvisible=retr('p8wasvisible');
if p8wasvisible==1
    switchui('multip08');
    sliderdisp
end
imgsavepath=retr('imgsavepath');
if isempty(imgsavepath)
    imgsavepath=retr('pathname');
end

if formattype==1
    [filename, pathname] = uiputfile({ '*.avi','movie (*.avi)'}, 'Save movie as',fullfile(imgsavepath, 'PIVlab_out'));
    if isequal(filename,0) || isequal(pathname,0)
        return
    end
    put('imgsavepath',pathname );
    compr=get(handles.usecompr,'value');
    
    
    if verLessThan('matlab','8.4')
        if compr==0
            compr='none';
        else
            compr='cinepak';
        end
        aviobj = avifile(fullfile(pathname,filename),'compression',compr,'quality', 100, 'fps', str2double(get(handles.fps_setting,'string')));
    else
        if compr==0
            compr='Uncompressed AVI';
        else
            compr='Motion JPEG AVI';
        end
        aviobj = VideoWriter(fullfile(pathname,filename),compr);
        aviobj.FrameRate = str2double(get(handles.fps_setting,'string'));
        open(aviobj);
    end
    
    for i=startframe:endframe
        set(handles.fileselector, 'value',i)
        sliderdisp
        hgca=gca;
        colo=get(gcf, 'colormap');
        axes_units = get(hgca,'Units');
        axes_pos = get(hgca,'Position');
        newFig=figure('visible', 'off');
        set(newFig,'visible', 'off');
        set(newFig,'Units',axes_units);
        set(newFig,'Position',[15 5 axes_pos(3)+30 axes_pos(4)+10]);
        axesObject2=copyobj(hgca,newFig);
        set(axesObject2,'Units',axes_units);
        set(axesObject2,'Position',[15 5 axes_pos(3) axes_pos(4)]);
        colormap(colo);
        if get(handles.displ_colorbar,'value')==1
            name=get(handles.derivchoice,'string');
            posichoice = get(handles.colorbarpos,'String');
            coloobj=colorbar (posichoice{get(handles.colorbarpos,'Value')},'FontWeight','bold','Fontsize',12);
            if strcmp(posichoice{get(handles.colorbarpos,'Value')},'East')==1 | strcmp(posichoice{get(handles.colorbarpos,'Value')},'West')==1
                set(coloobj,'YTickLabel',num2str(get(coloobj,'YTick')','%5.5g'))
                ylabel(coloobj,name{retr('displaywhat')},'fontsize',9,'fontweight','normal');
            end
            if strcmp(posichoice{get(handles.colorbarpos,'Value')},'North')==1 | strcmp(posichoice{get(handles.colorbarpos,'Value')},'South')==1
                set(coloobj,'XTickLabel',num2str(get(coloobj,'XTick')','%5.5g'))
                xlabel(coloobj,name{retr('displaywhat')},'fontsize',9,'fontweight','normal');
            end
        end
        delete(findobj('tag','smoothhint'));
        F=getframe(axesObject2);
        close(newFig)
        if verLessThan('matlab','8.4')
        aviobj = addframe(aviobj,F);
        else
            writeVideo(aviobj,F);
        end
    end
    if verLessThan('matlab','8.4')
    aviobj = close(aviobj);
    else
        close(aviobj);
    end
elseif formattype ==2 || formattype==3 || formattype==4 || formattype==5
    if formattype==2
        [filename, pathname] = uiputfile({ '*.jpg','images (*.jpg)'}, 'Save images as',fullfile(imgsavepath, 'PIVlab_out'));
    elseif formattype==3
        [filename, pathname] = uiputfile({ '*.bmp','images (*.bmp)'}, 'Save images as',fullfile(imgsavepath, 'PIVlab_out'));
    elseif formattype==4
        [filename, pathname] = uiputfile({ '*.eps','images (*.eps)'}, 'Save images as',fullfile(imgsavepath, 'PIVlab_out'));
    elseif formattype==5
        [filename, pathname] = uiputfile({ '*.pdf','PostScript (*.pdf)'}, 'Save images as',fullfile(imgsavepath, 'PIVlab_out'));
    end
    if isequal(filename,0) || isequal(pathname,0)
        return
    end
    put('imgsavepath',pathname );
    if formattype==2 || formattype==3
        reso=inputdlg(['Please enter scale factor' sprintf('\n') '(1 = render image at same size as currently displayed)'],'Specify resolution',1,{'1'});
        [reso status] = str2num(reso{1});  % Use curly bracket for subscript
        if ~status
            reso=1;
        end
    end
    
    for i=startframe:endframe
        set(handles.fileselector, 'value',i)
        sliderdisp
        hgca=gca;
        colo=get(gcf, 'colormap');
        axes_units = get(hgca,'Units');
        axes_pos = get(hgca,'Position');
        aspect=axes_pos(3)/axes_pos(4);
        newFig=figure;
        set(newFig,'visible', 'off');
        set(newFig,'Units',axes_units);
        set(newFig,'Position',[15 5 axes_pos(3)+30 axes_pos(4)+10]);
        axesObject2=copyobj(hgca,newFig);
        set(axesObject2,'Units',axes_units);
        set(axesObject2,'Position',[15 5 axes_pos(3) axes_pos(4)]);
        colormap(colo);
        if get(handles.displ_colorbar,'value')==1
            name=get(handles.derivchoice,'string');
            posichoice = get(handles.colorbarpos,'String');
            coloobj=colorbar (posichoice{get(handles.colorbarpos,'Value')},'FontWeight','bold','Fontsize',12);
            if strcmp(posichoice{get(handles.colorbarpos,'Value')},'East')==1 | strcmp(posichoice{get(handles.colorbarpos,'Value')},'West')==1
                set(coloobj,'YTickLabel',num2str(get(coloobj,'YTick')','%5.5g'))
                ylabel(coloobj,name{retr('displaywhat')},'fontsize',9,'fontweight','normal');
            end
            if strcmp(posichoice{get(handles.colorbarpos,'Value')},'North')==1 | strcmp(posichoice{get(handles.colorbarpos,'Value')},'South')==1
                set(coloobj,'XTickLabel',num2str(get(coloobj,'XTick')','%5.5g'))
                xlabel(coloobj,name{retr('displaywhat')},'fontsize',9,'fontweight','normal');
            end
        end
        delete(findobj('tag','smoothhint'));
        [Dir Name Ext] = fileparts(filename);
        newfilename=[Name sprintf('_%03d',i) Ext];
        drawnow
        if formattype==2 || formattype==3
            exportfig(newFig,fullfile(pathname,newfilename),'height',3,'color','rgb','format','bmp','resolution',96*reso,'FontMode','scaled','FontSizeMin',16);
        elseif formattype ==4
            exportfig(newFig,fullfile(pathname,newfilename),'preview','tiff','color','rgb','linemode','scaled','FontMode','scaled','FontSizeMin',16);
        elseif formattype==5
            exportfig(newFig,fullfile(pathname,newfilename),'format','pdf','preview','tiff','color','rgb','linemode','scaled','FontMode','scaled','FontSizeMin',16);
        end
        close(newFig)
        if formattype==2
            autocrop(fullfile(pathname,newfilename),1);
        elseif formattype==3
            autocrop(fullfile(pathname,newfilename),0);
        end
    end
end

% --- Executes on selection change in extraction_choice.
function extraction_choice_Callback(hObject, eventdata, handles)
if get(hObject, 'value') ~= 9
    handles=gethand;
    if get(handles.draw_what, 'value')==3
        set(handles.draw_what, 'value', 1)
    end
end

% --- Executes on selection change in draw_what.
function draw_what_Callback(hObject, eventdata, handles)
if get(hObject, 'value') == 3
    handles=gethand;
    set (handles.extraction_choice, 'value', 9);
    set (handles.extraction_choice, 'enable', 'off');
else
    set (handles.extraction_choice, 'enable', 'on');
end

function check_comma(who)
boxcontent=get(who,'String');% returns contents of time_inp as text
s = regexprep(boxcontent, ',', '.');
set(who,'String',s);


%__________________________________________________________________________
%unused callbacks:
function slider1_Callback(hObject, eventdata, handles)
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%__________________________________________________________________________

% --- Executes on button press in vorticity_roi.
function vorticity_roi_Callback(hObject, eventdata, handles)
resultslist=retr('resultslist');
%currentframe=2*floor(get(handles.fileselector, 'value'))-1;
currentframe=floor(get(handles.fileselector, 'value'));

if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    x=resultslist{1,currentframe};
    y=resultslist{2,currentframe};
    if size(resultslist,1)>6 %filtered exists
        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
            u=resultslist{10,currentframe};
            v=resultslist{11,currentframe};
        else
            u=resultslist{7,currentframe};
            if size(u,1)>1
                v=resultslist{8,currentframe};
            else
                u=resultslist{3,currentframe};
                v=resultslist{4,currentframe};
            end
        end
    else
        u=resultslist{3,currentframe};
        v=resultslist{4,currentframe};
    end
    caluv=retr('caluv');
    u=u*caluv-retr('subtr_u');
    v=v*caluv-retr('subtr_v');
    calxy=retr('calxy');
    derivative_calc(currentframe,2,1);
    derived=retr('derived');
    %currentimage=derived{1,(currentframe+1)/2};
    currentimage=derived{1,(currentframe)};
    delete(findobj('tag','vortarea'));
    h=figure;
    set (h,'DockControls', 'off', 'menubar', 'none', 'tag', 'vortarea')
    imagesc(currentimage);
    axis image
    hold on;
    quiver(u,v,'linewidth',str2double(get(handles.vecwidth,'string')))
    hold off;
    
    %draw ellipse
    for i=1:5
        [xellip(i),yellip(i),but] = ginput(1);
        if but~=1
            break
        end
        hold on;
        plot (xellip(i),yellip(i),'w*')
        hold off;
        if i==3
            line(xellip(2:3),yellip(2:3))
        end
        if i==5
            line(xellip(4:5),yellip(4:5))
        end
    end
    %click1=centre of vortical structure
    %click2=top of vortical structure
    %click3=bottom of vortical structure
    %click4=left of vortical structure
    %click5=right of vortical structure
    x0=(mean(xellip)+xellip(1))/2;
    y0=(mean(yellip)+yellip(1))/2;
    if xellip(2)<xellip(3)
        ang=acos((yellip(2)-yellip(3))/(sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)))-pi/2;
    else
        ang=asin((yellip(2)-yellip(3))/(sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)));
    end
    rb=sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)/2;
    ra=sqrt((xellip(4)-xellip(5))^2+(yellip(4)-yellip(5))^2)/2;
    text(xellip(1),yellip(1),int2str(ang/pi*180));
    ra=sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)/2;
    rb=sqrt((xellip(4)-xellip(5))^2+(yellip(4)-yellip(5))^2)/2;
    
    celllength=(x(1,2)-x(1,1))*calxy; %size of one cell
    cellarea=celllength^2; %area of one cell
    integralindex=0;
    for incr = -(ra+rb)/3 :0.5: (ra+rb)/2
        integralindex=integralindex+1;
        ra_new=ra+incr;
        if ra_new<0
            ra_new=0
        end
        if rb_new<0
            rb_new=0
        end
        rb_new=rb+incr;
        [outputx, outputy]=ellipse(ra_new,rb_new,ang,x0,y0);
        BW = roipoly(u,outputx,outputy);
        integral=0;
        for i=1:size(u,1)
            for j=1:size(u,2)
                if BW(i,j)==1
                    integral=integral+cellarea*currentimage(i,j);
                end
            end
        end
        integralseries(integralindex)=integral;
    end
    h2=figure;
    set(h2, 'tag', 'vortarea');
    plot(integralseries)
end

% --- Executes on button press in draw_area.
function draw_area_Callback(hObject, eventdata, handles)
%noch probleme wenn erster frame leer...
%dann geht er sofort zu datei asuwahl...
handles=gethand;
currentframe=floor(get(handles.fileselector, 'value'));
resultslist=retr('resultslist');

%NEU
if get(handles.extractareaall, 'value')==0
    startfr=currentframe;
    endfr=currentframe;
else
    %sollte erstes element sein mit inhalt...
    for findcontent=size(resultslist,2):-1:1
        if numel(resultslist{1,findcontent}) > 0
            startfr=findcontent;
        end
    end
    
    endfr=size(resultslist,2);
end
selected=0;
areaoperation=get(handles.areatype, 'value');
toolsavailable(0)
for i=startfr:endfr
    set(handles.fileselector, 'value',i)
    %sliderdisp
    currentframe=floor(get(handles.fileselector, 'value'));
    
    if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
        if areaoperation==1
            %area mean value
            sliderdisp
            filepath=retr('filepath');
            x=resultslist{1,currentframe};
            extractwhat=get(handles.area_para_select,'Value');
            derivative_calc(currentframe,extractwhat+1,0);
            derived=retr('derived');
            currentimage=imread(filepath{2*currentframe-1});
            sizeold=size(currentimage,1);
            sizenew=size(x,1);
            maptoget=derived{extractwhat,currentframe};
            maptoget=rescale_maps_nan(maptoget);
            if selected==0
                [BW,ximask,yimask]=roipoly;
            end
            if isempty(BW)
            else
                delete(findobj('tag','areaint'));
                delete(findobj('tag', 'extractline'))
                delete(findobj('tag', 'extractpoint'))
                numcells=0;
                summe=0;
                for i=1:size(BW,1)
                    for j=1:size(BW,2)
                        if BW(i,j)==1
                            if isnan(maptoget(i,j))==0
                                summe=summe+maptoget(i,j);
                                numcells=numcells+1;
                            end
                        end
                    end
                end
                average=summe/numcells;
                hold on;
                plot(ximask,yimask,'LineWidth',3, 'Color', [0,0.95,0],'tag','areaint');
                plot(ximask,yimask,'LineWidth',1, 'Color', [0.95,0.5,0.01],'tag','areaint');
                hold off;
                %get units
                unitpar=get(handles.area_para_select,'string');
                unitpar=unitpar{get(handles.area_para_select,'value')};
                unitpar=unitpar(strfind(unitpar,'[')+1:end-1);
                
                
                text(min(ximask),mean(yimask), ['area mean value = ' num2str(average) ' [' unitpar ']'], 'BackgroundColor', 'w','tag','areaint');
                areaoutput=average;
                varis='[mean]';
            end
        elseif areaoperation==2
            %area integral
            sliderdisp
            filepath=retr('filepath');
            x=resultslist{1,currentframe};
            extractwhat=get(handles.area_para_select,'Value');
            derivative_calc(currentframe,extractwhat+1,0);
            derived=retr('derived');
            maptoget=derived{extractwhat,currentframe};
            maptoget=rescale_maps_nan(maptoget);
            
            calxy=retr('calxy');
            currentimage=imread(filepath{2*currentframe-1});
            sizeold=size(currentimage,1);
            sizenew=size(x,1);
            if selected==0
                [BW,ximask,yimask]=roipoly; %select in currently displayed image
            end
            if isempty(BW)
            else
                delete(findobj('tag','areaint'));
                delete(findobj('tag', 'extractline'))
                delete(findobj('tag', 'extractpoint'))
                celllength=1*calxy; %size of one pixel
                cellarea=celllength^2; %area of one cell
                integral=0;
                for i=1:size(BW,1)
                    for j=1:size(BW,2)
                        if BW(i,j)==1
                            if isnan(maptoget(i,j))==0 %do not include nans and nan area in integral.
                                integral=integral+cellarea*maptoget(i,j);
                            end
                        end
                    end
                end
                hold on;
                plot(ximask,yimask,'LineWidth',3, 'Color', [0,0.95,0],'tag','areaint');
                plot(ximask,yimask,'LineWidth',1, 'Color', [0.95,0.5,0.01],'tag','areaint');
                hold off;
                
                %get units
                if retr('caluv')==1 && retr('calxy')==1
                    distunit='px^2';
                else
                    distunit='m^2';
                end
                
                unitpar=get(handles.area_para_select,'string');
                unitpar=unitpar{get(handles.area_para_select,'value')};
                unitpar=unitpar(strfind(unitpar,'[')+1:end-1);
                
                
                text(min(ximask),mean(yimask), ['area integral = ' num2str(integral) ' [' unitpar '*' distunit ']'], 'BackgroundColor', 'w','tag','areaint');
                areaoutput=integral;
                varis='[integral]';
            end
        elseif areaoperation==3
            % area only
            sliderdisp
            filepath=retr('filepath');
            currentimage=imread(filepath{2*currentframe-1});
            x=resultslist{1,currentframe};
            sizeold=size(currentimage,1);
            sizenew=size(x,1);
            if selected==0
                [BW,ximask,yimask]=roipoly;
            end
            if isempty(BW)
            else
                delete(findobj('tag','areaint'));
                delete(findobj('tag', 'extractline'))
                delete(findobj('tag', 'extractpoint'))
                calxy=retr('calxy');
                celllength=1*calxy;
                cellarea=celllength^2;
                summe=0;
                for i=1:size(BW,1)
                    for j=1:size(BW,2)
                        if BW(i,j)==1
                            summe=summe+cellarea;
                        end
                    end
                end
                hold on;
                plot(ximask,yimask,'LineWidth',3, 'Color', [0,0.95,0],'tag','areaint');
                plot(ximask,yimask,'LineWidth',1, 'Color', [0.95,0.5,0.01],'tag','areaint');
                hold off;
                
                %get units
                if retr('caluv')==1 && retr('calxy')==1
                    distunit='px^2';
                else
                    distunit='m^2';
                end
                
                
                text(min(ximask),mean(yimask), ['area = ' num2str(summe) ' [' distunit ']'], 'BackgroundColor', 'w','tag','areaint');
                areaoutput=summe;
                varis='[area]';
            end
        elseif areaoperation==4
            %area series
            if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
                x=resultslist{1,currentframe};
                y=resultslist{2,currentframe};
                if size(resultslist,1)>6 %filtered exists
                    if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
                        u=resultslist{10,currentframe};
                        v=resultslist{11,currentframe};
                    else
                        u=resultslist{7,currentframe};
                        if size(u,1)>1
                            v=resultslist{8,currentframe};
                        else
                            u=resultslist{3,currentframe};
                            v=resultslist{4,currentframe};
                        end
                    end
                else
                    u=resultslist{3,currentframe};
                    v=resultslist{4,currentframe};
                end
                caluv=retr('caluv');
                u=u*caluv-retr('subtr_u');
                v=v*caluv-retr('subtr_v');
                calxy=retr('calxy');
                
                extractwhat=get(handles.area_para_select,'Value');
                derivative_calc(currentframe,extractwhat+1,0);
                derived=retr('derived');
                currentimage=derived{extractwhat,currentframe};
                
                
                currentimage=rescale_maps_nan(currentimage);
                hgui=getappdata(0,'hgui');
                figure(hgui);
                sliderdisp
                
                delete(findobj('tag','vortarea'));
                
                %draw ellipse
                if selected==0
                    for i=1:5
                        [xellip(i),yellip(i),but] = ginput(1);
                        if but~=1
                            break
                        end
                        hold on;
                        plot (xellip(i),yellip(i),'w*')
                        hold off;
                        if i==3
                            line(xellip(2:3),yellip(2:3))
                        end
                        if i==5
                            line(xellip(4:5),yellip(4:5))
                        end
                    end
                end
                if size(xellip,2)==5
                    %click1=centre of vortical structure
                    %click2=top of vortical structure
                    %click3=bottom of vortical structure
                    %click4=left of vortical structure
                    %click5=right of vortical structure
                    x0=(mean(xellip)+xellip(1))/2;
                    y0=(mean(yellip)+yellip(1))/2;
                    if xellip(2)<xellip(3)
                        ang=acos((yellip(2)-yellip(3))/(sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)))-pi/2;
                    else
                        ang=asin((yellip(2)-yellip(3))/(sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)));
                    end
                    rb=sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)/2;
                    ra=sqrt((xellip(4)-xellip(5))^2+(yellip(4)-yellip(5))^2)/2;
                    ra=sqrt((xellip(2)-xellip(3))^2+(yellip(2)-yellip(3))^2)/2;
                    rb=sqrt((xellip(4)-xellip(5))^2+(yellip(4)-yellip(5))^2)/2;
                    
                    celllength=1*calxy;
                    %celllength=(x(1,2)-x(1,1))*calxy; %size of one cell
                    cellarea=celllength^2; %area of one cell
                    integralindex=0;
                    
                    if get(handles.usethreshold,'value')==1
                        %sign=currentimage(round(yellip(1)),round(xellip(1)));
                        condition=get(handles.smallerlarger, 'value'); %1 is larger, 2 is smaller
                        thresholdareavalue=str2num(get(handles.thresholdarea, 'string'));
                        
                        if condition==1
                            currentimage(currentimage>thresholdareavalue)=nan;
                        else
                            currentimage(currentimage<thresholdareavalue)=nan;
                        end
                        %{
                    %redraw map to show excluded areas
                    [xhelper,yhelper]=meshgrid(1:size(u,2),1:size(u,1));
                    areaincluded=ones(size(u));
                    areaincluded(isnan(currentimage)==1)=0;
                    imagesc(currentimage);
                    axis image
                    hold on;
                    quiver(xhelper(areaincluded==1),yhelper(areaincluded==1),u(areaincluded==1),v(areaincluded==1),'k','linewidth',str2double(get(handles.vecwidth,'string')))
                    scatter(xhelper(areaincluded==0),yhelper(areaincluded==0),'rx')
                    hold off;
                        %}
                    end
                    increasefactor=str2num(get(handles.radiusincrease,'string'))/100;
                    if ra<rb
                        minimumrad=ra;
                    else
                        minimumrad=rb;
                    end
                    %for incr = -(minimumrad)/1.5 :0.5: (ra+rb)/2*increasefactor
                    for incr = -(minimumrad)/1.5 :5: (ra+rb)/2*increasefactor
                        integralindex=integralindex+1;
                        [outputx, outputy]=ellipse(ra+incr,rb+incr,ang,x0,y0,'w');
                        %BW = roipoly(u,outputx,outputy);
                        BW = roipoly(currentimage,outputx,outputy);
                        ra_all(integralindex)=ra+incr;
                        rb_all(integralindex)=rb+incr;
                        
                        integral=0;
                        %for i=1:size(u,1)
                        for i=1:size(currentimage,1)
                            %for j=1:size(u,2)
                            for j=1:size(currentimage,2)
                                if BW(i,j)==1
                                    if isnan(currentimage(i,j))==0
                                        integral=integral+cellarea*currentimage(i,j);
                                    end
                                end
                            end
                        end
                        integralseries(integralindex)=integral;
                    end
                    put('ra',ra_all);
                    put('rb',rb_all)
                    put('ang',ang)
                    put('x0',x0)
                    put('y0',y0)
                    h2=figure;
                    %plot(integralseries)
                    set(h2, 'tag', 'vortarea');
                    
                    plot (1:size(integralseries,2), integralseries);
                    hold on;
                    scattergroup1=scatter(1:size(integralseries,2), integralseries, 80, 'ko');
                    hold off;
                    if verLessThan('matlab','8.4')
                        set(scattergroup1, 'ButtonDownFcn', @hitcircle2, 'hittestarea', 'off');
                    else
                        % >R2014a
                        set(scattergroup1, 'ButtonDownFcn', @hitcircle2, 'pickableparts', 'visible');
                    end
                    
                    title('Click the points of the graph to highlight it''s corresponding circle.')
                    put('integralseries',integralseries);
                    put ('hellipse',h2);
                    screensize=get( 0, 'ScreenSize' );
                    rect = [screensize(3)/4-300, screensize(4)/2-250, 600, 500];
                    set(h2,'position', rect);
                    
                    extractwhat=get(handles.area_para_select,'Value');
                    current=get(handles.area_para_select,'string');
                    current=current{extractwhat};
                    set(h2,'numbertitle','off','menubar','none','toolbar','figure','dockcontrols','off','name',[current ' area integral series, frame ' num2str(currentframe)]);
                    set (gca, 'xgrid', 'on', 'ygrid', 'on', 'TickDir', 'in')
                    xlabel('Ellipse series nr.');
                    
                    if retr('caluv')==1 && retr('calxy')==1
                        units='px^2';
                    else
                        units='m^2';
                    end
                    
                    current_2=current(1:strfind(current, '[')-1);
                    current_3=current(strfind(current, '[')+1:end-1);
                    
                    
                    ylabel([current_2 ' area integral [' current_3 '*' units ']']);
                    areaoutput=integralseries;
                    varis='[integral, starting at ellipse with smallest radius]';
                end
            end
        elseif areaoperation==5
            %weighted centroid
            if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
                x=resultslist{1,currentframe};
                y=resultslist{2,currentframe};
                if size(resultslist,1)>6 %filtered exists
                    if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
                        u=resultslist{10,currentframe};
                        v=resultslist{11,currentframe};
                    else
                        u=resultslist{7,currentframe};
                        if size(u,1)>1
                            v=resultslist{8,currentframe};
                        else
                            u=resultslist{3,currentframe};
                            v=resultslist{4,currentframe};
                        end
                    end
                else
                    u=resultslist{3,currentframe};
                    v=resultslist{4,currentframe};
                end
                caluv=retr('caluv');
                u=u*caluv-retr('subtr_u');
                v=v*caluv-retr('subtr_v');
                calxy=retr('calxy');
                extractwhat=get(handles.area_para_select,'Value');
                derivative_calc(currentframe,extractwhat+1,0);
                derived=retr('derived');
                currentimage=derived{extractwhat,currentframe};
                delete(findobj('tag','vortarea'));
                
                imagesc(currentimage);
                axis image
                hold on;
                quiver(u,v,'k','linewidth',str2double(get(handles.vecwidth,'string')))
                hold off;
                
                avail_maps=get(handles.colormap_choice,'string');
                selected_index=get(handles.colormap_choice,'value');
                if selected_index == 4 %HochschuleBremen map
                    load hsbmap.mat;
                    colormap(hsb);
                elseif selected_index== 1 %rainbow
                    load rainbow.mat;
                    colormap (rainbow);
                else
                    colormap(avail_maps{selected_index});
                end
                if selected==0
                    [BW,ximask,yimask]=roipoly;
                end
                if isempty(BW)
                else
                    
                    delete(findobj('tag', 'extractline'));
                    line(ximask,yimask,'tag', 'extractline');
                    [rows,cols] = size(currentimage);
                    
                    x = ones(rows,1)*[1:cols];
                    y = [1:rows]'*ones(1,cols);
                    area=0;
                    meanx=0;
                    meany=0;
                    for i=1:size(currentimage,1)
                        for j=1:size(currentimage,2)
                            if BW(i,j)==1
                                area=area+double(currentimage(i,j));%sum image intesity
                                meanx=meanx+x(i,j)*double(currentimage(i,j));%sum position*intensity
                                meany=meany+y(i,j)*double(currentimage(i,j));
                            end
                        end
                    end
                    meanx=meanx/area;%*(sizeold/sizenew)
                    meany=meany/area;%*(sizeold/sizenew)
                    hold on; plot(meanx,meany,'w*','markersize',20,'tag', 'extractline');hold off;
                    xecht=resultslist{1,currentframe};
                    yecht=resultslist{2,currentframe};
                    step=(xecht(1,2)-xecht(1,1))*calxy;
                    %+x(1,1)
                    areaoutput=[xecht(1,1)*calxy+(meanx-1)*step yecht(1,1)*calxy+(meany-1)*step];
                    
                    if retr('caluv')==1 && retr('calxy')==1
                        un=' px';
                    else
                        un=' m';
                    end
                    
                    text(x(1,1), y(1,1), ['x =' num2str(xecht(1,1)*calxy+(meanx-1)*step) un sprintf('\n') 'y =' num2str(yecht(1,1)*calxy+(meany-1)*step) un], 'margin', 0.01, 'fontsize', 10, 'color','w','fontweight','bold','BackgroundColor', [0 0 0],'verticalalignment','top','horizontalalignment','left');
                    
                    varis='[x coordinate, y coordinate]';
                end
            end
        elseif areaoperation==6
            %mean flow direction
            if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
                x=resultslist{1,currentframe};
                y=resultslist{2,currentframe};
                if size(resultslist,1)>6 %filtered exists
                    if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
                        u=resultslist{10,currentframe};
                        v=resultslist{11,currentframe};
                    else
                        u=resultslist{7,currentframe};
                        if size(u,1)>1
                            v=resultslist{8,currentframe};
                        else
                            u=resultslist{3,currentframe};
                            v=resultslist{4,currentframe};
                        end
                    end
                else
                    u=resultslist{3,currentframe};
                    v=resultslist{4,currentframe};
                end
                sliderdisp
                caluv=retr('caluv');
                u=u*caluv-retr('subtr_u');
                v=v*caluv-retr('subtr_v');
                calxy=retr('calxy');
                delete(findobj('tag','vortarea'));
                filepath=retr('filepath');
                x=resultslist{1,currentframe};
                y=resultslist{2,currentframe};
                if selected==0
                    [BW,ximask,yimask]=roipoly;
                end
                if isempty(BW)
                else
                    delete(findobj('tag', 'extractline'));
                    line(ximask,yimask,'tag', 'extractline');
                    umean=0;
                    vmean=0;
                    uamount=0;
                    u=rescale_maps_nan(u);
                    v=rescale_maps_nan(v);
                    for i=1:size(u,1)
                        for j=1:size(u,2)
                            if BW(i,j)==1
                                if isnan(u(i,j))==0 && isnan(v(i,j))==0
                                    umean=umean+u(i,j);
                                    vmean=vmean+v(i,j);
                                    uamount=uamount+1;
                                end
                            end
                        end
                    end
                    umean=umean/uamount;
                    vmean=vmean/uamount;
                    veclength=(x(1,2)-x(1,1))*6;
                    if vmean<=0
                        angle=-atan2(vmean,umean)*180/pi;
                    else
                        angle=360-atan2(vmean,umean)*180/pi;
                    end
                    magg=sqrt(umean^2+vmean^2);
                    areaoutput=[magg angle];
                    varis='[magnitude, angle in degrees, 0 = right, 90 = up, 180 = left, 270 = down, 360 = right]';
                    hold on;quiver(mean2(ximask), mean2(yimask), umean/sqrt(umean^2+vmean^2)*veclength,vmean/sqrt(umean^2+vmean^2)*veclength,'r','autoscale','off', 'autoscalefactor', 100, 'linewidth',2,'MaxHeadSize',3,'tag', 'extractline');hold off;
                    
                    if retr('caluv')==1 && retr('calxy')==1
                        un=' px/frame';
                    else
                        un=' m/s';
                    end
                    
                    
                    text(x(1,1), y(1,1), ['angle=' num2str(angle) '°' sprintf('\n') 'magnitude=' num2str(magg) un], 'margin', 0.01, 'fontsize', 10, 'color','w','fontweight','bold','BackgroundColor', [0 0 0],'verticalalignment','top','horizontalalignment','left');
                end
            end
        end %areaoperation
    end
    if get(handles.savearea,'Value')==1
        %nur wenn man es auch speichern will...
        if selected==0
            switch areaoperation
                case 1
                    whatoperation = 'mean_value';
                case 2
                    whatoperation = 'integral';
                case 3
                    whatoperation = 'area';
                case 4
                    whatoperation = 'integral_series';
                case 5
                    whatoperation = 'weighted centroid';
                case 6
                    whatoperation = 'mean_flow';
            end
            par = get(handles.area_para_select,'string');
            par=par{get(handles.area_para_select,'Value')};
            if areaoperation==3 || areaoperation==6
                par=[];
            end
            imgsavepath=retr('imgsavepath');
            if isempty(imgsavepath)
                imgsavepath=retr('pathname');
            end
            
            if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
                
                part1= par(1:strfind(par,'/')-1) ;
                part2= par(strfind(par,'/')+1:end);
                if isempty(part1)==1
                    parED=par;
                else
                    parED=[part1 ' per ' part2];
                end
                
                [FileName,PathName] = uiputfile('*.txt','Save extracted data as...',fullfile(imgsavepath,['PIVlab_Extr_' whatoperation '_' parED '.txt'])); %framenummer in dateiname
                selected=1;
                if isequal(FileName,0) | isequal(PathName,0)
                    break
                else
                    put ('imgsavepath',PathName);
                    fid = fopen(fullfile(PathName,FileName), 'w');
                    fprintf(fid, ['Frame Nr.,' par ': ' whatoperation ' ' varis '\r\n']);
                    fclose(fid);
                end
            end
        end
        if isequal(FileName,0) | isequal(PathName,0)
        else
            if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
                dlmwrite(fullfile(PathName,FileName), [currentframe areaoutput], '-append', 'delimiter', ',', 'precision', 10, 'newline', 'pc');
            end
        end
    end
    %areaoutput
end

toolsavailable(1)


function hitcircle2(src,eventdata)
posreal=get(gca,'CurrentPoint');
delete(findobj('tag','circstring'));
pos=round(posreal(1,1));
xposition=retr('xposition');
yposition=retr('yposition');
integralseries=retr('integralseries');
hgui=getappdata(0,'hgui');
h3plot=retr('hellipse');
figure(hgui);
delete(findobj('type', 'line', 'color', 'w')) %delete white ellipses
ra=retr('ra');
rb=retr('rb');
ang=retr('ang');
x0= retr('x0');
y0=retr('y0');

for m=1:size(ra,2)
    ellipse(ra(1,m),rb(1,m),ang,x0,y0,'w');
end
ellipse(ra(1,pos),rb(1,pos),ang,x0,y0,'b');
figure(h3plot);
marksize=linspace(80,80,size(ra,2))';
marksize(pos)=150;
set(gco, 'SizeData', marksize);
%units
handles=gethand;
extractwhat=get(handles.area_para_select,'Value');
current=get(handles.area_para_select,'string');
current=current{extractwhat};
if retr('caluv')==1 && retr('calxy')==1
    units='px^2';
else
    units='m^2';
end
current_3=current(strfind(current, '[')+1:end-1);
text(posreal(1,1)+0.25,posreal(2,2),['\leftarrow ' num2str(integralseries(pos)) ' ' current_3 '*' units],'tag','circstring','BackgroundColor', [1 1 1], 'margin', 0.01, 'fontsize', 7, 'HitTest', 'off')




% --- Executes on selection change in areatype.
function areatype_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'value')==4
    set(handles.text93, 'visible', 'on')
    set(handles.smallerlarger, 'visible', 'on')
    set(handles.text94, 'visible', 'on')
    set(handles.radiusincrease, 'visible', 'on')
    set(handles.thresholdarea, 'visible', 'on')
    set(handles.usethreshold, 'visible', 'on')
    set(handles.text95, 'visible', 'on')
else
    set(handles.text93, 'visible', 'off')
    set(handles.smallerlarger, 'visible', 'off')
    set(handles.text94, 'visible', 'off')
    set(handles.radiusincrease, 'visible', 'off')
    set(handles.thresholdarea, 'visible', 'off')
    set(handles.usethreshold, 'visible', 'off')
    set(handles.text95, 'visible', 'off')
end
if get(hObject,'value')==3 || get(hObject,'value')==6
    set(handles.area_para_select,'visible','off');
    set(handles.text89,'visible','off');
else
    set(handles.area_para_select,'visible','on');
    set(handles.text89,'visible','on');
end


% --- Executes on selection change in flow_sim.
function flow_sim_Callback(hObject, eventdata, handles)
handles=gethand;
contents = get(hObject,'value');
set(handles.rankinepanel,'visible','off');
set(handles.shiftpanel,'visible','off');
set(handles.rotationpanel,'visible','off');
set(handles.oseenpanel,'visible','off');
if contents==1
    set(handles.rankinepanel,'visible','on');
elseif contents==2
    set(handles.oseenpanel,'visible','on');
elseif contents==3
    set(handles.shiftpanel,'visible','on');
elseif contents==4
    set(handles.rotationpanel,'visible','on');
end


% --- Executes on selection change in singledoublerankine.
function singledoublerankine_Callback(hObject, eventdata, handles)
handles=gethand;
contents = get(hObject,'value');
set(handles.rankx1,'visible','off');
set(handles.rankx2,'visible','off');
set(handles.ranky1,'visible','off');
set(handles.ranky2,'visible','off');
set(handles.text102,'visible','off');
set(handles.text103,'visible','off');
set(handles.text104,'visible','off');
if contents==1
    set(handles.rankx1,'visible','on');
    set(handles.ranky1,'visible','on');
elseif contents==2
    set(handles.rankx1,'visible','on');
    set(handles.ranky1,'visible','on');
    set(handles.rankx2,'visible','on');
    set(handles.ranky2,'visible','on');
    set(handles.text102,'visible','on');
    set(handles.text103,'visible','on');
    set(handles.text104,'visible','on');
end

% --- Executes on selection change in singledoubleoseen.
function singledoubleoseen_Callback(hObject, eventdata, handles)
handles=gethand;
contents = get(hObject,'value');
set(handles.oseenx1,'visible','off');
set(handles.oseenx2,'visible','off');
set(handles.oseeny1,'visible','off');
set(handles.oseeny2,'visible','off');
set(handles.text110,'visible','off');
set(handles.text111,'visible','off');
set(handles.text112,'visible','off');
if contents==1
    set(handles.oseenx1,'visible','on');
    set(handles.oseeny1,'visible','on');
elseif contents==2
    set(handles.oseenx1,'visible','on');
    set(handles.oseeny1,'visible','on');
    set(handles.oseenx2,'visible','on');
    set(handles.oseeny2,'visible','on');
    set(handles.text110,'visible','on');
    set(handles.text111,'visible','on');
    set(handles.text112,'visible','on');
end

% --- Executes on button press in meanmaker.
function meanmaker_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
resultslist=retr('resultslist');
if isempty(resultslist)==0
    if size(filepath,1)>0
        sizeerror=0;
        typevectormittel=ones(size(resultslist{1,1}));
        ismean=retr('ismean');
        if isempty(ismean)==1
            ismean=zeros(size(resultslist,2),1);
        end
        
        %man könnte erstmal alle frames stapeln, und danach nur die verwenden deren
        %nummer gewählt war.
        %for count=1:size(resultslist,2)
        %eval('A(1,[1,6,3])')
        %in textbox eintragen 1,6,3
        %oder 1:5
        %oder 1:end
        %oder 1:3,8:10
        str = strrep(get(handles.selectedFramesMean,'string'),'-',':');
        endinside=findstr(str, 'end');
        if isempty(endinside)==0
            str = strrep(get(handles.selectedFramesMean,'string'),'end',num2str(max(find(ismean==0))));
        end
        selectionok=1;
        
        strnum=str2num(str);
        if isempty(strnum)==1 || isempty(findstr(str,'.'))==0 || isempty(findstr(str,';'))==0
            msgbox(['Error in frame selection syntax. Please use the following syntax (examples):' sprintf('\n') '1:3' sprintf('\n') '1,3,7,9' sprintf('\n') '1:3,7,8,9,11:13' ],'Error','error','modal')
            selectionok=0;
        end
        if selectionok==1
            mincount=(min(strnum));
            for count=mincount:size(resultslist,2)
                if size(resultslist,2)>=count && numel(resultslist{1,count})>0
                    x=resultslist{1,count};
                    y=resultslist{2,count};
                    if size(resultslist,1)>6 %filtered exists
                        if size(resultslist,1)>10 && numel(resultslist{10,count}) > 0 %smoothed exists
                            u=resultslist{10,count};
                            v=resultslist{11,count};
                            typevector=resultslist{9,count};
                            if numel(typevector)==0 %happens if user smoothes sth without NaN and without validation
                                typevector=resultslist{5,count};
                            end
                        else
                            u=resultslist{7,count};
                            if size(u,1)>1
                                v=resultslist{8,count};
                                typevector=resultslist{9,count};
                            else %filter was applied for other frames but not for this one
                                u=resultslist{3,count};
                                v=resultslist{4,count};
                                typevector=resultslist{5,count};
                            end
                        end
                    else
                        u=resultslist{3,count};
                        v=resultslist{4,count};
                        typevector=resultslist{5,count};
                    end
                    
                    %if count==mincount %besser: wenn orgsize nicht existiert
                    if exist('originalsizex')==0
                        originalsizex=size(u,2);
                        originalsizey=size(u,1);
                    else
                        
                        if size(u,2)~=originalsizex || size(u,1)~=originalsizey
                            sizeerror=1;
                        end
                    end
                    if ismean(count,1)==0 && sizeerror==0
                        umittel(:,:,count)=u;
                        vmittel(:,:,count)=v;
                    end
                    if sizeerror==0
                        typevectormittel(:,:,count)=typevector;
                    end
                end
                
            end
            if sizeerror==0
                for i=1:size(strnum,2)
                    if size(resultslist,2)>=strnum(i) %dann ok
                        x_tmp=resultslist{1,strnum(i)};
                        if isempty(x_tmp)==1 %dann nicht ok
                            msgbox('Your selected range includes non-analyzed frames.','Error','error','modal')
                            selectionok=0;
                            break
                        end
                    else
                        msgbox('Your selected range includes non-analyzed frames.','Error','error','modal')
                        selectionok=0;
                        break
                    end
                    if size(ismean,1)>=strnum(i)
                        if ismean(strnum(i))==1
                            msgbox('You must not include frames in your selection that already consist of mean vectors.','Error','error','modal')
                            selectionok=0;
                            break
                        end
                    else
                        msgbox('Your selected range exceeds the amount of analyzed frames.','Error','error','modal')
                        selectionok=0;
                        break
                    end
                end
                
                if selectionok==1
                    maskiererx=retr('maskiererx');
                    maskierery=retr('maskierery');
                    if isempty(maskiererx)==1
                        maskiererx=cell(1,1);
                        maskierery=cell(1,1);
                    end
                    maskiererx_temp=cell(1,1);
                    maskierery_temp=cell(1,1);
                    maskiererx_temp=maskiererx(:,1:2:end);
                    maskierery_temp=maskierery(:,1:2:end);
                    %kopieren in temp "originalmaske", dann alles löschen was nicht
                    %ausgewählt. (auf [] setzen)
                    % z.B.: maskiererxselected=maskiererx_temp(1,[1:6])
                    try
                        eval(['maskiererxselected=maskiererx_temp(:,[' str ']);']);
                        eval(['maskiereryselected=maskierery_temp(:,[' str ']);']);
                    catch
                        maskiererxselected=cell(1,1);
                        maskiereryselected=cell(1,1);
                    end
                    newmaskx=cell(0,0);
                    for i=1:size(maskiererxselected,1)
                        for j=1:size(maskiererxselected,2)
                            if numel(maskiererxselected{i,j})~=0
                                newmaskx{size(newmaskx,1)+1,1}=maskiererxselected{i,j};
                            end
                        end
                    end
                    for i=size(newmaskx,1):-1:2
                        if numel(newmaskx{i,1})==numel(newmaskx{i-1,1})
                            A=newmaskx{i-1,1};
                            B=newmaskx{i,1};
                            if mean(A-B)==0
                                newmaskx{i,1}={};
                            end
                        end
                    end
                    
                    try
                        newmaskx(cellfun(@isempty,newmaskx))=[];
                    catch
                        disp('Problems with old Matlab version... Please update Matlab or unexpected things might happen...')
                    end
                    newmasky=cell(0,0);
                    for i=1:size(maskiereryselected,1)
                        for j=1:size(maskiereryselected,2)
                            if numel(maskiereryselected{i,j})~=0
                                newmasky{size(newmasky,1)+1,1}=maskiereryselected{i,j};
                            end
                        end
                    end
                    for i=size(newmasky,1):-1:2
                        if numel(newmasky{i,1})==numel(newmasky{i-1,1})
                            A=newmasky{i-1,1};
                            B=newmasky{i,1};
                            if mean(A-B)==0
                                newmasky{i,1}={};
                            end
                        end
                    end
                    try
                        newmasky(cellfun(@isempty,newmasky))=[];
                    catch
                        disp('Problems with old Matlab version... Please update Matlab or unexpected things might happen...')
                    end
                    for i=1:size(newmaskx,1)
                        %ans Ende der originalmaske wird eine zusammengesetzte maske
                        %aus allen gewählten frames gehängt.
                        maskiererx{i,size(filepath,1)+1}=newmaskx{i,1};
                        maskiererx{i,size(filepath,1)+2}=newmaskx{i,1};
                        maskierery{i,size(filepath,1)+1}=newmasky{i,1};
                        maskierery{i,size(filepath,1)+2}=newmasky{i,1};
                    end
                    put('maskiererx',maskiererx);
                    put('maskierery',maskierery);
                    typevectoralle=ones(size(typevector));
                    
                    
                    
                    %Hier erst neue matrix erstellen mit ausgewählten frames
                    %typevectoralle ist ausgabe für gui
                    %typevectormean ist der mittelwert aller types
                    %typevectormittel ist der stapel aus allen typevectors
                    
                    eval(['typevectormittelselected=typevectormittel(:,:,[' str ']);']);
                    
                    typevectormean=mean(typevectormittelselected,3);
                    %for i=1:size(typevectormittelselected,3)
                    for i=1:size(typevectormittelselected,1)
                        for j=1:size(typevectormittelselected,2)
                            if mean(typevectormittelselected(i,j,:))==0
                                typevectoralle(i,j)=0;
                            end
                        end
                    end
                    %da wo ALLE null sidn auf null setzen.
                    %typevectoralle(typevectormittelselected(:,:,i)==0)=0; %maskierte vektoren sollen im Mean maskiert sein
                    % end
                    
                    typevectoralle(typevectormean>1.5)=2; %if more than 50% of vectors are interpolated, then mark vector in mean as interpolated too.
                    resultslist{5,size(filepath,1)/2+1}=typevectoralle;
                    resultslist{1,size(filepath,1)/2+1}=x;
                    resultslist{2,size(filepath,1)/2+1}=y;
                    
                    %hier neue matrix mit ausgewählten frames!
                    eval(['umittelselected=umittel(:,:,[' str ']);']);
                    eval(['vmittelselected=vmittel(:,:,[' str ']);']);
                    
                    resultslist{3,size(filepath,1)/2+1}=nanmean(umittelselected,3);
                    resultslist{4,size(filepath,1)/2+1}=nanmean(vmittelselected,3);
                    filepathselected=filepath(1:2:end);
                    eval(['filepathselected=filepathselected([' str '],:);']);
                    filepath{size(filepath,1)+1,1}=filepathselected{1,1};
                    filepath{size(filepath,1)+1,1}=filepathselected{1,1};
                    filename=retr('filename');
                    
                    filename{size(filename,1)+1,1}=['MEAN of frames ' str];
                    filename{size(filename,1)+1,1}=['MEAN of frames ' str];
                    
                    ismean(size(resultslist,2),1)=1;
                    put('ismean',ismean);
                    
                    put ('resultslist', resultslist);
                    put ('filepath', filepath);
                    put ('filename', filename);
                    put ('typevector', typevector);
                    sliderrange
                    try
                        set (handles.fileselector,'value',get (handles.fileselector,'max'));
                    catch
                    end
                    
                    sliderdisp
                end
            else %user tried to average analyses with different sizes
                errordlg('All analyses of one session have to be of the same size and have to be analyzed with identical PIV settings.','Averaging not possible...')
            end
        end
    end
end

function part_am_Callback(hObject, eventdata, handles)
check_comma(hObject)
function part_size_Callback(hObject, eventdata, handles)
check_comma(hObject)
function part_var_Callback(hObject, eventdata, handles)
check_comma(hObject)
function part_noise_Callback(hObject, eventdata, handles)
check_comma(hObject)
function oseenx1_Callback(hObject, eventdata, handles)
check_comma(hObject)
function rank_core_Callback(hObject, eventdata, handles)
check_comma(hObject)
function rank_displ_Callback(hObject, eventdata, handles)
check_comma(hObject)
function rotationdislacement_Callback(hObject, eventdata, handles)
check_comma(hObject)
function realdist_Callback(hObject, eventdata, handles)
check_comma(hObject)
function time_inp_Callback(hObject, eventdata, handles)
check_comma(hObject)
function subtr_u_Callback(hObject, eventdata, handles)
check_comma(hObject)
function subtr_v_Callback(hObject, eventdata, handles)
check_comma(hObject)
function mapscale_min_Callback(hObject, eventdata, handles)
check_comma(hObject)
function mapscale_max_Callback(hObject, eventdata, handles)
check_comma(hObject)
function stdev_thresh_Callback(hObject, eventdata, handles)
check_comma(hObject)
function loc_med_thresh_Callback(hObject, eventdata, handles)
check_comma(hObject)
function epsilon_Callback(hObject, eventdata, handles)
check_comma(hObject)
function thresholdarea_Callback(hObject, eventdata, handles)
check_comma(hObject)
function shiftdisplacement_Callback(hObject, eventdata, handles)
check_comma(hObject)
function sheetthick_Callback(hObject, eventdata, handles)
check_comma(hObject)
function ranky1_Callback(hObject, eventdata, handles)
check_comma(hObject)
function rankx1_Callback(hObject, eventdata, handles)
check_comma(hObject)
function rankx2_Callback(hObject, eventdata, handles)
check_comma(hObject)
function ranky2_Callback(hObject, eventdata, handles)
check_comma(hObject)
function oseen_displ_Callback(hObject, eventdata, handles)
check_comma(hObject)
function oseenx2_Callback(hObject, eventdata, handles)
check_comma(hObject)
function oseeny1_Callback(hObject, eventdata, handles)
check_comma(hObject)
function oseeny2_Callback(hObject, eventdata, handles)
check_comma(hObject)
function oseen_time_Callback(hObject, eventdata, handles)
check_comma(hObject)
function part_z_Callback(hObject, eventdata, handles)
check_comma(hObject)
function vecwidth_Callback(hObject, eventdata, handles)
check_comma(hObject)

% --- Executes on button press in savejpgseq.
function savejpgseq_Callback(hObject, eventdata, handles)

% --- Executes on button press in avifilesave.
function avifilesave_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value')==1
    set (handles.jpgfilesave, 'value', 0);
    set (handles.bmpfilesave, 'value', 0);
    set (handles.epsfilesave, 'value', 0);
    set (handles.pdffilesave, 'value', 0);
    set(handles.usecompr,'enable','on');
    set(handles.fps_setting,'enable','on');
else
    set (handles.avifilesave, 'value', 1);
end

% --- Executes on button press in jpgfilesave.
function jpgfilesave_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value')==1
    set (handles.avifilesave, 'value', 0);
    set (handles.bmpfilesave, 'value', 0);
    set (handles.epsfilesave, 'value', 0);
    set (handles.pdffilesave, 'value', 0);
    set(handles.usecompr,'value',0);
    set(handles.usecompr,'enable','off');
    set(handles.fps_setting,'enable','off');
    
else
    set (handles.jpgfilesave, 'value', 1);
end

% --- Executes on button press in bmpfilesave.
function bmpfilesave_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value')==1
    set (handles.avifilesave, 'value', 0);
    set (handles.jpgfilesave, 'value', 0);
    set (handles.epsfilesave, 'value', 0);
    set (handles.pdffilesave, 'value', 0);
    set(handles.usecompr,'value',0);
    set(handles.usecompr,'enable','off');
    set(handles.fps_setting,'enable','off');
    
else
    set (handles.bmpfilesave, 'value', 1);
end

function pdffilesave_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value')==1
    set (handles.avifilesave, 'value', 0);
    set (handles.jpgfilesave, 'value', 0);
    set (handles.epsfilesave, 'value', 0);
    set (handles.bmpfilesave, 'value', 0);
    set(handles.usecompr,'value',0);
    set(handles.usecompr,'enable','off');
    set(handles.fps_setting,'enable','off');
    
else
    set (handles.pdffilesave, 'value', 1);
end


% --- Executes on button press in drawstreamlines.
function drawstreamlines_Callback(hObject, eventdata, handles)
handles=gethand;
toggler=retr('toggler');
selected=2*floor(get(handles.fileselector, 'value'))-(1-toggler);
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    toolsavailable(0);
    x=resultslist{1,currentframe};
    y=resultslist{2,currentframe};
    typevector=resultslist{5,currentframe};
    if size(resultslist,1)>6 %filtered exists
        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
            u=resultslist{10,currentframe};
            v=resultslist{11,currentframe};
            typevector=resultslist{9,currentframe}; %von smoothed
        else
            u=resultslist{7,currentframe};
            if size(u,1)>1
                v=resultslist{8,currentframe};
                typevector=resultslist{9,currentframe}; %von smoothed
            else
                u=resultslist{3,currentframe};
                v=resultslist{4,currentframe};
                typevector=resultslist{5,currentframe};
            end
        end
    else
        u=resultslist{3,currentframe};
        v=resultslist{4,currentframe};
    end
    ismean=retr('ismean');
    if    numel(ismean)>0
        if ismean(currentframe)==1 %if current frame is a mean frame, typevector is stored at pos 5
            typevector=resultslist{5,currentframe};
        end
    end
    caluv=retr('caluv');
    u=u*caluv-retr('subtr_u');
    v=v*caluv-retr('subtr_v');
    u(typevector==0)=nan;
    v(typevector==0)=nan;
    calxy=retr('calxy');
    button=1;
    streamlinesX=retr('streamlinesX');
    streamlinesY=  retr('streamlinesY');
    if get(handles.holdstream,'value')==1
        if numel(streamlinesX)>0
            i=size(streamlinesX,2)+1;
            xposition=streamlinesX;
            yposition=streamlinesY;
        else
            i=1;
        end
    else
        i=1;
        put('streamlinesX',[]);
        put('streamlinesY',[]);
        xposition=[];
        yposition=[];
        delete(findobj('tag','streamline'));
    end
    while button == 1
        [rawx,rawy,button] = ginput(1);
        if button~=1
            break
        end
        xposition(i)=rawx;
        yposition(i)=rawy;
        h=streamline(mmstream2(x,y,u,v,xposition(i),yposition(i),'on'));
        set (h,'tag','streamline');
        i=i+1;
    end
    delete(findobj('tag','streamline'));
    if exist('xposition')==1
        h=streamline(mmstream2(x,y,u,v,xposition,yposition,'on'));
        set (h,'tag','streamline');
        contents = get(handles.streamlcolor,'String');
        set(h,'LineWidth',get(handles.streamlwidth,'value'),'Color', contents{get(handles.streamlcolor,'Value')})
        put('streamlinesX',xposition);
        put('streamlinesY',yposition);
    end
end
toolsavailable(1);

% --------------------------------------------------------------------
function Untitled_23_Callback(hObject, eventdata, handles)
switchui('multip18');

% --- Executes on button press in deletestreamlines.
function deletestreamlines_Callback(hObject, eventdata, handles)
put('streamlinesX',[]);
put('streamlinesY',[]);
delete(findobj('tag','streamline'));

% --- Executes on button press in streamrake.
function streamrake_Callback(hObject, eventdata, handles)
handles=gethand;
toggler=retr('toggler');
selected=2*floor(get(handles.fileselector, 'value'))-(1-toggler);
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    toolsavailable(0);
    x=resultslist{1,currentframe};
    y=resultslist{2,currentframe};
    typevector=resultslist{5,currentframe};
    if size(resultslist,1)>6 %filtered exists
        if size(resultslist,1)>10 && numel(resultslist{10,currentframe}) > 0 %smoothed exists
            u=resultslist{10,currentframe};
            v=resultslist{11,currentframe};
            typevector=resultslist{9,currentframe}; %von smoothed
        else
            u=resultslist{7,currentframe};
            if size(u,1)>1
                v=resultslist{8,currentframe};
                typevector=resultslist{9,currentframe}; %von smoothed
            else
                u=resultslist{3,currentframe};
                v=resultslist{4,currentframe};
                typevector=resultslist{5,currentframe};
            end
        end
    else
        u=resultslist{3,currentframe};
        v=resultslist{4,currentframe};
    end
    ismean=retr('ismean');
    if    numel(ismean)>0
        if ismean(currentframe)==1 %if current frame is a mean frame, typevector is stored at pos 5
            typevector=resultslist{5,currentframe};
        end
    end
    caluv=retr('caluv');
    u=u*caluv-retr('subtr_u');
    v=v*caluv-retr('subtr_v');
    u(typevector==0)=nan;
    v(typevector==0)=nan;
    calxy=retr('calxy');
    button=1;
    streamlinesX=retr('streamlinesX');
    streamlinesY=  retr('streamlinesY');
    if get(handles.holdstream,'value')==1
        if numel(streamlinesX)>0
            i=size(streamlinesX,2)+1;
            xposition=streamlinesX;
            yposition=streamlinesY;
        else
            i=1;
        end
    else
        i=1;
        put('streamlinesX',[]);
        put('streamlinesY',[]);
        xposition=[];
        yposition=[];
        delete(findobj('tag','streamline'));
    end
    [rawx,rawy,button] = ginput(1);
    hold on; scatter(rawx,rawy,'y*','tag','streammarker');hold off;
    [rawx(2),rawy(2),button] = ginput(1);
    delete(findobj('tag','streammarker'))
    rawx=linspace(rawx(1),rawx(2),str2num(get(handles.streamlamount,'string')));
    rawy=linspace(rawy(1),rawy(2),str2num(get(handles.streamlamount,'string')));
    
    xposition(i:i+str2num(get(handles.streamlamount,'string'))-1)=rawx;
    yposition(i:i+str2num(get(handles.streamlamount,'string'))-1)=rawy;
    h=streamline(mmstream2(x,y,u,v,xposition(i),yposition(i),'on'));
    set (h,'tag','streamline');
    i=i+1;
end
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    delete(findobj('tag','streamline'));
    h=streamline(mmstream2(x,y,u,v,xposition,yposition,'on'));
    contents = get(handles.streamlcolor,'String');
    set(h,'LineWidth',get(handles.streamlwidth,'value'),'Color', contents{get(handles.streamlcolor,'Value')});
    set (h,'tag','streamline');
    put('streamlinesX',xposition);
    put('streamlinesY',yposition);
end
toolsavailable(1);

% --- Executes on button press in applycolorwidth.
function applycolorwidth_Callback(hObject, eventdata, handles)
sliderdisp

% --- Executes on button press in putmarkers.
function putmarkers_Callback(hObject, eventdata, handles)
handles=gethand;
button=1;
manmarkersX=retr('manmarkersX');
manmarkersY=retr('manmarkersY');
if get(handles.holdmarkers,'value')==1
    
    if numel(manmarkersX)>0
        i=size(manmarkersX,2)+1;
        xposition=manmarkersX;
        yposition=manmarkersY;
    else
        i=1;
    end
else
    i=1;
    put('manmarkersX',[]);
    put('manmarkersY',[]);
    xposition=[];
    yposition=[];
    delete(findobj('tag','manualmarker'));
end
hold on;
toolsavailable(0)
while button == 1
    [rawx,rawy,button] = ginput(1);
    if button~=1
        break
    end
    xposition(i)=rawx;
    yposition(i)=rawy;
    plot(xposition(i),yposition(i), 'r*','Color', [0.55,0.75,0.9], 'tag', 'manualmarker');
    i=i+1;
end
toolsavailable(1)
delete(findobj('tag','manualmarker'));
plot(xposition,yposition, 'o','MarkerEdgeColor','k','MarkerFaceColor',[.2 .2 1], 'MarkerSize',9, 'tag', 'manualmarker');
plot(xposition,yposition, '*','MarkerEdgeColor','w', 'tag', 'manualmarker');
put('manmarkersX',xposition);
put('manmarkersY',yposition);
hold off

% --- Executes on button press in delmarkers.
function delmarkers_Callback(hObject, eventdata, handles)
put('manmarkersX',[]);
put('manmarkersY',[]);
delete(findobj('tag','manualmarker'));

% --- Executes on button press in holdmarkers.
function holdmarkers_Callback(hObject, eventdata, handles)

% --- Executes on button press in displmarker.
function displmarker_Callback(hObject, eventdata, handles)
sliderdisp;

% --- Executes on button press in dcc.
function dcc_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value')==1
    set(handles.fftmulti,'value',0)
    set(handles.uipanel42,'visible','off')
else
    set(handles.dcc,'value',1)
end
%{
if get(handles.enable_highpass,'value')==0
    button = questdlg(['For the DCC algorithm, it is recommended to enable the high-pass in image pre-processing.' sprintf('\n\n') 'Enable high-pass now?'],'High-pass recommended','Yes','No','Yes');
    if strmatch(button,'Yes')==1
        set(handles.enable_highpass,'value',1)
    end
end
%}
countparticles
set(handles.recommendation,'visible','on');
dispinterrog

% --- Executes on button press in fftmulti.
function fftmulti_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value') ==1
    set(handles.dcc,'value',0)
    set(handles.uipanel42,'visible','on')
else
    set(handles.fftmulti,'value',1)
end
%{
if get(handles.enable_highpass,'value')==1
    button = questdlg(['For the FFT multi-pass algorithm, it is recommended to disable the high-pass in image pre-processing.' sprintf('\n\n') 'Disable high-pass now?'],'High-pass not recommended','Yes','No','Yes');
    if strmatch(button,'Yes')==1
        set(handles.enable_highpass,'value',0)
    end
end
%}
set(handles.recommendation,'visible','off');
dispinterrog
% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value') == 0
    set(handles.edit50,'enable','off')
    set(handles.edit51,'enable','off')
    set(handles.edit52,'enable','off')
    set(handles.checkbox27,'value',0)
    set(handles.checkbox28,'value',0)
else
    set(handles.edit50,'enable','on')
end
dispinterrog

% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value') == 0
    set(handles.edit51,'enable','off')
    set(handles.edit52,'enable','off')
    set(handles.checkbox28,'value',0)
else
    set(handles.edit50,'enable','on')
    set(handles.edit51,'enable','on')
    set(handles.checkbox26,'value',1)
end
if get(handles.checkbox26,'value')==0
    set(handles.checkbox27,'value',0)
    set(handles.edit51,'enable','off')
end
dispinterrog

% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value') == 0
    set(handles.edit52,'enable','off')
else
    set(handles.edit52,'enable','on')
    set(handles.edit50,'enable','on')
    set(handles.edit51,'enable','on')
    set(handles.checkbox26,'value',1)
    set(handles.checkbox27,'value',1)
    
end
if get(handles.checkbox27,'value')==0
    set(handles.checkbox28,'value',0)
    set(handles.edit52,'enable','off')
end
dispinterrog

function edit50_Callback(hObject, eventdata, handles)
handles=gethand;
step=str2double(get(hObject,'String'));
set (handles.text126, 'string', int2str(step/2));
dispinterrog


function edit51_Callback(hObject, eventdata, handles)
handles=gethand;
step=str2double(get(hObject,'String'));
set (handles.text127, 'string', int2str(step/2));
dispinterrog

function edit52_Callback(hObject, eventdata, handles)
handles=gethand;
step=str2double(get(hObject,'String'));
set (handles.text128, 'string', int2str(step/2));
dispinterrog

function rangeRGB_Callback(hObject, eventdata, handles)
check_comma(hObject)
val=get(hObject,'string');
if str2double(val)>1
    set(hObject,'string',1);
end
if str2double(val)<0 || isempty(val)==1 || isnan(str2double(val))==1
    set(hObject,'string',0);
end

function vectorskip_Callback(hObject, eventdata, handles)
check_comma(hObject)
val=get(hObject,'string');
if str2double(val)<1 || isempty(val)==1|| isnan(str2double(val))==1
    set(hObject,'string',1);
end
resultslist=retr('resultslist');
handles=gethand;
currentframe=2*floor(get(handles.fileselector, 'value'))-1;
if size(resultslist,2)>=(currentframe+1)/2 && numel(resultslist{1,(currentframe+1)/2})>0
    if str2double(val) > size(resultslist{1,(currentframe+1)/2},1)/2-1 || str2double(val) > size(resultslist{2,(currentframe+1)/2},2)/2-1
        set(hObject,'string',min([size(resultslist{1,(currentframe+1)/2},1)/2-1 size(resultslist{2,(currentframe+1)/2},2)/2-1]));
    end
end


% --- Executes on button press in extractareaall.
function extractareaall_Callback(hObject, eventdata, handles)
% hObject    handle to extractareaall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=gethand;
if get(hObject,'Value')==1
    set(handles.savearea,'enable','off');
    set(handles.savearea,'value',1);
else
    set(handles.savearea,'enable','on');
end

function radiusincrease_Callback(hObject, eventdata, handles)
check_comma(hObject)
val=get(hObject,'string');
if str2double(val)>500
    set(hObject,'string',500);
end
if str2double(val)<0 || isempty(val)==1 || isnan(str2double(val))==1
    set(hObject,'string',0);
end
function derivdropdown(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'value')==10
    set(handles.LIChint1,'visible','on');
    set(handles.LIChint2,'visible','on');
    %set(handles.LIChint3,'visible','on');
    set(handles.licres,'visible','on');
else
    set(handles.LIChint1,'visible','off');
    set(handles.LIChint2,'visible','off');
    %set(handles.LIChint3,'visible','off');
    set(handles.licres,'visible','off');
end


function LICreschanged(hObject, eventdata, handles)
handles=gethand;
LICreso=round(get (handles.licres, 'value')*10)/10;
set(handles.LIChint2,'string',num2str(LICreso));


% --- Executes on selection change in derivchoice.
function derivchoice_Callback(hObject, eventdata, handles)
handles=gethand;
contents = get(hObject,'String');
currstring=contents{get(hObject,'Value')};
currstring=currstring(strfind(currstring,'['):end);
set(handles.text39,'String', ['min ' currstring ':']);
set(handles.text40,'String', ['max ' currstring ':']);
derivdropdown(hObject);




% --- Executes on button press in epsfilesave.
function epsfilesave_Callback(hObject, eventdata, handles)
handles=gethand;
if get(hObject,'Value')==1
    set (handles.avifilesave, 'value', 0);
    set (handles.jpgfilesave, 'value', 0);
    set (handles.bmpfilesave, 'value', 0);
    set (handles.pdffilesave, 'value', 0);
    set(handles.usecompr,'value',0);
    set(handles.usecompr,'enable','off');
    set(handles.fps_setting,'enable','off');
    
else
    set (handles.epsfilesave, 'value', 1);
end

function smallerlarger_Callback(hObject, eventdata, handles)
%do nothing


% --- Executes on button press in zoomon.
function zoomon_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    zoom(gca,'on')
    handles=gethand;
    set(handles.panon,'Value',0);
else
    zoom(gca,'off')
    put('xzoomlimit', get (gca, 'xlim'));
    put('yzoomlimit', get (gca, 'ylim'));
end


% --- Executes on button press in panon.
function panon_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    pan(gca,'on')
    handles=gethand;
    set(handles.zoomon,'Value',0);
else
    pan(gca,'off')
    put('xzoomlimit', get (gca, 'xlim'));
    put('yzoomlimit', get (gca, 'ylim'));
end



function selectedFramesMean_Callback(hObject, eventdata, handles)
% hObject    handle to selectedFramesMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectedFramesMean as text
%        str2double(get(hObject,'String')) returns contents of selectedFramesMean as a double


% --- Executes during object creation, after setting all properties.
function selectedFramesMean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectedFramesMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







function maskapplyselect_Callback(hObject, eventdata, handles)
% hObject    handle to maskapplyselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maskapplyselect as text
%        str2double(get(hObject,'String')) returns contents of maskapplyselect as a double


% --- Executes during object creation, after setting all properties.
function maskapplyselect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskapplyselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in paraview_current.
function paraview_current_Callback(hObject, eventdata, handles)
handles=gethand;
resultslist=retr('resultslist');
currentframe=floor(get(handles.fileselector, 'value'));
if size(resultslist,2)>=currentframe && numel(resultslist{1,currentframe})>0
    [FileName,PathName] = uiputfile('*.vtk','Save Paraview binary vtk as...','PIVlab.vtk'); %framenummer in dateiname
    if isequal(FileName,0) | isequal(PathName,0)
    else
        file_save(currentframe,FileName,PathName,3);
    end
end


% --- Executes on button press in paraview_all.
function paraview_all_Callback(hObject, eventdata, handles)
handles=gethand;
filepath=retr('filepath');
resultslist=retr('resultslist');
[FileName,PathName] = uiputfile('*.vtk','Save Paraview binary vtk as...','PIVlab.vtk'); %framenummer in dateiname
if isequal(FileName,0) | isequal(PathName,0)
else
    toolsavailable(0)
    for i=1:floor(size(filepath,1)/2)
        %if analysis exists
        if size(resultslist,2)>=i && numel(resultslist{1,i})>0
            [Dir Name Ext] = fileparts(FileName);
            FileName_nr=[Name sprintf('_%.4d', i) Ext];
            file_save(i,FileName_nr,PathName,3)
            set (handles.paraview_all, 'string', ['Please wait... (' int2str((i-1)/size(filepath,1)*200) '%)']);
            drawnow;
        end
    end
    toolsavailable(1)
    set (handles.paraview_all, 'string', 'Save all frames');
end


% --------------------------------------------------------------------
function paraview_Callback(hObject, eventdata, handles)
switchui('multip19')




% --------------------------------------------------------------------
function Website_Callback(hObject, eventdata, handles)
try
web('http://pivlab.blogspot.com','-browser')
catch
    %why does 'web' not work in v 7.1.0.246 ...?
    disp('Ooops, MATLAB couldn''t open the website.')
    disp('You''ll have to open the website manually:')
    disp('http://PIVlab.blogspot.de')
end

function Man_ROI_Callback
handles=gethand;
try
    x=round(str2num(get(handles.ROI_Man_x,'String')));
    y=round(str2num(get(handles.ROI_Man_y,'String')));
    w=round(str2num(get(handles.ROI_Man_w,'String')));
    h=round(str2num(get(handles.ROI_Man_h,'String')));
catch
end
if isempty(x)== 0 && isempty(y)== 0 && isempty(w)== 0 && isempty(h)== 0 && isnumeric(x) && isnumeric(y) && isnumeric(w) && isnumeric(h)
    roirect(1)=x;
    roirect(2)=y;
    roirect(3)=w;
    roirect(4)=h;

    toggler=retr('toggler');
    selected=2*floor(get(handles.fileselector, 'value'))-(1-toggler);
    filepath=retr('filepath');

    imagesize(1)=size(imread(filepath{selected}),1);
    imagesize(2)=size(imread(filepath{selected}),2);
    if roirect(1)<1
        roirect(1)=1;
    end
    if roirect(2)<1
        roirect(2)=1;
    end
    if roirect(3)>imagesize(2)-roirect(1)
        roirect(3)=imagesize(2)-roirect(1);
    end
    if roirect(4)>imagesize(1)-roirect(2)
        roirect(4)=imagesize(1)-roirect(2);
    end
    put ('roirect',roirect);
    dispROI
    set(handles.roi_hint, 'String', 'ROI active' , 'backgroundcolor', [0.5 1 0.5]);
end

function ROI_Man_x_Callback(hObject, eventdata, handles)
Man_ROI_Callback

function ROI_Man_y_Callback(hObject, eventdata, handles)
Man_ROI_Callback

function ROI_Man_w_Callback(hObject, eventdata, handles)
Man_ROI_Callback

function ROI_Man_h_Callback(hObject, eventdata, handles)
Man_ROI_Callback


% --------------------------------------------------------------------
function howtocite_Callback(hObject, eventdata, handles)
PIVlab_citing



% --------------------------------------------------------------------
function exitpivlab_Callback(hObject, eventdata, handles)
close(gcf)


% --------------------------------------------------------------------
function Forum_Callback(hObject, eventdata, handles)
try
web('http://pivlab.blogspot.de/p/forum.html','-browser')
catch
    %why does 'web' not work in v 7.1.0.246 ...?
    disp('Ooops, MATLAB couldn''t open the website.')
    disp('You''ll have to open the website manually:')
    disp('http://pivlab.blogspot.de/p/forum.html')
end
