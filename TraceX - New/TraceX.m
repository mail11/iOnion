function varargout = TraceX(varargin)
% TRACEX MATLAB code for TraceX.fig v1.05
%      TRACEX, by itself, creates a new TRACEX or raises the existing
%      singleton*.

%
%      H = TRACEX returns the handle to a new TRACEX or the handle to
%      the existing singleton*.
%
%      TRACEX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACEX.M with the given input arguments.
%
%      TRACEX('Property','Value',...) creates a new TRACEX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TraceX_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TraceX_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TraceX

% Last Modified by GUIDE v2.5 21-Oct-2015 09:59:32
% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @TraceX_OpeningFcn, ...
    'gui_OutputFcn',  @TraceX_OutputFcn, ...
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
% End initialization code - DO NOT EDIT

% --- Executes just before TraceX is made visible.
function TraceX_OpeningFcn(hObject, eventdata, handles, varargin)
% get(handles.Xlimit,'Position')
handles.ImagePlotPos=[.01 .09 .71 .71];%norm
handles.ImagePlotPosXlim=[.605 .0475 .05 .03];%norm
handles.ImagePlotPosYlim=[.041 .785 .047 .030];%norm

handles.GraphPlotPos=[.08 .09 .66 .71];%norm
handles.GraphPlotPosXlim=[.708 .05 .05 .033];%norm
handles.GraphPlotPosYlim=[.028 .782 .05 .033];%norm

handles.output = hObject;
% axes(handles.MainAxes)
%imshow('Electroceramics Group.png')
[A, map, alpha] = imread('TraceX_logo.png');
h = imshow(A, map);
set(h, 'AlphaData', alpha);
% set(handles.MainAxes, 'Visible','off');
% set(handles.MinorAxes,'Visible','off')
%%%Multiline tooltips
set(handles.Fit_Start, 'Tooltipstring', ...
    sprintf(['The distance from the surface before which the fitting algorithm will ignore. \n',...
    'If the box is green, then the value has been automatically set as the \n',...
    'pixel depth which first has 90%% of the mean counts.']));
set(handles.Fit_End, 'Tooltipstring', ...
    sprintf(['The distance from the surface after which the fitting algorithm will ignore. \n',...
    'If the box is green, then the value has been automatically set as the \n',...
    'pixel depth which first falls below the background.']));
set(handles.C_bg, 'Tooltipstring', ...
    sprintf(['Background 18O isotopic fraction \n',...
    'If the box is green, then the value has been automatically \n',...
    'set as the modal value of the profile']));
set(handles.AlignAngle, 'Tooltipstring', ...
    sprintf(['Select an alignment angle for the sample\n',...
    'If the box is green, then the value has been automatically \n',...
    'set using the method specified above']));
set(handles.ROI, 'Tooltipstring', ...
    sprintf(['Specify a region of interest using the mouse\n',...
    'to place vertices of polygone, then double-click \n',...
    'on first vertex to complete shape']));
set(handles.GenerateProfiles, 'Tooltipstring', ...
    sprintf(['Generate a range of profiles to illustrate the effect of masking and alignment.\n',...
    'Profile used for fitting is specified using dropdown box above.']));

guidata(hObject, handles);
% --- Outputs from this function are returned to the command line.
function varargout = TraceX_OutputFcn(hObject, eventdata, handles)
% set(gcf,'units','centimeters','outerposition',[0 3 25 20]);
varargout{1} = handles.output;
% set(gcf, 'Units' , 'Normalized');



%% --- Executes on button press in LoadProfileData.
function [handles]=LoadProfileData_Callback(hObject, eventdata, handles)
handles.ProLen_Or=500;
set(handles.Fit_End,'String',250);
[handles.FileName, handles.PathName] = uigetfile('*.*');%(FILTERSPEC, TITLE);
handles.DataFlag=1;
handles.PlotButtonFlag=0;
if ispc
    handles.ProfileData_Or=textread([handles.PathName,'\',handles.FileName],'%n')';
else
    handles.ProfileData_Or=textread([handles.PathName,handles.FileName],'%n')';
end
% handles.ProfileData=handles.ProfileData_Or;

set(handles.Xlimit,'String',get(handles.ProfileLength,'String'));
set(handles.CprimeY,'State','off')
set(handles.GBplot,'State','off')
set(handles.SurfPos,'String',0)

handles.SaveName_Or=[handles.PathName,'\',handles.FileName(1:end-4),...
    '_Norm18Prof.txt'];
set(handles.ProfDataSaveName, 'String',handles.FileName);
set(handles.C_bg,'String',...
    mode(roundsf(handles.ProfileData_Or,2,'round')));
set(handles.C_bg,'BackgroundColor',[0,1,0]);
set(handles.PixelNo,'String',length(handles.ProfileData_Or));
set(handles.PixelNo,'BackgroundColor',[0,1,0]);
set(handles.Fit_End,'String',get(handles.ProfileLength,'String'));
Fit_End_Callback(hObject, eventdata, handles)
set(handles.Fit_End,'BackgroundColor',[0,1,0]);

[handles]=DataPlotOnly(hObject, eventdata, handles);
guidata(hObject, handles);


function ProfileLength_Callback(hObject, eventdata, handles)
set(handles.PlaneSheet, 'Value',0)
if isfield(handles,'ProLen_Or')
    LenRat=str2double(get(handles.ProfileLength,'String'))/handles.ProLen_Or;
    handles.ProLen_Or=str2double(get(handles.ProfileLength,'String'));
    set(handles.SurfPos,'String',roundsf(LenRat*str2double(get(handles.SurfPos,'String')),3,'round'));
    set(handles.Fit_Start,'String',roundsf(LenRat*str2double(get(handles.Fit_Start,'String')),3,'round'));
    set(handles.Fit_End,'String',roundsf(LenRat*str2double(get(handles.Fit_End,'String')),3,'round'));
end
set(handles.Xlimit,'String',get(handles.ProfileLength,'String'));
Xlimit_Callback(hObject, eventdata, handles);
set(handles.ProfileLength,'BackgroundColor',[1,1,1]);
handles.ProLen_num=str2double(get(handles.ProfileLength,'String'));
if isfield(handles,'PlotType')
    if handles.PlotType==2;
        [handles]=PlotButton_Callback(hObject, eventdata, handles);
    elseif handles.PlotType==3;
        [handles]=DataPlotOnly(hObject, eventdata, handles);
    end
end

guidata(hObject, handles);
function ProfileLength_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function PixelNo_Callback(hObject, eventdata, handles)
set(handles.PixelNo,'BackgroundColor',[1,1,1]);
handles.PixelNo_num=str2double(get(handles.PixelNo, 'String'));
if isfield(handles,'PlotType')
    if handles.PlotType==2;
        [handles]=PlotButton_Callback(hObject, eventdata, handles);
    end
end
guidata(hObject, handles);
function PixelNo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function SurfPos_Callback(hObject, eventdata, handles)
set(handles.SurfPos,'BackgroundColor',[1,1,1]);
SurfPos=str2double(get(handles.SurfPos,'String'));
if SurfPos<0
    set(handles.SurfPos,'String',0)
elseif SurfPos>str2double(get(handles.ProfileLength,'String'))
    set(handles.SurfPos,'String',str2double(get(handles.ProfileLength,'String'))/2)
end

switch handles.CurrentPlot
    case 'Image'
    case 'Align'
        delete(handles.SurfLine2);
        SurfPos=str2double(get(handles.SurfPos,'String'));
        handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
            'LineWidth',2,'Color',[0 1 1],'LineStyle',':');
    case 'Mask'
        delete(handles.SurfLine2);
        SurfPos=str2double(get(handles.SurfPos,'String'));
        handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
            'LineWidth',2,'Color',[0 1 1],'LineStyle',':');
    case 'Generate'
        delete(handles.SurfLine2);
        SurfPos=str2double(get(handles.SurfPos,'String'));
        handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
            'LineWidth',2,'Color',[0 1 1],'LineStyle',':');
    case 'PlotButton'
        [handles]=PlotButton_Callback(hObject, eventdata, handles);
    case 'DataPlotOnly'
        [handles]=DataPlotOnly(hObject, eventdata, handles);
    case 'Fit'
        [handles]=PlotButton_Callback(hObject, eventdata, handles);
%         FitButton_Callback(hObject, eventdata, handles)
end
%
% if isfield(handles,'Masking');
%     if handles.Masking==1
%         if handles.PlotType==1
%             if isfield(handles,'SurfLine2')
%                 delete(handles.SurfLine2);
%             end
%             SurfPos=str2double(get(handles.SurfPos,'String'));
%             handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
%                 'LineWidth',2,'Color',[.1 1 .1],'LineStyle',':');
%         elseif isfield(handles,'ProfileData_Or')
%             [handles]=PlotButton_Callback(hObject, eventdata, handles);
%         end
%     else
%         if handles.PlotType==3
%             [handles]=DataPlotOnly(hObject, eventdata, handles);
%         end
%     end
% end
guidata(hObject, handles);
function SurfPos_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [handles]=DataPlotOnly(hObject, eventdata, handles)
handles.CurrentPlot='DataPlotOnly';
handles.PlotType=3;
handles.Masking=0;
handles.OverPlot=0;

set(handles.Ylimit,'Position',handles.GraphPlotPosYlim);
set(handles.Xlimit,'Position',handles.GraphPlotPosXlim);
set(handles.Ylimit,'Visible','on');
set(handles.Xlimit,'Visible','on');

ProfileLength=str2double(get(handles.ProfileLength, 'String'));
SurfPos=str2double(get(handles.SurfPos, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
DataStep=ProfileLength*1e-6/(PixelNo-1);
handles.Xdata_Or=0:DataStep:ProfileLength*1e-6;

[val_start handles.SurfPos_idx] = min(abs(handles.Xdata_Or-SurfPos*1e-6));
handles.Xdata=handles.Xdata_Or(1:end-handles.SurfPos_idx+1);
handles.ProfileData=handles.ProfileData_Or(handles.SurfPos_idx:end);

plot(handles.Xdata*1e6,handles.ProfileData);
Xlimit_Callback(hObject, eventdata, handles);
set(handles.Ylimit,'String',roundsf(max(handles.ProfileData_Or)*1.1,2,'ceil'));
if ~isfield(handles,'PlotFlag')
    Ylimit_Callback(hObject, eventdata, handles);
end
legend('Data');
xlabel('Depth /\mum');ylabel('Isotopic Fraction');
% set(gca,'position',[68 62 570 465]); %pixels
set(gca,'position',handles.GraphPlotPos); %norm
guidata(hObject, handles);

function C_bg_Callback(hObject, eventdata, handles)
%PlotButton_Callback(hObject, eventdata, handles);
if handles.PlotType==1
    if isfield(handles,'C_bgLine')
        delete(handles.C_bgLine);
    end
    C_bg=str2double(get(handles.C_bg,'String'));
    handles.C_bgLine=line([-5000 5000],[C_bg, C_bg],...
        'LineWidth',2,'Color',[0 1 0],'LineStyle',':');
end
set(handles.C_bg,'BackgroundColor',[1,1,1]);
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);
function C_bg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function C_gas_Callback(hObject, eventdata, handles)
%PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function C_gas_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function D1_Callback(hObject, eventdata, handles)
inc=get(handles.D1slider,'Value');
set(handles.D1slider,'Value',0);
D1=str2double(get(handles.D1, 'String'));
% D1=roundsf(D1*(1+inc*0.1),3,'round');
D1=roundsf(D1+inc*10^(floor(log10(D1)-1)),3,'round');
set(handles.D1,'String',num2str(abs(D1)));
PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function D1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function D2_Callback(hObject, eventdata, handles)
inc=get(handles.D2slider,'Value');
set(handles.D2slider,'Value',0);
D2=str2double(get(handles.D2, 'String'));
% D2=roundsf(D2*(1+inc*0.1),3,'round');
D2=roundsf(D2+inc*10^(floor(log10(D2)-1)),3,'round');
set(handles.D2,'String',num2str(abs(D2)));
PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function D2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function k1_Callback(hObject, eventdata, handles)
inc=get(handles.k1slider,'Value');
set(handles.k1slider,'Value',0);
k1=str2double(get(handles.k1, 'String'));
% k1=roundsf(k1*(1+inc*0.1),3,'round');
k1=roundsf(k1+inc*10^(floor(log10(k1)-1)),3,'round');
set(handles.k1,'String',num2str(abs(k1)));
PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

function k1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function k2_Callback(hObject, eventdata, handles)
inc=get(handles.k2slider,'Value');
set(handles.k2slider,'Value',0);
k2=str2double(get(handles.k2, 'String'));
% k2=roundsf(k2*(1+inc*0.1),3,'round');
k2=roundsf(k2+inc*10^(floor(log10(k2)-1)),3,'round');
set(handles.k2,'String',num2str(abs(k2)));
PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function k2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function t1_Callback(hObject, eventdata, handles)
%PlotButton_Callback(hObject, eventdata, handles);
set(handles.t1,'String',abs(str2double(get(handles.t1,'String'))));
guidata(hObject, handles);
function t1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function t2_Callback(hObject, eventdata, handles)
%PlotButton_Callback(hObject, eventdata, handles);
set(handles.t2,'String',abs(str2double(get(handles.t2,'String'))));
if str2num(get(handles.t2,'String')) == 0
    set(handles.k2FitCheck,'Value',0);
end
guidata(hObject, handles);
function t2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Xlimit_Callback(hObject, eventdata, handles)
set(handles.Xlimit,'Visible','on');
%  
if length(handles.CurrentPlot)==10;
    %     if length(get(handles.XprimeTog,'State'))==2
    %         handles.Xlimit_num=3;
    %     else
    handles.Xlimit_num=str2double(get(handles.Xlimit, 'String'));
    %     end
else
    handles.Xlimit_num=str2double(get(handles.Xlimit, 'String'));
end
xlim([0 handles.Xlimit_num]);
set(gca,'XTick',linspace(0,handles.Xlimit_num,11))
guidata(hObject, handles);
function Xlimit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Ylimit_Callback(hObject, eventdata, handles)
set(handles.Ylimit,'Visible','on');
handles.Ylimit_max=str2double(get(handles.Ylimit, 'String'));
if length(get(handles.GBplot,'State'))==2;
    ylim_min=handles.Ylimit_min;
else
    if get(handles.ErrorCheck, 'Value')==1;
        %     set(handles.ErrorCheck, 'Value',0); %%%possibly wrong
        ylim_min=handles.Ylimit_min;
    else
        ylim_min=0;
        %         set(gca,'YTick',linspace(0,handles.Ylimit_max,11))
    end
end
if ylim_min==0 && handles.Ylimit_max==0
    handles.Ylimit_max=1;
end
if ~isfinite(handles.Ylimit_max)
    handles.Ylimit_max=1;
end
% if isfield(handles,'PlotButtonFlag');
%     if handles.PlotButtonFlag==1
%         if isfinite(handles.Ylimit_min+handles.Ylimit_max)
%             %             ymin=handles.Ylimit_min;            ymin
%             %             ymax=handles.Ylimit_max;            ymax
%             ylim([handles.Ylimit_min handles.Ylimit_max]);
%         end
%         handles.PlotButtonFlag=0;
%     else
%     end;
% end
ylim([ylim_min handles.Ylimit_max]);

guidata(hObject, handles);
function Ylimit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ErrorCheck_Callback(hObject, eventdata, handles)
[handles]=PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

%%% PlotTown
%% --- Executes on button press in PlotButton.
function [handles]=PlotButton_Callback(hObject, eventdata, handles);%fun
handles.CurrentPlot='PlotButton';
handles.PlotType=2;
handles.OverPlot=0;
handles.PlotFlag=1;
if isfield(handles,'Masking')
    handles.Masking=0;
end
if isfield(handles, 'colo')
    delete(handles.colo);
end

set(handles.Ylimit,'Position',handles.GraphPlotPosYlim);
set(handles.Xlimit,'Position',handles.GraphPlotPosXlim);
set(handles.Ylimit,'Visible','on');
set(handles.Xlimit,'Visible','on');

C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
D1=str2double(get(handles.D1, 'String'));
D2=str2double(get(handles.D2, 'String'));
k1=str2double(get(handles.k1, 'String'));
k2=str2double(get(handles.k2, 'String'));
t1=str2double(get(handles.t1, 'String'));
t2=str2double(get(handles.t2, 'String'));
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));

dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain
PS=get(handles.PlaneSheet,'Value');
val_mirror=str2double(get(handles.MirrorPlane,'String'));
[val_mirror handles.MP_idx] = min(abs(X-val_mirror*1e-6));
MP=handles.MP_idx;
%Generate profiles
if PS==1;
    %     pro=zeros(size(X));
    if t2==0
        [X,pro]=InDiffsPS_inline(handles);
        ProTypeFlag=4;
    else
        [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP);
        ProTypeFlag=3;
    end
else
    if t2 >0 || t2 == inf
        if k1~=k2 || D1~=D2
            tic;
            [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP);
            %         if length(get(handles.Go,'Visible'))==2
            %             toc
            %             proCNop=pro;
            %             tic
            %             [X,pro,DepthFlag]=BackDiffsCN_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo);
            %             toc
            %             f1=figure(3);plot(X,pro-proCNop);
            %             f2=figure(2);plot(X,pro,X,proCNop);
            %             pause
            %             delete(f1);delete(f2);
            % %              
            %         end
            
            ProTypeFlag=3;
            if pro(end)>1.01*C_bg
                DepthFlag=1;
            end
            set(handles.WarningBox,'Visible','off')
            if DepthFlag==1
                %             set(handles.WarningBox,'Visible','on')
                %             set(handles.WarningBox,'String','Warning: Simulated profile exceeds domain')
            end
        else
            [X,pro]=BackDiffs1k_inline(C_gas,C_bg,D1,k1,t1,t2,ProfileLength,PixelNo,handles);
            ProTypeFlag=2;
            if ~isfinite(sum(pro))
                [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP);
                ProTypeFlag=3;
            end
        end
    else
        if D1==0
            D1=1e-15;
            pro=C_bg*ones(size(X));
        else
            [X,pro]=InDiffs_inline(C_gas,C_bg,D1,k1,t1,ProfileLength,PixelNo,handles);
        end
        ProTypeFlag=1;
        if ~isfinite(sum(pro))
            [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP);
            ProTypeFlag=3;
        end
    end %Various Conditions
end %Plane-Sheet
handles.pro=pro;

handles.X=X;
Fit_Start=str2double(get(handles.Fit_Start, 'String'));
tmp = abs(X-Fit_Start*1e-6);
[val_start handles.Fit_Start_idx] = min(tmp);
Fit_End=str2double(get(handles.Fit_End, 'String'));
tmp = abs(X-Fit_End*1e-6);
[val_end handles.Fit_End_idx] = min(tmp);

%%% Legend
switch ProTypeFlag
    case 1
        ProStr='Crank';
    case 2
        ProStr='Back-Crank';
    case 3
        ProStr='Crank-Nicolson FD Simulation';
    case 4
        ProStr='Plane Sheet';
end
% if t2==0
%     ProStr='Crank';
% else
%     if D1==D2 && k1==k2
%         ProStr='Back-Crank';
%     else
%         ProStr='Crank-Nicolson FD Simulation';
%     end
% end

%%%plots
if ~isfield(handles,'ProfileData') %% Just profile
    if length(get(handles.CprimeY,'State'))==2;
        pro=(pro-C_bg)/(C_gas-C_bg);
        if length(get(handles.GBplot,'State'))==2;
            handles.Ylimit_max=0;
            pro=log(abs(pro));
            handles.Ylimit_min=-20;
            %             handles.Ylimit_min=roundsf(min(pro(isreal(pro))),1,'floor');
            X=((X./ ((D1*t1*3600)^0.5) ).^(6/5))/1e6;
        else
            handles.Ylimit_max=ceil(20*max(pro))/20;
            handles.Ylimit_min=0;
        end
    else
        handles.Ylimit_max=ceil(20*max(pro))/20;
        handles.Ylimit_min=0;
    end
    if length(get(handles.XprimeTog,'State'))==3;
        plot(X*1e6,pro);
    else
        X=roundsf(X/(2*sqrt(D1*(t1+t2)*3600)),3,'round');
        plot(X,pro);
    end
    %     set(handles.Xlimit,'String',roundsf(max(X*1e6),2,'ceil'));
    
    axis normal;
    legend(ProStr);
else  %%%%%% Profile + Data
    [handles]=DataPlotOnly(hObject, eventdata, handles);
    Xdata=handles.Xdata;
    %%% Cprime or IF
    if length(get(handles.CprimeY,'State'))==2; %if is on
        pro=(pro-C_bg)/(C_gas-C_bg);
        data=(handles.ProfileData-C_bg)/(C_gas-C_bg);
    else
        data=handles.ProfileData;
    end
    %Grain boundary
    
    if length(get(handles.GBplot,'State'))==2;
        
        pro=log(pro);
        data=log(abs(data));
        X=((X/((D1*(t1+t2)*3600)^0.5)).^(6/5))/1e6;
        Xdata=((Xdata/((D1*(t1+t2)*3600)^0.5)).^(6/5))/1e6;
        handles.Ylimit_max=ceil(20*max(pro))/20;
        %         set(handles.Xlimit,'String',roundsf(max(X*1e6),2,'ceil'));
        %         handles.Ylimit_min=roundsf(min(min(pro),min(data)),1,'floor');
        %     else
        %         handles.Ylimit_max=ceil(20*max(max(pro),max(data)))/20;
        %                 set(handles.Xlimit,'String',...
        %                     roundsf(1e6*X(find(isfinite(data),1,'last')),2,'ceil')); %%%%changed
        
        
        %          
    end
    if length(get(handles.XprimeTog,'State'))==2;
        X=X/(2*sqrt(D1*(t1+t2)*3600))/1e6;
        Xdata=Xdata/(2*sqrt(D1*(t1+t2)*3600))/1e6;
        %         set(handles.Xlimit,'String',roundsf(1e6*X(find(1-(0<data<1),1,'last')),2,'ceil')); %%%%changed
    end
    prof_err=(pro(1:length(data))-data);
    plot(Xdata*1e6,data,'o',X*1e6,pro,'-r','linewidth',1,'MarkerSize',3);
    %%% Residuals
    if get(handles.ErrorCheck, 'Value')==1;
        if length(get(handles.GBplot,'State'))==2;
            handles.Ylimit_max=ceil(20*max(max(prof_err),-min(min(pro),min(data))/10))/20;
            handles.Ylimit_min=-20;
            %             handles.Ylimit_min=roundsf(min(min(pro),min(data)),1,'floor');
        else
            handles.Ylimit_max=ceil(20*max(max(pro),max(data)))/20;
            handles.Ylimit_min=roundsf(min(min(prof_err),handles.Ylimit_max/-10),1,'floor'); %make it minimum .1*max for equal step
        end
        hold on
        plot(Xdata*1e6,prof_err,'g','linewidth',1);
        legend('Data',ProStr,'Residual'); %add residual
        plot([0 1e4],[0 0],'-k'); %add x-axis
        hold off
    else
        if length(get(handles.GBplot,'State'))==2;
            handles.Ylimit_max=0;
            %             handles.Ylimit_min=roundsf(min(min(pro),min(data)),1,'floor');
            handles.Ylimit_min=-20;
        else
            handles.Ylimit_max=ceil(20*max(max(pro),max(data)))/20;
            handles.Ylimit_min=0;
        end
        legend('Data',ProStr);
    end
    % 
    if handles.Fit_Start_idx < find(isfinite(data),1,'first')
        handles.Fit_Start_idx = find(isfinite(data),1,'first');
        set(handles.Fit_Start,'String',...
            roundsf(1e6*X(handles.Fit_Start_idx),3,'round'));
    end
    
    if handles.Fit_End_idx > find(isfinite(data),1,'last')
        handles.Fit_End_idx = find(isfinite(data),1,'last');
        set(handles.Fit_End,'String',...
            roundsf(1e6*X(handles.Fit_End_idx),3,'round'));
    end
    
    [r2 rmse] = rsquare(handles.ProfileData(handles.Fit_Start_idx:handles.Fit_End_idx),...
        handles.pro(handles.Fit_Start_idx:handles.Fit_End_idx));
    handles.r2=r2;
    handles.rmse=rmse;
    set(handles.Rsq,'String',num2str(roundsf(r2,4,'round')));
    
    [r2All rmseAll] = rsquare(handles.ProfileData(1:end),...
        handles.pro(1:length(handles.ProfileData)));
    set(handles.RsqAll,'String',num2str(roundsf(r2All,4,'round')));
    %     handles.Ylimit_max=ceil(20*max(max(pro),max(data)))/20;
end
set(handles.Ylimit,'String',handles.Ylimit_max);

grid (get(handles.GridLines,'State'))
set(gca,'Xcolor',[0.2 0.2 0.2]); set(gca,'Ycolor',[0.2 0.2 0.2]);
if length(get(handles.GBplot,'State'))==2;
    xlim([0 10])
elseif length(get(handles.XprimeTog,'State'))==2;
    xlim([0 3])
    xlabel('Normalised depth, x''');
    set(handles.Xlimit,'Visible','off');
else
    xlabel('Depth, x /\mum');
    Xlimit_Callback(hObject, eventdata, handles);
end
Ylimit_Callback(hObject, eventdata, handles);
if length(get(handles.CprimeY,'State'))==2;
    if length(get(handles.GBplot,'State'))==2;
        ylabel('ln|C''|');
        xlabel('\eta^{6/5}');
        set(handles.Xlimit,'Visible','off')
        set(handles.Ylimit,'Visible','off')
    else
        ylabel('Normalised _{}^{18}O Fraction, C''');
    end
else
    ylabel('Isotopic Fraction');
end
handles.PlotButtonFlag=1;

% set(gca,'position',[68 62 570 465]); %pixels
set(gca,'position',handles.GraphPlotPos); %norm
% set(gca,'yscale','log');
Fit_End_Callback(hObject, eventdata, handles);%
handles.CurrentPlot='PlotButton';
PlotProperties(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Rsq_CreateFcn(hObject, eventdata, handles)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function MainAxes_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in FitButton.
function FitButton_Callback(hObject, eventdata, handles)
handles.CurrentPlot='Fit';
if ~isfield(handles,'ProfileData_Or')
    [handles]=LoadProfileData_Callback(hObject, eventdata, handles);
end
% Change appearance of Fit_Button

FitCheck(1)=get(handles.D1FitCheck, 'Value');
FitCheck(2)=get(handles.D2FitCheck, 'Value');
FitCheck(3)=get(handles.k1FitCheck, 'Value');
FitCheck(4)=get(handles.k2FitCheck, 'Value');

C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
D1=str2double(get(handles.D1, 'String'));
D2=str2double(get(handles.D2, 'String'));
k1=str2double(get(handles.k1, 'String'));
k2=str2double(get(handles.k2, 'String'));
t1=str2double(get(handles.t1, 'String'));
t2=str2double(get(handles.t2, 'String'));
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
PS=get(handles.PlaneSheet, 'Value');
[handles]=PlotButton_Callback(hObject, eventdata, handles);

Fit_Start=str2double(get(handles.Fit_Start, 'String'));
[val_start handles.Fit_Start_idx] = min(abs(handles.Xdata-Fit_Start*1e-6));
Fit_Start_idx=handles.Fit_Start_idx;

Fit_End=str2double(get(handles.Fit_End, 'String'));
[val_end handles.Fit_End_idx] = min(abs(handles.Xdata-Fit_End*1e-6));
Fit_End_idx=handles.Fit_End_idx;
MP=handles.MP_idx;
% FitRange=str2double(get(handles.FitRange, 'String'));
% FitPoints=str2double(get(handles.FitPoints, 'String'));
%%% Once all variables have been collected, begin sorting

% if FitCheck(2)==1 || FitCheck(4)==1 %Stop process from taking forever
%     if FitPoints>10
%         set(handles.FitPoints, 'String',5);
%         FitPoints=str2double(get(handles.FitPoints, 'String'));
%     end
% end
if ~isfield(handles,'ProfileData')
    [handles]=Data_Callback(hObject, eventdata, handles);
end

%Make the section of profile here
% DataStep=ProfileLength*1e-6/PixelNo;
% [m,n]=size(handles.ProfileData_Or);
% handles.Xdata_Or=0:DataStep:(max(m,n)-1)*DataStep;

DataStep=ProfileLength*1e-6/(PixelNo-1);
handles.Xdata_Or=0:DataStep:ProfileLength*1e-6;

SurfPos=str2double(get(handles.SurfPos, 'String'));
[val_start handles.SurfPos_idx] = min(abs(handles.Xdata_Or-SurfPos*1e-6));
handles.Xdata=handles.Xdata_Or(1:end-handles.SurfPos_idx+1);
handles.ProfileData=handles.ProfileData_Or(handles.SurfPos_idx:end);

if t2==0
    if PS==0 %Plane
        [D1,k1]=AutoFit_Crank(handles);
    else
        [D1,k1]=AutoFit_PlaneSheet(handles);
    end
else
%     if PS==1 %Plane
%         set(handles.FitButton,'ForegroundColor',[0.5,0.5,0.5]);
%         %         set(handles.FitButton,'Enable','off');
%         set(handles.FitButton,'String','Fitting...');pause(0.001);
%         [D1,k1,~,~]= AutoFit_BackDiffs(C_gas,C_bg,D1,k1,D2,k2,t1,t2,handles.ProfileData,...
%             handles.Xdata,Fit_Start_idx,Fit_End_idx,ProfileLength,PixelNo,handles,FitCheck);
%         [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP); %ImageLength
%         [r2 rmse] = rsquare(handles.ProfileData(Fit_Start_idx:Fit_End_idx),...
%             pro(Fit_Start_idx:Fit_End_idx));
%     end
    if FitCheck(2)==0 && FitCheck(4)==0
        [D1,k1]= AutoFit_BackCrank(C_gas,C_bg,D1,k1,t1,t2,handles.ProfileData,...
            handles.Xdata,Fit_Start_idx,Fit_End_idx);
        D2=D1;
        k2=k1;
        %     elseif FitCheck(2)==0 && FitCheck(4)==1
    else
        set(handles.FitButton,'ForegroundColor',[0.5,0.5,0.5]);
        %         set(handles.FitButton,'Enable','off');
        set(handles.FitButton,'String','Fitting...');pause(0.001);
        % Start by using the analytical solution to get a good/quick ball park...
        [D1,k1]= AutoFit_BackCrank(C_gas,C_bg,D1,k1,t1,t2,handles.ProfileData,...
            handles.Xdata,Fit_Start_idx,Fit_End_idx);
        D2=D1;
        k2=k1;
        
        %         [D1,D2,k1,k2,r2]=BackDiffs_fit_inline(...
        %             C_gas,C_bg,D1,D2,k1,k2,t1,t2,handles.ProfileData,FitRange,...
        %             FitPoints,FitCheck,ProfileLength,Fit_Start_idx,Fit_End_idx,handles);
        %         set(handles.Rsq,'String',num2str(r2));
        tic;
        r2=0;
        tim=0;
        i=1;
        r2check(1)=str2double(get(handles.Rsq,'String'));
        %          
        while r2<0.999 && tim<20;
            [D1,k1,D2,k2]= AutoFit_BackDiffs(C_gas,C_bg,D1,k1,D2,k2,t1,t2,handles.ProfileData,...
                handles.Xdata,Fit_Start_idx,Fit_End_idx,ProfileLength,PixelNo,handles,FitCheck);
            [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP); %ImageLength
            [r2 rmse] = rsquare(handles.ProfileData(Fit_Start_idx:Fit_End_idx),...
                pro(Fit_Start_idx:Fit_End_idx));
            
            tim=toc+tim;
            i=i+1;
            r2check(i)=r2;
            if r2check(i)==r2check(i-1);
                tim=20;
            end
            
        end
        disp(['Total time for fit ',num2str(round(toc)),' seconds.'])
        %         set(handles.FitButton,'Enable','on');
    end
end
D1=roundsf(D1,3,'round');
D2=roundsf(D2,3,'round');
k1=roundsf(k1,3,'round');
k2=roundsf(k2,3,'round');
%%% Set all the newly fitted values
set(handles.D1,'String',num2str(D1));
set(handles.D2,'String',num2str(D2));
set(handles.k1,'String',num2str(k1));
set(handles.k2,'String',num2str(k2));

FitPlot(hObject, eventdata, handles);

if isfield(handles,'O18pathname')
    handles.SaveName=[handles.PathName,handles.FileName(1:end-4)];
else
    handles.SaveName=[handles.PathName,handles.FileName(1:end-4)];
end
set(handles.FitButton,'ForegroundColor',[0,0,0]);
set(handles.FitButton,'String','Fit');
handles.CurrentPlot='Fit';
guidata(hObject, handles);

function FitPlot(hObject, eventdata, handles)
handles.CurrentPlot='Fit';
[handles]=PlotButton_Callback(hObject, eventdata, handles);
D1=str2double(get(handles.D1,'string'));
t1=str2double(get(handles.t1,'string'));
t2=str2double(get(handles.t2,'string'));
t_tot=(t1+t2)*3600;
% set(handles.RMSres,'String',num2str(handles.rmse));
hold on;
Fit_Start=str2double(get(handles.Fit_Start,'String'));
Fit_End=str2double(get(handles.Fit_End,'String'));
if length(get(handles.GBplot,'State'))==3;
    line([Fit_Start Fit_Start],[-40 1],'Color','k','LineStyle','--');
    line([Fit_End Fit_End],[-40 1],'Color','k','LineStyle','--');
else
    line(([Fit_Start Fit_Start]/ ((D1*t_tot*3600)^0.5)*1e-6 ).^(6/5),[-40 1],'Color','k','LineStyle','--');
    line(([Fit_End Fit_End]/ ((D1*t_tot*3600)^0.5)*1e-6 ).^(6/5),[-40 1],'Color','k','LineStyle','--');
end
hold off;
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);

function SaveFit_Callback(hObject, eventdata, handles)
[handles]=PlotButton_Callback(hObject, eventdata, handles);
SaveName=get(handles.ProfDataSaveName,'String');
% [handles.FileName, handles.PathName,~] = ...
%         uigetfile([handles.PathName,'*.txt'],'Oxygen 18 file');%(FILTERSPEC, TITLE);
if isfield(handles,'PathName')
    [FileName,PathName]=uiputfile([handles.PathName,SaveName(1:end-4),'_Profiles.txt']);
else
    [FileName,PathName]=uiputfile('Profiles.txt');
end

fid=fopen([PathName,FileName],'wt');

ProData=zeros(size(handles.X));
if isfield(handles,'ProfileData')
    ProData(1:length(handles.ProfileData))=handles.ProfileData;
end
fprintf(fid,'FileName = ;%s;\n',handles.SaveName);
fprintf(fid,'Image_Length = ;%s; [um]\n',get(handles.ProfileLength,'String'));
fprintf(fid,'Pixel_Number = ;%s; [pixels]\n',get(handles.PixelNo,'String'));
fprintf(fid,'Surface_Position = ;%s; [um]\n',get(handles.SurfPos,'String'));
fprintf(fid,'Align_Angle = ;%s; [deg]\n',get(handles.AlignAngle,'String'));
fprintf(fid,'Mask_Thresh = ;%s; []\n',get(handles.MaskThresh,'String'));

D1=str2double(get(handles.D1,'String'));
D2=str2double(get(handles.D2,'String'));
k1=str2double(get(handles.k1,'String'));
k2=str2double(get(handles.k2,'String'));
t1=str2double(get(handles.t1,'String'));
t2=str2double(get(handles.t2,'String'));
C_bg=str2double(get(handles.C_bg,'String'));
C_gas=str2double(get(handles.C_gas,'String'));

fprintf(fid,'D1 = ;%g; [m2/s]\n',D1);
fprintf(fid,'D2 = ;%g; [m2/s]\n',D2);
fprintf(fid,'k1 = ;%s; [m/s]\n',k1);
fprintf(fid,'k2 = ;%g; [m/s]\n',k2);
fprintf(fid,'t1 = ;%g; [hours]\n',t1);
fprintf(fid,'t2 = ;%g; [hours]\n',t2);
fprintf(fid,'C_bg = ;%g; []\n',C_bg);
fprintf(fid,'C_gas = ;%g; []\n',C_gas);
fprintf(fid,'R_Squared = ;%s; []\n',get(handles.RsqAll,'String'));
fprintf(fid,'R_Squared ROI= ;%s; []\n',get(handles.Rsq,'String'));
fprintf(fid,'Fit_Start = ;%s; [um]\n',get(handles.Fit_Start,'String'));
fprintf(fid,'Fit_End = ;%s; [um]\n',get(handles.Fit_End,'String'));

last = find(1-(0<ProData<1),1,'last')
Pros(:,1)=handles.X(1:last);
Pros(:,2)=ProData(1:last);
Pros(:,3)=handles.pro(1:last);
Pros(:,4)=ProData(1:last)-handles.pro(1:last);

Pros(:,5)=handles.X(1:last)./(4*D1*(t1+t2)*3600)^0.5;
Pros(:,6)=(ProData(1:last)-C_bg)/(C_gas-C_bg);
Pros(:,7)=(handles.pro(1:last)-C_bg)/(C_gas-C_bg);
Pros(:,8)=Pros(:,6)-Pros(:,7);
fprintf(fid,['Depth,x /m;Isotopic Fraction Data;Fitted Profile;Residual;;',...
    'Normalised Depth, x''=x/(4 D1 t)^1/2);Normalised Data;Normalised Fit;Residual\n']);
fprintf(fid,'%d;%d;%d;%d;;%d;%d;%d;%d\n',Pros');
fclose(fid);
% [FileName,PathName]=uiputfile([handles.PathName,SaveName(1:end-4),'_Fit.txt']);
% fid=fopen([PathName,FileName],'wt');
% fprintf(fid,'%d\n',handles.pro);
% fclose(fid);
% fid=fopen([PathName,FileName(1:end-4),'_ProfileData.txt'],'wt');
% fprintf(fid,'%d\n',handles.ProfileData);
% fclose(fid);
guidata(hObject, handles);

% function SaveVars_Callback(hObject, eventdata, handles)
% SaveName=get(handles.ProfDataSaveName,'String');
% [FileName,PathName]=uiputfile([handles.PathName,SaveName(1:end-4),'_Variables.txt']);
% fid=fopen([PathName,FileName],'wt');
% fprintf(fid,'FileName = %s\n',handles.SaveName);
% fprintf(fid,'Image_Length = %s [um]\n',get(handles.ProfileLength,'String'));
% fprintf(fid,'Pixel_Number = %s [pixels]\n',get(handles.PixelNo,'String'));
% fprintf(fid,'Surface_Position = %s [um]\n',get(handles.SurfPos,'String'));
%
% fprintf(fid,'D1 = %s [m2/s]\n',get(handles.D1,'String'));
% fprintf(fid,'D2 = %s [m2/s]\n',get(handles.D2,'String'));
% fprintf(fid,'k1 = %s [m/s]\n',get(handles.k1,'String'));
% fprintf(fid,'k2 = %s [m/s]\n',get(handles.k2,'String'));
% fprintf(fid,'t1 = %s [hours]\n',get(handles.t1,'String'));
% fprintf(fid,'t2 = %s [hours]\n',get(handles.t2,'String'));
% fprintf(fid,'C_bg = %s []\n',get(handles.C_bg,'String'));
% fprintf(fid,'C_gas = %s []\n',get(handles.C_gas,'String'));
% fprintf(fid,'R_Squared = %s []\n',get(handles.RsqAll,'String'));
% % fprintf(fid,'RMS residual = %s []\n',get(handles.rmse,'String'));
% fprintf(fid,'Fit_Start = %s [um]\n',get(handles.Fit_Start,'String'));
% fprintf(fid,'Fit_End = %s [um]\n',get(handles.Fit_End,'String'));
% fprintf(fid,'R_Squared ROI= %s []\n',get(handles.Rsq,'String'));
% fclose(fid);
% guidata(hObject, handles);

function D1FitCheck_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function D2FitCheck_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function k1FitCheck_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function k2FitCheck_Callback(hObject, eventdata, handles)
guidata(hObject, handles);


function Fit_Start_Callback(hObject, eventdata, handles)
set(handles.Fit_Start,'BackgroundColor',[1,1,1]);
% Fit_Start=str2double(get(handles.Fit_Start, 'String'));
% tmp = abs(X-Fit_Start*1e-6);
% [val_start handles.Fit_Start_idx] = min(tmp);
guidata(hObject, handles);
function Fit_Start_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Fit_End_Callback(hObject, eventdata, handles)
set(handles.Fit_End,'BackgroundColor',[1,1,1]);
t1=str2double(get(handles.t1,'String'));
t2=str2double(get(handles.t2,'String'));
t_tot=t1+t2;
D1=str2double(get(handles.D1,'String'));
set(handles.xStar_d,'String',...
    roundsf((str2double(get(handles.Fit_End,'String'))*...
    10^-6)/(2*sqrt(D1*(t_tot*3600))),3,'round'));
% Fit_End=str2double(get(handles.Fit_End, 'String'));
% tmp = abs(X-Fit_End*1e-6);
% [val_end handles.Fit_End_idx] = min(tmp);

guidata(hObject, handles);
function Fit_End_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FitRange_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function FitRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function FitPoints_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function FitPoints_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function D1slider_Callback(hObject, eventdata, handles);
D1_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function D1slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function D2slider_Callback(hObject, eventdata, handles)
D2_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function D2slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function k1slider_Callback(hObject, eventdata, handles)
k1_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function k1slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function k2slider_Callback(hObject, eventdata, handles)
k2_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function k2slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in LoadO16image.
function [handles]=LoadO16image_Callback(hObject, eventdata, handles)
handles.ProLen_Or=500;
handles.O16Flag=1;
if isfield(handles,'FileName')
    [handles.FileName, handles.PathName,~] = ...
        uigetfile([handles.PathName,'*.txt'],'Oxygen 16 file');%(FILTERSPEC, TITLE);
else
    [handles.FileName, handles.PathName,~] = ...
        uigetfile('*.txt','Oxygen 16 file');%(FILTERSPEC, TITLE);
end
if handles.PathName==0
    handles.O16Flag=0;
else
    [Image]=Loader2D(handles);
    handles.O16image=Image;
    handles.O16image_Or=Image;
    PlotO16image_Callback(hObject, eventdata, handles)
end
if isfield(handles,'O18Flag')
    if handles.O18Flag==1;
        if size(handles.O16image_Or)==size(handles.O18image_Or)
            handles.O16pO18image_Or=handles.O16image_Or+handles.O18image_Or;
            handles.O16pO18image=handles.O16pO18image_Or;
            NormO18image=(handles.O18image./handles.O16pO18image);
            NormO18image(isnan(NormO18image))=0;
            NormO18image(isinf(NormO18image))=0;
            handles.NormO18image_Or=NormO18image;
            handles.NormO18image=handles.NormO18image_Or;
        end
    end
end
guidata(hObject, handles);

% --- Executes on button press in LoadO18image.
function [handles]=LoadO18image_Callback(hObject, eventdata, handles)
handles.CurrentPlot='Image';
handles.O18Flag=1;
if isfield(handles,'FileName')
    [handles.FileName, handles.PathName,~] = ...
        uigetfile([handles.PathName,'*.txt'],'Oxygen 18 file');%(FILTERSPEC, TITLE);
else
    [handles.FileName, handles.PathName,~] = ...
        uigetfile('*.txt','Oxygen 18 file');%(FILTERSPEC, TITLE);
end
if handles.PathName==0
    handles.O18Flag=0;
else
    [Image]=Loader2D(handles);
    handles.O18image=Image;
    handles.O18image_Or=Image;
    PlotO18image_Callback(hObject, eventdata, handles)
end
if isfield(handles,'O16Flag')
    if handles.O16Flag==1;
        if size(handles.O16image_Or)==size(handles.O18image_Or)
            handles.O16pO18image_Or=handles.O16image_Or+handles.O18image_Or;
            handles.O16pO18image=handles.O16pO18image_Or;
            NormO18image=(handles.O18image./handles.O16pO18image);
            NormO18image(isnan(NormO18image))=0;
            NormO18image(isinf(NormO18image))=0;
            handles.NormO18image_Or=NormO18image;
            handles.NormO18image=handles.NormO18image_Or;
        end
    end
end
handles.SaveName_Or=[handles.FileName(1:end-4),...
    '_Norm18Prof.txt'];
set(handles.ProfDataSaveName, 'String',handles.FileName);
guidata(hObject, handles);

function [Image]=Loader2D(handles)
ext=handles.FileName(find(handles.FileName=='.',1,'last'):end);
if ext=='.txt'
    if ispc
        Image=dlmread([handles.PathName,'\',handles.FileName], ' ', 9, 2);
    else
        Image=dlmread([handles.PathName,handles.FileName], ' ', 9, 2);
    end
    PixelNo = sqrt(length(Image));
    set(handles.PixelNo,'String',num2str(PixelNo));
    Image=transpose(reshape(Image,PixelNo,PixelNo));
    
    %%%%%%%%%%%%%%%%%
    if ispc
        fileID = fopen([handles.PathName,'\',handles.FileName]);
    else
        fileID = fopen([handles.PathName,handles.FileName]);
    end
    LenStr=textscan(fileID,'%s',4,'Delimiter','\n');
    LenStr=LenStr{1,1}{3,1};
    LenStr=roundsf(str2double(LenStr(18:30)),3,'round');
    fclose(fileID);
    set(handles.ProfileLength,'String',LenStr);
    set(handles.Fit_End,'String',roundsf(...
        0.5*str2double(get(handles.ProfileLength,'String')),3,'floor'));
    set(handles.ProfileLength,'BackgroundColor',[0,1,0]);
    set(handles.Xlimit,'String',LenStr);
    %     ProfileLength=roundsf(dlmread([handles.PathName,'\',...
    %         handles.FileName], ':', [2 1 2 1]),4,'round');
    %     set(handles.ProfileLength,'String',ProfileLength);
    %     set(handles.ProfileLength,'BackgroundColor',[0,1,0]);
    %     set(handles.Xlimit,'String',ProfileLength);
    
    
    % elseif ext=='.tif' %otherwise people will think this is data
    %     Image=imread([handles.PathName,'\',handles.FileName]);
    %     Image=double(Image);
    %     PixelNo=length(Image);
elseif ext=='.ext'
    fileID = fopen([handles.PathName,handles.FileName]);
    tline = fgets(fileID);
    C = textscan(fileID,'%f','HeaderLines',1);
    H=textscan(tline, '%s %s %s %f %s %f %f','Delimiter',' ');
    rows=H{1,4};
    Image=vec2mat(C{1,1},rows);
    Image=double(Image);
    PixelNo=length(Image);
    fclose(fileID);
end
set(handles.PixelNo,'String',PixelNo);
set(handles.PixelNo,'BackgroundColor',[0,1,0]);



function PlotO16image_Callback(hObject, eventdata, handles)
set(handles.Xlimit,'Visible','off');
set(handles.Ylimit,'Visible','off');
if ~isfield(handles,'O16Flag')
    [handles]=LoadO16image_Callback(hObject, eventdata, handles);
end
% axes(handles.MainAxes)
imagesc(handles.O16image_Or);colormap(hot);
handles.colo=colorbar; ylabel(handles.colo,'Counts');
axis square; set(gca, 'XTick', [],'YTick', []);

% set(gca,'position',[90    60   465   465]); %pixels
PlotProperties(hObject, eventdata, handles)
set(gca,'position',handles.ImagePlotPos); %norm

function PlotO18image_Callback(hObject, eventdata, handles)
set(handles.Xlimit,'Visible','off');
set(handles.Ylimit,'Visible','off');
if ~isfield(handles,'O18Flag')
    [handles]=LoadO18image_Callback(hObject, eventdata, handles);
end
% axes(handles.MainAxes)
imagesc(handles.O18image_Or);colormap(hot);colo=colorbar;ylabel(colo,'Counts'); axis square; set(gca, 'XTick', [],'YTick', []);
% set(gca,'position',[90    60   465   465]);%pix
PlotProperties(hObject, eventdata, handles)
set(gca,'position',handles.ImagePlotPos); %norm


function [handles]=PlotO16pO18image_Callback(hObject, eventdata, handles)
set(handles.Xlimit,'Visible','off');
set(handles.Ylimit,'Visible','off');
if ~isfield(handles,'O16Flag')
    [handles]=LoadO16image_Callback(hObject, eventdata, handles);
else
    if handles.O16Flag==0
        [handles]=LoadO16image_Callback(hObject, eventdata, handles);
    end
end
if ~isfield(handles,'O18Flag')
    [handles]=LoadO18image_Callback(hObject, eventdata, handles);
else
    if handles.O18Flag==0
        [handles]=LoadO18image_Callback(hObject, eventdata, handles);
    end
end
% axes(handles.MainAxes)
imagesc(handles.O16pO18image_Or);colormap(hot);colo=colorbar;ylabel(colo,'Counts'); axis square; set(gca, 'XTick', [],'YTick', []);
%handles.N_O18image=handles.O18image./handles.O16pO18image;
%imagesc((handles.N_O18image));colormap(hot);colo=colorbar;ylabel(colo,'Counts'); axis square; set(gca, 'XTick', [],'YTick', []);
% set(gca,'position',[90    60   465   465]);%pix
set(gca,'position',handles.ImagePlotPos); %norm
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes on button press in PlotNormO18image.
function PlotNormO18image_Callback(hObject, eventdata, handles)
set(handles.Xlimit,'Visible','off');
set(handles.Ylimit,'Visible','off');
if ~isfield(handles,'O16Flag')
    [handles]=LoadO16image_Callback(hObject, eventdata, handles);
else
    if handles.O16Flag==0
        [handles]=LoadO16image_Callback(hObject, eventdata, handles);
    end
end
if ~isfield(handles,'O18Flag')
    [handles]=LoadO18image_Callback(hObject, eventdata, handles);
else
    if handles.O18Flag==0
        [handles]=LoadO18image_Callback(hObject, eventdata, handles);
    end
end
% axes(handles.MainAxes)
imagesc(handles.NormO18image_Or);colormap(hot);colo=colorbar;ylabel(colo,'Counts'); axis square; set(gca, 'XTick', [],'YTick', []);
%handles.N_O18image=handles.O18image./handles.O16pO18image;
%imagesc((handles.N_O18image));colormap(hot);colo=colorbar;ylabel(colo,'Counts'); axis square; set(gca, 'XTick', [],'YTick', []);

% set(gca,'position',[90    60   465   465]); %pixels
set(gca,'position',handles.ImagePlotPos); %norm
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes on button press in Align.
function [handles]=Align_Callback(hObject, eventdata, handles)
handles.CurrentPlot='Align';
handles.PlotType=1;
handles.Masking=0;
handles.Aligning=1;
AlignMode=get(handles.AlignMode, 'Value');
if length(get(handles.Advanced,'Visible'))==2
    Go=1;
else
    Go=0;
end
% set(handles.Xlimit,'Visible','off');
if ~isfield(handles,'O16Flag')
    [handles]=PlotO16pO18image_Callback(hObject, eventdata, handles);
end
handles.AlignFlag=1;
set(handles.Align,'ForegroundColor',[0.5,0.5,0.5]);
% set(handles.Align,'Enable','off'); %Don't turn it off, just make it grey
set(handles.Align,'String','Aligning...');pause(0.001);

[a,b]=size(handles.O16pO18image_Or);
Mask_Or=ones(a,b);
CoarseStep=10;
for theta=0:CoarseStep:360-CoarseStep
    i=1+theta/CoarseStep;
    O16pO18image_r=imrotate(handles.O16pO18image_Or,theta,'bilinear','crop'); %,'crop'
    mask_r=imrotate(Mask_Or,theta,'nearest','crop');
    sum_mod=sum(O16pO18image_r)./sum(mask_r);
    
    if Go==1
        sumlog_mod=sum(log(O16pO18image_r))./sum(mask_r);
        logsum_mod=log(sum(O16pO18image_r)./sum(mask_r));
        Max_LogGrad_C(i)=max(diff(logsum_mod));
        Sum_LogGrad_C(i)=sum(diff(log(sum_mod)));
        minsum(i)=min(sum_mod);
        Max_Grad_C(i)=max(diff(sum_mod));
        Sum_Grad_C(i)=sum(diff(sum_mod));
    elseif AlignMode==5
        sumlog_mod=sum(log(O16pO18image_r))./sum(mask_r);
        logsum_mod=log(sum(O16pO18image_r)./sum(mask_r));
        Max_LogGrad_C(i)=max(diff(logsum_mod));
        Sum_LogGrad_C(i)=sum(diff(log(sum_mod)));
    else
        minsum(i)=min(sum_mod);
        Max_Grad_C(i)=max(diff(sum_mod));
        Sum_Grad_C(i)=sum(diff(sum_mod));
    end
    
end
Theta_C=[0:CoarseStep:360-CoarseStep];

%% Go mode plotting
if length(get(handles.Advanced,'Visible'))==2
    set(handles.Xlimit,'Visible','off');
    %     a = annotation(gcf,'textbox',...
    %         [0.35 0.6 0.2 0.12],...
    %         'String','Press enter to continue',...
    %         'FontSize',20);
    %     plot(Theta_C,medfilt2(Max_Grad_C,[1 3])/max(medfilt2(Max_Grad_C,[1 3])),...
    %         Theta_C,Sum_Grad_C/max(Sum_Grad_C));
    %     legend('medfil(Max. Grad.)','Total Grad.');
    %     set(gca,'position',handles.GraphPlotPos);
    %     pause;
    plot(Theta_C,Max_Grad_C/max(Max_Grad_C),...
        Theta_C, minsum./max(minsum),...
        Theta_C,Sum_Grad_C/max(Sum_Grad_C),...
        Theta_C,Max_Grad_C.*Sum_Grad_C/(max(Max_Grad_C)*max(Sum_Grad_C)),...
        Theta_C,Max_LogGrad_C./max(Max_LogGrad_C(isfinite(Max_LogGrad_C))),...
        Theta_C,Sum_LogGrad_C./max(Sum_LogGrad_C(isfinite(Sum_LogGrad_C))));
    legend('Max. Grad.','MinSum','Total Grad.','Combo. Grad','Max LogGrad_C','Sum LogGrad_C');
    xlabel('Angle');set(gca,'position',handles.GraphPlotPos);
    %     f = warndlg('This is a warning.', 'A Warning Dialog');
    %     waitfor(f);
    pause(.1);
    %     delete(a);
    set(handles.Xlimit,'Visible','on');
end

%% Back to things
switch AlignMode
    case 1
        [~, CoarseMaxVal_idx]=max(Max_Grad_C.*Sum_Grad_C);
    case 2
        [CoarseMaxVal CoarseMaxVal_idx]=max(Max_Grad_C);
    case 3
        [CoarseMaxVal CoarseMaxVal_idx]=max(Sum_Grad_C);
    case 4
        [CoarseMaxVal CoarseMaxVal_idx]=min(minsum);
    case 5
        [CoarseMaxVal CoarseMaxVal_idx]=max(Max_LogGrad_C);
end
FA_init=Theta_C(CoarseMaxVal_idx);
% FA_init=(CoarseMaxVal_idx-1)*CoarseStep;
i=1;
AngleStep_Fine=0.2;
Wind=40;
Theta_F  =FA_init-Wind/2 : AngleStep_Fine : FA_init+Wind/2;
for theta=FA_init-Wind/2 : AngleStep_Fine : FA_init+Wind/2
    O16pO18image_r=imrotate(handles.O16pO18image_Or,theta,'bilinear','crop'); %,'crop'
    mask_r=imrotate(Mask_Or,theta,'nearest','crop');
    
    sum_mod=sum(O16pO18image_r)./sum(mask_r);
    minsum(i)=min(sum_mod);
    [Max_Grad_F(i,1) Max_Grad_F(i,2)]=max(diff(sum_mod));
    Sum_Grad_F(i)=sum(diff(sum_mod));
    
    
    sumlog_mod=sum(log(O16pO18image_r))./sum(mask_r);
    logsum_mod=log(sum(O16pO18image_r)./sum(mask_r));
    Max_LogGrad_F(i)=max(diff(logsum_mod));
    Sum_LogGrad_F(i)=sum(diff(log(sum_mod)));
    i=i+1;
end
%% Go mode plotting
if length(get(handles.Advanced,'Visible'))==2
    set(handles.Xlimit,'Visible','off');
    a = annotation(gcf,'textbox',...
        [0.35 0.6 0.2 0.12],...
        'String','Press enter to continue',...
        'FontSize',20);
    plot(Theta_F,Max_Grad_F(:,1)/max(Max_Grad_F(:,1)),...
        Theta_F, minsum./max(minsum),...
        Theta_F, Sum_Grad_F/max(Sum_Grad_F),...
        Theta_F,Sum_Grad_F'.*Max_Grad_F(:,1)/(max(Sum_Grad_F*max(Max_Grad_F(:,1)))),...
        Theta_F,Max_LogGrad_F./max(Max_LogGrad_F(isfinite(Max_LogGrad_F))),...
        Theta_F,Sum_LogGrad_F./max(Sum_LogGrad_F(isfinite(Sum_LogGrad_F))));
    legend('Max. Grad.','MinSum','Total Grad.','Combo. Grad','Max LogGrad','Sum LogGrad');
    xlabel('Angle');
    handles.GraphPlotPos=[.07 .09 .67 .71];%norm
    %     plot(Theta_C,Max_Grad_C/max(Max_Grad_C),...
    %         Theta_C, minsum./max(minsum),...
    %         Theta_C,Sum_Grad_C/max(Sum_Grad_C),...
    %         Theta_C,Max_Grad_C.*Sum_Grad_C/(max(Max_Grad_C)*max(Sum_Grad_C)),...
    %         Theta_C,Max_LogGrad_C./max(Max_LogGrad_C(isfinite(Max_LogGrad_C))),...
    %         Theta_C,Sum_LogGrad_C./max(Sum_LogGrad_C(isfinite(Sum_LogGrad_C))));
    %     legend('Max. Grad.','MinSum','Total Grad.','Combo. Grad','Max LogGrad_C','Sum LogGrad_C');
    %     xlabel('Angle');
    delete(a)
    set(handles.Xlimit,'Visible','on');
    %      
end

%% Back to align
switch AlignMode
    case 1;
        [Max_Grad_F_val Max_Grad_F_idx]=max(Sum_Grad_F'.*Max_Grad_F(:,1));
    case 2;
        [Max_Grad_F_val Max_Grad_F_idx]=max(Max_Grad_F(:,1));
    case 3;
        [Max_Grad_F_val Max_Grad_F_idx]=max(Sum_Grad_F);%Will also tell you where the surf is
    case 4;
        [Max_Grad_F_val Max_Grad_F_idx]=min(minsum);
    case 5
        [Max_Grad_F_val Max_Grad_F_idx]=max(Max_LogGrad_F);
end

AlignAngle=Theta_F(Max_Grad_F_idx);

handles.AlignAngle_FiAl_num=AlignAngle;
handles.AlignAngle_CoAl_num=round(AlignAngle/90)*90;
% handles.Mask_FiAl_Basic=imrotate(Mask_Or,AlignAngle,'nearest','crop');
% handles.Mask_CoAl_Basic=imrotate(Mask_Or,handles.AlignAngle_CoAl_num,'nearest','crop');

set(handles.AlignAngle,'String',num2str(AlignAngle));
[handles]=AlignAngle_Callback(hObject, eventdata, handles);
set(handles.AlignAngle,'BackgroundColor',[0 1 0]);
set(handles.Align,'ForegroundColor',[0,0,0]);
set(handles.Align,'String','Align');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch AlignMode
%     case 1
%         fun = @(theta)...
%             -sum(diff(...
%             sum(imrotate(handles.O16pO18image_Or,theta,'bilinear','crop'))./...
%             sum(imrotate(Mask_Or,theta,'nearest','crop'))))...
%             *max(diff(...
%             sum(imrotate(handles.O16pO18image_Or,theta,'bilinear','crop'))./...
%             sum(imrotate(Mask_Or,theta,'nearest','crop'))));
%         [theta,fminres] = fminbnd(fun,0,360);%,optimset('TolFun',1e-2,'TolX',1e-1));
%         theta
%     case 2
%         fun = @(theta)...
%             double(-max(diff(...
%             sum(imrotate(handles.O16pO18image_Or,theta,'bilinear','crop'))./...
%             sum(imrotate(Mask_Or,theta,'nearest','crop')))));
%         [theta,fminres] = fminbnd(fun,0,360);%,optimset('TolFun',1e-2,'TolX',1e-2));
%         theta
%
%         opts = optimoptions(@fmincon,'Algorithm','interior-point');
%         problem = createOptimProblem('fmincon','objective',...
%             @(theta) fun,'x0',3,'lb',0,'ub',360,'options',opts);
%         gs = GlobalSearch;
%         [x,f] = run(gs,problem)
%
%     case 3
%         fun = @(theta)...
%             -sum(diff(...
%             sum(imrotate(handles.O16pO18image_Or,theta,'bilinear','crop'))./...
%             sum(imrotate(Mask_Or,theta,'nearest','crop'))));
%         [theta,fminres] = fminbnd(fun,0,360);%,optimset('TolFun',1e-2,'TolX',1e-2));
%         theta
%     case 4
%     case 5
%         fun = @(theta)...
%             -max(diff(...
%             log(sum(imrotate(handles.O16pO18image_Or,theta,'bilinear','crop')))./...
%             sum(imrotate(Mask_Or,theta,'nearest','crop'))));
%         [theta,fminres] = fminbnd(fun,0,360);%,optimset('TolFun',1e-2,'TolX',1e-2));
%         theta
% end


guidata(hObject, handles);


% --- Executes on selection change in AlignMode.
function AlignMode_Callback(hObject, eventdata, handles)
[handles]=Align_Callback(hObject, eventdata, handles);
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);
function AlignMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [handles]=AlignAngle_Callback(hObject, eventdata, handles)
handles.CurrentPlot='Align';
if ~isfield(handles,'AlignAngle_CoAl_num')
    handles.AlignAngle_CoAl_num=0;
end
handles.Masking=0;
set(handles.Ylimit,'Visible','off');
handles.AlignAngle_FiAl_num=str2num(get(handles.AlignAngle,'String'));
set(handles.AlignAngle,'BackgroundColor',[1 1 1]);
set(gca,'position',handles.ImagePlotPos);
if isempty(handles.AlignAngle_FiAl_num)
    handles.AlignAngle_FiAl_num=0;
end
%%% Rotate all images
handles.O16image_FiAl=imrotate(handles.O16image_Or,...
    handles.AlignAngle_FiAl_num,'bilinear','crop');
handles.O18image_FiAl=imrotate(handles.O18image_Or,...
    handles.AlignAngle_FiAl_num,'bilinear','crop');
handles.O16pO18image_FiAl=imrotate(handles.O16pO18image_Or,...
    handles.AlignAngle_FiAl_num,'bilinear','crop');
handles.NormO18image_FiAl=imrotate(handles.NormO18image_Or,...
    handles.AlignAngle_FiAl_num,'bilinear','crop');
% handles.O16pO18image_FiAl=imrotate(handles.O16pO18image_Or,handles.AlignAngle_FiAl_num,'bicubic','crop');

%%% Coarse align all images if possible
if isfield(handles,'AlignAngle_CoAl_num')
    handles.O16image_CoAl=imrotate(handles.O16image_Or,...
        handles.AlignAngle_CoAl_num,'bilinear','crop');
    handles.O18image_CoAl=imrotate(handles.O18image_Or,...
        handles.AlignAngle_CoAl_num,'bilinear','crop');
    handles.O16pO18image_CoAl=imrotate(handles.O16pO18image_Or,...
        handles.AlignAngle_CoAl_num,'bilinear','crop');
    handles.NormO18image_CoAl=imrotate(handles.NormO18image_Or,...
        handles.AlignAngle_CoAl_num,'bilinear','crop');
end

%%% Rotate all masks
[a,b]=size(handles.O16pO18image_Or);
Mask_Or=ones(a,b);
handles.Mask_FiAl_Basic=imrotate(Mask_Or,...
    handles.AlignAngle_FiAl_num,'nearest','crop');
handles.Mask_CoAl_Basic=imrotate(Mask_Or,...
    handles.AlignAngle_CoAl_num,'nearest','crop');

handles.counts=sum(handles.O16pO18image_FiAl)./sum(handles.Mask_FiAl_Basic);
handles.counts(handles.counts==inf)=nan;
handles.counts(handles.counts==-inf)=nan;
%%% Plot fine align
% imagesc(handles.O16pO18image_FiAl);
% colormap(hot); axis square;set(gca, 'XTick', [],'YTick', []);
% set(gca,'position',[68    62   500   465]);
logCounts=log10(handles.counts);

ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
DataStep=ProfileLength*1e-6/(PixelNo-1); %spacial step
handles.X=0:DataStep:ProfileLength*1e-6; %domain

if get(handles.AlignMode, 'Value')==5
    imagesc([0 ProfileLength],...
        [floor(min(log10(handles.counts(isfinite(log10(handles.counts)))))),...
        ceil(max(log10(handles.counts(isfinite(log10(handles.counts))))))],...
        flipud(log10(handles.O16pO18image_FiAl))); axis square;
    hold on;
    
    plot(handles.X*1e6,log10(handles.counts),'c','linewidth',2.5);
    hold off
    handles.colo=colorbar; ylabel(handles.colo,'Log(Counts)');
    legend('Log(Mean Counts)');
    
    if isfinite(sum(logCounts))
        finitecountstart=2;
    else
        finitecountstart=2+find(~isfinite(logCounts),1,'last');
    end
    [SurfGrad SurfPos_idx]=max(diff(logCounts(finitecountstart:end/2)));
    %      
    SurfPos_idx=finitecountstart+SurfPos_idx;
    ylabel('Log(Mean Counts (^{16}O + ^{18}O))');
else
    imagesc([0 ProfileLength],[0 ceil(max(handles.counts))],...
        flipud(handles.O16pO18image_FiAl)); axis square;
    hold on;
    plot(handles.X*1e6,handles.counts,'c','linewidth',2.5);
    hold off
    handles.colo=colorbar; ylabel(handles.colo,'Counts');
    legend('Mean Counts');
    
    [SurfGrad SurfPos_idx]=max(diff(handles.counts(1:end/2)));
    SurfPos_idx=1+SurfPos_idx;
    ylabel('Mean Counts (^{16}O + ^{18}O)');
end

set(gca,'ydir','normal');
xlabel('Depth /\mum');

xlim([0 str2double(get(handles.Xlimit,'String'))]);
colormap(hot); axis square;
set(handles.Xlimit,'Visible','on');
set(handles.Xlimit,'Position',handles.ImagePlotPosXlim);


handles.SurfPos_idx=SurfPos_idx; %%% the gradient is between these two locations, we want to be on the high counts side
SurfPos=1e6*handles.X(SurfPos_idx);
set(handles.SurfPos,'String',...
    num2str(roundsf(SurfPos,3,'round')));
set(handles.SurfPos,'BackgroundColor',[0 1 0]);

handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
    'LineWidth',2,'Color',[.1 1 .1],'LineStyle',':');
% handles.SurfLine=imline(gca,[SurfPos, SurfPos], [-5000, 5000]);%,'Color','w','LineStyle','--');
% handles.h1 = get(handles.SurfLine,'Children');
% set(handles.h1,'LineStyle',':');
% set(handles.h1,'LineWidth',2);
% setColor(handles.SurfLine,[0.1 1 0.1]);
% pos = handles.SurfLine.getPosition();
PlotProperties(hObject, eventdata, handles)
pause(0.1)
guidata(hObject, handles);


function AlignAngle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in MaskMode.
function MaskMode_Callback(hObject, eventdata, handles)
switch get(handles.MaskMode, 'Value')
    case 1;
        set(handles.MaskThresh,'String','0.5')
    case 2;
        set(handles.MaskThresh,'String','1')
end
% MaskButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function MaskMode_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function MaskThresh_Callback(hObject, eventdata, handles)
MaskButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function MaskThresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [handles]=MaskButton_Callback(hObject, eventdata, handles)
handles.CurrentPlot='Mask';
handles.PlotType=1;
handles.Masking=1;
handles.ROIflag=0;
handles.Aligning=0;
if ~isfield(handles,'AlignAngle_FiAl_num')
    [handles]=Align_Callback(hObject, eventdata, handles);
end
set(handles.Xlimit,'Visible','on');
set(handles.Xlimit,'Position',handles.ImagePlotPosXlim);
set(handles.Ylimit,'Visible','off');
[a,b]=size(handles.O16pO18image_Or);
Mask_Or=ones(a,b);
Mask_FiAl_Basic=handles.Mask_FiAl_Basic;
Mask_CoAl_Basic=handles.Mask_CoAl_Basic;

MaskStd=std(handles.O16pO18image_Or(handles.O16pO18image_Or~=0));
MaskMean=mean(handles.O16pO18image_Or(handles.O16pO18image_Or~=0));

switch get(handles.MaskMode, 'Value')
    case 1;
        Threshold=MaskMean*str2double(get(handles.MaskThresh,'String'));
        Mask_Or(handles.O16pO18image_Or<Threshold)=0;
        handles.Mask_FiAl_Thresh=imrotate(Mask_Or,handles.AlignAngle_FiAl_num,...
            'nearest','crop');
        handles.Mask_CoAl_Thresh=imrotate(Mask_Or,handles.AlignAngle_CoAl_num,...
            'nearest','crop');
    case 2;
        Threshold_L=MaskMean-...
            str2double(get(handles.MaskThresh,'String'))*MaskStd;
        Threshold_H=MaskMean+...
            str2double(get(handles.MaskThresh,'String'))*MaskStd;
        % sum(sum(handles.O16pO18image~=0))
        Mask_Or(handles.O16pO18image_Or<Threshold_L | ...
            handles.O16pO18image_Or>Threshold_H)=0;
        handles.Mask_FiAl_Thresh=imrotate(Mask_Or,...
            handles.AlignAngle_FiAl_num,'nearest','crop');
        handles.Mask_CoAl_Thresh=imrotate(Mask_Or,...
            handles.AlignAngle_CoAl_num,'nearest','crop');
    case 3;
        %%% Where the norm18 in a column is more than 1 sd away from mean
        %%% of that column
end


ProLen=str2double(get(handles.ProfileLength,'String'));
x=[0,ProLen];
% set(gca,'position',[90    60   465   465]);%pix
set(gca,'position',handles.ImagePlotPos); %norm

%%% Plot mask
imagesc(x,x,handles.Mask_FiAl_Thresh);
set(gca, 'YTick', [])
xlim([0 str2double(get(handles.Xlimit,'String'))]);
colormap(hot); axis square;xlabel('Depth /\mum');
PlotProperties(hObject, eventdata, handles)
%%% Show surface line
SurfPos=str2double(get(handles.SurfPos,'String'));
handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
    'LineWidth',2,'Color',[1 0 0],'LineStyle',':');
% handles.SurfLine=imline(gca,[SurfPos, SurfPos], [-5000, 5000]);%,'Color','w','LineStyle','--');
% handles.h1 = get(handles.SurfLine,'Children');
% set(handles.h1,'LineStyle',':');
% set(handles.h1,'LineWidth',2);
% setColor(handles.SurfLine,[1 0 0]);
pause(0.4)

%% Plot masked image
imagesc(x,x,handles.O16pO18image_FiAl.*handles.Mask_FiAl_Thresh);
set(gca, 'YTick', [])
xlim([0 str2double(get(handles.Xlimit,'String'))]);
colormap(hot); axis square;xlabel('Depth /\mum');% set(gca, 'YTick', [])
if isfield(handles,'SurfLine2')
    delete(handles.SurfLine2)
end
handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
    'LineWidth',2,'Color',[.1 1 .1],'LineStyle',':');
% handles.SurfLine=imline(gca,[SurfPos, SurfPos], [-5000, 5000]);%,'Color','w','LineStyle','--');
% handles.h1 = get(handles.SurfLine,'Children');
% set(handles.h1,'LineStyle',':');
% set(handles.h1,'LineWidth',2);
% setColor(handles.SurfLine,[0.1 1 0.1]);

ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
handles.X=0:dx:ProfileLength*1e-6; %domain

% grid (get(handles.GridLines,'State'))
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);




function ROI_Callback(hObject, eventdata, handles)
[handles]=MaskButton_Callback(hObject, eventdata, handles);
handles.ROIflag=1;
%%% make variables
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
SurfPos=str2double(get(handles.SurfPos,'String'));

%%%spec ROI
[x,y,BW,xi,yi] = roipoly;
handles.ROImin=min(xi);%Know to change fit limits
handles.ROImax=max(xi);
handles.Mask_FiAl_Thresh=handles.Mask_FiAl_Thresh.*BW;
%%%plot ROI mask
imagesc(handles.X*1e6,handles.X*1e6,handles.Mask_FiAl_Thresh);
xlim([0 str2double(get(handles.Xlimit,'String'))]);
set(gca, 'YTick', []); colormap(hot);
axis square; xlabel('Depth /\mum');
% handles.SurfLine=imline(gca,[SurfPos, SurfPos], [-5000, 5000]);%,'Color','w','LineStyle','--');
% handles.h1 = get(handles.SurfLine,'Children');
% set(handles.h1,'LineStyle',':');
% set(handles.h1,'LineWidth',2);
% setColor(handles.SurfLine,[1 0 0]);
handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
    'LineWidth',2,'Color',[1 0 0],'LineStyle',':');

pause(0.5)

%%%plot thresholded region
imagesc([0 ProfileLength],[0 ProfileLength],handles.O16pO18image_FiAl.*handles.Mask_FiAl_Thresh);
xlim([0 str2double(get(handles.Xlimit,'String'))]);
set(gca, 'YTick', []); colormap(hot);
axis square;xlabel('Depth /\mum');% set(gca, 'YTick', [])
% handles.SurfLine=imline(gca,[SurfPos, SurfPos], [-5000, 5000]);%,'Color','w','LineStyle','--');
% handles.h1 = get(handles.SurfLine,'Children');
% set(handles.h1,'LineStyle',':');
% set(handles.h1,'LineWidth',2);
% setColor(handles.SurfLine,[0.1 1 0.1]);
if isfield(handles,'SurfLine2')
    delete(handles.SurfLine2)
end
handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
    'LineWidth',2,'Color',[.1 1 .1],'LineStyle',':');
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);

% --- Executes on button press in GenerateProfiles.
function GenerateProfiles_Callback(hObject, eventdata, handles)
handles.PlotType=1;
handles.CurrentPlot='Generate';
if ~isfield(handles,'Mask_FiAl_Thresh')
    [handles]=MaskButton_Callback(hObject, eventdata, handles);
else
    if handles.Masking==0;
        [handles]=MaskButton_Callback(hObject, eventdata, handles);
    end
end

handles.DataFlag=1;
Mask_FiAl_Thresh=handles.Mask_FiAl_Thresh; %Mask

%
NormO18image_FiAl_MaTh=handles.O18image_FiAl./handles.O16pO18image_FiAl.*Mask_FiAl_Thresh;
NormO18image_FiAl_MaTh(isnan(NormO18image_FiAl_MaTh))=0;
NormO18image_FiAl_MaTh(isinf(NormO18image_FiAl_MaTh))=0;
handles.NormO18image_FiAl_MaTh=NormO18image_FiAl_MaTh;

O16pO18prof_FiAl_MaTh=sum(handles.O16pO18image_FiAl.*handles.Mask_FiAl_Thresh)./sum(handles.Mask_FiAl_Thresh); %Masked O16pO18
% O16pO18prof_FiAl_MaTh(isnan(O16pO18prof_FiAl_MaTh))=0;
O16pO18prof_FiAl_MaTh(isinf(O16pO18prof_FiAl_MaTh))=0;
O18prof_FiAl_MaTh=sum(handles.O18image_FiAl.*handles.Mask_FiAl_Thresh)./sum(handles.Mask_FiAl_Thresh); %Masked O18
% O18prof_FiAl_MaTh(isnan(O18prof_FiAl_MaTh))=0;
O18prof_FiAl_MaTh(isinf(O18prof_FiAl_MaTh))=0;
% O16prof_FiAl_MaTh=sum(handles.O16image_FiAl.*handles.Mask_FiAl_Thresh)./sum(handles.Mask_FiAl_Thresh); %Masked O18
% O16prof_FiAl_MaTh(isnan(O16prof_FiAl_MaTh))=0;
% O16prof_FiAl_MaTh(isinf(O16prof_FiAl_MaTh))=0;

% UserProfSpec=1;
NormO18prof_FiAl_MaTh=O18prof_FiAl_MaTh./O16pO18prof_FiAl_MaTh; %Normalised profile
P(1,:)=NormO18prof_FiAl_MaTh;

%%% Crap profiles

%Profile based Aligned not masked
% O16pO18prof_AnM=sum(handles.O16pO18image);%.*handles.Mask_FiAl_Basic); %Masked O16pO18
% O18prof_AnM=sum(handles.O18image);%.*handles.Mask_FiAl_Basic);
% NormO18prof_AnM=O18prof_AnM./O16pO18prof_AnM;
% NormO18prof_AnM(isnan(NormO18prof_AnM))=nan;
% NormO18prof_AnM(isinf(NormO18prof_AnM))=nan;

O16pO18prof_AnM=sum(handles.O16pO18image_FiAl)./sum(handles.Mask_FiAl_Basic);
O18prof_AnM=sum(handles.O18image_FiAl)./sum(handles.Mask_FiAl_Basic);
NormO18prof_AnM=O18prof_AnM./O16pO18prof_AnM;
NormO18prof_AnM(isnan(NormO18prof_AnM))=nan;
NormO18prof_AnM(isinf(NormO18prof_AnM))=nan;

P(2,:)=NormO18prof_AnM;

%Profile based Masked not aligned
O16pO18prof_MnA=sum(handles.O16pO18image_CoAl.*handles.Mask_CoAl_Thresh);
O18prof_MnA=sum(handles.O18image_CoAl.*handles.Mask_CoAl_Thresh);
NormO18prof_MnA=O18prof_MnA./O16pO18prof_MnA;
NormO18prof_MnA(isnan(NormO18prof_MnA))=nan;
NormO18prof_MnA(isinf(NormO18prof_MnA))=nan;
P(3,:)=NormO18prof_MnA;

%Basic Profile Normailisation
O16pO18prof_Or=sum(handles.O16pO18image_CoAl.*handles.Mask_CoAl_Basic);
O18prof_Or=sum(handles.O18image_CoAl.*handles.Mask_CoAl_Basic);
NormO18prof_Or=O18prof_Or./O16pO18prof_Or;
NormO18prof_Or(isnan(NormO18prof_Or))=nan;
NormO18prof_Or(isinf(NormO18prof_Or))=nan;
NormO18prof_Or(NormO18prof_Or==1)=nan;
P(4,:)=NormO18prof_Or;

%Original image based normalisation
handles.NormO18prof_Or_image=handles.O18image_CoAl./handles.O16pO18image_CoAl;
handles.NormO18prof_Or_image(isnan(handles.NormO18prof_Or_image))=0;
P(5,:)=mean(handles.NormO18prof_Or_image);

handles.ProfileData_Or=NormO18prof_FiAl_MaTh;
[handles]=DataPlotOnly(hObject, eventdata, handles);
% set(handles.Ylimit,'Position',[26 525 40 20]);
%Setting evironment
ProfileSelection=get(handles.ProfileSelector, 'Value');
if ProfileSelection==6
    ProfileSelection=1;
end

% [SurfGrad SurfPos]=min(O16pO18prof_FiAl_MaTh(1:end-1)-O16pO18prof_FiAl_MaTh(2:end))
% [SurfGrad SurfPos]=min(O16pO18prof_AnM(1:end-1)-O16pO18prof_AnM(2:end))
%[SurfGrad SurfPos]=max(diff(O16pO18prof_AnM)); %Surface is at greatest gradient of rotated, but unmasked O16pO18image
% SurfPos=SurfPos+1

%Possibly start the fit where the surface roughness has stopped...
% set(handles.SurfPos,'String',...
%     num2str(roundsf(1e6*handles.Xdata_Or(SurfPos),2,'round')));

C_bg=median(roundsf(P(ProfileSelection,round(end/2):end),2,'round'));

% Pro_End=roundsf(max(1e6*handles.Xdata_Or)-...
%     str2double(get(handles.SurfPos,'String')),3,'floor');
% if C_bg<0.004 && C_bg>0.001
%     Fit_End=roundsf(1.2*(1e6*handles.Xdata_Or(find(...
%         roundsf(NormO18prof_FiAl_MaTh(handles.SurfPos_idx:end),2,'round' )<...
%         roundsf(C_bg,2,'round'),1,'first')...
%         )-str2double(get(handles.SurfPos,'String'))),3,'round');
%
%     %     find(NormO18prof_FiAl_MaTh
%
%     %     Fit_End=1.1*(num2str(roundsf(1e6*handles.Xdata_Or(find(...
%     %         roundsf(NormO18prof_FiAl_MaTh,2,'round' )<...
%     %         roundsf(C_bg,2,'round'),1,'first')),2,'round')-...
%     %         roundsf(str2double(get(handles.SurfPos,'String')),3,'floor')))
%     if Fit_End>Pro_End
%         Fit_End=Pro_End;
%     else
%     end
% else
%     C_bg=0.002;
%     Fit_End=Pro_End;
% end
handles.counts=sum(handles.O16pO18image_FiAl)./sum(handles.Mask_FiAl_Basic);
handles.counts(handles.counts==inf)=0;
% find(roundsf(handles.counts,2,'round')>0.99*mean(handles.counts),1,'first')
% [RatMax Fit_Start]=max(NormO18prof_FiAl_MaTh(SurfPos:end));
% set(handles.Fit_Start,'String',...
%     num2str(roundsf(1e6*handles.Xdata_Or(Fit_Start),2,'round')));
% pause(1)

Bmean_idx=find(roundsf(handles.counts,2,'round')>mean(handles.counts),1,'first');
ProStart_idx=find(P(ProfileSelection,:)>0,1);
if max(Bmean_idx,ProStart_idx)-handles.SurfPos_idx>0
    Fit_Start_idx=max(Bmean_idx,ProStart_idx)-handles.SurfPos_idx;
else
    Fit_Start_idx=1;
end
Fit_Start=1e6*handles.X(Fit_Start_idx);

% Bmean=1e6*handles.Xdata_Or(find(roundsf(handles.counts,2,'round')>...
%     mean(handles.counts),1,'first'));
% if Bmean<str2double(get(handles.SurfPos,'String'))
%     Fit_Start=0;
% else
%     Fit_Start=Bmean-str2double(get(handles.SurfPos,'String'));
% end
% % Fit_Start=(roundsf(1e6*handles.Xdata(find(...
% %     roundsf(handles.counts,2,'round')>0.99*mean(handles.counts)...
% %     ,1,'first')),2,'round'))-str2double(get(handles.SurfPos,'String'));
% Pro_Start=-roundsf(1e6*handles.Xdata_Or(find(P(ProfileSelection,:)>0,1))-...
%     str2double(get(handles.SurfPos,'String')),3,'floor')
% % find(handles.P(ProfileSelection,:)>0,1)
% % handles.SurfPos_idx
%
% if find(P(ProfileSelection,:)>0,1)>handles.SurfPos_idx
%     if 1e6*handles.X(find(P(ProfileSelection,:)>0,1)-handles.SurfPos_idx)>Fit_Start
%         Fit_Start=1e6*handles.X(find(P(ProfileSelection,:)>0,1)-handles.SurfPos_idx+1);
%     end
% end
%
% if find(P(ProfileSelection,:)>0,1,'last')>handles.SurfPos_idx
%     if 1e6*handles.X(find(P(ProfileSelection,:)>0,1,'last')-handles.SurfPos_idx)<Fit_End
%         Fit_End=1e6*handles.X(find(P(ProfileSelection,:)>0,1,'last')-handles.SurfPos_idx);
%     end
% end

% if handles.ROIflag==1
%     if handles.ROImin-str2double(get(handles.SurfPos,'String'))>Fit_Start
%         Fit_Start=handles.ROImin-str2double(get(handles.SurfPos,'String'))+2;
%     end
%     if handles.ROImax-str2double(get(handles.SurfPos,'String'))<Fit_End
%         Fit_End=handles.ROImax-str2double(get(handles.SurfPos,'String'))-2;
%     end
% % end
% if Fit_Start<Pro_Start
%     Fit_Start=Pro_Start;
% end
if ~isfinite(C_bg)
    C_bg=0.002;
end
set(handles.C_bg,'String',C_bg);
set(handles.C_bg,'BackgroundColor',[0,1,0]);
% if Fit_End<0 || ~isfinite(Fit_End)
%     Fit_End=0.5*str2double(get(handles.ProfileLength,'String'));
% else
%     set(handles.Fit_End,'BackgroundColor',[0,1,0]);
% end
% set(handles.Fit_End,'BackgroundColor',[0,1,0]);
if Fit_Start<0
    Fit_Start=0;
elseif Fit_Start>str2double(get(handles.ProfileLength, 'String'))/3;
    Fit_Start=0;
else
    set(handles.Fit_Start,'BackgroundColor',[0,1,0]);
end
set(handles.Fit_Start,'String',roundsf(Fit_Start,3,'ceil'));

handles.P=P;
%Plotting
[handles]=ProfileSelector_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
function GenerateProfiles_CreateFcn(hObject, eventdata, handles)

function [handles]=ProfileSelector_Callback(hObject, eventdata, handles)
%Plotting
handles.CurrentPlot='Generate';
handles.ProfileGenFlag=1;
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
SurfPos=str2double(get(handles.SurfPos, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
% DataStep=ProfileLength*1e-6/PixelNo;
% [m,n]=size(handles.ProfileData_Or);
% handles.Xdata_Or=0:DataStep:(max(m,n)-1)*DataStep;
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain
handles.Xdata_Or=X;
pX=handles.Xdata_Or*1e6;
P=handles.P;

%NormO18image_FiAl_MaTh=handles.O16pO18image.Aligned./(handles.O16image+handles.O18image);
if get(handles.ProfileSelector, 'Value')==6
    handles.ProfileData_Or=P(1,:);
    plot(pX,P(1,:),pX,P(2,:),pX,P(3,:),pX,P(4,:),pX,P(5,:));
    axis square;
    legend('Aligned and Masked','Aligned Only','Masked Only','Basic Norm. Prof.','Original Norm. Image')
    set(handles.Ylimit,'String',ceil(max(handles.ProfileData_Or)*110)/100);
    
    %     Ylimit_Callback(hObject, eventdata, handles);
    set(handles.Ylimit,'Visible','on');
    set(handles.Ylimit,'Position',handles.ImagePlotPosYlim);
    set(handles.Xlimit,'Position',handles.ImagePlotPosXlim);
    %     handles.SaveName=[handles.SaveName_Or(1:end-4),'(all).txt'];
    set(handles.ProfDataSaveName, 'String',[handles.SaveName_Or(1:end-4),'(all).txt']);
    % set(gca,'position',[90    60   465   465]);%pix
    set(gca,'position',handles.ImagePlotPos); %norm
    axis square;
else
    set(handles.Ylimit,'Visible','on');
    handles.OverPlot=1;
    switch get(handles.ProfileSelector, 'Value')
        case 1;
            handles.ProfileData_Or=P(1,:);
            %             handles.SaveName=[handles.SaveName_Or(1:end-4),'(A_M).txt'];
            set(handles.ProfDataSaveName, 'String',[handles.SaveName_Or(1:end-4),'(A_M).txt']);
        case 2;
            handles.ProfileData_Or=P(2,:);
            %             handles.SaveName=[handles.SaveName_Or(1:end-4),'(A).txt'];
            set(handles.ProfDataSaveName, 'String',[handles.SaveName_Or(1:end-4),'(A).txt']);
        case 3;
            handles.ProfileData_Or=P(3,:);
            %             handles.SaveName=[handles.SaveName_Or(1:end-4),'(M).txt'];
            set(handles.ProfDataSaveName, 'String',[handles.SaveName_Or(1:end-4),'(M).txt']);
        case 4;
            handles.ProfileData_Or=P(4,:);
            %             handles.SaveName=[handles.SaveName_Or(1:end-4),'(NormProf).txt'];
            set(handles.ProfDataSaveName, 'String',[handles.SaveName_Or(1:end-4),'(NormProf).txt']);
        case 5;
            handles.ProfileData_Or=P(5,:);
            %             handles.SaveName=[handles.SaveName_Or(1:end-4),'(Or).txt'];
            set(handles.ProfDataSaveName, 'String',[handles.SaveName_Or(1:end-4),'(Or).txt']);
    end
    %     set(handles.Ylimit,'String',roundsf(max(handles.ProfileData_Or)*1.1,2,'ceil'));
    set(handles.Ylimit,'String',ceil(max(handles.ProfileData_Or(X>SurfPos*1e-6))*105)/100);
    %     set(handles.Ylimit,'String',1);
    imagesc([0 str2num(get(handles.ProfileLength,'String'))],[0 str2num(get(handles.Ylimit,'String'))],...
        flipud(handles.NormO18image_FiAl_MaTh)); hold on; axis square;
    handles.colo=colorbar; ylabel(handles.colo,'Isotopic Fraction');
    set(handles.colo,'FontSize',10);
    set(handles.Ylimit,'Position',handles.ImagePlotPosYlim);
    set(handles.Xlimit,'Position',handles.ImagePlotPosXlim);
    % set(gca,'position',[90    60   465   465]);%pix
    set(gca,'position',handles.ImagePlotPos); %norm
    %     if length(get(handles.CprimeY,'State'))==2;
    %         C_bg=str2double(get(handles.C_bg, 'String'));
    %         C_gas=str2double(get(handles.C_gas, 'String'));
    %         pro=(handles.ProfileData_Or-C_bg)/(C_gas-C_bg);
    %         ylabel('Normalised _{}^{18}O Fraction, C''');
    %     else
    pro=handles.ProfileData_Or;
    %         ylabel('Isotopic Fraction');
    %     end
    hAx=plot(pX,pro,'w','linewidth',2.5);
    %     plot(pX,,pX,handles.ProfileData_Or,'w','linewidth',2.8);
    set(gca,'ydir','normal');
    %%%Show surface position
    
    %     handles.SurfLine=imline(gca,[SurfPos, SurfPos], [-5000, 5000]);%,'Color','w','LineStyle','--');
    %     handles.h1 = get(handles.SurfLine,'Children');
    %     set(handles.h1,'LineStyle',':');
    %     set(handles.h1,'LineWidth',2)
    %     setColor(handles.SurfLine,[0.1 1 0.1]);
    handles.SurfLine2=line([SurfPos, SurfPos],[-5000 5000],...
        'LineWidth',2,'Color',[0 1 1],'LineStyle',':');
    %%Show Cbg level
    C_bg=str2double(get(handles.C_bg,'String'));
    handles.C_bgLine=line([-5000 5000],[C_bg, C_bg],...
        'LineWidth',2,'Color',[.1 .1 1],'LineStyle',':');
    %     handles.C_bgLine=imline(gca, [-5000, 5000],[C_bg, C_bg]);%,'Color','w','LineStyle','--');
    %     h2 = get(handles.C_bgLine,'Children');
    %     set(h2,'LineStyle',':');
    %     set(h2,'LineWidth',2)
    %     legend(hAx,'Profile',handles.SurfLine,'Surface position',handles.C_bg,'Background concentration')
    %     setColor(handles.SurfLine,[0.1 1 0.1]);
end
handles.Ylimit_min=0;
Ylimit_Callback(hObject, eventdata, handles);
xlabel('Depth /\mum');ylabel('Isotopic Fraction');
Xlimit_Callback(hObject, eventdata, handles);
hold off;
grid (get(handles.GridLines,'State'))
set(gca,'Xcolor',[0.3 0.3 0.3]); set(gca,'Ycolor',[0.3 0.3 0.3]);
% handles.SaveName=[handles.FileName(1:end-4),...
%     '_Norm18Prof(',num2str(get(handles.ProfileSelector, 'Value')),').txt'];
% handles.SaveName=[handles.PathName,'\',handles.FileName(1:end-4),...
%     '_Norm18Prof(',num2str(get(handles.ProfileSelector, 'Value')),').txt'];
PlotProperties(hObject, eventdata, handles)
guidata(hObject, handles);
function ProfileSelector_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ProfDataSave_Callback(hObject, eventdata, handles)
%Saving profiles
% handles.Name=get(handles.ProfDataSaveName,'String');
SaveName=get(handles.ProfDataSaveName,'String');
[FileName,PathName]=uiputfile([handles.PathName,SaveName]);
fid=fopen([FileName,PathName],'wt');
if get(handles.ProfileSelector, 'Value')==6;
    fprintf(fid,'%d, %d, %d, %d \n',handles.P);
else
    fprintf(fid,'%d\n',handles.P(get(handles.ProfileSelector, 'Value'),:));
end
fclose(fid);
guidata(hObject, handles);


function ProfDataSaveName_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
function ProfDataSaveName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [y]=roundsf(number,sfs,method)
%opt = {'round','floor','ceil','fix'};
og = 10.^(floor(log10(abs(number)) - sfs + 1));
y = feval(method,number./og).*og;
y(find(number==0)) = 0;

function [r2 rmse] = rsquare(y,f,varargin)
if isempty(varargin); c = true;
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1};
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
        % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));

%BackDiffs_AutoCaller
function [pro]=BackDiffs_AutoCaller(p,b,ran,FitCheck,PS,MP);
if min(p)<=0
    pro=zeros(1,ran(2)-ran(1)+1);
else
    %     p=abs(p)
    if sum(FitCheck)==4
        D1=p(1);        D2=p(2);        k1=p(3);        k2=p(4);
    elseif FitCheck(2)==0
        D1=p(1);        D2=p(1);        k1=p(2);        k2=p(3);
    elseif FitCheck(4)==0
        D1=p(1);        D2=p(2);        k1=p(3);        k2=p(3);
    end
    
    %%%%%%%%%%%%
    % p(1:2)
    % p(3:4)
    %%%%%%%%%%%%
    C_gas=b(1);
    C_bg=b(2);
    t1=b(3);
    t2=b(4);
    ProfileLength=b(5);
    PixelNo=b(6);
    
    [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP);
    pro=pro(ran(1):ran(2));
end

function [X,pro,DepthFlag]=BackDiffsCN_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo)
DepthFlag=0;
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain
% dx=X(2)-X(1);
t_i=t1*3600; %s %initial exchange duration
t_b=t2*3600; %s %back exchange duration
t_tot=t_i+t_b; %s %total time
h1=k1/D1;
h2=k2/D2;
C_gas1=C_gas; % gas concentration for period 1
C_gas2=C_bg; % gas concentration for period 2

dt=5*round(dx^2/(2*max(D1,D2)),3,'significant'); %time step
% dt=round(dx^2/(2*max(D1,D2)),3,'significant');
Nt=round((t_tot)/dt); % Number of time steps
I=diag(ones(1,PixelNo));

sigma1=D1*dt/(2*dx^2);
sigma2=D2*dt/(2*dx^2);

%%% CN Diffusion matrix
A1=    full(gallery('tridiag',PixelNo, sigma1,1-2*sigma1, sigma1));
A2=    full(gallery('tridiag',PixelNo, sigma2,1-2*sigma2, sigma2));
A1_new=full(gallery('tridiag',PixelNo,-sigma1,1+2*sigma1,-sigma1));
A2_new=full(gallery('tridiag',PixelNo,-sigma2,1+2*sigma2,-sigma2));
%%% Exchange surface condition
A1(1,1:2)=    [1-2*dx*h1*sigma1-2*sigma1, 2*sigma1];
A2(1,1:2)=    [1-2*dx*h2*sigma2-2*sigma2, 2*sigma2];
A1_new(1,1:2)=[1+2*dx*h1*sigma1+2*sigma1,-2*sigma1];
A2_new(1,1:2)=[1+2*dx*h2*sigma2+2*sigma2,-2*sigma2];
beta1=-h1*C_gas1;
beta2=-h2*C_gas2;
G1=zeros(PixelNo,1); G1(1)=-4*dx*beta1*sigma1;
G2=zeros(PixelNo,1); G2(1)=-4*dx*beta2*sigma2;
%%% Mirror oundary condition
A1(end,end-2:end)=    [sigma1 -2*sigma1 1+sigma1]; %same second der
A2(end,end-2:end)=    [sigma2 -2*sigma2 1+sigma2];
A1_new(end,end-2:end)=[-sigma1 2*sigma1 1-sigma1];
A2_new(end,end-2:end)=[-sigma2 2*sigma2 1-sigma2];
A1_newI=A1_new\I;
A2_newI=A2_new\I;
% Time stepping
C=zeros(PixelNo,Nt+1);
C(:,1)=C_bg;
for t_idx=1:Nt
    if t_idx<=ceil(t_i/dt)
        C(:,t_idx+1)=A1_newI*(A1*C(:,t_idx)+G1);
    else
        C(:,t_idx+1)=A2_newI*(A2*C(:,t_idx)+G2);
    end
end
pro=C(:,end); pro=pro';
% pro
% combine with linear grad if profile extends outsde domain
if pro(end)>C_bg*1.1 %max(pro)-pro(end)>pro(end)*1.01
    C=zeros(PixelNo,Nt+1);
    C(:,1)=C_bg;
    %%% Linear grad.
    A1(end,end-2:end)=    [-sigma1 2*sigma1 1-sigma1]; %central dif
    A2(end,end-2:end)=    [-sigma2 2*sigma2 1-sigma2];
    A1_new(end,end-2:end)=[sigma1 -2*sigma1 1+sigma1];
    A2_new(end,end-2:end)=[sigma2 -2*sigma2 1+sigma2];
    A1_newI=A1_new\I;
    A2_newI=A2_new\I;
    % Time stepping
    for t_idx=1:Nt
        if t_idx<=ceil(t_i/dt)
            C(:,t_idx+1)=A1_newI*(A1*C(:,t_idx)+G1);
        else
            C(:,t_idx+1)=A2_newI*(A2*C(:,t_idx)+G2);
        end
    end
    
    pro=(pro'+C(:,end))./2; pro=pro';%pro=C(:,end); pro=pro';
    
end




%%%%%%%%%%%%% Optimisation trial
function [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,MP)
DepthFlag=0;
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain
% dx=X(2)-X(1);
t_i=t1*3600; %s %initial exchange duration
t_b=t2*3600; %s %back exchange duration
t_tot=t_i+t_b; %s %total time
h1=k1/D1;
h2=k2/D2;
C_gas1=C_gas; % gas concentration for period 1
C_gas2=C_bg; % gas concentration for period 2
Crank1=(C_bg+(C_gas-C_bg)*...
    (erfc(X./(2*sqrt(D1*t_i)))-...
    exp(h1.*X+t_i*D1*h1^2)...
    .*erfc(X./(2*sqrt(D1*t_i))+h1*sqrt(D1*t_i))))';

dt=5*round(dx^2/(2*max(D1,D2)),3,'significant'); %time step
% dt=round(dx^2/(2*max(D1,D2)),3,'significant');
Nt=round((t_tot)/dt); % Number of time steps

% Time stepping
if min([D1,D2,k1,k2])<0
    pro=2*ones(PixelNo,1);
else
    if PS==1
        PixelNo_Or=PixelNo;
        PixelNo=MP;
    end
    C=C_bg*ones(PixelNo,1);
    C1=C;
    I=diag(ones(1,PixelNo));
    %  
    if mean(~isfinite(Crank1))>0 || PS==1;
        %%% Step 1
        sigma1=D1*dt/(2*dx^2);
        %%% CN Diffusion matrix
        A1=    full(gallery('tridiag',PixelNo, sigma1,1-2*sigma1, sigma1));
        A1_new=full(gallery('tridiag',PixelNo,-sigma1,1+2*sigma1,-sigma1));
        %%% Exchange surface condition
        A1(1,1:2)=    [1-2*dx*h1*sigma1-2*sigma1, 2*sigma1];
        A1_new(1,1:2)=[1+2*dx*h1*sigma1+2*sigma1,-2*sigma1];
        beta1=-h1*C_gas1;
        G1=zeros(PixelNo,1); G1(1)=-4*dx*beta1*sigma1;
        %%% Mirror Boundary condition
        A1(end,end-2:end)=    [sigma1 -2*sigma1 1+sigma1];
        A1_new(end,end-2:end)=[-sigma1 2*sigma1 1-sigma1];
        A1_newI=A1_new\I;
        for t_idx=1:ceil(t_i/dt) %%%% double  check this
            C=A1_newI*(A1*C+G1);
        end
        C1=C;
    else
        C1=Crank1;
        t_idx=ceil(t_i/dt);
    end
    C=C1;
    
    %%% Step 2
    sigma2=D2*dt/(2*dx^2);
    %%% CN Diffusion matrix
    A2=    full(gallery('tridiag',PixelNo, sigma2,1-2*sigma2, sigma2));
    A2_new=full(gallery('tridiag',PixelNo,-sigma2,1+2*sigma2,-sigma2));
    %%% Exchange surface condition
    A2(1,1:2)=    [1-2*dx*h2*sigma2-2*sigma2, 2*sigma2];
    A2_new(1,1:2)=[1+2*dx*h2*sigma2+2*sigma2,-2*sigma2];
    beta2=-h2*C_gas2;
    G2=zeros(PixelNo,1); G2(1)=-4*dx*beta2*sigma2;
    %%% Mirror Boundary condition
    A2(end,end-2:end)=    [sigma2 -2*sigma2 1+sigma2];
    A2_new(end,end-2:end)=[-sigma2 2*sigma2 1-sigma2];
    A2_newI=A2_new\I;
    %  
    for t_idx=1:Nt-t_idx   %% double check
        C=A2_newI*(A2*C+G2);
    end
    if PS==1
        C_temp=C_bg*ones(PixelNo_Or,1);
        C_temp(1:length(C))=C;
        if MP<PixelNo_Or
            %          
            if PixelNo_Or<=2*MP
                C_temp(length(C)+1:end)=C(end:-1:1+end+MP-length(C_temp));
            else
                C_temp(length(C)+1:2*length(C))=C(end:-1:1);
                C_temp(2*length(C):end)=0;
            end
        end
        C=C_temp;
    end
    pro=C';
    % pro
    % combine with linear grad if profile extends outsde domain
    
    %%%% Maybe change this to if L*<2
    
    if PS==0;
        if pro(end)>C_bg*1.1 %max(pro)-pro(end)>pro(end)*1.01
            C=C_bg*zeros(PixelNo,1);
            %%% Linear grad.
            if ~isfinite(Crank1)
                A1(end,end-2:end)=    [-sigma1 2*sigma1 1-sigma1]; %central dif
                A1_new(end,end-2:end)=[sigma1 -2*sigma1 1+sigma1];
                A1_newI=A1_new\I;
                for t_idx=1:ceil(t_i/dt) %%%% double  check this
                    C1=A1_newI*(A1*C+G1);
                end
            else
                C1=Crank1;
                t_idx=ceil(t_i/dt);
            end
            C=C1;
            
            A2(end,end-2:end)=    [-sigma2 2*sigma2 1-sigma2];
            A2_new(end,end-2:end)=[sigma2 -2*sigma2 1+sigma2];
            A2_newI=A2_new\I;
            % Time stepping
            
            for t_idx=1:Nt-t_idx   %% double check
                C=A2_newI*(A2*C+G2);
            end
            pro=(pro+C')./2;%pro=C(:,end); pro=pro';
        end
    end %nagetive
end

% for i=1:2 % Build coefficient matricies
%     h(i)=k(i)/D(i);
%     beta(i)=-h(i)*C_gas(i);
%     sigma(i)=D(i)*dt/(2*dx^2);
%     G(:,i)=zeros(PixelNo,1);
%     G(1,i)=-4*dx*beta(i)*sigma(i);
%
%     A(:,:,i)=full(gallery('tridiag',PixelNo,sigma(i),1-2*sigma(i),sigma(i)));
%     A1(1,1:2,i)=[1-2*dx*h(i)*sigma(i)-2*sigma(i) , 2*sigma(i)];
%     A_new(:,:,i)=full(gallery('tridiag',PixelNo,-sigma(i),1+2*sigma(i),-sigma(i)));
%     A_new(:,:,i)=[1+2*dx*h(i)*sigma(i)+2*sigma(i) , -2*sigma(i)];
% end

function WarningBox_CreateFcn(hObject, eventdata, handles)

function [X,pro]=InDiffs_inline(C_gas,C_bg,D1,k1,t1,ProfileLength,PixelNo,handles)
set(handles.WarningBox,'Visible','off')
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain
h1=k1/D1;
t_tot=t1*3600;
if D1<0 || k1<0
    proPrime=ones(PixelNo);
else
    a=erfc(X./(2*sqrt(D1*t_tot)));
    b=exp(h1.*X+D1*t_tot*h1^2).*erfc(X./(2*sqrt(D1*t_tot))+h1*sqrt(D1*t_tot));
    %     if ~isfinite(b) % 1-exp((z).^2).*erfc(z)=1-0.0212=0.9788...
    %
    proPrime=a-b;
end
pro=proPrime*(C_gas-C_bg)+C_bg;

function [X,pro]=InDiffsPS_inline(handles)
set(handles.WarningBox,'Visible','off')
C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain

D1=str2double(get(handles.D1, 'String'));
k1=str2double(get(handles.k1, 'String'));
t1=str2double(get(handles.t1, 'String'));
t_tot=t1*3600;
MP=str2double(get(handles.MirrorPlane,'String'))*1e-6;
[val_mirror MP_idx] = min(abs(X-MP));
if D1<0 || k1<0
    C_prime=ones(PixelNo);
else
    L=MP*k1/D1
    RootNo=floor(3*(L+10)^.75);
    beta=zeros(1,RootNo);
    if L>30
        HuntStep=0.01;
    else
        HuntStep=0.01;
    end
    fun1 = @(beta) (beta*tan(beta)-L)^2;
    tic
    for k=1:RootNo
        beta(k)=fminbnd(fun1,pi*(k-1),pi*(k-0.5));
    end
    toc
    for i=1:length(beta)
        RootSum(i,:)=2*L*cos(beta(i)*(X-MP)/MP)*...
            exp(-beta(i)^2*D1*t_tot/MP^2)/((beta(i)^2+L^2+L)*cos(beta(i)));
    end
    C_prime=1-sum(RootSum,1);
end
pro=C_prime*(C_gas-C_bg)+C_bg;
if length(pro)>2*MP_idx
    pro(2*MP_idx:end)=0;
end

function [X,pro]=BackDiffs1k_inline(C_gas,C_bg,D1,k1,t1,t2,ProfileLength,PixelNo,handles)
set(handles.WarningBox,'Visible','off')
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain
h1=k1/D1;
t_tot=(t1+t2)*3600;
phi=t2/(t1+t2);
if D1<0 || k1<0
    proPrime=ones(PixelNo);
else
    proPrime=(erfc(X./(2*sqrt(D1*t_tot)))-exp(h1.*X+D1*t_tot*h1^2).*erfc(...
        X./(2*sqrt(D1*t_tot))+h1*sqrt(D1*t_tot))) - ...
        (erfc(X./(2*sqrt(D1*t_tot*phi)))-exp(h1.*X+D1*t_tot*phi*h1^2).*erfc(...
        X./(2*sqrt(D1*t_tot*phi))+h1*sqrt(D1*t_tot*phi)));
end
pro=proPrime*(C_gas-C_bg)+C_bg;

%% Fitter
% function [D1,D2,k1,k2,r2]=BackDiffs_fit_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,Data,...
%     range,points,FitCheck,ProfileLength,Fit_Start_idx,Fit_End_idx,handles)
%
% PixelNo=512;
% [m n]=size(Data);
% step=2*range/(points-1);
%
% R_D1=D1;
% R_D2=D2;
% R_k1=k1;
% R_k2=k2;
% if FitCheck(1)==1
%     R_D1=roundsf([(1-range):step:(1+range)]*D1,3,'round');
% end
% if FitCheck(2)==1
%     R_D2=roundsf([(1-range):step:(1+range)]*D2,3,'round');
% end
% if FitCheck(3)==1
%     R_k1=roundsf([(1-range):step:(1+range)]*k1,3,'round');
% end
% if FitCheck(4)==1
%     R_k2=roundsf([(1-range):step:(1+range)]*k2,3,'round');
% end
%
% %Make the matrix of values and then iterate through it
% results=0;
% counter=0;
% if t2==0
%     for D1=R_D1
%         for k1=R_k1
%             counter=counter+1;
%             [X,pro]=InDiffs_inline(C_gas,C_bg,D1,k1,t1,ProfileLength,PixelNo,handles);
%             [r2 rmse] = rsquare(Data(Fit_Start_idx:Fit_End_idx),...
%                 pro(Fit_Start_idx:Fit_End_idx)); %DataRange
%             results(counter,1)=D1;
%             results(counter,2)=D1;
%             results(counter,3)=k1;
%             results(counter,4)=k1;
%             results(counter,5)=r2;
%         end
%     end
% elseif t2>0 && FitCheck(2)==0 && FitCheck(4)==0
%     for D1=R_D1
%         for k1=R_k1
%             counter=counter+1;
%             [X,pro]=BackDiffs1k_inline(C_gas,C_bg,D1,k1,t1,t2,ProfileLength,PixelNo,handles);
%             [r2 rmse] = rsquare(Data(Fit_Start_idx:Fit_End_idx),...
%                 pro(Fit_Start_idx:Fit_End_idx)); %DataRange
%             results(counter,1)=D1;
%             results(counter,2)=D1;
%             results(counter,3)=k1;
%             results(counter,4)=k1;
%             results(counter,5)=r2;
%         end
%     end
% else
%     for D1=R_D1
%         for D2=R_D2
%             if FitCheck(2)==0
%                 D2=D1;
%             end
%             for k1=R_k1
%                 for k2=R_k2
%                     if FitCheck(4)==0
%                         k2=k1;
%                     end
%                     counter=counter+1;
%                     [X,pro,DepthFlag]=BackDiffsCN_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo); %ImageLength
%                     [r2 rmse] = rsquare(Data(Fit_Start_idx:Fit_End_idx),...
%                         pro(Fit_Start_idx:Fit_End_idx)); %DataRange
%                     results(counter,1)=D1;
%                     results(counter,2)=D2;
%                     results(counter,3)=k1;
%                     results(counter,4)=k2;
%                     results(counter,5)=r2;
%                 end
%             end
%         end
%     end
% end
% [r2,row]=max(results(:,5));
% D1=results(row,1);
% D2=results(row,2);
% k1=results(row,3);
% k2=results(row,4);
%% After thoughts

function [D1,k1]= AutoFit_Crank(handles)
C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
D1=str2double(get(handles.D1, 'String'));
k1=str2double(get(handles.k1, 'String'));
t1=str2double(get(handles.t1, 'String'));
X_temp = handles.Xdata(handles.Fit_Start_idx:handles.Fit_End_idx);
t_tot=t1*3600;
fun = @(p) sum((handles.ProfileData(handles.Fit_Start_idx:handles.Fit_End_idx) - ...
    (C_bg+(C_gas-C_bg)*...
    (erfc(X_temp./(2*sqrt(abs(p(1))*t_tot)))-...
    exp(abs(p(2)).*X_temp+abs(p(1))*t_tot*abs(p(2))^2)...
    .*erfc(...
    X_temp./(2*sqrt(abs(p(1))*t_tot))...
    +abs(p(2))*sqrt(abs(p(1))*t_tot)))...
    )).^2);
h1=abs(k1/D1); %5.5e3
%starting guess
pguess = [abs(D1),h1];
%optimise
[p,fminres] = fminsearch(fun,pguess);
% [p,fminres] = fminsearch(fun,pguess,optimset('Display','iter'));
D1=abs(p(1));
k1=abs(p(1)*p(2));

function [D1,k1]= AutoFit_PlaneSheet(handles)
C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
D1=str2double(get(handles.D1, 'String'));
k1=str2double(get(handles.k1, 'String'));
t1=str2double(get(handles.t1, 'String'));
X_temp = handles.Xdata(handles.Fit_Start_idx:handles.Fit_End_idx);
t_tot=t1*3600;
% MP=str2double(get(handles.MirrorPlane,'String'))*1e-6;
% [val_mirror MP_idx] = min(abs(X-MP));
fun = @(p) sum((...
    handles.ProfileData(handles.Fit_Start_idx:handles.Fit_End_idx)-...
    PlaneSheet_AutoCaller(p,handles)).^2);
%starting guess
pguess = [abs(D1),k1];
%optimise
[p,fminres] = fminsearch(fun,pguess);
% [p,fminres] = fminsearch(fun,pguess,optimset('Display','iter'));
D1=abs(p(1));
k1=abs(p(2));

function [pro]=PlaneSheet_AutoCaller(p,handles)
set(handles.WarningBox,'Visible','off')
C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
dx=ProfileLength*1e-6/(PixelNo-1); %spacial step
X=0:dx:ProfileLength*1e-6; %domain

D1=p(1);
k1=p(2);
t1=str2double(get(handles.t1, 'String'));
t_tot=t1*3600;
MP=str2double(get(handles.MirrorPlane,'String'))*1e-6;
[val_mirror MP_idx] = min(abs(X-MP));
if D1<0 || k1<0
    C_prime=ones(PixelNo);
else
    L=MP*k1/D1;
    RootNo=round(3*(L+10)^.75);
    beta=zeros(1,RootNo);
    if L>30
        HuntStep=0.01;
    else
        HuntStep=0.01;
    end
    fun1 = @(beta) (beta*tan(beta)-L)^2;
    for k=1:RootNo
        beta(k)=fminbnd(fun1,pi*(k-1),pi*(k-0.5));
    end
    for i=1:length(beta)
        RootSum(i,:)=2*L*cos(beta(i)*(X-MP)/MP)*...
            exp(-beta(i)^2*D1*t_tot/MP^2)/((beta(i)^2+L^2+L)*cos(beta(i)));
    end
    C_prime=1-sum(RootSum,1);
end
pro=C_prime(handles.Fit_Start_idx:handles.Fit_End_idx)*(C_gas-C_bg)+C_bg;
if length(pro)>2*MP_idx
    pro(2*MP_idx:end)=0;
end


function [D1,k1]= AutoFit_BackCrank(C_gas,C_bg,D1,k1,t1,t2,ProfileData,...
    Xdata,Fit_Start_idx,Fit_End_idx)
X_temp = Xdata(Fit_Start_idx:Fit_End_idx);
t_tot=(t1+t2)*3600;
phi=t2/(t1+t2);
% abs values used to ensure sqrt(p) is real
fun = @(p) sum((...
    ProfileData(Fit_Start_idx:Fit_End_idx)-...
    (C_bg+(C_gas-C_bg)*(...
    (erfc(X_temp./(2*sqrt(abs(p(1)*t_tot))))-exp(abs(p(2)).*X_temp+abs(p(1))*t_tot*abs(p(2))^2).*erfc(...
    X_temp./(2*sqrt(abs(p(1)*t_tot)))+abs(p(2))*sqrt(abs(abs(p(1))*t_tot)))) - ...
    (erfc(X_temp./(2*sqrt(abs(p(1)*t_tot*phi))))-exp(abs(p(2)).*X_temp+abs(p(1))*t_tot*phi*abs(p(2))^2).*erfc(...
    X_temp./(2*sqrt(abs(p(1)*t_tot*phi)))+abs(p(2))*sqrt(abs(p(1)*t_tot*phi))))...
    ))...
    ).^2);
h1=abs(k1/D1); %5.5e3
%starting guess
pguess = [abs(D1),h1];
%optimise
[p,fminres] = fminsearch(fun,pguess);
D1=p(1);
k1=p(1)*p(2);

%AutoFit_BackDiffs
function [D1,k1,D2,k2]= AutoFit_BackDiffs(C_gas,C_bg,D1,k1,D2,k2,t1,t2,ProfileData,...
    Xdata,Fit_Start_idx,Fit_End_idx,ProfileLength,PixelNo,handles,FitCheck)
ran(1)=Fit_Start_idx;
ran(2)=Fit_End_idx;

b(1)=C_gas;
b(2)=C_bg;
b(3)=t1;
b(4)=t2;
b(5)=ProfileLength;
b(6)=PixelNo;
PS=get(handles.PlaneSheet,'Value');
MP_idx=handles.MP_idx;
Mask=zeros(size(Xdata));
% Mask(Fit_Start_idx:Fit_End_idx)=1;
fun = @(p) sum((...
    ProfileData(ran(1):ran(2))-BackDiffs_AutoCaller(p,b,ran,FitCheck,PS,MP_idx)...
    ).^2);
%What about if you know D1 and k1 and are looking for D2 and k2?
if sum(FitCheck)==4
    pguess = [D1,D2,k1,k2];
elseif FitCheck(2)==0
    pguess = [D1,k1,k2];
elseif FitCheck(4)==0
    pguess = [D1,D2,k1];
end
[p,fminres] = fminsearch(fun,pguess,optimset('TolFun',1e-8));
% [p,fminres] = fseminf(fun,pguess,0);
% [p,resnorm] = lsqcurvefit(@fun,pguess,BackDiffs_AutoCaller(p,b,ran,FitCheck),ProfileData(ran(1):ran(2)));

if sum(FitCheck)==4
    D1=p(1);D2=p(2);k1=p(3);k2=p(4);
elseif FitCheck(2)==0
    D1=p(1);D2=p(1);k1=p(2);k2=p(3);
elseif FitCheck(4)==0
    D1=p(1);D2=p(2);k1=p(3);k2=p(3);
end

% --- Executes during object creation, after setting all properties.
function Align_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text30_CreateFcn(hObject, eventdata, handles)

function GridLines_ClickedCallback(hObject, eventdata, handles)
grid (get(handles.GridLines,'State'))
guidata(hObject, handles);

function D1_text_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in ExportFigure.
function ExportFigure_Callback(hObject, eventdata, handles)
%  [f p]=uiputfile('*.pdf','Save as PDF');
D1=str2double(get(handles.D1,'String'));
k1=str2double(get(handles.k1,'String'));
oldAxe=gca;
newFig=figure;
set(newFig,'Color',[1 1 1]);
set(newFig,'WindowStyle','normal')
set(newFig,'PaperPositionMode','auto');
set(newFig,'PaperOrientation','landscape');
set(newFig,'Position',[60 50 1200 800]);
set(gcf,'Position',[60 50 1200 800]);
set(gcf,'renderer','painters')

%
% hc  = get(handles.MainAxes,'children')
% hgc = get(hc, 'children')
axesObject2=copyobj(oldAxe,newFig);


hline=findobj(gcf, 'type', 'line');
% legend(hline([1,end-1,end]),{'Exchange Profiles','Max. Surface IF','Max. Peak Depth'}) % change order of legend entries
set(hline,'LineWidth',1.5);

set(legend,'Position',[0.67 0.82 0.10 0.03])
set(gca,'fontunits','points')
switch handles.CurrentPlot
    case 'Mask'
    case 'DataPlotOnly'
    case 'CoupledUncertainty'
        set(gca,'position',[0.14 0.12 0.8 0.82])
        set(hline(end),'LineWidth',2)
        names={'Data','Fits'};
        order=[1 length(hline)];
        legend(hline(order),names);
        set(legend,'Interpreter','latex');
        set(legend,'Position',[0.75 0.82 0.10 0.03])
        set(legend,'Box','off')
        for i=1:length(hline);
            set(hline(i),'Color',[.1 .6 1],'Linewidth',0.5);
        end
        set(hline(1),'Color',[1 .1 .1],'Linewidth',1.5);
        if length(get(handles.XprimeTog,'State'))==2;
            xlabel('Normalised depth, $x''$');
            set(gca,'xtick',[0:0.5:3]);
            %         set(gca,'xlables',[0:0.5:3]);
        else
            xlabel('Normalised Depth, $x''$');
        end
        if length(get(handles.CprimeY,'State'))==3;
            ylabel('Isotopic Fraction, $C$');
        else
            ylabel('Normalised Isotopic Fraction, $C''$');
        end
        set(gca,'YColor',[0.2 0.2 0.2],...
            'XColor',[0.2 0.2 0.2]);
        uistack(hline(end), 'top')
        xlim([0 3])
    case 'Contourf'
        set(gca,'position',[0.05 0.12 0.8 0.82])
        handles.colo=colorbar; ylabel(handles.colo,'Standard Deviation of Fit Residuals');
        set(handles.colo, 'XTick', [0, 0.002,0.004,0.006,0.008, 0.01])
        set(handles.colo, 'XTickLabel', {'0', '0.002','0.004','0.006','0.008', '0.01+'})
        set(handles.colo, 'TickLabelInterpreter','latex')
        set(handles.colo,'fontsize',24)
        colormap(flipud(jet))
        %          
        xlabel(['Self diffusivity, $D^*$ / 10$^{',num2str(round(log10(D1))),'}$m$^2$s$^{-1}$'])
        ylabel(['Effective exchange coefficient, $k^*$ / 10$^{',num2str(round(log10(k1))),'}$m s$^{-1}$'])
        axis square
        set(legend,'Interpreter','latex');
        set(get(gca,'xlabel'),'Interpreter','Latex');
        set(get(gca,'ylabel'),'Interpreter','Latex');
    case 'Generate'
        colormap(hot);
        set(gca,'position',[0.05 0.12 0.8 0.82])
        set(hline,'LineWidth',2.5)
        delete(hline(1))
        set(hline(2),'color',[0 1 0],...
            'linewidth',1.5,...
            'linestyle','-');
        set(gca,'xcolor',[0 0 0])
        set(gca,'ycolor',[0 0 0])
        xlabel('Depth, $x$ /$\mu$m');
        ylabel('Isotopic Fraction, $C$');
        legend('Isotopic Fraction','Surface Position')
        set(legend,'Interpreter','latex');
        set(legend,'Position',[0.49 0.83 0.21 0.08])
        set(legend,'color',[0.8 0.8 0.8]);
        set(legend,'edgecolor',[0.8 0.8 0.8]);
        set(legend,'fontsize',24);
        colo=colorbar; ylabel(colo,'Isotopic Fraction, $C$');
        set(colo,...
            'TickLabelInterpreter','latex',...
            'fontsize',24);
        %         set(gcf,'renderer','opengl')
        %         set(colo,'Position',[.82 .12 .03 .63]);
        %         set(legend,'Box','off')
        
    case 'PlotButton'
        set(gca,'position',[0.14 0.12 0.8 0.82])
        set(hline(1:2),'LineWidth',1)
        set(hline(end-1),'LineWidth',1.5)
        set(hline(end),'MarkerSize',5)
        set(hline(end),'LineWidth',1)
        
        if get(handles.ErrorCheck, 'Value')==0;
            legend('Data','Simulation')
        else
            legend('Data','Simulation','Residual')
        end
        set(legend,'Interpreter','latex');
        set(legend,'Position',[0.75 0.82 0.10 0.03])
        set(legend,'Box','off')
        xlabel('Depth /$\mu$m');
        if length(get(handles.CprimeY,'State'))==3;
            ylabel('Isotopic Fraction, C');
        else
            ylabel('Normalised Isotopic Fraction, $C''$');
        end
        set(get(gca,'xlabel'),'Interpreter','Latex');
        set(get(gca,'ylabel'),'Interpreter','Latex');
        if length(get(handles.XprimeTog,'State'))==2;
            xlabel('Normalised depth, $x''$');
            set(gca,'xtick',[0:0.5:3]);
            %         set(gca,'xlables',[0:0.5:3]);
        else
            xlabel('Depth, $x$ /$\mu$m');
        end
        
    case 'Fit'
        set(gca,'position',[0.14 0.12 0.8 0.82])
        set(hline(1:2),'LineWidth',1)
        set(hline(end-1),'LineWidth',1.5)
        set(hline(end),'MarkerSize',5)
        set(hline(end),'LineWidth',1)
        if get(handles.ErrorCheck, 'Value')==0;
            legend('Data','Simulation')
        else
            legend('Data','Simulation','Residual')
        end
        set(legend,'Interpreter','latex');
        set(legend,'Position',[0.75 0.82 0.10 0.03])
        set(legend,'Box','off')
        if length(get(handles.XprimeTog,'State'))==2;
            xlabel('Normalised depth, $x''$');
            set(gca,'xtick',[0:0.5:3]);
            %         set(gca,'xlables',[0:0.5:3]);
        else
            xlabel('Depth, $x$ /$\mu$m');
        end
        if length(get(handles.CprimeY,'State'))==3;
            ylabel('Isotopic Fraction, $C$');
        else
            ylabel('Normalised Isotopic Fraction, $C''$');
        end
        set(gca,'YColor',[0.2 0.2 0.2],...
            'XColor',[0.2 0.2 0.2]);
end
set(get(gca,'xlabel'),'Interpreter','Latex');
set(get(gca,'ylabel'),'Interpreter','Latex');
set(gca,'Parent',newFig,...
    'LineWidth',1.5,...
    'FontSize',24,...
    'TickLabelInterpreter','latex');
% set(gca,'FontUnits','Normalized','FontSize',0.05);
% disp('You are now in control of the figure and variables (type "return" to exit)')
%  
guidata(hObject, handles);
% xlim([0 10]);ylim([0 10]);

% h=gca;
% set(gca,'linewidth',2)
% hline=findobj(gcf, 'type', 'line');
% axes1 = axes('Parent',newFig,'LineWidth',2,'YColor',[0.2 0.2 0.2],...
%     'XColor',[0.2 0.2 0.2],...
%     'FontSize',16,...
%     'Position',[0.1 0.1 0.7 0.7]);
% printpreview(newFig);
% print(newFig,['-dpdf','-A5'],[p f]);
% close(newFig);
% hObject    handle to ExportFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% plot2svg

% --------------------------------------------------------------------
% function CprimeY_ClickedCallback(hObject, eventdata, handles)
%
% hObject    handle to CprimeY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function CprimeY_ClickedCallback(hObject, eventdata, handles)
if length(get(handles.CprimeY,'State'))==3;
    set(handles.GBplot,'State','off')
end
[handles]=PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% --------------------------------------------------------------------
function XprimeTog_ClickedCallback(hObject, eventdata, handles)
if length(get(handles.XprimeTog,'State'))==2;
    set(handles.GBplot,'State','off')
end
[handles]=PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);
% --------------------------------------------------------------------
function GBplot_ClickedCallback(hObject, eventdata, handles)
if length(get(handles.GBplot,'State'))==2;
    set(handles.CprimeY,'State','on')
    set(handles.XprimeTog,'State','off')
end
[handles]=PlotButton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes when uipanel2 is resized.
function uipanel2_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function MainAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MainAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function uipanel7_ButtonDownFcn(hObject, eventdata, handles)
if ~isfield(handles,'Secrets')
    handles.Secrets=1;
else
    handles.Secrets= handles.Secrets+1;
    if  handles.Secrets==4
        set(handles.Advanced,'Visible','on')
        %         set(handles.Build,'Visible','on')
    end
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)







% --- Executes on button press in Go.
function Go_Callback(hObject, eventdata, handles)
disp('Going...')
PixelNo=800;
set(handles.PixelNo,'String',PixelNo);
xStarLen=4;
ProfileLength=round(xStarLen*(1e6*(4*1e-13*3600)^0.5),3,'significant');
set(handles.ProfileLength,'String',ProfileLength);
[handles]=PlotButton_Callback(hObject, eventdata, handles);
tic
global SurfConc
global thetaFine
global G_pro
global G_pro3D
global k2
global G_x
global G_xStar
global G_ProData
global G_Peak
global G_xPeak_idx
global theta
D1=1e-13;
D2=1e-13;
k1=1e-5;
k2=1e-5;
set(handles.D1,'String',D1);
set(handles.D2,'String',D2);
set(handles.k1,'String',k1);
set(handles.k2,'String',k2);

t1=1;
t2=0;
set(handles.t1,'String',t1);
set(handles.t2,'String',t2);
C_bg=str2double(get(handles.C_bg,'String'));
C_gas=str2double(get(handles.C_gas,'String'));
t_tot=t1+t2;


set(handles.PixelNo,'String',PixelNo);
G_x=handles.X;
G_xStar=G_x/(4*D1*t_tot*3600)^0.5;

theta=[0:1/6:1];%[-3:0.1:3];
ls=[-20:1];
k2=k1*2.^ls;  %10.^ls;
k2(1)=0;k2(ls==0)=k2(ls==0)*1.01;
G_pro=zeros(length(k2),PixelNo);
G_pro3D=zeros(length(theta),length(k2),PixelNo);
for i=1:length(theta)
    if theta(i)==0 %%% All In Diffusion
        set(handles.t1,'String',1);
        set(handles.t2,'String',0);
        set(handles.k2,'String',k1*1.1);
        [handles]=PlotButton_Callback(hObject, eventdata, handles);
        G_pro(:,:)=meshgrid((handles.pro-C_bg)/(C_gas-C_bg),k2);
    elseif theta(i)==1 %%% All Back Diffusion
        G_pro=zeros(length(k2),PixelNo);
    else
        set(handles.t1,'String',t_tot*(1-theta(i)));
        set(handles.t2,'String',t_tot*(theta(i)));
        for j=1:length(k2)
            set(handles.k2,'String',k2(j))
            [handles]=PlotButton_Callback(hObject, eventdata, handles);
            G_pro(j,:)=(handles.pro-C_bg)/(C_gas-C_bg);
        end
    end
    G_pro3D(i,:,:)=G_pro;
end
toc

thetaFine=[0:0.005:0.195,0.2:0.04:0.76,0.8:0.005:1];
set(handles.k1,'String',1e-04);
set(handles.k2,'String',0);
for i=1:length(thetaFine)
    set(handles.t1,'String',t_tot*(1-thetaFine(i)));
    set(handles.t2,'String',t_tot*(thetaFine(i)));
    [handles]=PlotButton_Callback(hObject, eventdata, handles);
    SurfConc(i)=(handles.pro(1)-C_bg)/(C_gas-C_bg);
end
SurfConc(1)=1;SurfConc(end)=0; %correction for crappy erfc tables
toc
set(handles.k1,'String',1e-05);
set(handles.k2,'String',1e-04);
for i=1:length(thetaFine)
    set(handles.t1,'String',t_tot*(1-thetaFine(i)));
    set(handles.t2,'String',t_tot*(thetaFine(i)));
    [handles]=PlotButton_Callback(hObject, eventdata, handles);
    [G_xPeak_idx(i),G_xPeak_idx(i)]=max(handles.pro);
    %     k1=1e-05;k2=1e-04;
    %     t1=t_tot*(1-thetaFine(i));t2=t_tot*(thetaFine(i));
    %     [X,pro,DepthFlag]=BackDiffsCN_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo);
    %     [G_xPeak_idx(i),G_xPeak_idx(i)]=max(pro);
end

%%% Plotting
figure('Color',[1 1 1]);
plot31=plot3(G_xStar',zeros(length(G_xStar)),permute(G_pro3D(1,1,:),[3 1 2]));

G_xPeak_idx(1)=1;G_xPeak_idx(end)=G_xPeak_idx(end-1);
plot3(G_xStar(G_xPeak_idx)',thetaFine,meshgrid(0,thetaFine),'linewidth',2);hold on
plot3(meshgrid(0,thetaFine)',meshgrid(thetaFine,0),SurfConc,'linewidth',2)
for k=1:length(k2)
    if k==1 || k==length(k2)
        hold on; plot3(meshgrid(G_xStar,theta)',meshgrid(theta,G_xStar),...
            permute(G_pro3D(:,k,:),[3 1 2]),'k','linewidth',2);
    else
        hold on; plot3(meshgrid(G_xStar,theta)',meshgrid(theta,G_xStar),...
            permute(G_pro3D(:,k,:),[3 1 2]),'k');
    end
end
xlabel('Normalised Depth, x''','Rotation',-8);
ylabel('Exchange Time Ratio, \theta','Rotation',11);
zlabel('Normalised Concentration, C''');box on
set(gca,'YTick',theta);
set(gca,'YTickLabel',{'0','1/6','1/3','1/2','2/3','5/6','1'});
set(gca,'ZTick',[0 1]);
set(gca,'ZTickLabel',{'0','a_0'});
xlim([0 2.5]);ylim([0 1]);zlim([0 1]);
hold off
set(gca,'YDir','Reverse');

Leg={'Max. Peak Depth','Max. Surface Conc.','Exchange Profiles'};
legend(Leg,'Position',[.72 .67 .10 .085]);
set(gca,'FontSize',16);
ax=gca;
ax.XLabel.Position = [1.3 1 -.1];
ax.YLabel.Position = [2.6 .4 -.1];
view(45, 15);
annotation(gcf,'textarrow',[0.57 0.5],[0.48 0.510],...
    'String',{'\lambda=\infty'},'LineWidth',1,...
    'HeadWidth',6,'HeadStyle','cback1',...
    'HeadLength',6,'FontWeight','bold',...
    'FontSize',16);
annotation(gcf,'textarrow',[0.567 0.537],[0.549 0.563],...
    'String',{'\lambda=0'},'LineWidth',1,...
    'HeadWidth',6,'HeadStyle','cback1',...
    'HeadLength',6,'FontWeight','bold',...
    'FontSize',16);
set(gcf,'Position', [100, 100, 1049, 895]);
%export_fig('D:\PhD\BackDiffusion\MatLab Figures\OmniGraphHex_ef.pdf')
% set(gcf,'Legend','String',
% set(plot31(240),'DisplayName','Max Peak Depth, x*_{peak}');
% set(plot31(241),'DisplayName','Max. Surface Conc., C''_{max}');
% set(plot31(242),'DisplayName','Back Exchange Profiles_{}^{}','LineWidth',2);

% figure(11);plot(thetaFine,SurfConc(:,1)/max(SurfConc(:,1)));ylim([0 1]);
% plot(thetaFine,0.5*SurfConc(:,1)/SurfConc((1+end)/2,1),thetaFine,acos((thetaFine*2)-1)/pi());
% ylim([0 1]);set(gca,'xTick',([0:0.1:1]));grid;axis square
%%% divide everything everywhere by max(Crank(theta=0))
% k2_start=1e-11;
% ls=[-10:10];%[-3:0.1:3];
% k2=k1*2.^ls;%10.^ls;
% % theta=[0:0.01:1];
% % t_tot=100;
% % t1=[10:10:100];
% set(handles.k2,'String',k2_start)
% % for i=1:length(t1)
% % %     set(handles.t1,'String',t1(i));
% % for i=1:length(theta)
% %     set(handles.t1,'String',t_tot*(1-theta(i)))
% %     set(handles.t2,'String',t_tot*(theta(i)))
% for i=1:length(k2)
% %         k2(i)=k2_start*1.1^(i-1);
% %         k2(i)=50*i*k2_start;
%         set(handles.k2,'String',k2(i))
%     [handles]=PlotButton_Callback(hObject, eventdata, handles);
%     [G_Peak(i,1),G_Peak(i,2)]=max((handles.pro-C_bg)/(C_gas-C_bg));
%     SurfConc(i)=(handles.pro(1)-C_bg)/(C_gas-C_bg);
%     G_pro(i,:)=(handles.pro-C_bg)/(C_gas-C_bg);
%
%     %     set(handles.k2,'String',1e-12+(i-1)*1e-13)
%
% end
%
% % global G_x; global k2; global G_pro; global SurfConc;kRAT=k2/4.4e-10;
% % surf(G_x*1e6, kRAT, G_pro);xlabel('Depth /um');ylabel('Ratio k2/k1');zlabel('C''');
% figure; plot(theta,G_Peak(:,1),'k');xlabel('Theta');ylabel('C''_{peak}');
% figure; plot(1e6*G_x(1:308),G_pro(:,1:308),'k');xlabel('NormalisedDepth /um');ylabel('C''');
% hold on;plot(1e6*G_x(1:308),G_pro(61,1:308),'g:','LineWidth',1.5)
% hold on;plot(1e6*G_x(1:308),G_pro(41,1:308),'r:','LineWidth',1.5)

% figure; plot(G_xStar,G_pro,'k');xlabel('Normalised Depth, x*');ylabel('Normalised Concentration, C''');
% hold on;plot(G_xStar,G_pro(31,:),'r:','LineWidth',1.5)
% hold on;plot(G_xStar,G_pro(51,:),'g:','LineWidth',1.5)
toc
sound(1,1000);pause(0.2);sound(1,1000);
disp('Gone')
guidata(hObject, handles);



% --- Executes on button press in Build.
function Build_Callback(hObject, eventdata, handles)
[handles]=PlotButton_Callback(hObject, eventdata, handles);
% delete(gca)
% handles.AlignAngle_FiAl_num=0;%So that align isn't auto called

set(handles.AlignAngle,'String',0);
pro=zeros(size(handles.pro));
SurfPos=str2double(get(handles.SurfPos, 'String'));
ProLen=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));

SurfPos_idx=1+round(SurfPos/(ProLen/(PixelNo-1)));
pro(SurfPos_idx:end)=handles.pro(1:end-SurfPos_idx+1);
counts=100;
absErr=str2double(get(handles.absErr,'String'));
relErr=str2double(get(handles.relErr,'String'));

pro18=counts*pro;
% pro16=pro18/pro-pro18;
pro16=counts*ones(size(pro))-pro18;
pro16(pro==0)=0;

handles.O16Flag=1;handles.O18Flag=1;

handles.O18image=...
    (1-relErr*(1-2*rand(length(pro18)))).*...%
    (absErr*counts*(1-2*rand(length(pro18)))+...
    meshgrid(pro18,pro18));
handles.O18image(handles.O18image<0)=0;
handles.O18image_Or=handles.O18image;

handles.O16image=(1-relErr*2*(1-2*rand(length(pro16)))).*...
    (absErr*counts*2*(1-2*rand(length(pro16)))+...
    meshgrid(pro16,pro16));
handles.O16image(handles.O16image<0)=0;
handles.O16image_Or=handles.O16image;


handles.NormO18image_Or=handles.O18image./(handles.O18image+handles.O16image);
handles.O16pO18image_Or=handles.O18image+handles.O16image;

set(handles.ProfDataSaveName,'String','Dummy');
handles.PathName=[cd,'\'];
handles.FileName='Dummy';
handles.SaveName_Or='Dummy';
imagepro=meshgrid(pro,pro);

% imagesc(imagepro);
% colormap(hot);
% handles.colo=colorbar; ylabel(handles.colo,'Counts');
% axis square; set(gca, 'XTick', [],'YTick', []);
% set(handles.Xlimit,'Visible','off');set(handles.Ylimit,'Visible','off');
% set(handles.Fit_End,'String',roundsf(...
%     0.5*str2double(get(handles.ProfileLength,'String')),3,'floor'));
% pause(0.5)

si=[0.35 0.35];
sp(1)=subplot (3,2,3); imagesc(handles.O16image);
set(gca,'position',[0.03 0.45 si]);title('O16');
sp(2)=subplot (3,2,4); imagesc(handles.O18image);
set(gca,'position',[0.4 0.45 si]);title('O18');
sp(3)=subplot (3,2,5); imagesc(handles.O16pO18image_Or);
set(gca,'position',[0.03 0.03 si]);title('O16 + O18');
sp(4)=subplot (3,2,6); imagesc(handles.NormO18image_Or);
set(gca,'position',[0.4 0.03 si]);title('O18 Norm.');
% set(gcf,'position',handles.ImagePlotPos);
set(sp, 'XTick', [],'YTick', []);
set(handles.Xlimit,'Visible','off');set(handles.Ylimit,'Visible','off');
axis(sp,'square');colormap(hot);
disp('Built!')
waitforbuttonpress
delete(sp)
%  
% set(sp,'position',[.01  .1  .3 .3]);

% colormap(hot);
% handles.colo=colorbar; ylabel(handles.colo,'Counts');

% set(handles.Xlimit,'Visible','off');set(handles.Ylimit,'Visible','off');
%
% pause
% imagesc(imagepro-image18./(image16+image18))
% colormap(hot);
% handles.colo=colorbar; ylabel(handles.colo,'Counts');
% axis square; set(gca, 'XTick', [],'YTick', []);
% set(handles.Xlimit,'Visible','off');set(handles.Ylimit,'Visible','off');

guidata(hObject, handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over xStar_d.
function xStar_d_ButtonDownFcn(hObject, eventdata, handles)

function Xstar_ButtonDownFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
function uipanel2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D1=str2num(get(handles.D1,'String'));
t1=str2num(get(handles.t1,'String'));
t2=str2num(get(handles.t2,'String'));
t_tot=3600*(t1+t2);
set(handles.Fit_End,'String',roundsf(4e6*sqrt(D1*t_tot),3,'round'));
set(handles.xStar_d,'String',2);
guidata(hObject, handles);



function absErr_Callback(hObject, eventdata, handles)
% hObject    handle to absErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of absErr as text
%        str2double(get(hObject,'String')) returns contents of absErr as a double


% --- Executes during object creation, after setting all properties.
function absErr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to absErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function relErr_Callback(hObject, eventdata, handles)
% hObject    handle to relErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of relErr as text
%        str2double(get(hObject,'String')) returns contents of relErr as a double


% --- Executes during object creation, after setting all properties.
function relErr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to relErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Uncertainty.
function Uncertainty_Callback(hObject, eventdata, handles)
set(handles.Uncertainty,'String','Calculating...');pause(0.001);
% hObject    handle to Uncertainty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
C_bg=str2double(get(handles.C_bg, 'String'));
C_gas=str2double(get(handles.C_gas, 'String'));
D1=str2double(get(handles.D1, 'String'));
D2=str2double(get(handles.D2, 'String'));
k1=str2double(get(handles.k1, 'String'));
k2=str2double(get(handles.k2, 'String'));
t1=str2double(get(handles.t1, 'String'));
t2=str2double(get(handles.t2, 'String'));
ProfileLength=str2double(get(handles.ProfileLength, 'String'));
PixelNo=str2double(get(handles.PixelNo, 'String'));
t_tot=(t1+t2)*3600;
theta=t2/(t1+t2);
if get(handles.k2FitCheck, 'Value')==1
    %k2 and k1
    ls1=[-1:0.2:1];
    k1_n=k1*10.^ls1;
    ls2=[-2:0.4:2];
    k2_n=k2*10.^ls2;
    for i=1:length(k2_n)
        for j=1:length(k1_n)
            set(handles.k2,'String',k2_n(i));
            set(handles.k1,'String',k1_n(j));
            [handles]=PlotButton_Callback(hObject, eventdata, handles);
            RSqMap(i,j)=handles.r2;
            RMSeMap(i,j)=handles.rmse;
        end
        disp([num2str(i*length(ls1)),'/',num2str(length(ls1)*length(ls2))])
    end
    set(handles.k1,'String',roundsf(k1,3,'round'));
    set(handles.k2,'String',roundsf(k2,3,'round'));
    [handles]=PlotButton_Callback(hObject, eventdata, handles);
    set(handles.Xlimit,'Visible','off');
    set(handles.Ylimit,'Visible','off');
    k1_nI=linspace(min(k1_n),max(k1_n),200);
    k2_nI=linspace(min(k2_n),max(k2_n),200);
    contourf(k1_nI,k2_nI,interp2(k1_n,k2_n,RMSeMap,k1_nI,k2_nI'),64,'LineStyle','none')
    %     contourf(meshgrid(k1_n,k2_n),meshgrid(k2_n,k1_n)',RMSeMap)
    handles.colo=colorbar; ylabel(handles.colo,'Standard Deviation of Fit Residuals');
    set(handles.colo,'FontSize',10);
    caxis([0,2e-2]);
    set(handles.colo, 'XTick', [0, 0.004,0.008,0.012,0.016, 0.02])
    set(handles.colo, 'XTickLabel', {'0', '0.004','0.008','0.0012','0.0016', '0.02+'})
    set(gca,'position',handles.ImagePlotPos);
    colormap(flipud(parula))
    set(gca,'XScale','log');set(gca,'YScale','log');
    xlabel('Effective exchange coefficient, k1* / m s^{-1}')
    ylabel('Effective exchange coefficient, k2* / m s^{-1}')
    axis square
    
else
    %D1 and k1
    % should you change the fit region as D varies?????
    %  
    ls1=[0:0.2:2];
    D1_n=D1*ls1;
    ls2=[0:0.2:2];
    k1_n=k1*ls2;
    for i=1:length(k1_n)
        for j=1:length(D1_n)
            set(handles.k1,'String',k1_n(i));
            set(handles.D1,'String',D1_n(j));
            [handles]=PlotButton_Callback(hObject, eventdata, handles);
            RSqMap(i,j)=handles.r2;
            RMSeMap(i,j)=handles.rmse;
        end
        disp([num2str(i*length(ls1)),'/',num2str(length(ls1)*length(ls2)),' residuals quantified'])
    end
    set(handles.D1,'String',roundsf(D1,3,'round'));
    set(handles.k1,'String',roundsf(k1,3,'round'));
    [handles]=PlotButton_Callback(hObject, eventdata, handles);
    set(handles.Xlimit,'Visible','off');
    set(handles.Ylimit,'Visible','off');
    
    %% Calculate the coupled
    ProTemp_p=(handles.ProfileData(1:handles.Fit_End_idx)-C_bg)/(C_gas-C_bg);
    XTemp=handles.X(1:handles.Fit_End_idx);
%      
    if t2==0
        for i=1:length(handles.ProfileData(1:handles.Fit_End_idx))
            funCoupledDk = @(g) ((g-D1)/D1)^2+...
                ((Find_k_from_D_at_pCrank(g,t_tot,ProTemp_p(i),XTemp(i))-k1)/k1)^2;
            
            %% Optimise
            pguess = [D1];
            [g,fminres] = fminsearch(funCoupledDk,pguess,optimset('TolFun',1e-8));
            Dandk(i,1)=[g];
            Dandk(i,2)=Find_k_from_D_at_pCrank(g,t_tot,ProTemp_p(i),XTemp(i));
            NewCranks(i,:)=C_bg+(C_gas-C_bg)*...
                (erfc(handles.X./(2*sqrt(g(1)*t_tot)))-...
                exp(Dandk(i,2)/g(1).*handles.X+t_tot*Dandk(i,2)^2/g(1)).*erfc(handles.X./(2*sqrt(g(1)*t_tot))+Dandk(i,2)*sqrt(t_tot/g(1))));
            h=plot(handles.X/(2*sqrt(D1*t_tot)),NewCranks');
            hold on;plot(handles.Xdata(1:i)/(2*sqrt(D1*t_tot)),handles.ProfileData(1:i),'k','linewidth',1.5);hold off;
            pause(0.01)
        end
    else
        for i=1:length(handles.ProfileData(1:handles.Fit_End_idx))
            funCoupledDk = @(g) ((g-D1)/D1)^2+...
                ((Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))-k1)/k1)^2;
            %                         funCoupledDk = @(g) (log(g/D1))^2+...
            %                             (log(Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))/k1))^2;
            funCoupledDk = @(g) sum((ProTemp_p-(...
                (erfc(XTemp./(2*sqrt(g*t_tot)))-...
                exp(Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))/g.*XTemp+t_tot*Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))^2/g).*...
                erfc(XTemp./(2*sqrt(g*t_tot))+Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))*sqrt(t_tot/g)))-...
                (erfc(XTemp./(2*sqrt(g*t_tot*theta)))-...
                exp(Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))/g.*XTemp+t_tot*theta*Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))^2/g).*...
                erfc(XTemp./(2*sqrt(g*t_tot*theta))+Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i))*sqrt(t_tot*theta/g)))...
                )).^2);
            
            %% Optimise
            %             pguess = [D1];
            %             [g,fminres] = fminsearch(funCoupledDk,pguess,optimset('TolFun',1e-30,'TolX',1e-30,'MaxFunEvals',1000,'MaxIter',1000));%'TolX',1e-18,
            [g,fminres] = fminbnd(funCoupledDk,.1*D1,10*D1,optimset('TolFun',1e-18,'TolX',1e-18));%,'MaxFunEvals',1000,'MaxIter',1000));%'TolX',1e-18,
            Dandk(i,1)=[g];
            Dandk(i,2)=Find_k_from_D_at_pBackCrank(g,t_tot,theta,ProTemp_p(i),XTemp(i));
            NewCranks(i,:)=C_bg+(C_gas-C_bg)*(...
                (erfc(handles.X./(2*sqrt(g*t_tot)))-...
                exp(Dandk(i,2)/g.*handles.X+t_tot*Dandk(i,2)^2/g(1)).*erfc(handles.X./(2*sqrt(g*t_tot))+Dandk(i,2)*sqrt(t_tot/g)))-...
                (erfc(handles.X./(2*sqrt(g*t_tot*theta)))-...
                exp(Dandk(i,2)/g.*handles.X+t_tot*theta*Dandk(i,2)^2/g).*erfc(handles.X./(2*sqrt(g*t_tot*theta))+Dandk(i,2)*sqrt(t_tot*theta/g)))...
                );
            h=plot(handles.X/(2*sqrt(D1*t_tot)),NewCranks');xlim([0 3]);
            hold on;plot(handles.Xdata(1:i)/(2*sqrt(D1*t_tot)),handles.ProfileData(1:i),'k','linewidth',1.5);hold off;
            pause(0.01)
            %             Dandk0=meshgrid([D1,k1],Dandk(:,1));subplot(1,2,1);plot((Dandk-Dandk0)./Dandk0);subplot(1,2,2);plot(sum( ((Dandk-Dandk0)./Dandk0).^2 ,2))
        end
    end
    %      
    sorted_Dandk=sort(Dandk);
    RelevantPoints=length(handles.X(1:handles.Fit_End_idx));
    
    D_minus=D1-sorted_Dandk(ceil(RelevantPoints*0.025),1);
    D_plus=sorted_Dandk(floor(RelevantPoints*0.975),1)-D1;
    
    disp(['D1 = ',num2str(round(D1,2,'significant')),' - ',num2str(round(D_minus,2,'significant')),' + ',num2str(round(D_plus,2,'significant'))])
    
    k_minus=k1-sorted_Dandk(ceil(RelevantPoints*0.025),2);
    k_plus=sorted_Dandk(floor(RelevantPoints*0.975),2)-k1;
    
    disp(['k1 = ',num2str(round(k1,2,'significant')),' - ',num2str(round(k_minus,2,'significant')),' + ',num2str(round(k_plus,2,'significant'))])
    %% Plotting
    
    %     set(gca,'position',[0.03 0.03 0.4 0.4])
    handles.CurrentPlot='Contourf';
    scaleD1=round(log10(D1));
    scalek1=round(log10(k1));
    contourf(D1_n/10^scaleD1,k1_n/10^scalek1,RMSeMap,32,'LineStyle',':')
    handles.colo=colorbar; ylabel(handles.colo,'Standard Deviation of Fit Residuals');
    set(handles.colo,'FontSize',10);
    caxis([0,1e-2]);
    set(handles.colo, 'XTick', [0, 0.002,0.004,0.006,0.008, 0.01])
    set(handles.colo, 'XTickLabel', {'0', '0.002','0.004','0.006','0.008', '0.01+'})
    set(gca,'position',handles.ImagePlotPos);
    colormap(flipud(jet))
    xlabel(['Self Diffusivity, D* / 1e',num2str(round(log10(D1))),'m^2s^{-1}'])
    ylabel(['Effective exchange coefficient, k* / 1e',num2str(round(log10(k1))),'m s^{-1}'])
    axis square
    guidata(hObject, handles);
    %      
    waitforbuttonpress
    pause(0.1)
    waitforbuttonpress
    h=plot(handles.X/(2*sqrt(D1*t_tot)),NewCranks',handles.Xdata/(2*sqrt(D1*t_tot)),handles.ProfileData,'k');
    %     hold on;plot();hold off
    set(gca,'position',handles.GraphPlotPos);
    %     title(['Curves required to fit each point up to x*=',get(handles.xStar_d,'String')]);
    xlabel('Normalised Depth, x''');ylabel('Isotopic Fraction');
    handles.CurrentPlot='CoupledUncertainty';
    guidata(hObject, handles);
%      
    guidata(hObject, handles);
    waitforbuttonpress
    pause(0.1)
    waitforbuttonpress
    PlotButton_Callback(hObject, eventdata, handles);
    %     set(gca,'position',[0.03 0.45 0.4 0.4])
end
sound(1,1000);pause(0.2);sound(1,1000);
set(handles.Uncertainty,'String','Uncertainty')
%  


function [k]=Find_k_from_D_at_pCrank(D,t_tot,pointC,pointX);

% pointC=(pointC-0.002)/(0.95-0.002);
if D<=0
    k=10;
else
    A=erfc(pointX/(2*sqrt(D*t_tot)));
    fun_k_from_D = @(g) (pointC-...
        (A-exp(g/D.*pointX+t_tot*g^2/D).*erfc(pointX./(2*sqrt(D*t_tot))+g*sqrt(t_tot/D)))...
        )^2;
    [k,fminres] = fminbnd(fun_k_from_D,0,sqrt(26^2*D/t_tot),optimset('TolX',1e-18,'TolFun',1e-18));
end

function [k]=Find_k_from_D_at_pBackCrank(D,t_tot,theta,pointC,pointX);

% pointC=(pointC-0.002)/(0.95-0.002);
if D<=0
    k=10;
else
    A=erfc(pointX/(2*sqrt(D*t_tot)));
    B=erfc(pointX/(2*sqrt(D*t_tot*theta)));
    t_222=t_tot*theta;
    fun_k_from_D = @(g) (pointC-(...
        (A-exp(g/D.*pointX+t_tot*g^2/D).*erfc(pointX./(2*sqrt(D*t_tot))+g*sqrt(t_tot/D)))-...
        (B-exp(g/D.*pointX+t_222*g^2/D).*erfc(pointX./(2*sqrt(D*t_222))+g*sqrt(t_222/D)))...
        ))^2;
    
    [k,fminres] = fminbnd(fun_k_from_D,0,sqrt(26^2*D/t_tot),optimset('TolX',1e-18,'TolFun',1e-18));
end

% --- Executes during object creation, after setting all properties.
function PlotO16image_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotO16image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in PlaneSheet.
function PlaneSheet_Callback(hObject, eventdata, handles)
if get(handles.PlaneSheet, 'Value')==1;
%     set(handles.MirrorPlane,'String',get(handles.Xlimit,'String'))
    set(handles.MirrorPlane,'Enable','on');
else
    set(handles.MirrorPlane,'Enable','off');
end
[handles]=PlotButton_Callback(hObject, eventdata, handles);%fun

function MirrorPlane_Callback(hObject, eventdata, handles)
[handles]=PlotButton_Callback(hObject, eventdata, handles);
% hObject    handle to MirrorPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MirrorPlane as text
%        str2double(get(hObject,'String')) returns contents of MirrorPlane as a double


% --- Executes during object creation, after setting all properties.
function MirrorPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MirrorPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function PlotProperties(hObject, eventdata, handles)
set(gca,...
    'LineWidth',1,...
    'FontUnits','Normalized',...
    'FontSize',0.034);
% set(get(gca,'xlabel'),'FontSize','Latex');
% set(get(gca,'ylabel'),'Interpreter','Latex');

guidata(hObject, handles);
