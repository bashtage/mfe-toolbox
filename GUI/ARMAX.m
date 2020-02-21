function varargout = ARMAX(varargin)
% ARMAX M-file for ARMAX.fig
%      ARMAX, by itself, creates a new ARMAX or raises the existing
%      singleton*.
%
%      H = ARMAX returns the handle to a new ARMAX or the handle to
%      the existing singleton*.
%
%      ARMAX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARMAX.M with the given input arguments.
%
%      ARMAX('Property','Value',...) creates a new ARMAX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ARMAX_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ARMAX_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ARMAX

% Last Modified by GUIDE v2.5 29-Jan-2007 16:09:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ARMAX_OpeningFcn, ...
    'gui_OutputFcn',  @ARMAX_OutputFcn, ...
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


% --- Executes just before ARMAX is made visible.
function ARMAX_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ARMAX (see VARARGIN)

% Choose default command line output for ARMAX
handles.output = hObject;

handles.model_ARMAX_AR_order=[];
handles.model_ARMAX_MA_order=[];
handles.model_ARMAX_AR_order_string='';
handles.model_ARMAX_MA_order_string='';
handles.y = varargin{1};
handles.model_ARMAX_residuals=[];
handles.ACF_lags=min(24,floor(length(handles.y)/8));
handles.model_InferenceMethod = 1;
handles.model_IncludeConstant=1;
handles.models={};
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes ARMAX wait for user response (see UIRESUME)
% uiwait(handles.figure1);

axes(handles.plot_handle);
clear_axes(handles)
h=plot(handles.y);
set(h,'LineWidth',2,'Color',[0 0 .6])
axis tight;
title('Plot of original data')
AX=axis;
range = AX(4)-AX(3);
AX(3)=AX(3)-.1*range;
AX(4)=AX(4)+.1*range;
axis(AX);

set(handles.heterorobust_menu,'Checked','on')
set(handles.homoerror_menu,'Checked','off')
set(handles.ACF_lags_edit,'String',num2str(handles.ACF_lags))

% --- Outputs from this function are returned to the command line.
function varargout = ARMAX_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function EXOG_edit_Callback(hObject, eventdata, handles)
% hObject    handle to EXOG_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EXOG_edit as text
%        str2double(get(hObject,'String')) returns contents of EXOG_edit as a double


% --- Executes during object creation, after setting all properties.
function EXOG_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EXOG_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in estimate_pushbutton.
function estimate_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to estimate_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Thsi is the main function.  Call ARMACFILTER and then call teh display
%routine to display the results.  Make sure there are no errorr, if there
%are display a message that something went wrrong.
y=handles.y;
ar = handles.model_ARMAX_AR_order;
ma = handles.model_ARMAX_MA_order;
constant = handles.model_IncludeConstant;

% A couple of quick sanity checks
readytorun = true;
notalreadyrun = true;
if (isempty(ar) || all(ar==0)) && (isempty(ma) || all(ma==0))
    readytorun=false;
end

% Check to make sure the this model hasn't been estimated;  if it has,
% do not reestimate but load the results of the previous estimation.
% Probably print a message in th3 window.

% To check, compare teh constant, and the other things,  if a match is
% founf, set readyto run to false, update a message and display the model
% results by changing the current model number.



models = handles.models;
matchedModel=0;
for i=1:length(handles.models)
    match=true;
    thisModel = models{i};
    % conditions for maths
    % constant == modelcontant
    % length(ar)==length(modelar)
    %   if true than check that they are all the same
    % length(ma)==length(modelma)
    %   if true than check that they are all the same

    if constant ~= thisModel.constant
        match=false;
    end

    if length(ar)~=length(thisModel.ARlags)
        match=false;
    elseif ~all(unique(ar)==unique(thisModel.ARlags))
        match=false;
    end

    if length(ma)~=length(thisModel.MAlags)
        match=false;
    elseif ~all(unique(ma)==unique(thisModel.MAlags))
        match=false;
    end

    if match
        matchedModel=i;
        notalreadyrun = false;
    end
end

if matchedModel
    handles.currentModelNumber = matchedModel;
    errors = models{matchedModel}.errors;
end





if readytorun && notalreadyrun
    try
        [parameters, ll, errors, seregression, diagnostics, vcvrobust, vcv, likelihoods, scores]=armaxfilter(y,constant,ar,ma);
        handles.model_ARMAX_residuals=errors;
        modelno=length(handles.models)+1;
        handles.models{modelno}.parameters=parameters;
        handles.models{modelno}.errors=errors;
        handles.models{modelno}.ll=ll;
        handles.models{modelno}.seregression=seregression;
        handles.models{modelno}.diagnostics=diagnostics;
        if handles.model_InferenceMethod
            handles.models{modelno}.covariance=vcvrobust;
        else
            handles.models{modelno}.covariance=vcv;
        end
        handles.models{modelno}.likelihoods=likelihoods;
        handles.models{modelno}.scores=scores;
        handles.models{modelno}.ARlags=ar;
        handles.models{modelno}.MAlags=ma;
        handles.models{modelno}.constant=constant;
        handles.models{modelno}.y=y;
        handles.models{modelno}.yhat=y-errors;
        [aic,bic]=aicsbic(errors,constant,ar,ma);
        handles.models{modelno}.AIC=aic;
        handles.models{modelno}.BIC=bic;
        handles.models{modelno}.K=length(parameters);
        handles.models{modelno}.StdErr=sqrt(diag(handles.models{modelno}.covariance));
        handles.models{modelno}.Tstat=parameters./handles.models{modelno}.StdErr;
        handles.models{modelno}.Pval=2-2*normcdf(abs(handles.models{modelno}.Tstat));

        % Update the string in the popup
        ar = handles.models{modelno}.ARlags;
        ma = handles.models{modelno}.MAlags;
        irregular = false;
        if ~isempty(ar)
            if length(ar)<max(ar)
                irregular = true;
            end
        end
        if ~isempty(ma)
            if length(ma)<max(ma)
                irregular = true;
            end
        end
        this_str = 'ARMA(';
        if ~isempty(max(ar))
            this_str=[this_str num2str(max(ar)) ','];
        else
            this_str=[this_str '0,'];
        end

        if ~isempty(max(ma))
            this_str=[this_str num2str(max(ma)) ')'];
        else
            this_str=[this_str '0)'];
        end

        if constant
            this_str = [this_str ' w/ constant'];
        else
            this_str = [this_str ' w/o constant'];
        end

        if irregular
            this_str = [this_str ' (Irregular)'];
        end
        handles.models{modelno}.StringID = this_str;

        str =[];
        for i=1:length(handles.models)
            str=strvcat(str,handles.models{i}.StringID);
        end
        set(handles.model_selector_popup,'String',str,'Value',size(str,1));

        set(handles.message_text,'String','The model has been sucessfully estimated.','ForegroundColor',[0 0 .8])
        handles.currentModelNumber=length(models);
    catch
        err=lasterror;
        errmsg = err.message;
        set(handles.message_text,'String',['There was an error: ' errmsg],'ForegroundColor',[.9 0 0])
    end
elseif ~readytorun
    set(handles.message_text,'String','The model is not ready to run.  Please check the inputs.','ForegroundColor',[.9 0 0])
elseif ~notalreadyrun
    set(handles.message_text,'String','Model previously estimated; results loaded.','ForegroundColor',[0 0 .8])
end

if readytorun || ~notalreadyrun
    % Plot the residuals and update the radio button
    set(handles.resid_raw_radio,'Value',1)
    axes(handles.plot_handle);
    clear_axes(handles)
    h=plot(errors);
    set(h,'LineWidth',2,'Color',[0 0 .6])
    axis tight;
    title('Plot of model errors')
    AX=axis;
    range = AX(4)-AX(3);
    AX(3)=AX(3)-.1*range;
    AX(4)=AX(4)+.1*range;
    axis(AX);
end

% Store the results in the display model structure
guidata(hObject, handles);
% Call the display function with the display model structure and the model
% number
if readytorun
    ARMAX_viewer(handles.models)
elseif ~notalreadyrun
    ARMAX_viewer(handles.models,matchedModel)
end




function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in close_pushbutton.
function close_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to close_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos_size = get(handles.figure1,'Position');
% Call modaldlg with the argument 'Position'.
%user_response ='Yes';
user_response = ARMAX_close_dialog('Title','Confirm Close');
switch user_response
    case {'No'}
        %     take no action
    case 'Yes'
        delete(handles.figure1)
end

% --- Executes during object creation, after setting all properties.
function MA_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MA_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function reset_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.model_ARMAX_MA_order= [];
handles.model_ARMAX_AR_order = [];
% Update handles structure
guidata(hObject, handles);

set(findobj('Tag','MA_lags_edit'),'String','');
set(findobj('Tag','AR_lags_edit'),'String','');
set(findobj('Tag','EXOG_edit'),'String',' ');
set(findobj('Tag','orig_raw_radio'),'Value',1)
set(findobj('Tag','message_text'),'String','')
set(findobj('Tag','hold_back_edit'),'String','0')

y=handles.y;
axes(handles.plot_handle);
clear_axes(handles)
hold off;
% Need to produce a plot
h=plot(y);
set(h,'LineWidth',2,'Color',[0 0 .6])
axis tight;
title('Plot of original data')
AX=axis;
range = AX(4)-AX(3);
AX(3)=AX(3)-.1*range;
AX(4)=AX(4)+.1*range;
axis(AX);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function residuals_save_residuals_Callback(hObject, eventdata, handles)
% hObject    handle to residuals_save_residuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function residuals_save_std_residuals_Callback(hObject, eventdata, handles)
% hObject    handle to residuals_save_std_residuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_tiff_Callback(hObject, eventdata, handles)
% hObject    handle to export_tiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

this_handle = handles.plot_handle;
HandleNewFigure = figure;
copyobj(handles.plot_handle, HandleNewFigure);
HandleAxes = get(HandleNewFigure,'Children');
set(HandleAxes,'Units', 'Pixels');
HandleAxesPosition = get(HandleAxes, 'Position');
HandleNewFigurePosition = get(HandleNewFigure, 'Position');
set(HandleAxes, 'Position', [70 70 HandleAxesPosition(3)+50 HandleAxesPosition(4)+50]);
set(HandleNewFigure,'Position', [HandleNewFigurePosition(1) HandleNewFigurePosition(2) HandleAxesPosition(3)+110+50 HandleAxesPosition(4)+110+50]);

% --------------------------------------------------------------------
function export_png_Callback(hObject, eventdata, handles)
% hObject    handle to export_png (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_eps_Callback(hObject, eventdata, handles)
% hObject    handle to export_eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function residual_menu_Callback(hObject, eventdata, handles)
% hObject    handle to residual_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_menu_Callback(hObject, eventdata, handles)
% hObject    handle to export_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in EXOG_clear_pushbutton.
function EXOG_clear_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to EXOG_clear_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(findobj('Tag','EXOG_edit'),'String',' ');


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function plot_handle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot_handle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate plot_handle

function AR_lags_edit_Callback(hObject, eventdata, handles)
% hObject    handle to AR_lags_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AR_lags_edit as text
%        str2double(get(hObject,'String')) returns contents of AR_lags_edit as a double
[lags, str, OK] = validate_AR_MA_lags(get(hObject,'String'),handles.model_ARMAX_AR_order,handles.model_ARMAX_AR_order_string);
if OK && ~isempty(lags)
    handles.model_ARMAX_AR_order = lags;
    % Update handles structure
    guidata(hObject, handles);
    set(handles.message_text,'String',['AR lags set sucessfully to ' num2str(lags) '.'],'ForegroundColor',[0 0 0])
elseif OK && isempty(lags)
    set(handles.message_text,'String','No AR lags included.','ForegroundColor',[0 0 0])
else
    set(handles.message_text,'String','There was a problem setting the AR lags.  Please try again.','ForegroundColor',[.8 0 0])
end

% --- Executes during object creation, after setting all properties.
function AR_lags_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AR_lags_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MA_lags_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MA_lags_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MA_lags_edit as text
%        str2double(get(hObject,'String')) returns contents of MA_lags_edit as a double
[lags, str, OK] = validate_AR_MA_lags(get(hObject,'String'),handles.model_ARMAX_MA_order,handles.model_ARMAX_MA_order_string);
if OK && ~isempty(lags)
    handles.model_ARMAX_MA_order = lags;
    % Update handles structure
    guidata(hObject, handles);
    set(handles.message_text,'String',['MA lags set sucessfully to ' num2str(lags) '.'],'ForegroundColor',[0 0 0])
elseif OK && isempty(lags)
    set(handles.message_text,'String',['No MA lags included.'],'ForegroundColor',[0 0 0])
else
    set(handles.message_text,'String','There was a problem setting the MA lags.  Please try again.','ForegroundColor',[.8 0 0])
end

% --- Executes during object creation, after setting all properties.
function MA_lags_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MA_lags_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in last_model_pushbutton.
function last_model_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to last_model_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function copy_figure_copy_Callback(hObject, eventdata, handles)
% hObject    handle to copy_figure_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% FIXME  I don't work with two axes, need to copy them seperately!
f=figure;
set(f,'PaperSize',[11 8.5],'PaperPosition',[0.25 0.25 10.5 8],...
    'PaperOrientation','landscape','Visible','on');
% Copy the current axes to the new figure

copyobj(handles.plot_axes2,f);
if get(handles.radiobutton18,'Value') || get(handles.radiobutton19,'Value')
    hold on;
    copyobj(handles.plot_handle,f);
    hold off;
end
% Get the axes handle
ax = get(f,'Children');
% Make the figure standard
set(ax,'ActivePositionProperty','outerposition','Units','normalized',...
    'Position',[0.1300 0.1100 0.7750 0.8150]);




function hold_back_edit_Callback(hObject, eventdata, handles)
% hObject    handle to hold_back_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hold_back_edit as text
%        str2double(get(hObject,'String')) returns contents of hold_back_edit as a double


% --- Executes during object creation, after setting all properties.
function hold_back_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hold_back_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function initialization_Callback(hObject, eventdata, handles)
% hObject    handle to initialization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the initialization popup window and return the correct info

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos_size = get(handles.figure1,'Position');
% Call modaldlg with the argument 'Position'.
%user_response  = 'Yes';
user_response = ARMAX_close_dialog('Title','Confirm Close');
switch user_response
    case {'No'}
        %     take no action
    case 'Yes'
        delete(handles.figure1)
end


function ACF_lags_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ACF_lags_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ACF_lags_edit as text
%        str2double(get(hObject,'String')) returns contents of ACF_lags_edit as a double
ACF_lags=handles.ACF_lags;
y=handles.y;
val = get(hObject,'String');
[val, OK] = bounded_integer_validate(val,ACF_lags,1,length(y)-1);
if ~OK
    set(hObject,'String',num2str(val));
end

handles.ACF_lags = val;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ACF_lags_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ACF_lags_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function [val, OK] = bounded_integer_validate(val,orig_val,LB,UB)
val = str2double(val);
% Test cases
%   Nan: return old value and  OK=0
%   <0 return 0 and OK=0
%   + but not integer, floor and OK=0
%   +, integer, OK =1
OK = false;
if isnan(val)
    val = orig_val;
elseif val<LB
    val = LB;
elseif val>UB;
    val = UB;
elseif floor(val)~=val
    val=floor(val);
else
    OK = true;
end




% --------------------------------------------------------------------
function about_menu_Callback(hObject, eventdata, handles)
% hObject    handle to about_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value = ARMAX_about('Title','About');



function no_model_run_plot(ax,handles)
axes(ax);
clear_axes(handles)
h=plot([0 1],[0 1]);
axis tight
set(h,'Visible','off')
t=text(.1,.5,'You must estimate a model before residual plots can be produced');
set(t,'FontWeight','bold','FontSize',8,'FontName','Tahoma')
set(get(h,'Parent'),'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1])


% --------------------------------------------------------------------
function stderror_menu_Callback(hObject, eventdata, handles)
% hObject    handle to stderror_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function heterorobust_menu_Callback(hObject, eventdata, handles)
% hObject    handle to heterorobust_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.model_InferenceMethod = 1;
set(hObject,'Checked','on')
set(handles.homoerror_menu,'Checked','off')
guidata(hObject, handles);
set(handles.message_text,'String',['Inference will be made using a heteroskedasticity robust covariace estimator.'],'ForegroundColor',[0 0 0])





% --------------------------------------------------------------------
function homoerror_menu_Callback(hObject, eventdata, handles)
% hObject    handle to homoerror_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.model_InferenceMethod = 0;
set(hObject,'Checked','on')
set(handles.heterorobust_menu,'Checked','off')
guidata(hObject, handles);
set(handles.message_text,'String',['Inference will be made assuming homoskedastic errors.'],'ForegroundColor',[0 0 0])




function [lags, str, OK] = validate_AR_MA_lags(str,old_lags,old_str)
% This code will do the input validation for the strings used for the AR
% and MA lags

% Lags must be less than T-1
% Must be > 0

% first things to do is to find the spaces.  From here I can break up the
% string into units which have to be of one type, either a number (#) or a
% number-colon-number combo (#:#).  Any failure and the OK = 0 and the
% program should terminate.

if isempty(str) || all(ismember(str',' '))
    OK = true;
    lags = [];
    str = '';
    return
end

% The fisr check is to verify that all the elements are in {1,2,...,9,0,:}

valid=all(ismember(str',{'1','2','3','4','5','6','7','8','9','0',':',' '}));
lags=[];
OK = false;
if valid
    try
        eval(['lags=[' str '];']);
    catch
        OK = false;
    end
end
if ~isempty(lags)
    lags = unique(lags);
    OK=true;
else
    str = old_str;
    OK = false;
end






% --- Executes on selection change in model_selector_popup.
function model_selector_popup_Callback(hObject, eventdata, handles)
% hObject    handle to model_selector_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns model_selector_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_selector_popup

% 1. Update AR Lags, MA Lags and Constant
model = handles.models{get(hObject,'Value')};
set(handles.AR_lags_edit,'String',num2str(model.ARlags))
set(handles.MA_lags_edit,'String',num2str(model.MAlags))
% 2. Change residuals
handles.model_ARMAX_residuals=model.errors;
% 3. Graph raw residuals
y=handles.model_ARMAX_residuals;
axes(handles.plot_handle);
clear_axes(handles)
hold off;
% Need to produce a plot
h=plot(y);
set(h,'LineWidth',2,'Color',[0 0 .6])
axis tight;
title('Plot of estimated errors')
AX=axis;
range = AX(4)-AX(3);
AX(3)=AX(3)-.1*range;
AX(4)=AX(4)+.1*range;
axis(AX);
% 4. Change Radio button
set(handles.resid_raw_radio,'Value',1)
% 2. Update Message
set(handles.message_text,'String',['Results from model ' model.StringID ' loaded.'],'ForegroundColor',[0 0 .9])


% --- Executes during object creation, after setting all properties.
function model_selector_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_selector_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function constant_include_CreateFcn(hObject, eventdata, handles)
% hObject    handle to constant_include (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in constant_include.
function constant_include_Callback(hObject, eventdata, handles)
% hObject    handle to constant_include (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constant_include
handles.model_IncludeConstant=1;
guidata(hObject, handles);
set(handles.message_text,'String','Model will include a constant.','ForegroundColor',[0 0 0])



% --- Executes on button press in constant_exclude.
function constant_exclude_Callback(hObject, eventdata, handles)
% hObject    handle to constant_exclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of constant_exclude
handles.model_IncludeConstant=0;
guidata(hObject, handles);
set(handles.message_text,'String','Model will not include a constant.','ForegroundColor',[0 0 0])




% --- Executes on button press in orig_raw_radio.
function orig_raw_radio_Callback(hObject, eventdata, handles)
% hObject    handle to orig_raw_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orig_raw_radio
y=handles.y;
axes(handles.plot_handle);
clear_axes(handles)
hold off;
% Need to produce a plot
h=plot(y);
set(h,'LineWidth',2,'Color',[0 0 .6])
axis tight;
title('Plot of original data')
AX=axis;
range = AX(4)-AX(3);
AX(3)=AX(3)-.1*range;
AX(4)=AX(4)+.1*range;
axis(AX);




% --- Executes on button press in orig_acf_radio.
function orig_acf_radio_Callback(hObject, eventdata, handles)
% hObject    handle to orig_acf_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orig_acf_radio

y=handles.y;
ACF_lags=handles.ACF_lags;
robust=handles.model_InferenceMethod;
axes(handles.plot_handle);
clear_axes(handles);
hold off;
[ac,acstd]=sacf(y,ACF_lags,robust,0);

h = bar(ac);
set(h,'FaceColor',[.5 .5 1]);
axis tight;
ax = axis;
spread = .2*(ax(4)-ax(3));
ax(1) = 0;
ax(2) = ACF_lags+1;
ax(3) = ax(3) - spread;
ax(4) = ax(4) + spread;
axis(ax);
hold on;
h2 = plot((0:ACF_lags+1)',[[acstd(1);acstd;acstd(ACF_lags)] -[acstd(1);acstd;acstd(ACF_lags)]]);

if robust
    title('Sample Autocorrelations and Robust Standard Errors');
else
    title('Sample Autocorrelations and Non-robust Standard Errors');
end

set(h2(1),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
set(h2(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);


% --- Executes on button press in orig_pacf_radio.
function orig_pacf_radio_Callback(hObject, eventdata, handles)
% hObject    handle to orig_pacf_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of orig_pacf_radio
y=handles.y;
ACF_lags=handles.ACF_lags;
robust=handles.model_InferenceMethod;
axes(handles.plot_handle);
clear_axes(handles)
hold off;
[pac,pacstd]=spacf(y,ACF_lags,robust,0);

h = bar(pac);
set(h,'FaceColor',[.5 .5 1]);
axis tight;
ax = axis;
spread = .2*(ax(4)-ax(3));
ax(1) = 0;
ax(2) = ACF_lags+1;
ax(3) = ax(3) - spread;
ax(4) = ax(4) + spread;
axis(ax);
hold on;
h2 = plot((0:ACF_lags+1)',[2*[pacstd(1);pacstd;pacstd(ACF_lags)] -2*[pacstd(1);pacstd;pacstd(ACF_lags)]]);

if robust
    title('Sample Partial Autocorrelations and Robust Standard Errors');
else
    title('Sample Partial Autocorrelations and Non-robust Standard Errors');
end

set(h2(1),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
set(h2(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);



% --- Executes on button press in resid_raw_radio.
function resid_raw_radio_Callback(hObject, eventdata, handles)
% hObject    handle to resid_raw_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resid_raw_radio

if isempty(handles.model_ARMAX_residuals)
    no_model_run_plot(handles.plot_handle,handles);
else
    y=handles.model_ARMAX_residuals;
    axes(handles.plot_handle);
    clear_axes(handles)
    hold off;
    % Need to produce a plot
    h=plot(y);
    set(h,'LineWidth',2,'Color',[0 0 .6])
    axis tight;
    title('Plot of estimated errors')
    AX=axis;
    range = AX(4)-AX(3);
    AX(3)=AX(3)-.1*range;
    AX(4)=AX(4)+.1*range;
    axis(AX);
end



% --- Executes on button press in resid_fit_radio.
function resid_fit_radio_Callback(hObject, eventdata, handles)
% hObject    handle to resid_fit_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resid_fit_radio
if isempty(handles.model_ARMAX_residuals)
    no_model_run_plot(handles.plot_handle,handles);
else
    yhat=handles.y-handles.model_ARMAX_residuals;
    y=handles.y;
    axes(handles.plot_handle);
    clear_axes(handles)
    hold off;
    % Need to produce a plot
    h=plot([y yhat]);
    set(h(1),'LineWidth',2,'Color',[0 0 .6])
    set(h(2),'LineWidth',2,'Color',[.16 .38 .27])
    axis tight;
    title('Plot of estimated errors')
    legend('Original Data','Fit Data')
    AX=axis;
    range = AX(4)-AX(3);
    AX(3)=AX(3)-.1*range;
    AX(4)=AX(4)+.1*range;
    axis(AX);
end




% --- Executes on button press in resid_acf_radio.
function resid_acf_radio_Callback(hObject, eventdata, handles)
% hObject    handle to resid_acf_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resid_acf_radio

if isempty(handles.model_ARMAX_residuals)
    no_model_run_plot(handles.plot_handle,handles);
else
    y=handles.model_ARMAX_residuals;
    ACF_lags=handles.ACF_lags;
    robust=handles.model_InferenceMethod;
    axes(handles.plot_handle);
    clear_axes(handles);
    hold off;
    [ac,acstd]=sacf(y,ACF_lags,robust,0);

    h = bar(ac);
    set(h,'FaceColor',[.5 .5 1]);
    axis tight;
    ax = axis;
    spread = .2*(ax(4)-ax(3));
    ax(1) = 0;
    ax(2) = ACF_lags+1;
    ax(3) = ax(3) - spread;
    ax(4) = ax(4) + spread;
    axis(ax);
    hold on;
    h2 = plot((0:ACF_lags+1)',[[acstd(1);acstd;acstd(ACF_lags)] -[acstd(1);acstd;acstd(ACF_lags)]]);

    if robust
        title('Sample Autocorrelations and Robust Standard Errors');
    else
        title('Sample Autocorrelations and Non-robust Standard Errors');
    end

    set(h2(1),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
    set(h2(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
end



% --- Executes on button press in resid_pacf_radio.
function resid_pacf_radio_Callback(hObject, eventdata, handles)
% hObject    handle to resid_pacf_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resid_pacf_radio

if isempty(handles.model_ARMAX_residuals)
    no_model_run_plot(handles.plot_handle,handles);
else
    y=handles.model_ARMAX_residuals;
    ACF_lags=handles.ACF_lags;
    robust=handles.model_InferenceMethod;
    axes(handles.plot_handle);
    clear_axes(handles)

    hold off;
    [pac,pacstd]=spacf(y,ACF_lags,robust,0);

    h = bar(pac);
    set(h,'FaceColor',[.5 .5 1]);
    axis tight;
    ax = axis;
    spread = .2*(ax(4)-ax(3));
    ax(1) = 0;
    ax(2) = ACF_lags+1;
    ax(3) = ax(3) - spread;
    ax(4) = ax(4) + spread;
    axis(ax);
    hold on;
    h2 = plot((0:ACF_lags+1)',[2*[pacstd(1);pacstd;pacstd(ACF_lags)] -2*[pacstd(1);pacstd;pacstd(ACF_lags)]]);

    if robust
        title('Sample Partial Autocorrelations and Robust Standard Errors');
    else
        title('Sample Partial Autocorrelations and Non-robust Standard Errors');
    end

    set(h2(1),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);
    set(h2(2),'LineWidth',2,'LineStyle',':','Color',[0 0 0]);

end



% --- Executes on button press in radiobutton18.
function radiobutton18_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton18


lags = handles.ACF_lags;
y = handles.y;
if handles.model_InferenceMethod
    [stat,pval] = lmtest1(y,lags);
    titleString = 'LM Test Statistic Values (left) and P-values (right)';
else
    [stat,pval] = ljungbox(y,lags);
    titleString = 'Ljung-Box Q Statistics (left) and P-values (right)';
end

clear_axes(handles)
% 1. Set up the axes and their limits
barUL = max(stat)*1.1;
plotUL = 1.1;
plotLL = 0;
barLL = -max(stat)*.1;
axes(handles.plot_handle);
ax=gca;
axis([0 lags+1 barLL barUL])
set(ax,'XTickMode', 'auto','YTickMode', 'auto','Units','pixels','Position',[32 37 585 392]);
axes(handles.plot_axes2);
axis([0 lags+1 -.1 1.1])
ax=gca;
set(ax,'YTickMode','auto','Units','pixels','Position',[32 37 585 392],'YAxisLocation','right');

% 2. Draw teh bar graph and establish it's limits
axes(handles.plot_handle);
hold on;
b = bar(stat);
set(b,'EdgeColor',[1 1 1],'FaceColor',[.5 .5 1])
set(get(b(1),'BaseLine'),'LineWidth',1,'LineStyle','none')
hold on;
h = plot(max(stat)*pval);
set(h,'LineWidth',1,'LineStyle','none','Color',[0 0 0],'Marker','s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
legend('Statistic Value','P-value','Location','NorthWest')

% 3. Labels and titles
title(titleString)





% --- Executes on button press in radiobutton19.
function radiobutton19_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton19
if isempty(handles.model_ARMAX_residuals)
    no_model_run_plot(handles.plot_handle,handles);
else
    lags = handles.ACF_lags;
    y=handles.model_ARMAX_residuals;
    if handles.model_InferenceMethod
        [stat,pval] = lmtest1(y,lags);
        titleString = 'Error LM Test Statistic Values (left) and P-values (right)';
    else
        [stat,pval] = ljungbox(y,lags);
        titleString = 'Error Ljung-Box Q Statistics (left) and P-values (right)';
    end

    clear_axes(handles)
    % 1. Set up the axes and their limits
    barUL = max(stat)*1.1;
    plotUL = 1.1;
    plotLL = 0;
    barLL = -max(stat)*.1;
    axes(handles.plot_handle);
    ax=gca;
    axis([0 lags+1 barLL barUL])
    set(ax,'XTickMode', 'auto','YTickMode', 'auto','Units','pixels','Position',[32 37 585 392]);
    axes(handles.plot_axes2);
    axis([0 lags+1 -.1 1.1])
    ax=gca;
    set(ax,'YTickMode','auto','Units','pixels','Position',[32 37 585 392],'YAxisLocation','right');

    % 2. Draw teh bar graph and establish it's limits
    axes(handles.plot_handle);
    hold on;
    b = bar(stat);
    set(b,'EdgeColor',[1 1 1],'FaceColor',[.5 .5 1])
    set(get(b(1),'BaseLine'),'LineWidth',1,'LineStyle','none')
    hold on;
    h = plot(max(stat)*pval);
    set(h,'LineWidth',1,'LineStyle','none','Color',[0 0 0],'Marker','s','MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
    legend('Statistic Value','P-value','Location','NorthWest')

    % 3. Labels and titles
    title(titleString)
end


function clear_axes(handles)
axes(handles.plot_handle);
cla
set(handles.plot_handle,'XTick',[],'YTick',[],'Units','pixels','Position',[32 37 585 392])
title('')
axes(handles.plot_axes2);
cla
set(handles.plot_axes2,'XTick',[],'YTick',[],'Units','pixels','Position',[32 37 585 392])
title('')
