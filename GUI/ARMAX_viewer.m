function varargout = ARMAX_viewer(varargin)
% ARMAX_viewer M-file for ARMAX_viewer.fig
%      ARMAX_viewer, by itself, creates a new ARMAX_viewer or raises the existing
%      singleton*.
%
%      H = ARMAX_viewer returns the handle to a new ARMAX_viewer or the handle to
%      the existing singleton*.
%
%      ARMAX_viewer('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARMAX_viewer.M with the given input arguments.
%
%      ARMAX_viewer('Property','Value',...) creates a new ARMAX_viewer or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ARMAX_viewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ARMAX_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% TODO : Add parameter display for hyperparapmeter for distributions
% FIXME: Fix the behavior for very large models

% Edit the above text to modify the response to help ARMAX_viewer

% Last Modified by GUIDE v2.5 09-Sep-2006 14:45:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ARMAX_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @ARMAX_viewer_OutputFcn, ...
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

% --- Executes just before ARMAX_viewer is made visible.
function ARMAX_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ARMAX_viewer (see VARARGIN)

% Choose default command line output for ARMAX_viewer
handles.output  = hObject;
handles.Results = varargin{1};
if length(varargin)==2
    handles.DisplayNumber = varargin{2};
else
    handles.DisplayNumber = length(handles.Results);
end

% Update handles structure
guidata(hObject, handles);
popup_cell=[];
for i=1:length(handles.Results);
    popup_cell{i}=handles.Results{i}.StringID;
end
set(handles.model_popupmenu,'String',popup_cell)
set(handles.model_popupmenu,'Value',length(handles.Results))
set(gcf,'Position',[10 10 1000 706])
%
display_model(handles,hObject);



% UIWAIT makes ARMAX_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ARMAX_viewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over axes background.
function model_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to model_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function generate_text(handles,hObject);
% Get rid of results not needed
number = handles.DisplayNumber;
results=handles.Results{number};
% Compute the needed strings
[ParameterString,VariableString]= generate_strings(results);
% Layout the model
[model_layout_info]=model_layout(ParameterString,VariableString,handles);
% Layout the model statistics
[model_stats_layout_info]=model_stats_layout(results,handles);
% Resize the Parameter Estimates Window
[estimate_window]=resize_parameter_estimate_uipanel(handles);
% Layout the column titles and return their position
[title_positions]=title_layout(handles);
% Compute the text elements for the main display and their positions
[text_positions]=compute_text_positions(handles,ParameterString,title_positions);
% Update the results with a flag that the text layout is done and with all
% of the relevant information, also set the handles
handles.Results{number}.ParameterString=ParameterString;
handles.Results{number}.VariableString=VariableString;
handles.Results{number}.model_layout_info=model_layout_info;
handles.Results{number}.model_stats_layout_info=model_stats_layout_info;
handles.Results{number}.estimate_window=estimate_window;
handles.Results{number}.title_position=title_positions;
handles.Results{number}.text_position=text_positions;
handles.Results{number}.TextFinished=true;
guidata(hObject, handles);



function [ParameterString,VariableString]= generate_strings(results);
% Not finished
% TO DO: ARARCH, EGARCH and AGARCH
index=1;

ParameterString{index}='y_t =';
c{index}='';
index=index+1;

if results.constant
    ParameterString{index}='\phi_0';
    VariableString{index}='';
    index=index+1;
end

if ~isempty(results.ARlags)
    if max(results.ARlags)>0
        for i=1:length(results.ARlags)
            ParameterString{index}=['\phi_{' num2str(results.ARlags(i)) '}'];
            VariableString{index}=['y_{t-' num2str(results.ARlags(i)) '}'];
            index=index+1;
        end
    end
end


if ~isempty(results.MAlags)
    if max(results.MAlags)>0
        for i=1:length(results.MAlags)
            ParameterString{index}=['\theta_{' num2str(results.MAlags(i)) '}'];
            VariableString{index}=['\epsilon_{t-' num2str(results.MAlags(i)) '}'];
            index=index+1;
        end
    end
end

ParameterString{index}='';
VariableString{index}='\epsilon_t';




function [layout_info]=model_layout(ParameterString,VariableString,handles)

axes(handles.model_axes);
%set(gcf,'Position',[ -3          35        1024         664],'Units','pixels')
%set(ax,'Units','pixels','Position',[134.1200   500.0400  793.6000  120])
line=0.75;
count=1;

result = handles.Results{handles.DisplayNumber};
ModelString = result.StringID;
% Add on the model specific string
final_str{count}=ModelString;
count=count+1;

%FIX ME, need to worry about no constant!
if result.constant
    final_str{count}=[ParameterString{1} ParameterString{2}];
else
    final_str{count}=[ParameterString{1} ParameterString{2} VariableString{2}];
end


xpos=.05;
for i=3:length(ParameterString)
    t=text(xpos,line,['$' final_str{count} '+' ParameterString{i} VariableString{i} '$' ],'Interpreter','latex','FontSize',12);
    ex=get(t,'Extent');
    delete(t);
    if ex(1)+ex(3)<1
        % No problems, add and try next
        final_str{count}=[final_str{count} '+' ParameterString{i} VariableString{i}];
    else
        %Problem, no not add and start a new line
        xpos=.1;
        count=count+1;
        final_str{count}=['+' ParameterString{i} VariableString{i}];
    end
end

t=[];
line=0.75;
for i=1:count;
    if i<3
        xpos=.02;
    else
        xpos=.07;
    end
    if i>1
        t{i}=text(xpos,line,['$' final_str{i} '$'] ,'Interpreter','latex','FontSize',12);
    else
        t{i}=text(xpos,line,['$' final_str{i} '$'] ,'FontSize',12);
    end
    line=line-0.25;
end

% Now that the text is layed out, the axes and uipanel need to be resized
% accordingly and the text strings need to be vertically repositioned
% axes height is 33* no lines, y position is always 9
% uipanel height is 33*no lines+24, y position is 650-uiheight
no_lines=length(final_str);
uiHeight = no_lines*33+24;
uiY = 670 - uiHeight;
uiP=get(handles.model_uipanel,'Position');
uiP(2)=uiY;
uiP(4)=uiHeight;
set(handles.model_uipanel,'Position',uiP);

axHeight = 33*no_lines;
axY = 9;
axP=get(handles.model_axes,'Position');
axP(2)=axY;
axP(4)=axHeight;
set(handles.model_axes,'Position',axP);

% Finally uniformly distribute the t{i} between at
% [1:no_lines]/(no_lines+1)
positions = (no_lines:-1:1)/(no_lines+1);
for i=1:no_lines;
    tP=get(t{i},'Position');
    tP(2)=positions(i);
    set(t{i},'Position',tP);
end

% Need to return everything needed to make relayout easy:
%   uipanle position
%   axes position
%   text strings
%   text string position
layout_info.model_uipanel_position = get(handles.model_uipanel,'Position');
layout_info.model_axes_position = get(handles.model_axes,'Position');
layout_info.model_text_strings = final_str;
for i=1:length(final_str);
    layout_info.model_text_strings_position{i}=get(t{i},'Position');
end
% Clear the axis
cla






function [layout_info]=model_stats_layout(results,handles)
% 1. Fix the position of and size of the UI and axes
% 2. Layout the important stats and their labels

axes(handles.stats_axes);
set(handles.stats_axes,'Position',[15     9   920   60])
model_uipanel_position=get(handles.model_uipanel ,'Position');
set(handles.model_stats_uipanel,'Position',[25  model_uipanel_position(2)-84    950   84])

titles={'Log-likelihood','No. of Params',' Akaike IC','Schwartz/Bayesian IC','$\hat{\sigma}^2$'};
positions=linspace(.02,.9,length(titles));
total_space=0;
for i=1:length(titles)
    if i==5
        t(i)=text(positions(i),.66,titles{i},'FontName','Tahoma','FontSize',10,'Interpreter','latex');
    else
        t(i)=text(positions(i),.66,titles{i},'FontName','Tahoma','FontSize',10);
    end
    ex(i,:)=get(t(i),'Extent');
    if i>1
        total_space=total_space+ex(i,1)-(ex(i-1,3)+ex(i-1,1));
    end
end
total_space=total_space+.98-(ex(5,3)+ex(5,1));
% Compute the total amount of space between the objects and uniformly
% respace them
spaces=[4 4 4 4]*total_space/22;
for i=2:length(titles)
    ex=get(t(i-1),'Extent');
    set(t(i),'Position',[ex(1)+ex(3)+spaces(i-1) .66 0]);
end

% Set up the text strings.  I need to limit the width
data{1}=sprintf('%0.7g',results.ll);
data{2}=num2str(results.K);
data{3}=sprintf('%0.7g',results.AIC);
data{4}=sprintf('%0.7g',results.BIC);
data{5}=sprintf('%0.7g',results.seregression);


for i=1:5
    d(i)=text(positions(i),.33,data{i},'FontName','Tahoma','FontSize',10);
    % Now move it so that it is right aligned
    p1 = get(d(i),'Extent');
    p2 = get(t(i),'Extent');
    set(d(i),'Position',[p2(1)+p2(3)-p1(3) .33 0]);
end
AX=axis;
hold on;
for i=1:5
    % Now move it so that it is right aligned
    p2 = get(t(i),'Extent');
    axis(AX);
    x_line_data{i}=[p2(1) p2(1)+p2(3)]';
    y_line_data{i}=[p2(2) p2(2)];
    plot([p2(1) p2(1)+p2(3)],[p2(2) p2(2)],'k');
end
hold off;

% Need to return all of teh strings along with all of their positions
layout_info.Titles=titles;
layout_info.Data=data;
layout_info.TitlePosition=[];
layout_info.DataPosition=[];

for i=1:5;
    layout_info.TitleHandles{i}=t(i);
    layout_info.TitlePosition{i}=get(t(i),'Position');
    layout_info.DataHandle{i}=d(i);
    layout_info.DataPosition{i}=get(d(i),'Position');
end
layout_info.model_stats_uipanel_position = get(handles.model_stats_uipanel,'Position');
layout_info.stats_axes_position = get(handles.stats_axes,'Position');
layout_info.x_line_data=x_line_data;
layout_info.y_line_data=y_line_data;
% Clear the axis
cla


% --- Executes on selection change in model_popupmenu.
function model_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to model_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns model_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from model_popupmenu
handles.DisplayNumber = get(hObject,'Value');
guidata(hObject, handles);
display_model(handles,hObject);

% --- Executes during object creation, after setting all properties.
function model_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function page_number_edit_Callback(hObject, eventdata, handles)
% hObject    handle to page_number_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of page_number_edit as text
%        str2double(get(hObject,'String')) returns contents of page_number_edit as a double
val = get(hObject,'String');
result = handles.Results{handles.DisplayNumber};
LastPage = handles.LastPage;
NumberOfPages = size(result.text_position.x_pos,3);
[CurrentPage, OK] = bounded_integer_validate(val,LastPage,1,NumberOfPages);
if ~OK
    set(hObject,'String',num2str(CurrentPage));
end
draw_parameter_table(result,CurrentPage,handles)
handles.LastPage=CurrentPage;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function page_number_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to page_number_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function num_pages_text_Callback(hObject, eventdata, handles)
% hObject    handle to num_pages_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_pages_text as text
%        str2double(get(hObject,'String')) returns contents of num_pages_text as a double


% --- Executes during object creation, after setting all properties.
function num_pages_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_pages_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in increase_page_pushbutton.
function increase_page_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to increase_page_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.Results{handles.DisplayNumber};
NumberOfPages = size(result.text_position.x_pos,3);
CurrentPage = str2double(get(handles.page_number_edit,'String'));
if CurrentPage < NumberOfPages
    CurrentPage = CurrentPage +1;
    set(handles.page_number_edit,'String',num2str(CurrentPage));

    draw_parameter_table(result,CurrentPage,handles)
    handles.LastPage=CurrentPage;
    % Update handles structure
    guidata(hObject, handles);

end

% --- Executes on button press in decrease_page_pushbutton.
function decrease_page_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to decrease_page_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.Results{handles.DisplayNumber};
CurrentPage = str2double(get(handles.page_number_edit,'String'));
if CurrentPage >1
    CurrentPage = CurrentPage - 1;
    set(handles.page_number_edit,'String',num2str(CurrentPage));

    draw_parameter_table(result,CurrentPage,handles)
    handles.LastPage=CurrentPage;
    % Update handles structure
    guidata(hObject, handles);
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)




function [window_stats]=resize_parameter_estimate_uipanel(handles)
model_uipanel_pos=get(handles.model_stats_uipanel,'Position');
estimates_uipanel_pos=get(handles.estimates_uipanel ,'Position');
estimates_uipanel_pos(4)=estimates_uipanel_pos(4)+model_uipanel_pos(2)-sum(estimates_uipanel_pos([2 4]));
set(handles.estimates_uipanel ,'Position',estimates_uipanel_pos);

estimates_axes_pos=get(handles.estimates_axes,'Position');
estimates_axes_pos(2)=12;
set(handles.estimates_axes,'Position',estimates_axes_pos)

estimates_uipanel_pos=get(handles.estimates_uipanel ,'Position');
estimates_axes_pos=get(handles.estimates_axes,'Position');
estimates_axes_pos(4)=estimates_uipanel_pos(4)-24;
set(handles.estimates_axes,'Position',estimates_axes_pos);
estimates_axes_pos=get(handles.estimates_axes,'Position');
% Note, I want ot make sure that the estimates window is at least 330
% pixels high

% FIXME: Fix the behavior for very large models
% ScreenSize = get(0,'ScreenSize')
% set(handles.figure1,'Position',[1 ScreenSize-HeightofWindow-33 WidthofWindow HeightofWindow])
% Return the position of the uipanel, the axes and the figure!
window_stats.axes_pos=estimates_axes_pos;
window_stats.uipanel_pos=estimates_uipanel_pos;
window_stats.figure_pos=get(handles.figure1,'Position');





function [title_position]=title_layout(handles)
titles={'Parameter','Estimate','Std. Error','T-stat','P-val'};
axes(handles.estimates_axes);
t(1)=text(.02,.5,titles{1},'FontName','Tahoma','FontSize',10);
axes_pos=get(handles.estimates_axes,'Position');
set(t(1),'Units','pixels');
t_pos=get(t(1),'Position');
t_ex = get(t(1),'Extent');
y=floor(axes_pos(4)-t_ex(4)/2-10);
t_pos(2)=y;
set(t(1),'Position',t_pos)
position{1}=t_pos;
for i=2:length(titles);
    t(i)=text((i-1)/6+.02,.5,titles{i},'FontName','Tahoma','FontSize',10);
    set(t(i),'Units','pixels') ;
    t_pos=get(t(i),'Position');
    t_pos(2)=y;
    set(t(i),'Position',t_pos);
    position{i}=t_pos;
end
AX=axis;
hold on;
for i=1:length(titles)
    % Now move it so that it is right aligned
    set(t(i),'Units','normalized');
    p2 = get(t(i),'Extent');
    % Need to get the right points of each label to do right alignment
    right_points{i}=p2(1)+p2(3);
    axis(AX);
    x_plot{i}=[p2(1) p2(1)+p2(3)];
    y_plot{i}=[p2(2) p2(2)]';
    plot(x_plot{i},y_plot{i},'k');
end

hold off;
tt=text(.5,.5,'$\alpha \beta \gamma \lambda \omega \nu \mu$','Interpreter','latex','FontSize',14);
row_height=get(tt,'Extent');
delete(tt);
row_height=row_height(4);
title_position.Titles=titles;
title_position.TitlePosition=position;
title_position.x_line_data=x_plot;
title_position.y_line_data=y_plot;
title_position.RowHeight=row_height;
title_position.NumRows=floor((1-0.5*row_height)/row_height);
title_position.RightPoints=right_points;
% Clear the axis
cla


function [text_positions]=compute_text_positions(handles,str,title_position)
% First set up the texts

results = handles.Results{handles.DisplayNumber};
right_points = title_position.RightPoints;
actual_rows=title_position.NumRows-1;
for i=2:length(str)-1
    % first compute the page
    page = ceil((i-1)/actual_rows);

    if mod(i-1,actual_rows)==0
        table_position = actual_rows;
    else
        table_position = mod(i-1,actual_rows);
    end
    for j=1:5
        switch j
            case 1
                % Latex string
                window_text{table_position,j,page} = ['$' str{i} '$'];
            case 2
                % Estimate
                temp = results.parameters(i-1);
            case 3
                % Estimate
                temp = results.StdErr(i-1);
            case 4
                % Estimate
                temp = results.Tstat(i-1);
            case 5
                % Estimate
                temp = results.Pval(i-1);
        end
        if j~=1
            if temp<.000001
                window_text{table_position,j,page} = sprintf('%0.3g',temp);
            else
                window_text{table_position,j,page} = sprintf('%0.5f',temp);
            end
        end
    end
    if mod(i-1,actual_rows)==0
        actual_position = actual_rows+1;
    else
        actual_position = mod(i-1,actual_rows)+1;
    end
    row_height=title_position.RowHeight;
    for j=1:5
        y_pos{table_position,j,page} = 1-(actual_position)*row_height;
        if j==1
            temp_text =text(.5,y_pos{table_position,j,page},window_text{table_position,j,page},'Interpreter','latex','FontSize',12);
        else
            temp_text =text(.5,y_pos{table_position,j,page},window_text{table_position,j,page},'FontName','Tahoma','FontSize',10);
        end
        ex=get(temp_text,'Extent');
        delete(temp_text);
        x_pos{table_position,j,page}=right_points{j}-ex(3);
    end
end
% Set up the return data
text_positions.x_pos = x_pos;
text_positions.y_pos = y_pos;
text_positions.window_text = window_text;

function draw_data(result,handles,hObject)
% Resize the figure
set(handles.figure1,'Position',result.estimate_window.figure_pos);
% Resize the uipanels
set(handles.model_uipanel,'Position',result.model_layout_info.model_uipanel_position);
set(handles.model_stats_uipanel,'Position',result.model_stats_layout_info.model_stats_uipanel_position);
set(handles.estimates_uipanel,'Position',result.estimate_window.uipanel_pos)
% Resize the axes
set(handles.model_axes,'Position',result.model_layout_info.model_axes_position);
set(handles.stats_axes ,'Position',result.model_stats_layout_info.stats_axes_position);
set(handles.estimates_axes,'Position',result.estimate_window.axes_pos);
% Draw the model
axes(handles.model_axes)
% Make sure itis clear
cla
str = result.model_layout_info.model_text_strings;
pos = result.model_layout_info.model_text_strings_position;
N = length(str);
for i=1:N
    x=pos{i}(1);
    y=pos{i}(2);
    if i>1
        text(x,y,['$' str{i} '$'],'Unit','normalized','interpreter','latex','FontSize',12)
    else
        text(x,y,[str{i}],'Unit','normalized','FontSize',12)
    end
end
% Draw the model statistic titles, lines, and values
axes(handles.stats_axes)
% Make sure itis clear
cla
str = result.model_stats_layout_info.Titles;
pos = result.model_stats_layout_info.TitlePosition;
N = length(str);
for i=1:N
    x=pos{i}(1);
    y=pos{i}(2);
    if i~=5
        text(x,y,str{i},'Unit','normalized','FontSize',10,'FontName','Tahoma')
    else
        text(x,y,str{i},'Unit','normalized','FontSize',10,'FontName','Tahoma','Interpreter','latex')
    end
end
hold on;
for i=1:N
    % Now move it so that it is right aligned
    x=result.model_stats_layout_info.x_line_data{i};
    y=result.model_stats_layout_info.y_line_data{i};
    plot(x,y,'k');
end
hold off;
str = result.model_stats_layout_info.Data;
pos = result.model_stats_layout_info.DataPosition;
N = length(str);
for i=1:N
    x=pos{i}(1);
    y=pos{i}(2);
    text(x,y,str{i},'Unit','normalized','FontSize',10,'FontName','Tahoma')
end
% Call the function to draw the correct page, presumably one
draw_parameter_table(result,1,handles) % In the initial draw we will always draw page 1

handles.LastPage=1;
% Update handles structure
guidata(hObject, handles);


% Finally update the UI elements that have to do with the pages
NumberOfPages = size(result.text_position.x_pos,3);
if NumberOfPages==1
    % Easy case
    set(handles.page_number_edit,'BackgroundColor',[0.5 0.5 0.5],'String','1')
    set(handles.increase_page_pushbutton,'Enable','off')
    set(handles.decrease_page_pushbutton,'Enable','off')
    set(handles.number_of_pages_text,'String',num2str(NumberOfPages));
else
    set(handles.page_number_edit,'BackgroundColor',[1 1 1],'String','1')
    set(handles.increase_page_pushbutton,'Enable','on')
    set(handles.decrease_page_pushbutton,'Enable','on')
    set(handles.number_of_pages_text,'String',num2str(NumberOfPages));
end






function draw_parameter_table(result,page,handles)
axes(handles.estimates_axes)
% Make sure it is clear
cla
% Draw the column titles
str = result.title_position.Titles;
pos = result.title_position.TitlePosition;
N = length(str);
for i=1:N
    x=pos{i}(1);
    y=pos{i}(2);
    text(x,y,str{i},'Unit','pixel','FontSize',10,'FontName','Tahoma')
end
% Draw the lines
hold on;
for i=1:N
    % Now move it so that it is right aligned
    x=result.title_position.x_line_data{i};
    y=result.title_position.y_line_data{i};
    plot(x,y,'k');
end
hold off;
x = result.text_position.x_pos(:,:,page);
y = result.text_position.y_pos(:,:,page);
str = result.text_position.window_text(:,:,page);

N = size(x,1);

for i=1:N
    for j=1:5
        % Draw the data on the correct page
        if j==1
            text(x{i,j},y{i,j},str{i,j},'Interpreter','latex','FontSize',12);
        else
            text(x{i,j},y{i,j},str{i,j},'FontName','Tahoma','FontSize',10);
        end
    end
end



function display_model(handles,hObject)
% Check to see if the moodel has been processed.  If it has, go striaght to
% the draw model, else compute te strings then draw the model
results = handles.Results{handles.DisplayNumber};
if isfield(results,'TextFinished')
    draw_data(results,handles,hObject)
else
    % Generate the text labels
    generate_text(handles,hObject);
    % update handles
    handles = guidata(hObject);
    results = handles.Results{handles.DisplayNumber};
    % Draw the labels
    draw_data(results,handles,hObject)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to do integer vallidation for ACF Lags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
