function varargout = regression_scrambler(varargin)
% REGRESSION_SCRAMBLER MATLAB code for regression_scrambler.fig
%      REGRESSION_SCRAMBLER, by itself, creates a new REGRESSION_SCRAMBLER or raises the existing
%      singleton*.
%
%      H = REGRESSION_SCRAMBLER returns the handle to a new REGRESSION_SCRAMBLER or the handle to
%      the existing singleton*.
%
%      REGRESSION_SCRAMBLER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGRESSION_SCRAMBLER.M with the given input arguments.
%
%      REGRESSION_SCRAMBLER('Property','Value',...) creates a new REGRESSION_SCRAMBLER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before regression_scrambler_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to regression_scrambler_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help regression_scrambler

% Last Modified by GUIDE v2.5 03-Aug-2016 16:17:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @regression_scrambler_OpeningFcn, ...
    'gui_OutputFcn',  @regression_scrambler_OutputFcn, ...
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


% --- Executes just before regression_scrambler is made visible.
function regression_scrambler_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to regression_scrambler (see VARARGIN)

% Choose default command line output for regression_scrambler
handles.output = hObject;

%Initilizing required variables :
set(handles.parameter_selected, 'String', 'None Selected');
set(handles.reg_method, 'String', 'None Selected');

handles.X = 0;
handles.Y = 0;
handles.Y_complex = 0;
handles.wavelength = 0;
handles.var_wavelength = 0;
handles.range = 0;

handles.data.x_train = 0;
handles.data.x_train_var_sel = 0;
handles.data.y_train = 0;
handles.data.x_test = 0;
handles.data.x_test_var_sel = 0;
handles.data.y_test = 0;

handles.data.reference = 0;
handles.continue_flag = 0;

addpath Regress_Files;
addpath supporting_files;
addpath variable_main_files;
addpath Variable_Sel_supporting;
addpath results_regSelection;
addpath results_varSelection;

%For Variable selection and Regression Mtd
handles.var_sel_processing = 'None';
handles.reg_sel_processing = 'None';
handles.func_prop_sel = 'None';


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes regression_scrambler wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = regression_scrambler_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.continue_flag == 0)
    [X, Y_complex, wavelength] = load_data();
    handles.X = X;
    handles.Y_complex = Y_complex;
    handles.wavelength = wavelength;
    handles.Y = property_selector(Y_complex,handles.func_prop_sel);
    data_size= strcat(num2str(size(X,1)),'','X','', num2str(size(X,2)));
    set(handles.nir_size,'string',data_size);
end
set(handles.preprocess,'visible','on');
set(handles.load_data,'visible','off');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in preprocess.
function preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of preprocess

set(handles.preprocess,'visible','off');
X = handles.X ;
Y = handles.Y ;
if (get(handles.dependency_data,'value') == 1)
interim1 = cov(handles.Y);
interim2 = cov(handles.Y) * tapering(size(handles.Y,1),1);

Y = inv(sqrt(interim2)) * handles.Y;
X = inv(sqrt(interim2)) * handles.X;

end
    

[handles.data.x_train,handles.data.y_train...
    ,handles.data.x_test,handles.data.y_test, handles.range] = preprocess_nir_data( X , Y , handles.wavelength, handles.continue_flag, handles.func_prop_sel);

set(handles.reg_sel_mtd,'visible','on');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in var_select.
function var_select_Callback(hObject, eventdata, handles)
% hObject    handle to var_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of var_select

if(get(handles.var_select,'value') == 1)
    set(handles.var_sel_mtd,'visible','on');
else
    set(handles.var_sel_mtd,'visible','off');
end

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on selection change in var_sel_mtd.
function var_sel_mtd_Callback(hObject, eventdata, handles)
% hObject    handle to var_sel_mtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns var_sel_mtd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from var_sel_mtd

full_list2 = get(handles.var_sel_mtd,'string');
sel_val=get(handles.var_sel_mtd,'value');
sel_items=full_list2(sel_val);
set(handles.parameter_selected,'string',sel_items);
handles.var_sel_processing = sel_items ;

if (isequal(sel_items{1},'None')== 0)
    set(handles.reg_sel_mtd,'visible','off');
    [ wavelength_final ] = variable_selection( handles.data.x_train,handles.data.y_train,handles.data.x_test,handles.data.y_test,handles.var_sel_processing,handles.range );
    set(handles.reg_sel_mtd,'visible','on');
    handles.var_wavelength = wavelength_final;
    handles.data.x_train_var_sel = handles.data.x_train(:,wavelength_final);
    handles.data.x_test_var_sel = handles.data.x_test(:,wavelength_final);
end

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function var_sel_mtd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to var_sel_mtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in reg_sel_mtd.
function reg_sel_mtd_Callback(hObject, eventdata, handles)
% hObject    handle to reg_sel_mtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns reg_sel_mtd contents as cell array
%        contents{get(hObject,'Value')} returns selected item from reg_sel_mtd

full_list = get(handles.reg_sel_mtd,'string');
sel_val=get(handles.reg_sel_mtd,'value');
sel_items=full_list(sel_val);
set(handles.reg_method,'string',sel_items);
set(handles.regress,'visible','on');
handles.reg_sel_processing = sel_items ;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function reg_sel_mtd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reg_sel_mtd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in regress.
function regress_Callback(hObject, eventdata, handles)
% hObject    handle to regress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.regress,'visible','off');

if(get(handles.var_select,'value') == 1)
    variable_selected_method = handles.var_sel_processing;
else
    variable_selected_method = 'None';
end

regression_selected_method = handles.reg_sel_processing;

%x_train = handles.data.x_train;
y_train = handles.data.y_train;
%x_test = handles.data.x_test;
y_test = handles.data.y_test;

if (isequal(handles.var_wavelength,0) == 1)
    x_train = handles.data.x_train;
    x_test = handles.data.x_test;
else
    x_train = handles.data.x_train_var_sel;
    x_test = handles.data.x_test_var_sel;
end
range = handles.range ;

switch regression_selected_method{1}
    case 'PLS'
        [ r2train,r2test,rmsecvC,rmsepC  ] = pls_regression( x_train, y_train, x_test , y_test, handles.func_prop_sel );
    case 'PCA'
        [ r2train,r2test,rmsecvC,rmsepC  ] = pca_regression( x_train, y_train, x_test , y_test , range,handles.func_prop_sel);
    case 'ICA'
        [ r2train,r2test,rmsecvC,rmsepC  ] = ica_regression( x_train, y_train, x_test , y_test ,range);
    case 'Ridge'
        [ r2train,r2test,rmsecvC,rmsepC  ] = ridge_regression( x_train, y_train, x_test , y_test , handles.func_prop_sel);
    case 'LASSO'
        [ r2train,r2test,rmsecvC,rmsepC  ] = lasso_regression( x_train, y_train, x_test , y_test, handles.func_prop_sel );
    case 'LS-SVM'
        [ r2train,r2test,rmsecvC,rmsepC  ] = lssvm_regression( x_train, y_train, x_test , y_test);
    case 'RegressAll'
        regressAll( x_train, y_train, x_test , y_test , range, handles.func_prop_sel);
    otherwise
        warning('no Regression selected')
end


% --- Executes on key press with focus on reg_sel_mtd and none of its controls.
function reg_sel_mtd_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to reg_sel_mtd (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
clear all;
clc;
regression_scrambler


% --- Executes on selection change in functional_property.
function functional_property_Callback(hObject, eventdata, handles)
% hObject    handle to functional_property (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns functional_property contents as cell array
%        contents{get(hObject,'Value')} returns selected item from functional_property

full_list3 = get(handles.functional_property,'string');
sel_val=get(handles.functional_property,'value');
sel_items=full_list3(sel_val);
set(handles.load_data,'visible','on');
handles.func_prop_sel = sel_items;
set(handles.func_selected,'string',sel_items);
set(handles.functional_property,'visible','off');
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function functional_property_CreateFcn(hObject, eventdata, handles)
% hObject    handle to functional_property (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in next_round.
function next_round_Callback(hObject, eventdata, handles)
% hObject    handle to next_round (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.continue_flag = 1 ;
set(handles.functional_property,'visible','on');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in dependency_data.
function dependency_data_Callback(hObject, eventdata, handles)
% hObject    handle to dependency_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dependency_data


% --- Executes on button press in coll_varSel.
function coll_varSel_Callback(hObject, eventdata, handles)
% hObject    handle to coll_varSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
collate_var ;


% --- Executes on button press in coll_regSel.
function coll_regSel_Callback(hObject, eventdata, handles)
% hObject    handle to coll_regSel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
collate_reg ;


% --- Executes on button press in clean_var.
function clean_var_Callback(hObject, eventdata, handles)
% hObject    handle to clean_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete('D:\NIR Gui Project\results_varSelection\*.mat');


% --- Executes on button press in clean_reg.
function clean_reg_Callback(hObject, eventdata, handles)
% hObject    handle to clean_reg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete('D:\NIR Gui Project\results_regSelection\*.mat');
