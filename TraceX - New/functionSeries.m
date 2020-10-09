function FitButton_Callback(hObject, eventdata, hand)
hand.CurrentPlot='Fit';
if ~isfield(hand,'ProfileData_Or') % is there hand.ProfileData_Or (i.e. is there any user data) - if not, it will prompt the user to load profile data
    [hand]=LoadProfileData_Callback(hObject, eventdata, hand);
    if hand.FileName==0 %if they just press esc and don't load any data, stops the whole fitting process
        return
    end
end
% Change appearance of Fit_Button

%which things do the user want to fit
FitCheck(1)=get(hand.D1FitCheck, 'Value'); 
FitCheck(2)=get(hand.D2FitCheck, 'Value');
FitCheck(3)=get(hand.k1FitCheck, 'Value');
FitCheck(4)=get(hand.k2FitCheck, 'Value');
 
C_bg=str2double(get(hand.C_bg, 'String'));
C_gas=str2double(get(hand.C_gas, 'String'));
D1=str2double(get(hand.D1, 'String'));
D2=str2double(get(hand.D2, 'String'));
k1=str2double(get(hand.k1, 'String'));
k2=str2double(get(hand.k2, 'String'));
t1=str2double(get(hand.t1, 'String'));
t2=str2double(get(hand.t2, 'String'));
ProfileLength=str2double(get(hand.ProfileLength, 'String'));
PixelNo=str2double(get(hand.PixelNo, 'String'));
PS=get(hand.PlaneSheet, 'Value');
[hand] =PlotButton_Callback(hObject, eventdata, hand); % Press plot button

Fit_Start=str2double(get(hand.Fit_Start, 'String'));
[val_start hand.Fit_Start_idx] = min(abs(hand.Xdata-Fit_Start*1e-6)); % Selects data point
Fit_Start_idx=hand.Fit_Start_idx;

Fit_End=str2double(get(hand.Fit_End, 'String'));
[val_end hand.Fit_End_idx] = min(abs(hand.Xdata-Fit_End*1e-6));
Fit_End_idx=hand.Fit_End_idx;
MP=hand.MP_idx;

if ~isfield(hand,'ProfileData')
    [hand]=Data_Callback(hObject, eventdata, hand);
end


DataStep=ProfileLength*1e-6/(PixelNo-1);
hand.Xdata_Or=0:DataStep:ProfileLength*1e-6;

SurfPos=str2double(get(hand.SurfPos, 'String'));
[val_start hand.SurfPos_idx] = min(abs(hand.Xdata_Or-SurfPos*1e-6));
hand.Xdata=hand.Xdata_Or(1:end-hand.SurfPos_idx+1);
hand.ProfileData=hand.ProfileData_Or(hand.SurfPos_idx:end);

if t2==0
    if PS==0 %Plane
        [D1,k1]=AutoFit_Crank(hand);
        if get(hand.Check_LeClaire,'Value')==1
            [D1,k1,A_gb,Z_gb]=AutoFit_Crank_LC(hand);
        end
    else
        [D1,k1]=AutoFit_PlaneSheet(hand);
    end
else

    if FitCheck(2)==0 && FitCheck(4)==0
        [D1,k1]= AutoFit_BackCrank(C_gas,C_bg,D1,k1,t1,t2,hand.ProfileData,...
            hand.Xdata,Fit_Start_idx,Fit_End_idx);
        D2=D1;
        k2=k1;
        %     elseif FitCheck(2)==0 && FitCheck(4)==1
    else
        set(hand.FitButton,'ForegroundColor',[0.5,0.5,0.5]);
        %         set(hand.FitButton,'Enable','off');
        set(hand.FitButton,'String','Fitting...');pause(0.001);
        % Start by using the analytical solution to get a good/quick ball park...
        [D1,k1]= AutoFit_BackCrank(C_gas,C_bg,D1,k1,t1,t2,hand.ProfileData,...
            hand.Xdata,Fit_Start_idx,Fit_End_idx);
        D2=D1;
        k2=k1;
        
        tic;
        r2=0;
        tim=0;
        i=1;
        r2check(1)=str2double(get(hand.Rsq,'String'));
        %
        while r2<0.999 && tim<20;
            
            [D1,k1,D2,k2]= AutoFit_BackDiffs(C_gas,C_bg,D1,k1,D2,k2,t1,t2,hand.ProfileData,...
                hand.Xdata,Fit_Start_idx,Fit_End_idx,ProfileLength,PixelNo,hand,FitCheck);
            [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,hand); %ImageLength
            [r2 rmse] = rsquare(hand.ProfileData(Fit_Start_idx:Fit_End_idx),...
                pro(Fit_Start_idx:Fit_End_idx));
            
            tim=toc+tim;
            i=i+1;
            r2check(i)=r2;
            if r2check(i)==r2check(i-1);
                tim=20;
            end
            
        end
        disp(['Total time for fit ',num2str(round(toc)),' seconds.'])
        %         set(hand.FitButton,'Enable','on');
    end
end
D1=roundsf(D1,3,'round');
D2=roundsf(D2,3,'round');
k1=roundsf(k1,3,'round');
k2=roundsf(k2,3,'round');
%%% Set all the newly fitted values
set(hand.D1,'String',num2str(D1));
set(hand.D2,'String',num2str(D2));
set(hand.k1,'String',num2str(k1));
set(hand.k2,'String',num2str(k2));

FitPlot(hObject, eventdata, hand);

if isfield(hand,'O18pathname')
    hand.SaveName=[hand.PathName,hand.FileName(1:end-4)];
else
    hand.SaveName=[hand.PathName,hand.FileName(1:end-4)];
end
set(hand.FitButton,'ForegroundColor',[0,0,0]);
set(hand.FitButton,'String','Fit');
hand.CurrentPlot='Fit';
guidata(hObject, hand);


function [D1,k1,D2,k2]= AutoFit_BackDiffs(C_gas,C_bg,D1,k1,D2,k2,t1,t2,ProfileData,...
    Xdata,Fit_Start_idx,Fit_End_idx,ProfileLength,PixelNo,hand,FitCheck)
ran(1)=Fit_Start_idx;
ran(2)=Fit_End_idx;

b(1)=C_gas;
b(2)=C_bg;
b(3)=t1;
b(4)=t2;
b(5)=ProfileLength;
b(6)=PixelNo;
PS=get(hand.PlaneSheet,'Value');
MP_idx=hand.MP_idx;
Mask=zeros(size(Xdata));
% Mask(Fit_Start_idx:Fit_End_idx)=1;
fun = @(p) sum((...
    ProfileData(ran(1):ran(2))-BackDiffs_AutoCaller(p,b,ran,FitCheck,PS,MP_idx,hand)...
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

% [p] = fmincon(fun,pguess,[],[],[],[],0.001*pguess,1000*pguess)
% [p,fminres] = fseminf(fun,pguess,0);
% [p,resnorm] = lsqcurvefit(@fun,pguess,BackDiffs_AutoCaller(p,b,ran,FitCheck),ProfileData(ran(1):ran(2)));
%

if sum(FitCheck)==4
    D1=p(1);D2=p(2);k1=p(3);k2=p(4);
elseif FitCheck(2)==0
    D1=p(1);D2=p(1);k1=p(2);k2=p(3);
elseif FitCheck(4)==0
    D1=p(1);D2=p(2);k1=p(3);k2=p(3);
end

function [pro]=BackDiffs_AutoCaller(p,b,ran,FitCheck,PS,MP,hand);
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
    
    [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,hand);
    pro=pro(ran(1):ran(2));
end