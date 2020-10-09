function [X,pro,DepthFlag]=BackDiffsCNop_inline(C_gas,C_bg,D1,D2,k1,k2,t1,t2,ProfileLength,PixelNo,PS,hand)
DepthFlag=0;
PixelNo_Data=PixelNo;
PixelNo=300;
dx=ProfileLength*1e-6/(PixelNo-1); %spatial step
X=0:dx:ProfileLength*1e-6; %domain
MP=str2double(get(hand.MirrorPlane,'String'))*1e-6;
[val_mirror MP] = min(abs(X-MP));
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

dt=5*round(dx^2/(2*min(D1,D2)),3,'significant'); %time step
% dt=round(dx^2/(2*max(D1,D2)),3,'significant');
Nt=round((t_tot)/dt); % Number of time steps
%
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
    if mean(~isfinite(Crank1))>0 || PS==1 %% If Crank1 is infinite or it isn't, but the PS is activated
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
        AA1=A1_newI*A1;
        AG1=A2_newI*G1;
        for t_idx=1:ceil(t_i/dt) %%%% double  check this
            C=AA1*C+AG1;
        end
        C1=C;
    else
        C1=Crank1;
        t_idx=ceil(t_i/dt);
    end
    C=C1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    AA2=A2_newI*A2;
    AG2=A2_newI*G2;
    for t_idx=1:Nt-t_idx   %% double check
        C=AA2*C+AG2;
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
                AA1=A1_newI*A1;
                AG1=A1_newI*G1;
                for t_idx=1:ceil(t_i/dt) %%%% double  check this
                    C1=AA1*C+AG1;
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
            AA2=A2_newI*A2;
            AG2=A2_newI*G2;
            for t_idx=1:Nt-t_idx   %% double check
                C=AA2*C+AG2;
            end
            pro=(pro+C')./2;%pro=C(:,end); pro=pro';
        end
    end %nagetive
end

% if pro(1)>pro(2) || pro(1)<0 || ~isfinite(pro(1))
%     pro(1:4)=interp1(X(5:10),pro(5:10),X(1:4),'pchip');
% end
PixelNo=PixelNo_Data;
dx=ProfileLength*1e-6/(PixelNo-1); %spatial step
X_sim=X;
X=0:dx:ProfileLength*1e-6; %domain
pro=interp1(X_sim,pro,X,'spline');
