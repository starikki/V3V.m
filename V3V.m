%V3V DATA takes the gridded data exported by V3V, calculates Velocity
%Magnitude, Vorticity Components, and Q-Criterion and exports into a
%tecplot format.

% Version:  EXP4 P005->EXP3 P001
% 


clear
% --- Outputs from this function are returned to the command line.
% function V3VData()
VolNum=1;
Blanking = 0.8; %blank vort by vortx value
Wing_blanking = 0; % switch on and off wing blanking
Trim = 0; % Matrix Box trim
Xbskip = 2;
Ybskip = 2;
Zbskip = 1;

%%%%%%%%%% sAR = 3 %%%%%%%%%%%%%%%
%FreestreamVel=[0.1065; 0.1051; 0.1051; 0.1045; 0.1045;];
%FLEXIBLE (not necessarily best!!!)
% Xshift=[26.8; 21.3; 27.8; 30.8; 21.8; -18; -19; -11];
% Yshift=[33.5; 32; 32; 30; 40; 705; 708; 704];
% Zshift=[725; 790; 846; 897; 937; 79; 192; 280];

%% Base Info
%RIGID
Xshift=[36; 18.3]; %x-axis + to shift fwd
Yshift=[48; 30];  % Z-axis, + to shift up
Zshift=[640; 790]; % Y-axis + to shift right

FreestreamVel=[0.295; 0.295];
chord = 0.06;
span = 0.24;

XVolSize=[-120 60; -200 200];        %-96 30
YVolSize=[-72 24; -100 100];
ZVolSize=[-12 42; 0 360];

Rotate=[0 0 0 0 0 1 1 1]; %Rotate volumes


Xshift=Xshift(VolNum);
Yshift=Yshift(VolNum);
Zshift=Zshift(VolNum);
FreestreamVel=FreestreamVel(VolNum);

XMin = XVolSize(VolNum, 1);
XMax = XVolSize(VolNum, 2);
YMin = YVolSize(VolNum, 1);
YMax = YVolSize(VolNum, 2);
ZMin = ZVolSize(VolNum, 1);
ZMax = ZVolSize(VolNum, 2);

Ucut = 1.5; %velocity Cut off proportion
Vcut  = 0.6; % reduction coefficient of celocity judgement on VW
Wcut = 0.6;

HomeCD=cd;

%% Load Files
[filename,pathname] = uigetfile('*.EV3D','select Ensembled', 'MultiSelect','on');



    %fid(f)=fopen(Fullname);
Nfile = size(filename,2);   
if  ischar(filename) ==1
    Nfile = 1;
end
for f  = 1:Nfile
tic;    
    
if Nfile == 1
        Fullname = [pathname,filename];
        Current_name = ['Grp_',filename];
else
        Fullname = [pathname,filename{1,f}];
        Current_name = ['Grp_',filename{1,f}];
end

disp ('-------------------------------');
disp (['loaded:' Current_name ]);
toc
    Matrix1 = csvread(Fullname,1,0); %load matrix from file
    
    X = Matrix1(:,1); 
    Y = Matrix1(:,2); 
    Z = Matrix1(:,3); 
    U = Matrix1(:,4);
    V = Matrix1(:,5); 
    W = Matrix1(:,6); 
    CHC = Matrix1(:,7);

    m = zeros(1,2);
%% %%%%%%%%%%%%%%%%%%%Change for Columns into 2D Matrix%%%%%%%%%%%%%%%%%%%%%%%%%
    m = find(X==X(1,1),2, 'first');     %diff of m = Number of rows
    n=length(X)/(m(2)-m(1));            %Number of Columns
   
    %RESHAPE FROM LIST INTO NEW MATRIX
    Xnew = (reshape(X,m(2)-m(1),n));
    Ynew = (reshape(Y,m(2)-m(1),n));
    Znew = (reshape(Z,m(2)-m(1),n));
    Unew = (reshape(U,m(2)-m(1),n));
    Vnew = (reshape(V,m(2)-m(1),n));
    Wnew = (reshape(W,m(2)-m(1),n));
    CHCnew = (reshape(CHC,m(2)-m(1),n));
    
%     Xnew(CHC==-1)=NaN;
%     Ynew(CHC==-1)=NaN;
%     Znew(CHC==-1)=NaN;    
%     Unew(CHC==-1)=NaN;
%     Vnew(CHC==-1)=NaN;
%     Wnew(CHC==-1)=NaN;
    
    Unew(CHC<=500)=NaN;
    Vnew(CHC<=500)=NaN;
    Wnew(CHC<=500)=NaN;

% %   DELETE ROWS AND COLUMS WITH FULL NAN
%     Unew = Unew(:,any(~isnan(Unew)));  % for columns
%     Unew = Unew(any(~isnan(Unew),2),:);   %for rows     
%     Vnew = Vnew(:,any(~isnan(Vnew)));  % for columns
%     Vnew = Vnew(any(~isnan(Vnew),2),:);   %for rows      
%     Wnew = Wnew(:,any(~isnan(Wnew)));  % for columns
%     Wnew = Wnew(any(~isnan(Wnew),2),:);   %for rows   

 

Xnew = Xnew + Xshift;
Xnew = Xnew * -1;   % Flip
Ynew = Ynew + Yshift;
Ynew = Ynew * -1;   % Flip
Znew = Znew + Zshift;
Znew = Znew * -1;
Unew = Unew*1;
Vnew = Vnew*1;
Wnew = Wnew*1;



if Trim == 1

XnewO = Xnew; % Xnew matrix for deletion index

YnewO(Xnew<XMin | Xnew>XMax)=NaN;
XnewO(Ynew<YMin | Ynew>YMax)=NaN;
XnewO(Znew<ZMin | Znew>ZMax)=NaN;


Xnew = Xnew(:,any(~isnan(XnewO)));  % for columns
Xnew = Xnew(any(~isnan(XnewO),2),:);   %for rows     
Ynew = Ynew(:,any(~isnan(XnewO)));  % for columns
Ynew = Ynew(any(~isnan(XnewO),2),:);   %for rows  
Znew = Znew(:,any(~isnan(XnewO)));  % for columns
Znew = Znew(any(~isnan(XnewO),2),:);   %for rows  
Unew = Unew(:,any(~isnan(XnewO)));  % for columns
Unew = Unew(any(~isnan(XnewO),2),:);   %for rows  
Vnew = Vnew(:,any(~isnan(XnewO)));  % for columns
Vnew = Vnew(any(~isnan(XnewO),2),:);   %for rows  
Wnew = Wnew(:,any(~isnan(XnewO)));  % for columns
Wnew = Wnew(any(~isnan(XnewO),2),:);   %for rows  

else
Xnew = Xnew(:,any(~isnan(Xnew)));  % for columns
Xnew = Xnew(any(~isnan(Xnew),2),:);   %for rows     
Ynew = Ynew(:,any(~isnan(Xnew)));  % for columns
Ynew = Ynew(any(~isnan(Xnew),2),:);   %for rows  
Znew = Znew(:,any(~isnan(Xnew)));  % for columns
Znew = Znew(any(~isnan(Xnew),2),:);   %for rows  
Unew = Unew(:,any(~isnan(Xnew)));  % for columns
Unew = Unew(any(~isnan(Xnew),2),:);   %for rows  
Vnew = Vnew(:,any(~isnan(Xnew)));  % for columns
Vnew = Vnew(any(~isnan(Xnew),2),:);   %for rows  
Wnew = Wnew(:,any(~isnan(Xnew)));  % for columns
Wnew = Wnew(any(~isnan(Xnew),2),:);   %for rows   
end


% Unew(Unew<0 | Unew>Ucut*FreestreamVel )=NaN;
% Vnew(Unew<0 | Unew>Ucut*FreestreamVel )=NaN;
% Wnew(Unew<0 | Unew>Ucut*FreestreamVel )=NaN;
% 
% Unew(Vnew<-Vcut*FreestreamVel | Vnew>Vcut*FreestreamVel )=NaN;
% Vnew(Vnew<-Vcut*FreestreamVel | Vnew>Vcut*FreestreamVel )=NaN;
% Wnew(Vnew<-Vcut*FreestreamVel | Vnew>Vcut*FreestreamVel )=NaN;
% 
% Unew(Wnew<-Wcut*FreestreamVel | Wnew>Wcut*FreestreamVel )=NaN;
% Vnew(Wnew<-Wcut*FreestreamVel | Wnew>Wcut*FreestreamVel )=NaN;
% Wnew(Wnew<-Wcut*FreestreamVel | Wnew>Wcut*FreestreamVel )=NaN;
%% Convert to 3D matrix
    %%%%%%%%%%%%%%%%%%%Change for Columns into 2D Matrix%%%%%%%%%%%%%%%%%%%%%%%%%
    m = find(Ynew(1,:)==Ynew(1,1),2, 'first'); %num of columns per z plane
    n=length(Ynew(1,:))/(m(2)-m(1)); % num of Z planes
    
    XI = zeros(size(Ynew,1),m(2)-1,n-1);
    YI = zeros(size(Ynew,1),m(2)-1,n-1);
    ZI = zeros(size(Ynew,1),m(2)-1,n-1);
    UI = zeros(size(Ynew,1),m(2)-1,n-1);
    VI = zeros(size(Ynew,1),m(2)-1,n-1);
    WI = zeros(size(Ynew,1),m(2)-1,n-1);
    CHCI= zeros(size(Ynew,1),m(2)-1,n-1);
    
for Zplane = 1:n-1
% [ (Zplane-1)*(m(2)-1)+1,Zplane*(m(2)-1)]
% end
    XI(:,:,Zplane) = Xnew( :,[(Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);
    YI(:,:,Zplane) = Ynew( :,[(Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);
    ZI(:,:,Zplane) = Znew( :,[(Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);
    UI(:,:,Zplane) = Unew( :,[(Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);
    VI(:,:,Zplane) = Vnew( :,[(Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);
    WI(:,:,Zplane) = Wnew( :,[(Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);  
    CHCI(:,:,Zplane) = CHCnew(:,[ (Zplane-1)*(m(2)-1)+1:Zplane*(m(2)-1)]);   
end
% 

% % % % % Correction of Matrix and Axis

% % in = size(XI,1);
% % jn = size(XI,3);
% % kn = size(XI,2);
% % 
% % for i = 1:in
% %     for j = 1:jn
% %         for k = 1:kn
% % XI_new(i,j,k) = XI(i,k,j);
% % YI_new(i,j,k) = ZI(i,k,j);
% % ZI_new(i,j,k) = YI(i,k,j);
% % CHCI_new(i,j,k) = CHCI(i,k,j);
% % UI_new(i,j,k) = UI(i,k,j);
% % VI_new(i,j,k) = WI(i,k,j);
% % WI_new(i,j,k) = VI(i,k,j);
% %         end
% %     end
% % end

XI_new = permute(XI,[1,3,2]);
YI_new = permute(ZI,[1,3,2]);
ZI_new = permute(YI,[1,3,2]);
UI_new = permute(UI,[1,3,2]);
VI_new = permute(WI,[1,3,2]);
WI_new = permute(VI,[1,3,2]);
CHCI_new = permute(CHCI,[1,3,2]);

XI = XI_new;
YI = YI_new;
ZI = ZI_new;
UI = UI_new;
VI = VI_new;
WI = WI_new;
CHCI = CHCI_new;




%% Data Blanking under wing
if Wing_blanking == 1
    j_min = find(YI(1,:,1)<=0,1, 'first') ;  
    j_max = find(YI(1,:,1)<=0,1, 'last') ;

    i_min = find(XI(:,1,1)<=0,1, 'first') +1;  
    i_max = find(XI(:,1,1)>=-59.77,1, 'last') -1;

    for j = (j_min+Ybskip):j_max
        for i = i_min+Xbskip:i_max-Xbskip
            XX = 60+XI(i,1,1);
            z_lim = 5*0.12*60*(0.2969*sqrt(XX/60)+(-0.126)*(XX/60)+(-0.3516)*(XX/60)^2+0.2843*(XX/60)^3+(-0.1015)*(XX/60)^4)+((60-XX)*sin(5/180*pi())); %Blanking under wing profile
            z_max = find(ZI(1,1,:)<z_lim,1,'last');

            for k = 1:(z_max-Zbskip)  %(z_max-n)  
                UI(i,j,k) = NaN;
                VI(i,j,k) = NaN;
                WI(i,j,k) = NaN;
            end
        end
    end

end

%% generate new columns    
Xcol = reshape(XI, size(XI,1)*size(XI,2)*size(XI,3),1);
Ycol = reshape(YI, size(YI,1)*size(YI,2)*size(YI,3),1);
Zcol = reshape(ZI, size(ZI,1)*size(ZI,2)*size(ZI,3),1);

% Xcol = Xcol*-1;

Ycol = Ycol*-1;

Zcol = Zcol*-1;


Ucol = reshape(UI, size(XI,1)*size(XI,2)*size(XI,3),1);
Vcol = reshape(VI, size(YI,1)*size(YI,2)*size(YI,3),1);
Wcol = reshape(WI, size(ZI,1)*size(ZI,2)*size(ZI,3),1);

CHCcol= reshape(CHCI, size(CHCI,1)*size(CHCI,2)*size(CHCI,3),1);

FirstX = min(Xcol);
% LastX = max(Xcol);

FirstY = min(Ycol);
% LastY = max(Ycol);

FirstZ = min(Zcol);
% LastZ = max(Zcol);
 %Spacial resolution           
    dx=abs(Xnew(2,1)-Xnew(1,1));
    dy=abs(Ynew(1,2)-Ynew(1,1));    
       dz = abs(ZI(1,1,2)-ZI(1,1,1));
%     dz=abs(Znew(1, m(2))-Znew(1, m(2)-1));      
%% %%%%%%%%%%%%%%%%%%%%   CALCULATE VORTICITY  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VorticityXI=zeros(size(XI,1),size(XI,2),size(XI,3));
VorticityYI=zeros(size(XI,1),size(XI,2),size(XI,3));
VorticityZI=zeros(size(XI,1),size(XI,2),size(XI,3));
VorticityMagI=zeros(size(XI,1),size(XI,2),size(XI,3));
% w = waitbar(0,'vorticity');
for i = 2:size(XI,1)-1
    for j = 2:size(XI,2)-1
        for k = 2:size(XI,3)-1
            VorticityXI(i,j,k) = 1000*((WI(i,j+1,k)-WI(i,j-1, k))/(YI(i,j+1, k)-YI(i, j-1, k)) - (VI(i,j, k+1) - VI(i, j, k-1))/(ZI(i, j, k+1)-ZI(i, j, k-1)));
            VorticityYI(i,j,k) = 1000*((UI(i, j,k+1)-UI(i,j, k-1))/(ZI(i,j, k+1)-ZI(i, j, k-1)) - (WI(i+1,j, k) - WI(i-1, j, k))/(XI(i+1, j, k)-XI(i-1, j, k)));
            VorticityZI(i,j,k) = 1000*((VI(i+1,j, k)-VI(i-1, j, k))/(XI(i+1,j, k)-XI(i-1, j, k)) - (UI(i,j+1, k) - UI(i, j-1, k))/(YI(i, j+1, k)-YI(i, j-1, k)));
            VorticityMagI(i,j,k) = sqrt(VorticityXI(i,j,k)^2 + VorticityYI(i,j,k)^2 + VorticityZI(i,j,k)^2);
        end
    end
    progress  = floor(i/(size(XI,1)-1)*100);
%   waitbar(f/Nfile,w,sprintf('%d',progress))
end
%  close(w)
VorticityXcol=reshape(VorticityXI, size(VorticityXI,1)*size(VorticityXI,2)*size(VorticityXI,3),1);
VorticityYcol=reshape(VorticityYI, size(VorticityYI,1)*size(VorticityYI,2)*size(VorticityYI,3),1);
VorticityZcol=reshape(VorticityZI, size(VorticityZI,1)*size(VorticityZI,2)*size(VorticityZI,3),1);
VorticityMagcol=reshape(VorticityMagI, size(VorticityMagI,1)*size(VorticityMagI,2)*size(VorticityMagI,3),1);


%%%%%%%%%%%%%%%%%%%%%   CALCULATE Q-CRITERION   %%%%%%%%%%%%%%%%%%%%%%%%%%%
[UX, UY, UZ] = gradient(UI, dx/1000, dy/1000, dz/1000);
[VX, VY, VZ] = gradient(VI, dx/1000, dy/1000, dz/1000);
[WX, WY, WZ] = gradient(WI, dx/1000, dy/1000, dz/1000);

QI=zeros(size(XI,1),size(XI,2),size(XI,3));
VelMagI=zeros(size(XI,1),size(XI,2),size(XI,3));
% for i = 1:size(XI,1)
%     for j = 1:size(XI,2)
%         for k = 1:size(XI,3)
%             DIV = [UX(i,j,k) UY(i,j,k) UZ(i,j,k); VX(i,j,k) VY(i,j,k) VZ(i,j,k); WX(i,j,k) WY(i,j,k) WZ(i,j,k)];
%             QI(i,j,k)=0.5*((norm(0.5*(DIV-transpose(DIV))))^2-(norm(0.5*(DIV+transpose(DIV))))^2);
%             
%             VelMagI(i,j,k) = sqrt(UI(i,j,k)^2+VI(i,j,k)^2+WI(i,j,k)^2);
%         end
%     end
% end
Qcol=reshape(QI, size(QI,1)*size(QI,2)*size(QI,3),1);
VelMagcol=reshape(VelMagI, size(VelMagI,1)*size(VelMagI,2)*size(VelMagI,3),1);


% % %%%%%%%%%%%%%%%%%%%    Calculate GAMMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Gsize = 1;
% % GammaX = zeros(size(XI,1),size(XI,2),size(XI,3));
% % 
% % for i = 1:size(XI,1)
% % 
% %         for jc = (Gsize+1):1:((size(YI,2))- Gsize)            %define global ID of processing box
% %             for kc = (Gsize+1):1:((size(ZI,3)) - Gsize)
% %         %matrix for gamma box
% %         UM =    zeros(3,2*Gsize+1,2*Gsize+1);
% %         PM =    zeros(3, 2*Gsize+1,2*Gsize+1);
% %         Udiff = zeros(3,2*Gsize+1,2*Gsize+1);
% %         gamma = zeros(3,2*Gsize+1,2*Gsize+1);
% %         gamma_sum = zeros(2*Gsize+1,2*Gsize+1);
% % 
% %         %fill value for each location of the box
% %         UM (1, :,:) = VI(i,[jc-Gsize:jc+Gsize],[kc-Gsize:kc+Gsize]);
% %         UM (2, :,:) = WI(i,[jc-Gsize:jc+Gsize],[kc-Gsize:kc+Gsize]);
% %         PM (1, :,:) = YI(i,[jc-Gsize:jc+Gsize],[kc-Gsize:kc+Gsize])- YI(i,jc,kc);
% %         PM (2, :,:) = ZI(i,[jc-Gsize:jc+Gsize],[kc-Gsize:kc+Gsize])- ZI(i,jc,kc);
% % 
% %         %empty core
% %         UM(1,Gsize+1,Gsize+1) = 0;
% %         UM(2,Gsize+1,Gsize+1) = 0;
% %         PM(1,Gsize+1,Gsize+1) = 0;
% %         PM(2,Gsize+1,Gsize+1) = 0;
% % 
% %         % Cross / dot for each loc
% %                     for j = 1: 1: 2*Gsize+1
% %                         for k = 1: 1: 2*Gsize+1
% % %                                 gamma(:,i,j) = cross(PM(:,i,j),Udiff(:,i,j))./ (   (norm(PM(:,i,j)) .*  norm(Udiff(:,i,j))));
% %                         gamma(:,j,k) = cross(PM(:,j,k),UM(:,j,k))./ (   (norm(PM(:,j,k)) .*  norm(UM(:,j,k))));
% %                         end 
% %                     end
% %                     gamma_sum (:,:) = gamma(3,:,:);
% %                     GammaX(i,jc,kc) = (1/Gsize)*norm(nansum(gamma_sum(:)));
% %             end
% %         end
% % end
% % 
% % GammaXcol=reshape(GammaX, size(GammaX,1)*size(GammaX,2)*size(GammaX,3),1);

%% %%%%%%%%%%%%%%%   NON-DIMENSIONALIZE AND SHIFT   %%%%%%%%%%%%%%%%%%%%%%%%
VorticityXcol=VorticityXcol*chord/FreestreamVel;
VorticityYcol=VorticityYcol*chord/FreestreamVel;
VorticityZcol=VorticityZcol*chord/FreestreamVel;
VorticityMagcol=VorticityMagcol*chord/FreestreamVel;
VelMagcol=VelMagcol/FreestreamVel;
Qcol=Qcol*chord/FreestreamVel;

%% TRIM TOP AND BOTTOM OF DATA - if required
%  GammaXcol(Xcol<XMin | Xcol>XMax)=NaN;
% VorticityXcol(Xcol<XMin | Xcol>XMax)=NaN;
% VorticityYcol(Xcol<XMin | Xcol>XMax)=NaN;
% VorticityZcol(Xcol<XMin | Xcol>XMax)=NaN;
% VorticityMagcol(Xcol<XMin | Xcol>XMax)=NaN;
% VelMagcol(Xcol<XMin | Xcol>XMax)=NaN;
% Qcol(Xcol<XMin | Xcol>XMax)=NaN;
% 
% %  GammaXcol(Ycol<YMin | Ycol>YMax)=NaN;
% VorticityXcol(Ycol<YMin | Ycol>YMax)=NaN;
% VorticityYcol(Ycol<YMin | Ycol>YMax)=NaN;
% VorticityZcol(Ycol<YMin | Ycol>YMax)=NaN;
% VorticityMagcol(Ycol<YMin | Ycol>YMax)=NaN;
% VelMagcol(Ycol<YMin | Ycol>YMax)=NaN;
% Qcol(Ycol<YMin | Ycol>YMax)=NaN;
% 
% %  GammaXcol(Zcol<ZMin | Zcol>ZMax)=NaN;
% VorticityXcol(Zcol<ZMin | Zcol>ZMax)=NaN;
% VorticityYcol(Zcol<ZMin | Zcol>ZMax)=NaN;
% VorticityZcol(Zcol<ZMin | Zcol>ZMax)=NaN;
% VorticityMagcol(Zcol<ZMin | Zcol>ZMax)=NaN;
% VelMagcol(Zcol<ZMin | Zcol>ZMax)=NaN;
% Qcol(Zcol<ZMin | Zcol>ZMax)=NaN;
%% Conditional Value blaking!
 VorticityMagcol(VorticityXcol>-Blanking & VorticityXcol<Blanking)=NaN; %Trim off 
% VorticityXcol(VorticityXcol>-0.5 & VorticityXcol<0.5)=NaN; %Clean the mess
%REMOVE NANs
% Xcol(isnan(Xcol)==1)=0;
% Ycol(isnan(Ycol)==1)=0;
% Zcol(isnan(Zcol)==1)=0;
Ucol(isnan(Ucol)==1)=0;
Vcol(isnan(Vcol)==1)=0;
Wcol(isnan(Wcol)==1)=0;
VelMagcol(isnan(VelMagcol)==1)=0;
VorticityXcol(isnan(VorticityXcol)==1)=0;
VorticityYcol(isnan(VorticityYcol)==1)=0;
VorticityZcol(isnan(VorticityZcol)==1)=0;
VorticityMagcol(isnan(VorticityMagcol)==1)=0;
Qcol(isnan(Qcol)==1)=0;
%  GammaXcol(isnan(GammaXcol)==1)=0;
 
 Xcol = Xcol/chord/1000;
 Ycol = Ycol/chord/1000;
 Zcol = Zcol/chord/1000;
%  

VorticityXcol(isinf(VorticityXcol)==1)=0;
VorticityYcol(isinf(VorticityYcol)==1)=0;
VorticityZcol(isinf(VorticityZcol)==1)=0;
VorticityMagcol(isinf(VorticityMagcol)==1)=0;
Qcol(isinf(Qcol)==1)=0;
%  GammaX(isinf(GammaX)==1)=0;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%   WRITE TO FILE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(HomeCD)
HeaderText = ['Title="' Current_name(1:10) '" '...
    'VARIABLES="X","Y","Z","U/U<math>%</math>","V/U<math>%</math>","W/U<math>%</math>","CHC","Velocity Magnitude","<greek>w</greek><sub>X</sub>c /U<math>%</math>","<greek>w</greek><sub>Y</sub>c /U<math>%</math>","<greek>w</greek><sub>Z</sub>c /U<math>%</math>","Vorticity Magnitude","QCriterion",DATASETAUXDATA Position_Unit="mm",DATASETAUXDATA Velocity_Unit="m/s",DATASETAUXDATA Vorticity_Unit="1/s",DATASETAUXDATA DataType="V",DATASETAUXDATA Dimension="3",DATASETAUXDATA HasVelocity="Y",DATASETAUXDATA ExtraDataNumber="6",DATASETAUXDATA GridSpacingX="'...
    num2str(dx) '",DATASETAUXDATA GridSpacingY="' num2str(dy) '",DATASETAUXDATA GridSpacingZ="' num2str(dz) '",DATASETAUXDATA FirstNodeX="' num2str(FirstX) '",DATASETAUXDATA FirstNodeY="' num2str(FirstY) '",DATASETAUXDATA FirstNodeZ="' num2str(FirstZ)...
    '",ZONE T="' Current_name(1:15) '",I=' num2str(size(XI, 1)) ',J=' num2str(size(XI, 2)) ',K=' num2str(size(XI, 3)) ',F=POINT,'];

fid=fopen(strcat(Fullname, '.plt'), 'w');
fprintf(fid, '%s \n', HeaderText);
fclose(fid);
dlmwrite(strcat(Fullname, '.plt'), [Xcol Ycol Zcol Ucol Vcol Wcol CHCcol VelMagcol VorticityXcol VorticityYcol VorticityZcol VorticityMagcol Qcol], '-append')
toc
end
% for VV = min(VorticityXI(:)):1:max(VorticityXI(:))
%     isosurface(XI,YI,ZI,VorticityXI,VV)
%     refreshdata
%     drawnow
%     VV
% end
Sys_time = datestr(now);
disp('done')
disp( Sys_time )