%{
    The MATLAB code for the VF-MC model is available here. 
    You can run this script program now.

Function:
    Calculate the feedback coefficient, the low-frequency vibration level (LVL) factor, αi, and the low-frequency vibration divergence (LVD) factor, βi, for robotic milling.

Remarks:
    The symbols used in this code have been kept as similar as possible to those in the paper.
%}

%% Initialization
clear; clf; clc
%% Input parameters
% ！！！！ Tool parameters
global N DD r beta lambda rou
N = 2;                  % Number of teeth
DD = 20e-3;             % Diameter (m)
r = 0e-3;               % Fillet radius (m)
beta = 30*pi/180;       % Helical angle
rou = 0.0125e-3;        % Radius of tool runout(m)
lambda = 55/180*pi;     % Angle of tool runout

% ！！！！ Milling parameters
global up_or_down ae Kt Kr Ka w nn fe theta
up_or_down = -1;        % Milling mode. 1: up milling; -1: down milling
ae = 5e-3;              % Cutting width (m)
w = 1e-3;           	% Cutting depth (m)
nn = 7200;              % Spindle speed (rpm)
fe = 0.1e-3;            % Feed per tooth (m)
theta = 90*pi/180;      % Feed Angle

% ！！！！ Cutting force coefficients (N/m^2)
Kt = 632.6e6;           % Tangential
Kr = -205.8e6;          % Radial 
Ka = -43.5e6;           % Axial 

% ！！！！ Simulation parameters
global mm dz
mm = 100;               % Discrete number of spindle rotation period
dz = 0.1e-3;            % Axial discrete height of the tool
L = 4;                  % Maximum multiple of the harmonics of spindle rotation
J = 4;                  % Maximum number of feedbacks considered

% ！！！！ Robotic modal parameters
xb = [- 0.34, 0.90, 0.23     % Modal direction
- 0.47, 0.21, - 0.85
- 0.92, - 0.17, 0.33
- 0.32, 0.92, 0.20
- 0.76, 0.62, 0.14
- 0.91, 0.04, 0.40 
-0.53, 0.84, - 0.14 
- 0.72, 0.35, 0.58
- 0.79, - 0.12, 0.59
0.75, 0.35, 0.55 
- 0.64, - 0.55, 0.51
- 0.83, - 0.39, 0.38
0.40, 0.81, - 0.42];
ModeNum = size(xb,1);               % Number of modes
for iM = 1:ModeNum                  % Normalize xb
    xb(iM,:) = xb(iM,:)/norm(xb(iM,:));
end
xb = xb';
fre = [7.07	10.29	16.02	18.49	34.39	37.81	72.67	107.81	116.76	121.93	141.74	193.88	253.25]';       % Modal frequency
damp = [3.74	2.21	2.3	3.32	4	3.84	5.05	4.84	7.6	3.44	2.07	3.67	3.88]'*1e-2;            % Damping ratio
residue = [- 1.43e-6 + 4.38e-5j                                                                                     % Residues in the modal direction
(- 1.55e-6 - 3.45e-5j) * 0.9
(- 5.08e-6 - 3.62e-5j) * 0.65
(- 4.94e-6 + 2.15e-5j) * 1.25
- 2.45e-6 - 1.28e-5j
(+ 5.74e-7 - 1.78e-5j)
(+ 1.23e-5 + 4.86e-5j) *0.8
(- 6.21e-6 - 2.69e-5j) *0.5
(- 8.63e-6 - 5.84e-5j) *0.5
(+ 9.92e-7 - 5.05e-6j) *0.5
+ 5.22e-7 - 2.53e-6j
- 2.44e-6 - 7.04e-6j
+ 5.86e-7 + 3.72e-6j];
[~,~,K] = modal(fre,damp,imag(residue));    % Modal stiffness
wn = fre*2*pi;          % Modal angular frequency
xb = Rz(-theta)*xb;     % Rotate the modal direction in reverse to consider the change in feed direction

%% Some secondary parameters
T_tooth = 60/nn/N;    % Tooth passing period
f_spin = nn/60;       % Frequency of spindle rotation
w_spin = 2*pi*f_spin; % Angular frequency of spindle rotation
dt = N*T_tooth/mm;    % Sampling Interval
Fs = round(1/dt);     % Sampling frequency

%% Parameter setting of the initial vibration, ＠ri(t)
Ai = 1;     % Amplitude; no effect on feedback coefficients
fii = 0;    % Phase; no effect on feedback coefficients
i = 1;      % The mode to be calculated for its LFC
Amp_Gamma = [Ai/2,Ai/2]'; Pha_Gamma = [-fii,fii]';  % Amplitudes and phases of the impulses of Δγi(ω)

%% Calculate the coefficient matrix h(t) and its amplitude and phase in the frequency domain
z_Num = int32(w/dz);z_Num = double(z_Num);          % Number of tool discretizations in the axial direction

% ！！！！ Calculate the matrix h(t)
ht = ht_3D(z_Num);
ht(:,:,end) = [];
htt = repmat(ht, [1 1 10]);                         % Repeat for 10 periods to get more accurate spectra

% ！！！！ Calculate the frequency-domain amplitude and phase of h(t)
A_mnl = zeros( size(ht,1), size(ht,2), L+1 ); fi_mnl = A_mnl;   % Preset amplitude and phase
for m = 1:3
    for n = 1:3
        [Y,f,P1] = FourTran( reshape(htt(m,n,:),[],1) , Fs );
        AY = angle(Y);
        for l = 0:L 
            P1_lf = P1( f > l*f_spin-1/4*f_spin & f < l*f_spin+1/4*f_spin );
            f_lf  =  f( f > l*f_spin-1/4*f_spin & f < l*f_spin+1/4*f_spin );
            lf_index = find ( P1==max(P1_lf) );
            lf_index = lf_index(1);     % Prevents errors with multiple numbers in lf_index
            A_mnl(m,n,l+1)  = P1( lf_index );
            fi_mnl(m,n,l+1) = AY( lf_index );
        end
    end
end
% ！！！！ Assign values to Hmn(ω) (i.e., the spectrum of h(t))
Amp_Hmn = zeros( size(ht,1), size(ht,2), 2*L+1 ); Pha_Hmn = Amp_Hmn;  % Preset amplitude and phase
for m = 1:3 
    for n = 1:3
        for l = -L:L
            if l == 0
                Amp_Hmn(m,n,L+1) = A_mnl(m,n,1);
                Pha_Hmn(m,n,L+1) = fi_mnl(m,n,1);
            else
                Amp_Hmn(m,n,l+L+1) = A_mnl(m,n,abs(l)+1)/2;
                Pha_Hmn(m,n,l+L+1) = sign(l) * fi_mnl(m,n,abs(l)+1);
            end
            
        end
    end
end
%% Calculate 1st feedback
% ！！！！ Influence coefficients of the time delay term on the amplitude and phase
[DA1,DP1] = TimeDelay(-wn(i),T_tooth); 
[DA2,DP2] = TimeDelay(wn(i),T_tooth);

% ！！！！ Calculate 1st feedback
Amp_Gamma_iTOi1 = zeros( ModeNum, 4*L+2); Pha_Gamma_iTOi1 = Amp_Gamma_iTOi1; 
for i1 = 1:ModeNum
    for l = -L:L
        [H_i1_lwomwi,Fi_i1_lwomwi] = FRF(wn(i1),damp(i1),K(i1),abs(l*w_spin-wn(i)));
        Fi_i1_lwomwi = sign(l*w_spin-wn(i)) * Fi_i1_lwomwi;   % Phase change sign when frequency is less than zero
        [H_i1_lwopwi,Fi_i1_lwopwi] = FRF(wn(i1),damp(i1),K(i1),abs(l*w_spin+wn(i)));
        Fi_i1_lwopwi = sign(l*w_spin+wn(i)) * Fi_i1_lwopwi;
        
        X_lwomwi = 0; Y_lwomwi = 0;
        X_lwopwi = 0; Y_lwopwi = 0;
        for m = 1:3
            for n = 1:3
                Amp_Gamma_iTOi1mn_lwomwi = Amp_Gamma(1)*DA1 * Amp_Hmn(m,n,l+L+1) * H_i1_lwomwi;
                Amp_Gamma_iTOi1mn_lwopwi = Amp_Gamma(2)*DA2 * Amp_Hmn(m,n,l+L+1) * H_i1_lwopwi;
                Pha_Gamma_iTOi1mn_lwomwi = Pha_Gamma(1)+DP1 + Pha_Hmn(m,n,l+L+1) + Fi_i1_lwomwi;
                Pha_Gamma_iTOi1mn_lwopwi = Pha_Gamma(2)+DP2 + Pha_Hmn(m,n,l+L+1) + Fi_i1_lwopwi;
                
                X_lwomwi = X_lwomwi + xb(m,i1) * xb(n,i) * Amp_Gamma_iTOi1mn_lwomwi * cos(Pha_Gamma_iTOi1mn_lwomwi);
                X_lwopwi = X_lwopwi + xb(m,i1) * xb(n,i) * Amp_Gamma_iTOi1mn_lwopwi * cos(Pha_Gamma_iTOi1mn_lwopwi);
                Y_lwomwi = Y_lwomwi + xb(m,i1) * xb(n,i) * Amp_Gamma_iTOi1mn_lwomwi * sin(Pha_Gamma_iTOi1mn_lwomwi);
                Y_lwopwi = Y_lwopwi + xb(m,i1) * xb(n,i) * Amp_Gamma_iTOi1mn_lwopwi * sin(Pha_Gamma_iTOi1mn_lwopwi);
            end
        end
        Amp_Gamma_iTOi1(i1,2*l+2*L+1) = norm([X_lwomwi Y_lwomwi]); 
        Amp_Gamma_iTOi1(i1,2*l+2*L+2) = norm([X_lwopwi Y_lwopwi]);
        Pha_Gamma_iTOi1(i1,2*l+2*L+1) = angle(X_lwomwi+1i*Y_lwomwi);
        Pha_Gamma_iTOi1(i1,2*l+2*L+2) = angle(X_lwopwi+1i*Y_lwopwi);
    end
end
%% Calculate j-th feedback
Amp_Gamma_j_ij = zeros( ModeNum, 4*L+2, J); Pha_Gamma_j_ij = Amp_Gamma_j_ij;
Amp_Gamma_j_ij(:,:,1) = Amp_Gamma_iTOi1;
Pha_Gamma_j_ij(:,:,1) = Pha_Gamma_iTOi1;

for j = 2:J 
    for ij = 1:ModeNum
        X_j_ij = zeros(4*L+2,1); Y_j_ij = zeros(4*L+2,1); 
        for ijm1 = 1:ModeNum 
            Amp_Gamma_ijm1TOij_ijm1 = zeros(4*L+2,1);  Pha_Gamma_ijm1TOij_ijm1 = zeros(4*L+2,1); 
            
            for l = -L:L
                kmin = max(-L,l-L);
                kmax = min(L,l+L);
                
                [H_ij_lwomwi,Fi_ij_lwomwi] = FRF(wn(ij),damp(ij),K(ij),abs(l*w_spin-wn(i)));
                Fi_ij_lwomwi = sign(l*w_spin-wn(i)) * Fi_ij_lwomwi;   % Phase change sign when frequency is less than zero
                [H_ij_lwopwi,Fi_ij_lwopwi] = FRF(wn(ij),damp(ij),K(ij),abs(l*w_spin+wn(i))); 
                Fi_ij_lwopwi = sign(l*w_spin+wn(i)) * Fi_ij_lwopwi;
                
                X_lwomwi = 0; Y_lwomwi = 0; 
                X_lwopwi = 0; Y_lwopwi = 0;
                
                for m = 1:3
                    for n = 1:3
                        Xk_mwi = 0; Yk_mwi = 0;
                        Xk_pwi = 0; Yk_pwi = 0;
                        for k = kmin:kmax
                            [DA1,DP1] = TimeDelay( k*w_spin-wn(i),T_tooth );
                            [DA2,DP2] = TimeDelay( k*w_spin+wn(i),T_tooth );
                            Amp_Gammajm1_ijm1_kwomwi = Amp_Gamma_j_ij(ijm1,2*k+2*L+1,j-1) * DA1; 
                            Amp_Gammajm1_ijm1_kwopwi = Amp_Gamma_j_ij(ijm1,2*k+2*L+2,j-1) * DA2;
                            Pha_Gammajm1_ijm1_kwomwi = Pha_Gamma_j_ij(ijm1,2*k+2*L+1,j-1) + DP1;
                            Pha_Gammajm1_ijm1_kwopwi = Pha_Gamma_j_ij(ijm1,2*k+2*L+2,j-1) + DP2;
                            
                            Ak_mwi = Amp_Gammajm1_ijm1_kwomwi * Amp_Hmn(m,n,l-k+L+1) * H_ij_lwomwi;
                            Ak_pwi = Amp_Gammajm1_ijm1_kwopwi * Amp_Hmn(m,n,l-k+L+1) * H_ij_lwopwi;
                            Pha_mwi = Pha_Gammajm1_ijm1_kwomwi + Pha_Hmn(m,n,l-k+L+1) + Fi_ij_lwomwi;
                            Pha_pwi = Pha_Gammajm1_ijm1_kwopwi + Pha_Hmn(m,n,l-k+L+1) + Fi_ij_lwopwi;
                            
                            Xk_mwi = Xk_mwi + Ak_mwi*cos(Pha_mwi);
                            Xk_pwi = Xk_pwi + Ak_pwi*cos(Pha_pwi);
                            Yk_mwi = Yk_mwi + Ak_mwi*sin(Pha_mwi);
                            Yk_pwi = Yk_pwi + Ak_pwi*sin(Pha_pwi);
                        end
                        Amp_Gammaj_ijm1TOijmn_lwomwi = norm([Xk_mwi Yk_mwi]); 
                        Amp_Gammaj_ijm1TOijmn_lwopwi = norm([Xk_pwi Yk_pwi]);
                        Pha_Gammaj_ijm1TOijmn_lwomwi = angle(Xk_mwi+1i*Yk_mwi);
                        Pha_Gammaj_ijm1TOijmn_lwopwi = angle(Xk_pwi+1i*Yk_pwi);
                        
                        X_lwomwi = X_lwomwi + xb(m,ij) * xb(n,ijm1) * Amp_Gammaj_ijm1TOijmn_lwomwi * cos(Pha_Gammaj_ijm1TOijmn_lwomwi);
                        X_lwopwi = X_lwopwi + xb(m,ij) * xb(n,ijm1) * Amp_Gammaj_ijm1TOijmn_lwopwi * cos(Pha_Gammaj_ijm1TOijmn_lwopwi);
                        Y_lwomwi = Y_lwomwi + xb(m,ij) * xb(n,ijm1) * Amp_Gammaj_ijm1TOijmn_lwomwi * sin(Pha_Gammaj_ijm1TOijmn_lwomwi);
                        Y_lwopwi = Y_lwopwi + xb(m,ij) * xb(n,ijm1) * Amp_Gammaj_ijm1TOijmn_lwopwi * sin(Pha_Gammaj_ijm1TOijmn_lwopwi);
                    end
                end
                Amp_Gamma_ijm1TOij_ijm1(2*l+2*L+1) = norm([X_lwomwi Y_lwomwi]); 
                Amp_Gamma_ijm1TOij_ijm1(2*l+2*L+2) = norm([X_lwopwi Y_lwopwi]);
                Pha_Gamma_ijm1TOij_ijm1(2*l+2*L+1) = angle(X_lwomwi+1i*Y_lwomwi);
                Pha_Gamma_ijm1TOij_ijm1(2*l+2*L+2) = angle(X_lwopwi+1i*Y_lwopwi);
            end
            X_j_ij = X_j_ij + Amp_Gamma_ijm1TOij_ijm1 .* cos(Pha_Gamma_ijm1TOij_ijm1);
            Y_j_ij = Y_j_ij + Amp_Gamma_ijm1TOij_ijm1 .* sin(Pha_Gamma_ijm1TOij_ijm1);
        end
        Amp_Gamma_j_ij(ij,:,j) = sqrt(X_j_ij.^2 + Y_j_ij.^2);
        Pha_Gamma_j_ij(ij,:,j) = angle( X_j_ij + 1i*Y_j_ij );
    end
end
%% Calculate the feedback coefficients
% ！！！！ Calculate the feedback coefficient of each mode
Alpha_j_i = zeros(J,1);     % Amplitude of the feedback coefficient, i.e., the feedback coefficient α(j)i in the paper
Alpha_Phi_j_i = zeros(J,1); % Phase of the feedback coefficient
Xall = 0; Yall = 0;
for j = 1:J
    X = Amp_Gamma_j_ij(i,2*L+2,j) * cos(Pha_Gamma_j_ij(i,2*L+2,j)) + Amp_Gamma_j_ij(i,2*L+1,j) * cos(-Pha_Gamma_j_ij(i,2*L+1,j));
    Y = Amp_Gamma_j_ij(i,2*L+2,j) * sin(Pha_Gamma_j_ij(i,2*L+2,j)) + Amp_Gamma_j_ij(i,2*L+1,j) * sin(-Pha_Gamma_j_ij(i,2*L+1,j));
    Alpha_j_i(j) = norm([X Y]);
    Alpha_Phi_j_i(j) = angle(X + 1i*Y);
    
    Xall = Xall + X;
    Yall = Yall + Y;
end

figure(1); plot(1:J,Alpha_j_i)
xlabel('Number of feedbacks'); ylabel('Feedback coefficients')

disp(['Feedback coefficient:    αji =    ',num2str(Alpha_j_i')])

% ！！！！ Calculate the general feedback coefficient, i.e., the low-frequency vibration level (LVL) factor, αi
Xall = Xall + cos(fii); Yall = Yall + sin(fii); 
Alpha_i = norm([Xall Yall]); % LVL factor, αi
disp(['LVL factor    αi =    ',num2str(Alpha_i)])

%% Calculate the low-frequency vibration divergence (LVD) factor, βi, to judge the LFC
Betai = 0; 
for j = 1:J-2 
    Betai = Betai + log( mean([Alpha_j_i(j+1),Alpha_j_i(j+2)]) / mean([Alpha_j_i(j),Alpha_j_i(j+1)]) );
end
Betai = exp( Betai / (J-2) ); % LVD factor, βi
disp(['LVD factor    βi =    ',num2str(Betai)])

