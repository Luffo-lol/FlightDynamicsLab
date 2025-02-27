%% Stability Derivative Calculator - Flight Dynamics
% Suhail Halaby - GROUP 7
% All derivatives calculated with regards to the methods introduced in the
% lecture notes by Errikos Levis

clear
clc
close all

% note: partial derivative parameters are given as D_D_
% note: Tx is thrust in the x direction, Tz is in the z directin
% note: V is the freestream velocity (equillibrium) (V = U_e)
% note: _e denotes equillibrium conditions

%% Part 1: Symmetric Stability Derivatives 

% U-derivatives: 
    % X-wise velocity perturbation impact on aerodynamic forces
    % The primary mechanisms at play are compressibility effects and
    % changes to induced aerodynamic forces
   
X_u = @(DTxDV,rho,V,S_ref,Cd_e,DCdDV) (DTxDV) - rho*V*S_ref*Cd_e - 0.5*rho*(V^2)*S_ref*(DCdDV);
Z_u = @(DTzDV,rho,V,S_ref,Cl_e,DClDV) (DTzDV) - rho*V*S_ref*Cl_e - 0.5*rho*(V^2)*S_ref*(DClDV);
M_u = @(l_t,DTDV,rho,V,S_ref,MAC,DCmDV,Cl_e,DKnDV,Kn,DClDV) l_t*(DTDV) + rho*(V^2)*S_ref*MAC*( DCmDV - Cl_e * DKnDV - Kn*(DClDV));

% W-derivatives:
    % Z-wise velocity perturbation impact on aerodynamic forces
    % The primary mechanisms at play involve changes to the angle of attack

X_w = @(V,DTxDalpha,rho,S_ref,DCdDalpha,Cl_e) (1/V)*(DTxDalpha)+0.5*rho*V*S_ref*(-(DCdDalpha)+Cl_e);
Z_w = @(V,DTzDalpha,rho,S_ref,DClDalpha,Cd_e) (1/V)*(DTzDalpha)-0.5*rho*V*S_ref*((DClDalpha)+Cd_e);
M_w = @(V,l_t,DTDalpha,rho,S_ref,MAC,Kn,DClDalpha) (l_t/V)*(DTDalpha) - 0.5*rho*V*S_ref*MAC*Kn*DClDalpha ;

    % note: DCdDalpha_e can be derived from the aircraft's polars
    % note: DCmDalpha_e = -Kn*DClDalpha_e

% Tailplane Only W-derivatives:

X_w_H = @(V,rho,S_H,DCd_HDalpha,Cl_H) 0.5*rho*V*S_H*(-(DCd_HDalpha)+Cl_H);
Z_w_H = @(V,rho,S_H,DCl_HDalpha,Cd_H) -0.5*rho*V*S_H*((DCl_HDalpha)+Cd_H);
M_w_H = @(V,rho,S_H,MAC_H,DCm_HDalpha) 0.5*rho*V*S_H*MAC_H*(DCm_HDalpha);

% q-Derivatives
    % pitch rate perturbation impact on aerodynamic forces
    % The primary mechanisms at play involve changes to the local angle of
    % attack, behaviour is H-stab dominated

X_q = @(rho,V,S_H,l_H,DCd_HDalpha,Cl_H) 0.5*V*S_H*l_H*(-(DCd_HDalpha)+Cl_H);
Z_q = @(rho,V,S_H,h_H,DCl_HDalpha,Cd_H) 0.5*V*S_H*h_H*((DCl_HDalpha)+Cd_H);
M_q = @(l_H,Z_w_H,h_H,X_w_H,M_w_h) l_H^2*Z_w_H - h_H*l_H*X_w_H + l_H*M_w_h;

% W_dot-Derivatives
    % sink acceleration perturbation impact on aerodynamic forces
    % The primary mechanisms at play involve changes to the downwash lag
    % charachteristics

X_w_dot = @(l_H,V,DepsilonDalpha,X_w_H) (l_H/V)*(DepsilonDalpha)*X_w_H;
Z_w_dot = @(l_H,V,DepsilonDalpha,Z_w_H) (l_H/V)*(DepsilonDalpha)*Z_w_H;
M_w_dot = @(l_H,V,DepsilonDalpha,h_H,X_w_H,Z_w_H,M_w_H) (l_H/V)*(DepsilonDalpha)*(-h_H*X_w_H+l_H*Z_w_H+M_w_H);


%% Step 2: Asymmetric Stability Derivatives

info = readmatrix('derivatives.csv');
asym_d = @(rho,V,S_ref,b) infor(:,3)./(0.5*rho*V*S_ref.*(b^info(:,2)));
