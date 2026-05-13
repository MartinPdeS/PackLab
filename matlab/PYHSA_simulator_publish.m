%% Function description

% this function computes the scattering properties of a suspension of
% spherical scatterers considering both independent scattering (valid at low 
% volume fraction) and dependent scattering (always valid but relevant at 
% high volume fraction). Dependent scattering is modeled using the
% Percus-Yevick hard-sphere approximation (PYHSA).

%input: cfg struct containing:

%cfg.n_medium:      (1)     refractive index of medium          [-]
%cfg.r_vec:         (Nx1)   size vector for PSD (bin centers)   [m]
%cfg.n_vec:         (Nx1)   number density vector for PSD       [1/m3]
%cfg.lam_test:      (1)     wavelength (vacuum)                 [m]
%cfg.Nt:            (1)     number of points for theta in Mie   [-]
%cfg.Mie            {Nx2}   cell array containing S_matrices & sigma_s

%output: output struct containing:

%output.mu_s_indep  (1)     scat. coeff. under IS               [1/m]
%output.mu_s_dep    (1)     scat. coeff. under DS               [1/m]
%output.g_indep     (1)     anisotropy under IS                 [-]
%output.g_dep       (1)     anisotropy under DS                 [-]
%output.P_indep     (Ntx1)  phase function under IS             [-]
%output.P_dep       (Ntx1)  phase function under DS             [-]

function [output] = PYHSA_simulator_publish(n_medium,r_vec,n_vec,lam_test,Nt,Mie_sols)

%wavenumber in medium
k0 = 2*pi*n_medium/lam_test;
%linear sampling of scattering vector
p_vec   = linspace(0,2*k0,Nt);
%defining angular values in radians
ang_rad = 2*asin(p_vec/2/k0);

%computing partial correlations
[H_save] = H_mat_calc_fast(p_vec,r_vec,n_vec);

%init. tmp Mie vectors
S1_tmp  = zeros(length(r_vec),length(p_vec));
S2_tmp  = zeros(length(r_vec),length(p_vec));
C_tmp   = zeros(length(r_vec),1);

for i = 1:length(r_vec)
    %storing Mie results
    S1_tmp(i,:) = Mie_sols{i,1}(1,1,:);
    S2_tmp(i,:) = Mie_sols{i,1}(2,2,:);
    C_tmp(i)    = Mie_sols{i,2}.sca;
end

%independent scattering component
mu_s_indep  = sum(n_vec.*C_tmp);
P_indep     = sum(n_vec/2/k0^2.*(abs(S1_tmp).^2+abs(S2_tmp).^2),1);
%expanding Mie form factor
F_Mie_vec(1,:,:) =  1i/k0*S1_tmp;
F_Mie_vec(2,:,:) = -1i/k0*S2_tmp;
%dependent scatetring
F_Mie_vec_alpha = repmat(F_Mie_vec,[1 length(r_vec) 1]);
F_Mie_vec_beta  = repelem(F_Mie_vec,1,length(r_vec),1);
F_Mie_prod      = 1/2*squeeze(sum(F_Mie_vec_alpha.*conj(F_Mie_vec_beta),1));
F_Mie_prod      = F_Mie_prod.*reshape(H_save,prod(size(H_save,[1 2])),size(H_save,3));
n1n2_vec        = sqrt(n_vec.*n_vec');
sum_tot         = sum(F_Mie_prod.*n1n2_vec(:),1);

%dependent phase function includes indep and dep components
P_dep               = real(P_indep + sum_tot);
%normalizing
P_indep             = P_indep/(2*pi*(trapz(ang_rad,P_indep.*sin(ang_rad))));
P_dep               = P_dep/(2*pi*(trapz(ang_rad,P_dep.*sin(ang_rad))));
%computing anisotropy
g_indep             = 2*pi*trapz(ang_rad,P_indep.*cos(ang_rad).*sin(ang_rad));
g_dep               = 2*pi*trapz(ang_rad,P_dep.*cos(ang_rad).*sin(ang_rad));

%storing results
output.mu_s_indep   = mu_s_indep;
output.mu_s_dep     = mu_s_indep + 2*pi*trapz(ang_rad,real(sum_tot).*sin(ang_rad));
output.P_indep      = P_indep';
output.P_dep        = P_dep';
output.g_indep      = g_indep;
output.g_dep        = g_dep;


end

%% Function description

%this function computes the Fourier Transform of the total correlation for
%a polydisperse sample (H-matrix). Each row/column corresponds to one
%particle size and each element is the FT of the total correlation between
%the particle size (row) and the particle size (column).

%in this function this matrix is computed for every value of the scattering
%vector p_vec. The H-matrix is computed by solving Eq. 8.2.7 in chapter 8,
%section 2 of the reference below. Equation numbering may change for
%different editions.

%Tsang, L., Kong, J. A., Ding, K. H., & Ao, C. O. (2004). Scattering of 
%electromagnetic waves: numerical simulations. John Wiley & Sons.

%input: 

%p_vec:     	[Ntx1]  scattering vector
%r_vec:         [Nx1]   radius vector for particle sizes
%n_vec:         [Nx1]   number density vector

%output:

%H_save:        [NxNxNt] 3D H-matrix (one slice per scattering vector element)

function [H_save] = H_mat_calc_fast(p_vec,r_vec,n_vec)

%initializing H_matrix storage
H_save = zeros(length(r_vec),length(r_vec),length(p_vec));

%comuting H-matrix for all values of p
for i = 1:length(p_vec)
    %computing C-matrix
    tmp_C   = C_mat_calc_fast(p_vec(i),r_vec,n_vec);
    %computing H-matrix (8.2.7)
    tmp_H   = (eye(size(tmp_C))-tmp_C)\tmp_C;
    %storing result
    H_save(:,:,i) = tmp_H;
end

end

%% Function description

%this function computes the Fourier transform of the direct correlation
%matrix (C-matrix). Each row/column corresponds to one
%particle size and each element is the FT of the direct correlation between
%the particle size (row) and the particle size (column).

%The C-matrix is computed following Eqs. 8.2.8-8.2.15 of chapter 8, section
%2 in the reference below. Equation numbering may change for different eds.

%Tsang, L., Kong, J. A., Ding, K. H., & Ao, C. O. (2004). Scattering of 
%electromagnetic waves: numerical simulations. John Wiley & Sons.

%input: 

%p_val:     	[1]     scattering vector
%r_vec:         [Nx1]   radius vector for particle sizes
%n_vec:         [Nx1]   number density vector

%output:

%C_mat:         [NxN]   2D C-matrix

function [C_mat] = C_mat_calc_fast(p_val,r_vec,n_vec)

%computing xi values (8.2.15)
xi_0    = pi/6*sum(n_vec.*(2*r_vec).^0);
xi_1    = pi/6*sum(n_vec.*(2*r_vec).^1);
xi_2    = pi/6*sum(n_vec.*(2*r_vec).^2);
xi_3    = pi/6*sum(n_vec.*(2*r_vec).^3);

%defining all vectors
Ri      = 2*r_vec;
Rj      = Ri';
Xi      = p_val*Ri/2;
Xj      = Xi';
Ni      = Ri.^2.*sinc(Xi/pi);
Nj      = Ni';
Mi      = 3*Ri.^3.*(sin(Xi)./Xi.^3-cos(Xi)./Xi.^2);

%correcting for numerical errors
Mi(Xi<1e-4) = 3*Ri(Xi<1e-4).^3*1/3;
Mj      = Mi';

%(8.2.16) line 1
tmp1 = Mj.*(cos(Xi) + Xi.*sin(Xi) + 3.*xi_2.*Ri.*cos(Xi)/(1-xi_3) + 3.*xi_1.*Ni/(1-xi_3) + 9.*xi_2^2.*Ni/(1-xi_3)^2 );
%(8.2.16) line 2
tmp2 = Mi.*(cos(Xj) + Xj.*sin(Xj) + 3.*xi_2.*Rj.*cos(Xj)/(1-xi_3) + 3.*xi_1.*Nj/(1-xi_3) + 9.*xi_2^2.*Nj/(1-xi_3)^2 );
%(8.2.16) line 3
tmp3 = Mi.*Mj.*(xi_0/(1-xi_3) + p_val^2.*xi_2/4/(1-xi_3) + 6.*xi_1.*xi_2/(1-xi_3)^2 + 9.*xi_2^3/(1-xi_3)^3 );
%(8.2.16) line 4
tmp4 = 3.*Ni.*Rj.*cos(Xj) + 3.*Nj.*Ri.*cos(Xi) + 9.*xi_2.*Ni.*Nj/(1-xi_3);
%(8.2.16) total
C_mat = -pi/6*sqrt(n_vec'.*n_vec)./(1-xi_3).*(tmp1 + tmp2 + tmp3 + tmp4);


end
