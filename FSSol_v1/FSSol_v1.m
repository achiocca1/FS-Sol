%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% FSSol_v1.0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directoryRES = 'Results'; % Directory where all the ANSYS files are located
blocklength = 2;                     % How many load step are present in the .csv file
ref_times = [1,2];                   % Reference load steps in the .csv file
kFS  = 0.4;                          % Material constant of Fatemi-Socie
Sy = 355;                            % Material yield strength (MPa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% START DATA READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Files = dir(directoryRES);
% Read Stress Strain from ANSYS files
FileNames = [];
for K=1:length(Files)

    if Files(K).bytes == 0
        continue % Eliminate the files with zero bytes extension
    else
        file = Files(K).name;
        FileNames = [FileNames, string(file)];
    end
end
% Import data from files
index = 1; % Necessary for preallocation of matrix seizes
f = waitbar(0,'Simulation in progress...'); % PROGRESS BAR Option

FS_MAX_vect = zeros( length(FileNames), 1); % Preallocation of FS values
Angles_vect = zeros( length(FileNames), 4); % Preallocation of Angles values
NODENUMB_vect = zeros( length(FileNames), 1); % Preallocation of node number values
for ii = 1 : length(FileNames)
    file = csvread(strcat(directoryRES,"\",FileNames(ii)),2,1);
    nodenumber = sscanf(FileNames(ii), 'risultati_%d.csv'); % get the node number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% FINISH DATA READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stress and Strain matrix creation trough cell array
    % Preallocate cell memory
    E = cell(1, blocklength);
    S = cell(1, blocklength);
    for j = 1 : blocklength
        E0 = [file(j,1) + file(j,7), (file(j,4) + file(j,10))/2, (file(j,6) + file(j,12))/2;
            (file(j,4) + file(j,10))/2, file(j,2) + file(j,8), (file(j,5) + file(j,11))/2;
            (file(j,6) + file(j,12))/2, (file(j,5) + file(j,11))/2, file(j,3) + file(j,9)];
        E{j} = E0; % Strain Tensor
    end
    for j = 1 : blocklength
        S0 = [file(j,13), file(j,16), file(j,18);
            file(j,16), file(j,14), file(j,17);
            file(j,18), file(j,17), file(j,15)];
        S{j} = S0; % Strain Tensor
    end

index = index + 1;
waitbar(ii/length(FileNames))
%%
%Delta eps
DeltaE = E{ref_times(1)} - E{ref_times(2)};
% Calculate and sort principal strains
[V, D] = eig(DeltaE, 'vector');
[D0, ind] = sort(D, 'descend');
V0 = V(:, ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  ATTENTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This filter is necessary because there is a bug in the "eig" 
% function that does not always produce right eigenvectors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(abs(D0(1)) - abs(D0(3))) > 1e-6 && det(V0) > 0
elseif abs(abs(D0(1)) - abs(D0(3))) < 1e-6 && det(V0) < 0 % FOR PURE TORSION
else
    DeltaE = E{ref_times(2)} - E{ref_times(1)};
    [V, D] = eig(DeltaE, 'vector');
    [D0, ind] = sort(D, 'descend');
    V0 = V(:, ind);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate and sort principal stresses at load step 1
[W01, G1] = eig(S{ref_times(1)});
G01 = diag(G1);
% Calculate and sort principal stresses at load step 2
[W02, G2] = eig(S{ref_times(2)});
G02 = diag(G2);

EigenvalE = sort(D0, 'descend');
EigenvalS1 = sort(G01, 'descend'); % Stress Tensor 1
EigenvalS2 = sort(G02, 'descend'); % Stress Tensor 2

tolerance = 0.01*max(abs(EigenvalE));
%% Warning and error messages %%
if abs(EigenvalE(1) - EigenvalE(2)) < tolerance && abs(EigenvalE(2) - EigenvalE(3)) < tolerance && abs(EigenvalE(1) - EigenvalE(3)) < tolerance
    error('In this case FS = 0! No critical plane exists')
elseif abs(EigenvalE(1) - EigenvalE(2)) < tolerance || abs(EigenvalE(2) - EigenvalE(3)) < tolerance || abs(EigenvalE(1) - EigenvalE(3)) < tolerance
    warning('The value of FS found is correct! However, infinite critical plane orientations may exist in addition to those found. Use the standard scanning plane method to look for all the possible plane orientations.')
end


a = (EigenvalE(1) - EigenvalE(3))/2;
%% Load step 1 (S1)
b = 0.5*(EigenvalS1(1) + EigenvalS1(3))/Sy;
c = 0.5*(EigenvalS1(1) - EigenvalS1(3))/Sy;
k = kFS;

if b >= -1/k
    FS_MAX_1 = (a*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)+3*b*k+3)*sqrt((k*(b*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)-2)+b^2*(-k)+4*c^2*k)+sqrt(k^2*(b^2+8*c^2)+2*b*k+1)-1)/(c*k)))/(8*sqrt(2)*sqrt(c)*sqrt(k));
    x = (sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)-b.*k-1)./(c.*k);
    y = (sqrt(2).*sqrt((k.*(b.*(sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)-2)+b.^2.*(-k)+4.*c.^2.*k)+sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)-1)./(c.*k)))./(sqrt(c).*sqrt(k));
    OMEGA_MAX = 0.5*atan(y/x);
else
    FS_MAX_1 = (a*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)-3*b*k-3)*sqrt(-(k*(b*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)+2)+b^2*k-4*c^2*k)+sqrt(k^2*(b^2+8*c^2)+2*b*k+1)+1)/(c*k)))/(8*sqrt(2)*sqrt(c)*sqrt(k));
    x = -(sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)+b.*k+1)./(c.*k);
    y = -(sqrt(2).*sqrt(-(k.*(b.*(sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)+2)+b.^2.*k-4.*c.^2.*k)+sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)+1)./(c.*k)))./(sqrt(c).*sqrt(k));
    OMEGA_MAX = 0.5*atan(y/x);
end
% + OMEGA max angle
% Calculate the angles based on the rotation R = Rotz(Psi)*Roty(Theta)
RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
    0        1        0;
    -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

Matr = V0*RotY;

Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
idx = ( Psi>= 0) & ( Psi<= 2*pi);
Psi = Psi(idx);

for i = 1:length(Psi)
    SINTHETA = Matr(2,3)/sin(Psi(i));
    Theta = atan2(SINTHETA, Matr(3,3));
    if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
        Theta_1 = Theta;
        Psi_1 = Psi(i);
        break
    else
    end
    Theta_1 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
    Psi_1 = atan2(Matr(2,3), Matr(1,3));
end


% - OMEGA max angle
OMEGA_MAX = - OMEGA_MAX;
% Calculate the angles based on the rotation R = Rotz(Psi)*Roty(Theta)
RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
    0        1        0;
    -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

Matr = V0*RotY;

Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
idx = ( Psi>= 0) & ( Psi<= 2*pi);
Psi = Psi(idx);

for i = 1:length(Psi)
    SINTHETA = Matr(2,3)/sin(Psi(i));
    Theta = atan2(SINTHETA, Matr(3,3));
    if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
        Theta_2 = Theta;
        Psi_2 = Psi(i);
        break
    else
    end
    Theta_2 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
    Psi_2 = atan2(Matr(2,3), Matr(1,3));
end

%% Load step 2 (S2)
b = 0.5*(EigenvalS2(1) + EigenvalS2(3))/Sy;
c = 0.5*(EigenvalS2(1) - EigenvalS2(3))/Sy;

if b >= -1/k
    FS_MAX_2 = (a*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)+3*b*k+3)*sqrt((k*(b*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)-2)+b^2*(-k)+4*c^2*k)+sqrt(k^2*(b^2+8*c^2)+2*b*k+1)-1)/(c*k)))/(8*sqrt(2)*sqrt(c)*sqrt(k));
    x = (sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)-b.*k-1)./(c.*k);
    y = (sqrt(2).*sqrt((k.*(b.*(sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)-2)+b.^2.*(-k)+4.*c.^2.*k)+sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)-1)./(c.*k)))./(sqrt(c).*sqrt(k));
    OMEGA_MAX = 0.5*atan(y/x);
else
    FS_MAX_2 = (a*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)-3*b*k-3)*sqrt(-(k*(b*(sqrt(k^2*(b^2+8*c^2)+2*b*k+1)+2)+b^2*k-4*c^2*k)+sqrt(k^2*(b^2+8*c^2)+2*b*k+1)+1)/(c*k)))/(8*sqrt(2)*sqrt(c)*sqrt(k));
    x = -(sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)+b.*k+1)./(c.*k);
    y = -(sqrt(2).*sqrt(-(k.*(b.*(sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)+2)+b.^2.*k-4.*c.^2.*k)+sqrt(k.^2.*(b.^2+8.*c.^2)+2.*b.*k+1)+1)./(c.*k)))./(sqrt(c).*sqrt(k));
    OMEGA_MAX = 0.5*atan(y/x);
end
% + OMEGA max angle
% Calculate the angles based on the rotation R = Rotz(Psi)*Roty(Theta)
RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
    0        1        0;
    -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

Matr = V0*RotY;

Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
idx = ( Psi>= 0) & ( Psi<= 2*pi);
Psi = Psi(idx);

for i = 1:length(Psi)
    SINTHETA = Matr(2,3)/sin(Psi(i));
    Theta = atan2(SINTHETA, Matr(3,3));
    if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
        Theta_3 = Theta;
        Psi_3 = Psi(i);
        break
    else
    end
    Theta_3 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
    Psi_3 = atan2(Matr(2,3), Matr(1,3));
end

% - OMEGA max angle
OMEGA_MAX = - OMEGA_MAX;
% Calculate the angles based on the rotation R = Rotz(Psi)*Roty(Theta)
RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
    0        1        0;
    -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

Matr = V0*RotY;

Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
idx = ( Psi>= 0) & ( Psi<= 2*pi);
Psi = Psi(idx);

for i = 1:length(Psi)
    SINTHETA = Matr(2,3)/sin(Psi(i));
    Theta = atan2(SINTHETA, Matr(3,3));
    if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
        Theta_4 = Theta;
        Psi_4 = Psi(i);
        break
    else
    end
    Theta_4 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
    Psi_4 = atan2(Matr(2,3), Matr(1,3));
end

FS_MAX_vect(ii,:) = max(FS_MAX_1, FS_MAX_2);
if FS_MAX_1 > FS_MAX_2
    X = [Theta_1,Psi_1];
    Y = [Theta_2,Psi_2];
else
    X = [Theta_3,Psi_3];
    Y = [Theta_4,Psi_4];
end
Angles_vect(ii,:) = [X, Y];
NODENUMB_vect(ii,:) =  nodenumber;
end
delete(f)

[FS_MAX, ind] = max(FS_MAX_vect);
Theta_1 = Angles_vect(ind,1);
Psi_1 = Angles_vect(ind,2);
Theta_2 = Angles_vect(ind,3);
Psi_2 = Angles_vect(ind,4);
nodenumber_max = NODENUMB_vect(ind,1);
Q = ['The critical plane factor is maximum for the node numbered ',  num2str(nodenumber_max), ' and has a value of FS = ', num2str(FS_MAX)];
X = ['  - the first critical plane orientation is given by the angles θ = ', num2str(Theta_1), ' and Ψ = ', num2str(Psi_1)];
Y = ['  - the second critical plane orientation is given by the angles θ = ', num2str(Theta_2), ' and Ψ = ', num2str(Psi_2)];
Z = 'Considering a rotation RotZ(Ψ)RotY(θ):';
disp(Q)
disp(Z)
disp(X)
disp(Y)