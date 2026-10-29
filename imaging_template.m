% Template script for radio astronomy image formation
% Millad Sardarabadi, Sept 2015
%
clear all
close all
clc
%% Loading Data
% The data from LOFAR core will be loaded in this cell.
%
% Rh: A p x p sample covariance matrix, where p = number of antennas
% poslocal: A p x 3 matrix containing the earth-bound coordinates of the
%           receiving elements.
% freq: Central frequency for the measured channel.
% t_obs_matlab: Observation time.
%               (datestr(t_obs_matlab): 11-Jul-2012 22:22:20) 
% lon: Longitude of the LOFAR core
% lat: Latitude of the LOFAR core

load lofar_DSP_data_1.mat

%% Initialization
c = 299792458; % speed of light
p = size(Rh,1);% number of receiving elements
lambda = c / freq; % wavelength at the observation frequency

% define a distance matrix of all baselines (=vectors between pairs of antennas)
dist = sqrt(  (meshgrid(poslocal(:,1)) - meshgrid(poslocal(:,1)).').^2 ...
            + (meshgrid(poslocal(:,2)) - meshgrid(poslocal(:,2)).').^2);
D = max(dist(:)); % Maximum Baseline


%% Coordinate system
% In astronomy it is common to use the (l,m,n) coordinates which 
% are the components of the more familiar direction vector in spherical 
% coordinates:
%    s = [cos(theta) * cos(phi)
%        cos(theta) * sin(phi)
%        sin(theta)]


% define a grid of image coordinates (l,m) at the right resolution
FracAngle = 0.25;% fraction of theoretical limit for angle resolution
                 % 0.25 means ~4x4 pixels in the main lobe
dl = FracAngle * lambda / D; % width/height of a pixel
l = 0:dl:1; % First image coordinate := cos(theta)cos(phi)
l = [-fliplr(l(2:end)),l]';
m = l;  % Second image coordinate := cos(theta)sin(phi)
        % m = l means an aspect ratio 1:1


%% Your Imaging Algorithm should come hereafter...
Z = 2*pi*freq/c * poslocal
len = length(l)
for a = 1:len
    for b = 1:len
        if l(a)^2 + m(b)^2 <= 1
            directionCosines = [l(a), m(b), sqrt(1 - l(a)^2 - m(b)^2)];
            
            steeringVector = exp(-1i * Z * directionCosines.');
            image(a,b) = 1 / real(steeringVector' * (Rh \ steeringVector));
            image2(a,b) =  real(steeringVector' * (Rh \ steeringVector)) / real(steeringVector' * (Rh \ (Rh \ steeringVector)))^2;
            image3(a,b) =  real(steeringVector' * (Rh \ steeringVector)) / real(steeringVector' * (Rh \ (Rh \ steeringVector)));

        else
            image(a,b) = NaN;
            image2(a,b) = NaN;
            image3(a,b) = NaN;
        end
    end
end

figure;
h = imagesc(image);

% 2. Personalizar la barra de colores y los ejes
colorbar;
title('Mapa de valores con NaN');
xlabel('Columna');
ylabel('Fila');

figure;
h2 = imagesc(image2);

% 2. Personalizar la barra de colores y los ejes
colorbar;
title('Mapa de valores con NaN');
xlabel('Columna');
ylabel('Fila');

figure;
h3 = imagesc(image3);

% 2. Personalizar la barra de colores y los ejes
colorbar;
title('Mapa de valores con NaN');
xlabel('Columna');
ylabel('Fila');