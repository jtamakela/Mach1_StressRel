
% For analyse of mechanical data from Mach-1
% Stress-relaxation
% 4 steps
% Indentation

clear all, close all, clc

% % %% %% %% %% %% %%
velocity = 100; %STRAIN VELOCITY IN PERCENTS
stepsize = 5; %STEP SIZE IN PERCENTS
% % %% %% %% %% %% %%
r = 0.001/2; %RADIUS
% % %% %% %% %% %% %%
Sampling_freq = 100;

% File name
measurement = 'Stress_relaxation_test.txt';
% Containing folder
datafolder = '/media/janne/Data/UEF/Measurements/Ali/Ponies/Mechanics/Preliminary_ponies';


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

kolors = 'brgk'; % For plotting
scrsz = get(0,'ScreenSize'); % For plotting


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Importing stress-relaxation data % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Has to be done with textread 
no_of_headerlines = 25;
fid = fopen([datafolder, '/', measurement]); %Opens file1 for read access
temp=fread(fid, [1 inf]); % luetaan tiedosto
fclose(fid);

%These are used to find steps in the data. 
%<divider> divides stress-relaxation steps
alut = findstr('<DATA>', temp);
dividers = findstr('<divider>', temp);
loput = findstr('<END DATA>', temp);

temp_Stress = char(temp(alut(1)+75:dividers(1)));
Stressrelaxation{1} = textscan(temp_Stress, '%f %f %f %f %f'); %Muutetaan matriisiksi
for i = 1:3;
    temp_Stress = char(temp(dividers(i)+10: dividers(i+1)));
    Stressrelaxation{i+1} = textscan(temp_Stress, '%f %f %f %f %f'); %Muutetaan matriisiksi
end

% Plotting the imported data
figure(1);
for i = 1:4
    yyaxis left
    plot(Stressrelaxation{1,i}{1,1}, Stressrelaxation{1,i}{1,5}, '-', 'color', kolors(i))
    ylabel('F (g)');

    hold on;
    yyaxis right
    plot(Stressrelaxation{1,i}{1,1}, Stressrelaxation{1,i}{1,2}, '-')
    ylabel('Displacement (mm)');
end
xlim([-50 3650]);
xlabel('t (s)');


% Combining the imported data
t = [];
h = [];
F = [];

for i = 1:4
    t = [t; Stressrelaxation{1,i}{1,1}];
    h = [h; Stressrelaxation{1,i}{1,2}];
    F = [F; Stressrelaxation{1,i}{1,5}];
end

%Conversion from g to N
F = F/1000*9.81;

% ANALYSIS


% Area
A = pi*r^2;


%Stress (Pressure) P=F/A
P = F./A;


figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/3 scrsz(4)/3])

subplot(2,1,1)
ax = plotyy(t, P,t, h,'plot'); %Plotting time in seconds

hold on;

set(get(ax(1),'Ylabel'),'String','Stress (Pa)') 
set(get(ax(2),'Ylabel'),'String','Displacement (mm)') 

xlabel('Time (s)');


max_ind = length(P); %Size of measurement
max_windows = [1 round(max_ind*1/4) round(max_ind*1/2) round(max_ind*3/4) max_ind]; %Windows for peak values

% Peak values and indexes
[P_dyn(1) ind_dyn(1)] = max(P(1:max_windows(2)));
[P_dyn(2) ind_dyn(2)] = max(P(max_windows(2):max_windows(3)));
[P_dyn(3) ind_dyn(3)] = max(P(max_windows(3):max_windows(4)));
[P_dyn(4) ind_dyn(4)] = max(P);
ind_dyn = [ind_dyn(1) ind_dyn(2)+max_windows(2) ind_dyn(3)+max_windows(3) ind_dyn(4)];

plot(t(ind_dyn), P_dyn, 'bo','markersize', 10, 'LineWidth', 3);

ind_eq_end = [round(ind_dyn - 7) max_ind]; %step time (~16) + 1 second
ind_eq_Start = ind_eq_end-20; %2s window for mean P calculation
ind_eq = [mean([ind_eq_Start; ind_eq_end])]; %For plotting. Center point


for k = 1:5;
    P_eq(k) = mean(P(ind_eq_Start(k):ind_eq_end(k))); %Calculate mean P values at equilibrium instead of taking one value
end


plot(t(ind_eq), P_eq, 'ro','markersize', 10, 'LineWidth', 3);


% %%
% Taring and zeroing
% Calculating thickness based on displacement

%Tare 1 F_eq out
P = P-P_eq(1);

for k = 1:5;
    h_eq(k) = mean(h(ind_eq_Start(k):ind_eq_end(k))); %Calculate mean h from 100 seconds to equilibrium
end

h_zeroed = h_eq-h_eq(1); %Subtracting 0-point before actual measurements


% %% Assumption made -> strain 
strain = [0 0.05 0.1 0.15 0.2];

%Thickness based on the assumption that strains are 0-0.2. 
h_plug = h_zeroed/strain;

% %% 
% Least-squares fitting for Pressure. Young's modulus calculated from that then

S = [strain', ones(length(strain),1)];
E_eq = (S\P_eq');

plot(t(ind_eq), strain.*E_eq(1)+E_eq(2), 'g', 'LineWidth', 2);

% And the same without 1 point %%%%%%%%%%%%%%%%%%%%%%%%

S = [strain(2:end)', ones(length(strain(2:end)),1)];
E_eq_w1 = (S\P_eq(2:end)');

plot(t(ind_eq(2:end)), strain(2:end).*E_eq_w1(1)+E_eq_w1(2), 'r--', 'LineWidth', 2);


% Dynamic Moduli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P_diff = [P_dyn - P_eq(1:end-1)];
E_dyn = P_diff./0.05; %Strain 5%


%Plotting %% %% %%
subplot(2,1,2)
plot(strain*100, P_eq./1e3, 'ro','markersize', 10, 'LineWidth', 3);
title('Equilibrium stiffnesses')
ylabel('Stress [kPa]');
xlabel('Strain [%]');
hold on;
plot(strain*100, (E_eq(1)./1e3)*strain+(E_eq(2)./1e3), 'g', 'LineWidth', 2);
plot(strain*100, (E_eq_w1(1)./1e3)*strain+(E_eq_w1(2)./1e3), 'r', 'LineWidth', 2);
legend('equilibrium stiffness', 'LS-Fit', 'Fit without 1. value', 'Location', 'NW')


% %% 
% Displaying values
disp(' ')
% disp(['Sample ', fn(i).name]);
disp(['Calculated thickness = ', num2str(h_plug), ' mm'])
disp(' ')
disp(['Equilibrium Modulus = ', num2str(E_eq(1)/1e6), ' MPa'])
disp(['Equilibrium Modulus without 1. point = ', num2str(E_eq_w1(1)/1e6), ' MPa'])
disp(' ')
disp(['Dynamic moduli = ', num2str(E_dyn/1e6), ' MPa'])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('');

%Results %% %% %% %%
% Final_name{i,:} = fn(i).name;
Final_thickness(i,:) = h_plug;
Final_E_eq(i,:) = E_eq(1);
Final_E_eq_w1(i,:) = E_eq_w1(1);
Final_E_dyn(i,1:4) = E_dyn;

Final_all = [Final_thickness Final_E_eq Final_E_eq_w1 Final_E_dyn];


%close(h_wait); %Close waitbar


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps
% For modeling purposes

%Strain rate is 5um/s

strain_time = h_plug*0.5/0.05;

for i = 1:4;
    step_beginning(i) = max(find(t<t(ind_dyn(i))-strain_time));
end

steptimes = [t(step_beginning); t(end)];

figure; 
for i = 1:4;
subplot(2,2,i);
plot(t,h);
hold on;
plot(t(ind_dyn(i)), h(ind_dyn(i)), 'ro','markersize', 10, 'LineWidth', 3);

plot(t(step_beginning(i)), h(step_beginning(i)), 'go','markersize', 10, 'LineWidth', 3);
xlim([steptimes(i)-5 t(ind_dyn(i))+5]);
title(['step ', num2str(i)]);
end


%plot(t(ind_dyn), P_dyn, 'bo','markersize', 10, 'LineWidth', 3);

disp('');
disp(['Displacements start at = ', num2str(steptimes')])
disp(['Step time = ', num2str(strain_time)])
disp(['Step size = ', num2str(0.05*h_plug), 'mm']) 
disp('');

StepsForInp = [steptimes(1) steptimes(2)-steptimes(1)-strain_time...
    steptimes(3)-steptimes(2)-strain_time steptimes(4)-steptimes(3)-strain_time...
    steptimes(5)-steptimes(4)-strain_time];




% % % % % % % % % % % % % % % %% Dynamic Sinusoidal
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % %Dynamic data time restarts after every step
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % for i = 2:length(alut)-1 %Removing first (stress-relaxation) and the last (lift).
% % % % % % % % % % % % % % %     te = char(temp(alut(i+1)+75:loput(i+1)));
% % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % %     Dynamic{i} = textscan(te, '%f %f %f %f %f'); %Muutetaan matriisiksi
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % 1 Time, s; 2 Position (z), mm; 3 Position (x), mm; 4 Position (y), mm; 5 Fz, gf
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % figure(1);
% % % % % % % % % % % % % % % for i = 2:length(alut)-1 
% % % % % % % % % % % % % % %     plot(Dynamic{1,i}{1,1}, Dynamic{1,i}{1,5})
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % 













