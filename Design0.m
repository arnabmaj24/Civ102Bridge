clear; close all;

%% 0. Initialize Parameters
L = 1200; % Length of bridge
n = 1200; % Discretize into 1 mm seg.
P = 400; % Total weight of train

x = linspace(0, L, n+1);

x_train = [52 228 392 568 732 908] * 10^-3;
x_trainCentered = [52+120 228+120 392+120 568+120 732+120 908+120];
x_train_start = [1 177 341 517 681 857];% Train Load Locations
P_train_withloc = [66.65 66.65 66.65 66.65 90 90];
x_train_start0 = [1 177 341 517 681 857]-856;% Train Load Locations
%by load cases
loc_wieght = 420/1.35;
P_trainLC1 = [1 1 1 1 1 1] * P/6;
P_trainLC2 = [loc_wieght/2 loc_wieght/2 loc_wieght/2 loc_wieght/2 1.35*loc_wieght/2 1.35*loc_wieght/2];
sum(P_trainLC2)
P_trainLC2
n_t = 344; % num of locations
SFDlist = zeros(n_t, n+1); %n_t by n+1 array, each row is an sfd
BMDlist = zeros(n_t, n+1); %same


%UDL Because quirky
UDL = 3.01 * 10^-3; % N / mm
UseUDL = false;




%for all locations with local motive
Postemp = x_train_start0;
Ptemp = P_trainLC1;
domain = n+856;

SFE = zeros(1, domain);
BME = zeros(1, domain);

SFDlist = zeros(domain, n+1); %n_t by n+1 array, each row is an sfd
BMDlist = zeros(domain, n+1); %same


for i = 1:domain
    
    pos = zeros(1, length(Postemp));
    P = zeros(1, length(Ptemp));

    for u = 1:length(pos)
        if Postemp(u) > 0 && Postemp(u) < n+2
            pos(u) = Postemp(u);
            P(u) = Ptemp(u);
        end

    end

    %reaction forces
    R_b = (dot((pos* 10^-3), P))/(1.2);
    R_a = sum(P) - R_b;


    %SFD
    SFD(1:n+1) = R_a; %start
    for g = 1:length(pos)
        if pos(g) > 0
            %cases for when train is partially on bridge
            if pos(g) == 1
                SFD(pos(g):n+1) = SFD(pos(g):n+1) - P(g);
            end
            if pos(g) > 1
                SFD(pos(g)+1:n+1) = SFD(pos(g)+1:n+1) - P(g);
            end
        end
        
    end
    %setting it to start at 0
    SFD(1) = 0;
    SFD(end) = SFD(end) + R_b; %final reac force
    
    %BMD calc
    BMD = zeros(1, n+1);
    for j = 1:n+1
        BMD(j) = sum(SFD(1:j))*10^-3; %changing from M to mm
    end 
    BMD;

    %for next iteration
    Postemp(1:length(pos)) = Postemp(1:length(pos)) + 1;

    %storage
    SFDlist(i,1:n+1) = SFD;
    BMDlist(i,1:n+1) = BMD;

end

SFDlist(:,n+1) = 0;

%Shear force envelope
SFE = zeros(1,n+1);
for v = 1:n+1
    if abs(max(SFDlist(1:end,v))) > abs(min(SFDlist(1:end,v)))
    SFE(v) = (max(SFDlist(1:end,v)));
    end
    if abs(max(SFDlist(1:end,v))) < abs(min(SFDlist(1:end,v)))
    SFE(v) = (min(SFDlist(1:end,v)));
    end

end

%Bending moment envelope
BME = zeros(1, n+1);
for w = 1:n+1
    BME(w) = max(BMDlist(1:end, w));
end

% figure
% hold on; grid on; grid minor;
% plot(1:n+1, SFE, "g-")
% % plot(1:1400,zeros(1,n_t), "black-")
% xlabel('Length (mm)')
% ylabel('Max Shear Force (N)')
% % 
% figure
% hold on; grid on; grid minor;
% plot(1:n+1,BME,"b-")
% % plot(1:1400,zeros(1,1400), "black-")
% xlabel('Length (mm)')
% ylabel('Moment (Nm)')
% 



%Y bar and I
numsecs = 1;

cross_section_0 = [80, 1.27, 1.27/2;
                            1.27, 75-1.27, (75-1.27)/2+1.27;
                            1.27, 75.127, (75-1.27)/2+1.27;
                            5, 1.27, 75-1.27/2;
                            5, 1.27, 75-1.27/2;
                            100, 1.27, 75+1.27/2];
                   
cross_sections = {cross_section_0};
Is = zeros(1,numsecs);
Ys = zeros(1, numsecs);
for y = 1:numsecs
    num_shapes = length(cross_sections{y});
    cross_section = cross_sections{y};
    cross_section =cross_section';    
    sum_num = 0;
    sum_den = 0;
    for i = 1:num_shapes
        sum_num = sum_num + cross_section(1,i)*cross_section(2,i)*cross_section(3,i);
        sum_den = sum_den + cross_section(1,i)*cross_section(2,i);
    end
    
    y_bar = sum_num / sum_den
    I = 0;
    for j = 1:num_shapes
        I = I + ((cross_section(1,j)*cross_section(2,j)^3)/12) + (cross_section(1,j)*cross_section(2,j)*(y_bar-cross_section(3,j))^2);
    end
    Is(y) = I
    Ys(y) = y_bar
end





%flextural stresses and more shear stuff

Maxflexprofile = zeros(2, n+1);
Iprofile = zeros(1, n+1);
Iprofile(1:n+1) = Is(1);
% Iprofile(301:900) = Is(2);
Yprofile = zeros(2, n+1);
Yprofile(1, 1:n+1) = 76.127;
% Yprofile(1, 301:900) = 200.738;
Yprofile(2, 1:n+1) = Ys(1);
% Yprofile(2, 301:900) = Ys(2);

for q = 1:n+1
    Maxflexprofile(1, q) = -BME(q) * (abs(Yprofile(1,q)-Yprofile(2,q))) / Iprofile(q) * 10^3;
    Maxflexprofile(2, q) = BME(q) * (Yprofile(2, q)) / Iprofile(q) * 10^3;

%     sigma_top = max_moment * (y_bar-dist_ground_to_top) / I * 10^3;
% sigma_bottom = max_moment * y_bar / I * 10^3;
end
Maxflexprofile;

%Shear
sametimeglue = 1; %how many glue points are simultaneously present
Shearprofile = zeros(sametimeglue, n+1);
%manually input distance to bottom of glue points and their range
g1 = 75;
g2 = 1.27;
g3 = 250;
tabtotal = 2*5;
findybar(Asabove(cross_section_0, Yprofile(2,1)))


% Q = Ad for all glues
Q1 = Aabove(cross_section_0, g1) * abs((Yprofile(2,1) - findybar(Asabove(cross_section_0, g1))));

Qcenareas = Asabove(cross_section_0, Yprofile(2,1));
Qcen1 = (dot(Qcenareas(:,1),Qcenareas(:,2)) * (findybar(Qcenareas)- Yprofile(2,1)));

%add to profile
Shearprofile(1, 1:end) = (abs(SFE(1:end)) * Q1) ./ (Iprofile(1,1) * 2*5);
% Shearprofile(1, 901:n+1) = (abs(SFE(901:n+1)) * Q1) ./ (Iprofile(1,1) * 2*6.27);
% Shearprofile(2, 901:n+1) = (abs(SFE(901:n+1)) * Q2) ./ (Iprofile(1,1) * 2*6.27);
% Shearprofile(2, 1:300) = (abs(SFE(1:300)) * Q2) ./ (Iprofile(1,1) * 2*6.27);
% Shearprofile(1, 301:900) = (abs(SFE(301:900)) * Q3) ./ (Iprofile(1,500) * 2*6.27);

Shearcenprofile(1:n+1) = (abs(SFE(1:n+1)) * Qcen1) ./ (Iprofile(1,1) * 2*1.27);
% Shearcenprofile(301:900) = (abs(SFE(301:900)) * Qcen2) ./ (Iprofile(1,500) * 2*1.27);




%thin plate buckling
fish = 0.2;
E = 4000 %N/mm
Diaphragmlos = [1,401,801,1201];

%top middle flange (case 1)
top_middle_flange = zeros(1, n+1);
thicknesstop = 1.27;
basetopmiddle = 80;
top_middle_flange(1:n+1) = ( ...
    (4*pi^2*E)/(12*(1-fish^2)))*(thicknesstop/basetopmiddle)^2;

%top sides (case 2)
top_side_flange = zeros(1,n+1);
baseside = 10;
top_side_flange(1:n+1) = ( ...
    (0.4254*pi^2*E)/(12*(1-fish^2)))*(thicknesstop/baseside)^2; 


%sides (case 3)
Sides = zeros(1, n+1);
Sides(1:n+1) = ((6*pi^2*E)/(12*(1-fish^2)))*(1.27./(Yprofile(1,q)-Yprofile(2,q)-1.27)).^2;

thinplate = [top_middle_flange;
    top_side_flange;
    Sides]

%shear buckling
b_web = 1.27;
h_web = 75-1.27;
shearbuckling = zeros(1,n+1);
for i = 1:length(Diaphragmlos)-1
    shearbuckling(Diaphragmlos(i):Diaphragmlos(i+1)) = ((5*pi^2*E)/(12*(1-fish^2)))*((b_web/h_web)^2+(b_web/(Diaphragmlos(i+1)-Diaphragmlos(i))^2));

end

% Graphing SFE
figure
subplot(1,2,1)
hold on; grid on; grid minor;
plot(1:n+1, SFE, "g-")
xlabel('Length (mm)')
ylabel('Max Shear Force (N)')
title("Shear Force Envelope")
% 
subplot(1,2,2)
hold on; grid on; grid minor;
plot(1:n+1,BME,"b-")
xlabel('Length (mm)')
ylabel('Moment (Nm)')
title("Bending Moment Envelope")

%graphing failrue flexture
sixs(1:1400) = -6;
figure
hold on; grid on; grid minor;
set(gca, 'YDir', 'reverse')
plot(1:n+1,Maxflexprofile,"black-")
plot(1:n+1,-top_middle_flange,"red-")
plot(sixs, "magenta-")
plot(zeros(1,1400), 'black-')
xlabel('Length (mm)')
ylabel('Flextural Stress (Mpa)')
title("Flextural Stress")
legend('Bridge Flexture Compression (Top)', 'Bridge FLexture Tension (Bottom)','Case 1 Buckling', 'Compression Failure')
hold off

g_fail(1:n+1) = 2;
figure
subplot(1,2,1)
hold on; grid on; grid minor;
plot(1:n+1,Shearprofile,"black-")
plot(g_fail, "blue-")
xlabel('Length (mm)')
ylabel('Shear Stress at glue joints(Mpa)')
title("Glue points shear")
legend('shear at glue', 'failure point')
hold off

shrfail(1:1400) = 4;
subplot(1,2,2)
hold on; grid on; grid minor;
plot(1:n+1,shearbuckling,"magenta-")
plot(1:n+1,Shearcenprofile, "black-")
plot(shrfail, 'red-')
xlabel('Length (mm)')
ylabel('Shear Stress at Centroid(Mpa)')
title("Shear Centroid")
legend( 'failure shear buckling of web', 'shear at centroid', 'failure Matboard shear at centroid')
hold off
% 



%final Display
figure
subplot(1,2,1)
hold on; grid on; grid minor;
plot(1:n+1,SFDlist(1:domain,1:end),"r-")
xlabel('Length (mm)')
ylabel('Shear Force (N)')
subplot(1,2,2)
hold on; grid on; grid minor;
set(gca, 'YDir', 'reverse')
plot(1:n+1,BMDlist(1:domain,1:end),"b-")
xlabel('Length (mm)')
ylabel('Moment (Nm)')
hold off


%Failure Table
Tensile = 30;
Compression = 6;
Shearmatboard = 4;
Shearglue = 2;

%display everything
FosComp = Compression/abs(min(Maxflexprofile(1,1:end)));
FosTen = Tensile/max(abs(Maxflexprofile(2,1:end)));
FosGlueShear = Shearglue/max(Shearprofile, [], "all");
FosMatboardShear = Shearmatboard/max(Shearcenprofile, [], "all");
FosFlangeBucklingTopMiddle = min(thinplate(1), [], "all")/abs(min(Maxflexprofile(1,1:end)));
FosFlangeBucklingSides = min(thinplate(3), [], "all")/abs(min(Maxflexprofile(1,1:end)));
FosFlangeBucklingTopside = min(thinplate(2), [], "all")/abs(min(Maxflexprofile(1,1:end)));
FosShearBuckling = min(shearbuckling, [], "all")/max(Shearcenprofile, [], "all");
Fails = [FosComp, FosTen, FosGlueShear, FosMatboardShear, FosFlangeBucklingTopside, FosFlangeBucklingTopMiddle, FosFlangeBucklingSides, FosShearBuckling];

FosDiagram = [Tensile./abs(Maxflexprofile(2,1:end));
              Compression./abs(Maxflexprofile(1,1:end));
              Shearmatboard./Shearcenprofile
              thinplate(1)./abs(Maxflexprofile(1,1:end));
              thinplate(2)./abs(Maxflexprofile(1,1:end));
              shearbuckling./Shearcenprofile];
%remove infinities
FosDiagram(1:end,1) = 1000000;
FosDiagram(1:end,end) = 1000000;

Names = ["Tensile", "Compression", "Thinplate middle top", "Thinplate sides/webs", "Shear matboard", "Shear Buckling", "Glue Shearing", "Thinplate side flange"]

FailsImproved = [min(FosDiagram(1,:)),min(FosDiagram(2,:)),min(FosDiagram(3,:)),FosGlueShear, min(FosDiagram(4,:)), min(FosDiagram(5,:)), FosFlangeBucklingSides, min(FosDiagram(6,:))]

fprintf(['Factors of Safety\n' ...
         '  Matboard Tensile: %.2f\n' ...
         '  Matboard Compressive: %.2f\n' ...
         '  Matboard Shear: %.2f\n' ...
         '  Glue Shear: %.2f\n' ...
         '  Buckling Case 1: %.2f\n' ...
         '  Buckling Case 2: %.2f\n' ...
         '  Buckling Case 3: %.2f\n' ...
         '  Shear Buckling: %.2f\n'], FailsImproved)

function ybar = findybar(cross_section)
    num_shapes = length(cross_section);
    cross_section = cross_section';
    sum_num = 0;
    sum_den = 0;
    for i = 1:num_shapes
        sum_num = sum_num + cross_section(1,i)*cross_section(2,i)*cross_section(3,i);
        sum_den = sum_den + cross_section(1,i)*cross_section(2,i);
    end
    
    ybar = sum_num / sum_den;
end



function area = Aabove(cross_section, limit)
    area = 0;
    for i = 1:length(cross_section)
        if cross_section(i, 3) > limit
            area = area + cross_section(i, 1) * cross_section(i, 2);
        end
    end
end

function areas = Asabove(cross_section, limit)
    areas = zeros(length(cross_section), 3);
    for i = 1:length(cross_section)
        if cross_section(i, 3) > limit
            areas(i, 1:3) = cross_section(i, 1:3);
        
        else 
            if cross_section(i,3) + cross_section(i,2)*0.5 > limit
            areas(i,1) = cross_section(i,1);
            areas(i,2) = cross_section(i,2) - (limit - (cross_section(i,3)- cross_section(i,2)*0.5));
            areas(i,3) = areas(i,2)/2 + limit;
            end
        end

    end
end

