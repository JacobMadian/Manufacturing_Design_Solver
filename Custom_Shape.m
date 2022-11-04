clear all; close all; clc;

% Given Materials from Design Assignment
steel = [28 75 0.29 2.00];
aluminum = [10 40 0.1 4.50];
GFRP = [4 45 0.065 7.00];

%%%%%%%%% To Do: Find and add two Additional Materials %%%%%%%%
titanium = [17 35 0.163 110];
cast_iron_g1800 = [26 14 0.272 3.6];

% Form arrays for material data access in below sections and names for
% later reference
materials = [steel; aluminum; GFRP; titanium; cast_iron_g1800];
names = ["steel"; "aluminum"; "GFRP";"titanium";"cast iron"];

width = 4;
base_cross_section = 0.25;

P = 500*width;
P_1 = 500*width;
L = 25;
D = 0.75;

[rows,cols] = size(materials);

% Empty array to hold maximum height per given material and geometry
data_matrix = [];

for n = 1 : rows
    E = materials(n,1);
    S = materials(n,2);

    stiff_requirement = (P*L^3)/(48*D*E*10^6);
    strength_requirement = (P*L)/(8*S*10^3);

    for height = 0:0.001:10
    % Max height of two equations used
        I_stiff_Calc = ((base_cross_section) * height ^3)/3;
        I_stren_Calc = ((base_cross_section) * height ^2)/3;
    
        if I_stiff_Calc > stiff_requirement && I_stren_Calc > strength_requirement
            disp(strcat("Successfull height found for ",string(names(n))," at height ",string(height)," in"))
            break
        end
    end
    data_matrix(1,n) = height;

    Total_Volume = 25 * width * height;
    Net_Volume = Total_Volume - (width*(23.83 * (2/3) * height + pi * height * (1/3)^2));
    data_matrix(2,n) = Net_Volume;

    Weight = Net_Volume * materials(n,3);
    data_matrix(3,n) = Weight;

    penalty_total = (Weight - width*2) * 0.5;
    disp(penalty_total)
    if penalty_total > 0
        cost_tread = Weight * materials(n,4) + penalty_total;
    else
        cost_tread = Weight * materials(n,4);
    end
    
    cost_total = cost_tread * 1/width * 12 * 60;
    data_matrix(4,n) = cost_total;
    disp(strcat("Total Cost for ",string(names(n))," is ",string(cost_total)))
end

