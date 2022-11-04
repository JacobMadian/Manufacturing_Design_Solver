%% ME 501 Design Project - Z. Bacon, C. Neher, S. Park, J. Madian
%% Material Data
clear all; close all; clc;

% Given Materials from Design Assignment
steel = [28 75 0.29 2.00];
aluminum = [10 40 0.1 4.50];
GFRP = [4 45 0.065 7.00];

%%%%%%%%% To Do: Find and add two Additional Materials %%%%%%%%
titanium = [16 20.3 0.163 17];
cast_iron_g1800 = [21.5 62 0.261 3.6];

% Form arrays for material data access in below sections and names for
% later reference
materials = [steel; aluminum; GFRP; titanium; cast_iron_g1800];
names = ["steel"; "aluminum"; "GFRP";"titanium";"cast iron"];

%% Height Calculation Per Material
clc;

% Request whether user wants to customize beam shape
x = input("Would you like to customize tread shape? enter 'Y' for yes, 'N' for no " + newline);
if x == 'N'
    % Typical Parameters for Size and Force Requirements
    base = 3;
    P = base*500;
    P_1 = base*500;
    L = 25;

    D = 0.75;

% Inputs for custom beam shape
elseif x == 'Y'
    base = input("Input base length of tread (units: inches)"+ newline);
    P = input("Input load per tread (units: lb/inch)"+ newline);
    L = input("Input length of tread (units: inches)"+ newline);
    D = input("Input deflection limit of tread (units: inches)"+ newline);
end

% rows are used for iterating for loops
[rows,cols] = size(materials);

% Empty array to hold maximum height per given material and geometry
height_per_mat = [];

for n = 1 : rows
    %height_1 = sqrt(18.75/materials(n,2));
    %height_2 = nthroot(((2.604)/materials(n,1)),3);

    % Height 1 - dependent on stress
    % Height 2 - dependent on deflection
    %height_1 = sqrt((3*P*l)/(2*base*materials(n,2)*1000));
    %height_2 = nthroot((P_1*l^3)/(4*base*deflection_limit*materials(n,1)*1000000),3);
    E = materials(n,1);
    S = materials(n,2);
    
    stiff_requirement = (P*L^3)/(48*D*E*10^6);
    strength_requirement = (P*L)/(8*S*10^3);
    
    for height = 0:0.001:10
    % Max height of two equations used
        I_stiff_Calc = (base * height ^3)/12;
        I_stren_Calc = (base * height ^2)/12;
    
        if I_stiff_Calc > stiff_requirement && I_stren_Calc > strength_requirement
            disp(strcat("Successfull height found for ",string(names(n))," at height ",string(height)," in"))
            break
        end
    end

    height_per_mat(n) = height;
    disp("Height per Material (max h)")
    disp(strcat(names(n), " : ", string(height_per_mat(n))," in"))
end

%% Checking Materials Against Strength/Stiffness Limitations
clc;
moment_of_inertia = @(b,h) (1/12)*b*h^3;
c = @(h) h/2;
stress = @(M,c,I) ((M*c)/I)*10^-6;
stiffness = @(P,L,E,I) ((P*L^3)/(48*E*I))*10^-6;

stress_results = [];
stiff_results = [];

for n = 1:rows
    I_calc = moment_of_inertia(base, height_per_mat(n));
    c_calc = c(height_per_mat(n));
    stress_results(length(stress_results) + 1) = stress(((P*l^3)/4),c_calc,I_calc);
    stiff_results(length(stiff_results) + 1) = stiffness(P,l,materials(n,1),I_calc);
    disp(strcat(string(stress_results(length(stress_results))), " stress"))
    disp(strcat(string(stiff_results(length(stiff_results))), " stiff"))

    if stress_results(length(stress_results)) >= materials(n,2)
        disp(strcat(names(n)," failed stress limit with: ", string(stress_results(end))))
    end

    if ge(stiff_results(length(stiff_results)),0.75)
        disp(strcat(names(n)," failed deflection limit with: ",string(stiff_results(end))))
    end
end


%% Weight Calculation Per Material (based on max h)
clc;

% Empty array to hold the calculated weight per material and geometry
weight_per_mat = [];
volume_per_mat = [];

for n = 1: rows
    % Weight calcualted by H,L,W (aka base), and Density
    volume = height_per_mat(n) * L * base;
    volume_per_mat(n) = volume;
    weight = height_per_mat(n) * L * base * materials(n,3);
    weight_per_mat(n) = weight;

    disp("Weight per Material (max h)")
    disp(strcat(names(n), " : ", string(weight_per_mat(n)), " lb"))
end
%% Volume Add for Part 4) Grooves. SKIP SECTION IF DOING TASK 1
clc;

% Empty array to hold the calculated weight per material and geometry
weight_per_mat = [];
volume_per_mat = [];

for n = 1: rows
    % Weight calcualted by H,L,W (aka base), and Density
    volume = height_per_mat(n) * l * base;
    Volume_with_Grooves = (1/8) * (1/8) * 25 * base + volume;
    volume_per_mat(n) = Volume_with_Grooves;

    weight = Volume_with_Grooves * materials(n,3);
    weight_per_mat(n) = weight;

    disp("Weight per Material (max h)")
    disp(strcat(names(n), " : ", string(weight_per_mat(n)), " lb"))
end

%% Cost Calculation Per Material (based on max h)
clc;

% Empty arrays for the given cost per material, and the total cost of the
% entire bridge
cost_per_mat = [];
total_cost = [];

% Empty complete matrix for holding pre-penalty, post-penalty, and TC
cost_matrix = zeros(3,3);

for n = 1: rows

    % Cost calculated by lb * cost/lb
    cost = weight_per_mat(n) * materials(n,4);
    cost_per_mat(n) = cost;
    cost_matrix(n,1) = cost;

    %disp(' ');disp("Cost (pre penalty) per Material & Tread")
    %disp(strcat(names(n), " : $", string(cost_per_mat(n)), "/Tread"))

    %disp(' ');disp("Cost (post penalty) per Material")

    % Penalty calculated per given equation in assignment, honestly unclear
    penalty = (weight_per_mat(n) - (base * 2)) * 0.5;
    if penalty >= 0
        cost_matrix(n,2) = cost_per_mat(n) + penalty;
        cost_per_mat(n) = cost_per_mat(n) + penalty;
    else
        cost_matrix(n,2) = cost_per_mat(n);
    end
    %disp(strcat(names(n), " : $", string(cost_per_mat(n)), "/Tread"))

    % Total bridge cost calculated
    %disp(' ');disp("Cost per Span")
    total_cost(n) = cost_per_mat(n) * 1/base * 12 * 60;
    cost_matrix(n,3) = cost_per_mat(n) * 1/base * 12 * 60;
    %disp(strcat(names(n), " : $", string(total_cost(n))))
    
    disp(strcat(names(n), " pre-penalty: ", string(cost_matrix(n,1)* 1/base * 12 * 60)))
    disp(strcat(names(n), " post-penalty: ", string(cost_matrix(n,2)* 1/base * 12 * 60)))


    
end

[cheapest_cost, cheapest_index] = min(total_cost);
[highest_cost, highest_index] = max(total_cost);

disp(strcat("Cheapest: ", names(cheapest_index),"   Cost: ", string(cheapest_cost)))
disp(strcat("Most Expensive: ", names(highest_index),"   Cost: ", string(highest_cost)))




