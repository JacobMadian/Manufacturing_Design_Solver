%% Question 5) Part 1

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

box_outer = 3;
box_inner = @(t,box_outer) box_outer - 2*t;

height_outer = @(h) h;
height_inner = @(t,height_outer) height_outer - 2*t;


square_solve = [];

P = 2000;
L = 25;
D = 0.75;

for material = 1:length(materials)

    E = materials(material,1);
    S = materials(material,2);
    
    stiff_requirement = (P*L^3)/(48*D*E*10^6);
    strength_requirement = (P*L)/(8*S*10^3);

    for thickness = 0:0.001:10
        
    b_outer = box_outer;
    b_inner = box_inner(thickness,b_outer);

    h_outer = b_outer;
    h_inner = b_inner;

    I_stiff_calc = (1/12)*((b_outer * h_outer^3) - (b_inner * h_inner^3));
    I_strength_calc = (1/12)*((b_outer * h_outer^2) - (b_inner * h_inner^2));
    
    if I_stiff_calc > stiff_requirement && I_strength_calc > strength_requirement
        
        disp(strcat(string(names(material))," thickness requirement met at: ", string(thickness)))
        square_solve(1,material) = thickness;

        Volume_Outer = 25 * b_outer * h_outer;
        Volume_Inner = 25 * b_inner * h_inner;
        Volume_Net = Volume_Outer - Volume_Inner;
        square_solve(2,material) = Volume_Net;

        mass_net = Volume_Net * materials(material,3);
        square_solve(3,material) = mass_net;
        
        cost_per_tread = mass_net * materials(material,4);
        penalty_total = (mass_net - b_outer*2) * 0.5;

        if penalty_total > 0
            cost_post_penalty = cost_per_tread + penalty_total;
        else
            cost_post_penalty = cost_per_tread;
        end
        square_solve(4,material) = cost_post_penalty;
        break
    end

    end
end

%% Question 5) Part 2 & 3
clc; close all;

rectangle_solve = [];
processing_charge = [1 1.5 2 3 4];

cost_matrix = [];

for thickness = [0.05 0.1]
    disp(thickness)
    for material = 1:length(materials)
    
        E = materials(material,1);
        S = materials(material,2);
        
        stiff_requirement = (P*L^3)/(48*D*E*10^6);
        strength_requirement = (P*L)/(8*S*10^3);
    
        for height = 0:0.001:10
            
        b_outer = box_outer;
        b_inner = box_inner(thickness,b_outer);
    
        h_outer = height;
        h_inner = height_inner(thickness, h_outer);
    
        I_stiff_calc = (1/12)*((b_outer * h_outer^3) - (b_inner * h_inner^3));
        I_strength_calc = (1/12)*((b_outer * h_outer^2) - (b_inner * h_inner^2));
        
        if I_stiff_calc > stiff_requirement && I_strength_calc > strength_requirement
            
            disp(strcat(string(names(material))," height requirement met at: ", string(height)))
            rectangle_solve(1,material) = height;
    
            Volume_Outer = 25 * b_outer * h_outer;
            Volume_Inner = 25 * b_inner * h_inner;
            Volume_Net = Volume_Outer - Volume_Inner;
            rectangle_solve(2,material) = Volume_Net;
    
            mass_net = Volume_Net * materials(material,3);
            rectangle_solve(3,material) = mass_net;
            
            for cost_adjustment = 1:length(processing_charge)
                cost_per_tread = mass_net * materials(material,4) * processing_charge(cost_adjustment);
                penalty_total = (mass_net - b_outer*2) * 0.5;
        
                if penalty_total > 0
                    cost_post_penalty = cost_per_tread + penalty_total;
                else
                    cost_post_penalty = cost_per_tread;
                end

                if thickness == 0.05
                    cost_matrix(material,cost_adjustment) = cost_post_penalty;
                elseif thickness == 0.1
                    cost_matrix(material + 5,cost_adjustment) = cost_post_penalty;
                end
                rectangle_solve(4,material) = cost_post_penalty;
            end
            break
        end
    
        end
    end
end

index = 1:5;
hold on

colors = ["r","k","b","m","c"];
colorcount = 1;

for material = 1:length(materials)
    if material == 4
    else
        plot(processing_charge(index), cost_matrix(material + 5,index),'color',colors(colorcount),'marker','o')
        plot(processing_charge(index), cost_matrix(material,index),'color',colors(colorcount),'marker','*')
        colorcount = colorcount + 1;
    end
end

xlabel("Processing Charge (cost * scale factor)")
ylabel("Span Cost ($)")

legend

%% Question 5) Part 4
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

box_inner = @(t,box_outer) box_outer - 2*t;

height_outer = 2;
height_inner = @(t,height_outer) height_outer - 2*t;


rectange_h_solve = [];

P = 2000;
L = 25;
D = 0.75;

processing_charge = [1 1.5 2 3 4];

for box_width = [3 4 6]
    disp(strcat("For Box Width: ",string(box_width)))
    for material = 1:length(materials)
    
        E = materials(material,1);
        S = materials(material,2);
        
        stiff_requirement = (P*L^3)/(48*D*E*10^6);
        strength_requirement = (P*L)/(8*S*10^3);
    
        for thickness = 0:0.001:10
    
        h_outer = height_outer;
        h_inner = height_inner(thickness,h_outer);
    
        b_outer = box_width;
        b_inner = box_inner(thickness,b_outer);
    
        I_stiff_calc = (1/12)*((b_outer * h_outer^3) - (b_inner * h_inner^3));
        I_strength_calc = (1/12)*((b_outer * h_outer^2) - (b_inner * h_inner^2));
        
        if I_stiff_calc > stiff_requirement && I_strength_calc > strength_requirement
            
            disp(strcat(string(names(material))," thickness requirement met at: ", string(thickness)))
            rectangle_h_solve(1,material) = thickness;
    
            Volume_Outer = 25 * b_outer * h_outer;
            Volume_Inner = 25 * b_inner * h_inner;
            Volume_Net = Volume_Outer - Volume_Inner;
            rectangle_h_solve(2,material) = Volume_Net;
    
            mass_net = Volume_Net * materials(material,3);
            rectangle_h_solve(3,material) = mass_net;
            
            for cost_adjustment = 1:length(processing_charge)
                cost_per_tread = mass_net * materials(material,4) * processing_charge(cost_adjustment);
                penalty_total = (mass_net - b_outer*2) * 0.5;
        
                if penalty_total > 0
                    cost_post_penalty = cost_per_tread + penalty_total;
                else
                    cost_post_penalty = cost_per_tread;
                end
    
                if box_width == 3
                    cost_matrix(material,cost_adjustment) = cost_post_penalty;
                elseif box_width == 4
                    cost_matrix(material + 5,cost_adjustment) = cost_post_penalty;
                elseif box_width == 6
                    cost_matrix(material + 10,cost_adjustment) = cost_post_penalty;
                end
                rectangle_h_solve(4,material) = cost_post_penalty;
            end
            break
        end
    
        end
    end
end
