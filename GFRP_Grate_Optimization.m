%% 4-Inch 3 Hole Optimization
clear all; close all; clc;

colors = ["r","g","b","c","m","y","k"];
colorcount = 1;

% Base variables (CHANGE HERE)
r = 0.5;
Y = 0.25;
skel_width = 0.1;

x_offset = @(r,skel_width,override)(((2*r)+ skel_width)*cosd(60))*override;
y_offset = @(r,skel_width,override)(2*r+skel_width)*sind(60)*override;

end_double_cir = @(r,skel_width)((2*r)+ skel_width)*cosd(60) + r;

% Circular Functions
y_center = @(t,r) sqrt(r^2 - t.^2);
y_top = @(t,r,skel_width,override_x,override_y) sqrt(r^2 - (t - x_offset(r,skel_width,override_x)).^2) + y_offset(r,skel_width,override_y);
y_bottom = @(t,r,skel_width,override_x,override_y,vert_flip) vert_flip*(sqrt(r^2 - (t - x_offset(r,skel_width,override_x)).^2) - y_offset(r,skel_width,override_y));

width = 3;

hold on

if width == 4
    ylim([0 3])
    yline(2)
    r_end = 2;
    y_border = 2;

elseif width == 6
    ylim([0 4])
    yline(3)
    r_end = 3;
    y_border = 3;
elseif width == 3
    ylim([0 2.5])
    yline(1.5)
    r_end = 1.5;
    y_border = 1.5;

end



base_data = [];
data_count = 1;
r_step = 1/16;

for d = 3/16:r_step:r_end
    for skel_width = 0.1:0.01:4
        t = 0:0.01:end_double_cir(d/2,skel_width)+0.01;

        y_c = [];
        y_t = [];
        y_b = [];
        y_tot = [];
        
       
        for i = 1:length(t)
            y_c(i) = real(y_center(t(i),d/2));
            y_t(i) = real(y_top(t(i),d/2,skel_width,1,0));
            y_b(i) = real(y_bottom(t(i),d/2,skel_width,1,0,1));
        
            y_tot(i) = y_c(i) + y_t(i) + y_b(i);
        end
        [value, index] = max(y_tot);
        
        %for j = 1:length(t)
        %    y_c_1(j) = real(y_center(t(j),r));
        %    y_t_1(j) = real(y_top(t(j),r,skel_width,1,1));
        %    y_b_1(j) = real(y_bottom(t(j),r,skel_width,1,1,-1));
        %end
        y_c_1 = real(y_center(t,d/2));
        y_t_1 = real(y_top(t,d/2,skel_width,1,1));
        y_b_1 = real(y_bottom(t,d/2,skel_width,1,1,-1));

        if y_border - max(y_t_1) < Y
           %xline(t(index),colors(colorcount))
            %disp(strcat("skel: "),string(skel_width))
           if skel_width ~= 0.1
                %disp(strcat("r: ",string(d)))
                %disp(strcat("skel: ",string(skel_width)))
                disp(strcat("SUCCESS: DIAMETER ", string(d), " SKEL: ",string(skel_width)))
                plot(t,y_c_1,colors(colorcount), 'DisplayName',strcat(string(d)," in"))
                plot(t,y_t_1,colors(colorcount), 'HandleVisibility','off')
                plot(t,y_b_1,colors(colorcount), 'HandleVisibility','off')
                title("Maximum Skeleton Width per. Diameter", strcat("(Tread Width = ", string(width),")"))
                xlabel("Distance (in)")
                ylabel("Distance (in)")
                %legend(d,colors(colorcount))
                %xline(t(index),colors(colorcount))
                

                base_data(1,data_count) = tl;
                base_data(2,data_count) = r_top;
                base_data(3,data_count) = ml;
                base_data(4,data_count) = r_bottom;
                base_data(5,data_count) = value;
                base_data(6,data_count) = d;
                base_data(7,data_count) = t(index);
                base_data(8,data_count) = skel_width;

                data_count = data_count + 1;
           else
                disp(strcat("ERROR: DIAMETER ", string(d), " SKEL: ",string(skel_width)))
           end
           break
        end

        tl = width/2 - y_t_1(index);
        ml = y_b_1(index) - y_c_1(index);
        
        if abs(width/2 - tl - ml - value) > 0.000001
            delta = abs(width/2 - tl - ml - value);
            disp(strcat("ERROR: TOTAL WIDTH ~= 0. DELTA = ",string(delta)))
        end
        
        r_top = (width/2 - y_t_1(index))/2 + y_t_1(index);
        r_bottom = (y_b_1(index) - y_c_1(index))/2 + y_c_1(index);

    end
    if colorcount ~= 7
        colorcount = colorcount + 1;
    else
        colorcount = 1;
    end
end
legend
legend('NumColumns',2)

%% 4-Inch 2 Hole Optimization
clear all; close all; clc;

colors = ["r","g","b","c","m","y","k"];
colorcount = 1;

% Base variables (CHANGE HERE)
r = 0.5;
Y = 0.25;
skel_width = 0.1;

x_offset = @(r,skel_width,override)(((2*r)+ skel_width)*cosd(60))*override;
y_offset = @(r,skel_width,override)(2*r+skel_width)*sind(60)*override;

end_double_cir = @(r,skel_width)((2*r)+ skel_width)*cosd(60) + r;

% Circular Functions
y_far = @(t,r) sqrt(r^2 - (t).^2) + 1;
y_close = @(t,r) sqrt(r^2 - (t).^2) - 1;
t = -3:0.01:3;


r_step = 1/8;
hold on
ylim([-2 2])
xlim([ -2 2])
yline(2)
yline(2 - 0.25)
for d = 0.5:r_step:2
    y_t_1 = real(y_far(t,d/2));
    y_b_1 = real(y_close(t,d/2));
    %plot(t,y_t_1)
    if 2 - max(y_t_1) < Y
        disp(strcat("ERROR: ",string(d)))
    else
        plot(t,y_t_1,colors(colorcount))
        plot(t,-1*y_t_1,colors(colorcount))
        plot(t,y_b_1,colors(colorcount))
        plot(t,-1*y_b_1,colors(colorcount))
        disp(strcat("SUCCESS: ",string(d)))
    end

    if colorcount ~= 7
        colorcount = colorcount + 1;
    else
        colorcount = 1;
    end
end

%% post Processing
clc
t_base = 1:length(base_data);

figure
plot(t_base,base_data)
legend("tl","r_top","ml","r_bottom","value")

for i = 1:length(base_data)
    delta = abs(2 - base_data(1,i) - base_data(3,i) - base_data(5,i));
    if abs(delta) > 0.001
        disp(strcat("ERROR POST PROCESSING, DELTA: ",string(delta)))
    end
end
%% I CALCULATION 10/30
clc;

P = 500*width;
L = 25;
D = 0.75;
E = 4*10^6;
S = 45*10^3;

I_stiff = @(top_base, mid_base,h) 2*(top_base * h^3)/12 + 2*(mid_base * h^3)/12;
I_strength = @(top_base, mid_base,h) 2*(top_base * h^2)/12 + 2*(mid_base * h^2)/12;

stiff_requirement = (P*L^3)/(48*E*D);
strength_requirement = (P*L)/(8*S);

for data_counter = 1:length(base_data)
    for h = 0:0.01:20
        I_stiff_calc = I_stiff(base_data(1,data_counter), base_data(3,data_counter),h);
        I_strength_calc = I_strength(base_data(1,data_counter), base_data(3,data_counter),h);

    
        if I_stiff_calc > stiff_requirement && I_strength_calc > strength_requirement
            disp(strcat("SUCCESSFUL HEIGHT FOUND AT ",string(h), " FOR DIAMETER: ",string(base_data(6,data_counter)), " FOR I VALUE: ",string(I_stiff_calc)));
            base_data(9,data_counter) = h;
            disp(string(I_strength_calc))
            disp(string(I_stiff_calc))
            break
        elseif I_stiff_calc > stiff_requirement
            if h == 20
                %disp(strcat("FAILURE TO FIND HEIGHT FOR: ",string(base_data(6,data_counter), " - STRENGTH")))
            end
        elseif I_strength_calc > strength_requirement
            if h == 20
                %disp(strcat("FAILURE TO FIND HEIGHT FOR: ",string(base_data(6,data_counter), " - STIFFNESS")))
            end
        else
            if h == 20
                %disp(strcat("DOUBLE FAILURE TO FIND HEIGHT FOR: ",string(base_data(6,data_counter)), " - STIFFNESS + STRENGTH"))
            else
                %disp(strcat("DOUBLE FAILURE TO FIND HEIGHT FOR: ",string(base_data(6,data_counter)), " - STIFFNESS + STRENGTH ",string(h)))
            end        
        end
    end
end

%%
clc

for data_counter = 1:length(base_data)
   
    xoff = x_offset(base_data(6,data_counter)/2,base_data(8,data_counter),1);
    edge_length = 0.25 + base_data(6,data_counter)/2;
    set_length = xoff*2;

    lower_count = floor((25 - edge_length)/set_length);
    outside_rows = 2*lower_count;

    full_dist = 0.25 + base_data(6,data_counter)/2 + lower_count*set_length;
    
    if full_dist + base_data(6,data_counter)/2 > 25
        full_dist = full_dist - set_length;
        %disp(strcat("Last lower circle is cutoff for: ",string(base_data(6,data_counter))))
        if full_dist + xoff + base_data(6,data_counter)/2 > 25
            %disp(strcat("cannot fit additional center at: ", string(base_data(6,data_counter))))
            middle_count = lower_count - 1;
        else
            %disp("CENTER IS SAME AS LOWER")
            middle_count = lower_count;
        end
    else
        lower_count = lower_count + 1;
        %disp(strcat("Last lower cirlce is NOT cutoff for: ",string(base_data(6,data_counter))))
        if full_dist + xoff + base_data(6,data_counter)/2 > 25
            %disp("CENTER IS ONE LESS THAN LOWER")
            middle_count = lower_count - 1;
        else
            %disp("CENTER IS SAME AS LOWER")
            middle_count = lower_count;
        end
    end
    disp(strcat(string(base_data(6,data_counter)), " MIDDLE ROW: ",string(middle_count), " LOWER ROW: ",string(lower_count)))
    base_data(10,data_counter) = 2*lower_count + middle_count;

    solid_volume = 25*width*base_data(9,data_counter);
    circle_volume = (2*lower_count + middle_count) * (base_data(9,data_counter) * pi * (base_data(6,data_counter)/2)^2) ;

    total_volume = solid_volume - circle_volume;
    base_data(10,data_counter) = total_volume;

    mass_total = total_volume * 0.065;
    base_data(11,data_counter) = mass_total;

    penalty_total = (mass_total - width*2) * 0.5;
    if penalty_total > 0
        cost_tread = mass_total * 7.2 + penalty_total;
    else
        cost_tread = mass_total * 7.2;
    end
    
    cost_total = cost_tread * 1/width * 12 * 60;
    base_data(12,data_counter) = cost_total;
end

figure
plot(base_data(6,:),base_data(12,:))
title(strcat("Diameter vs. Bridge Cost ", "(",string(width), '" Width)'))
xlabel("Diameter (in)")
ylabel("Total Cost ($)")