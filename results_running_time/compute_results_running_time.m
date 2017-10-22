total_time_U_SR = zeros(2,1);
total_time_S_U = zeros(4,1);
total_time_L_4 = zeros(4,1);
total_time_L_6 = zeros(4,1);
total_time_L_8 = zeros(4,1);

for i = 1:10
    load(['G1_' num2str(i) '.mat']);
    total_time_U_SR(1) = total_time_U_SR(1) + time_U_SR(end);
    total_time_S_U(1) = total_time_S_U(1) + time_S_U(end);
    total_time_L_4(1) = total_time_L_4(1) + time_L_4(end);
    total_time_L_6(1) = total_time_L_6(1) + time_L_6(end);
    total_time_L_8(1) = total_time_L_8(1) + time_L_8(end);
    
    load(['G2_' num2str(i) '.mat']);
    total_time_U_SR(2) = total_time_U_SR(2) + time_U_SR(end);
    total_time_S_U(2) = total_time_S_U(2) + time_S_U(end);
    total_time_L_4(2) = total_time_L_4(2) + time_L_4(end);
    total_time_L_6(2) = total_time_L_6(2) + time_L_6(end);
    total_time_L_8(2) = total_time_L_8(2) + time_L_8(end);
    
    load(['G3_' num2str(i) '.mat']);
    total_time_S_U(3) = total_time_S_U(3) + time_S_U(end);
    total_time_L_4(3) = total_time_L_4(3) + time_L_4(end);
    total_time_L_6(3) = total_time_L_6(3) + time_L_6(end);
    total_time_L_8(3) = total_time_L_8(3) + time_L_8(end);
    
    load(['G4_' num2str(i) '.mat']);
    total_time_S_U(4) = total_time_S_U(4) + time_S_U(end);
    total_time_L_4(4) = total_time_L_4(4) + time_L_4(end);
    total_time_L_6(4) = total_time_L_6(4) + time_L_6(end);
    total_time_L_8(4) = total_time_L_8(4) + time_L_8(end);
end

total_time_U_SR = total_time_U_SR / 10;
total_time_S_U = total_time_S_U / 10;
total_time_L_4 = total_time_L_4 / 10;
total_time_L_6 = total_time_L_6 / 10;
total_time_L_8 = total_time_L_8 / 10;