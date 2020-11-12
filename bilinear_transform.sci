// take as input the normalized discrete time specs
// give out the analog frequency specifications



function [analog_specs]=bilinear_transform_wtoW(critical_w)
    analog_specs = tan(0.5*critical_w);
endfunction

// for testing the function
//criticalw1 =[1.3099489, 1.3861088, 1.7669079, 1.8430677];
//criticalw2 = [1.3388018,  1.4354662,  2.2426138,  2.3392782];
//criticalW1 = bilinear_transform_wtoW(criticalw1);
//criticalW2 = bilinear_transform_wtoW(criticalw2);
//disp("Analog W for filter 1");
//disp(criticalW1);
//disp("Analog W for filter 2");
//disp(criticalW2);

//step = 0.01;
//w = 0:step:5;
//filter_specs = zeros(1,floor(criticalW(2)/step));
//L_pass = floor(criticalW(3)/step)-floor(criticalW(2)/step);
//filter_specs = [filter_specs, ones(1, L_pass)];
//filter_specs = [filter_specs, zeros(1, length(w)-length(filter_specs))]; // subtracting old length
//plot(w, filter_specs); // the ideal low pass filter in the analog domain

//plot(criticalw, criticalw);

//test = rand(200,1);
//test_i = 3*gsort(test, 'g', 'i')
//analog = transform_ztos(test_i);
//plot(test_i, analog, '*')

