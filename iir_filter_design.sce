// use functions defined to get all results

// first, common parameters for both filters
M = 37; // my filter design number
delta = 0.15; // tolerance in passband and stopband
B_tran = 4e3; // transition bandwidth in Hz
//...................................................
// parameters specific to the two filters
// Filter 1: A Bandpass filter
B_signal_analog_1 = 160e3; // bw of analog signal, inconsequential
F_sampling_1 = 330e3; // sampling frequency in Hz for filter 1
filter_type_1 = 'bpf';
filter_nature_1 = 'bu';
//...................................................

// Filter 2: Band Stop filter

B_signal_analog_2 = 120e3; //analog signal BW in Hz
F_sampling_2 = 260e3; // Sampling frequency in Hz
filter_type_2 = 'bsf';
filter_nature_2 = 'ch';
//..................................................


// Normalized filter specifications

criticalf_1 = un_norm_filter_edges(M, B_tran, filter_type_1);
criticalf_2 = un_norm_filter_edges(M, B_tran, filter_type_2);

//disp("Un-normalized Filter_1 (BPF) Specifications [fs1, fp1, fp2, fs2]:");
//disp(criticalf_1);
//disp("Un-normalized Filter_2 (BSF) Specifications [fp1, fs1, fs2, fp2]:");
//disp(criticalf_2);


criticalw_1 = get_critical_w(M, filter_type_1, B_tran, F_sampling_1);
criticalw_2 = get_critical_w(M, filter_type_2, B_tran, F_sampling_2);

//disp("Normalized Filter_1 (BPF) Specifications [ws1, wp1, wp2, ws2]:");
//disp(criticalw_1);
//disp("Normalized Filter_2 (BSF) Specifications [wp1, ws1, ws2, wp2]:");
//disp(criticalw_2);

//....................................................

// Bilinear transform to convert to analog frequency W

criticalW_1 = bilinear_transform_wtoW(criticalw_1);
criticalW_2 = bilinear_transform_wtoW(criticalw_2);

//disp("Analog Filter_1 (BPF) Specifications [Ws1, Wp1, Wp2, Ws2]");
//disp(criticalW_1);
//disp("Analog Filter_2 (BSF) Specifications [Wp1, Ws1, Ws2, Wp2]");
//disp(criticalW_2);

 //....................................................

//Frequency transform to convert to LPF

temp_W_lpf_1 = freq_trans_init(criticalW_1, filter_type_1);
temp_W_lpf_2 = freq_trans_init(criticalW_2, filter_type_2); // temp since need to choose to get stopband, passband

//disp("Frequency transformed LPF values Filter_1 ");
//disp(temp_W_lpf_1);
//disp("Frequency transformed LPF values Filter_2");
//disp(temp_W_lpf_2);

criticalW_lpf_1 = analog_freq_trans(criticalW_1, filter_type_1);
criticalW_lpf_2 = analog_freq_trans(criticalW_2, filter_type_2);

//disp("Analog LPF specifications Filter_1:");
//disp(criticalW_lpf_1);
//disp("Analog LPF specifications Filter_2:");
//disp(criticalW_lpf_2);

//.........................................................
// get useful parameters 

parameters_1=lpf_parameters(delta, criticalW_lpf_1, filter_nature_1); //[N, Wc]
parameters_2=lpf_parameters(delta, criticalW_lpf_2, filter_nature_2); //[N, eps, Wp]

//disp("Analog LPF parameters Filter 1 (Butterworth) [N, Wc]");
//disp(parameters_1);
//disp("Analog LPF parameters Filter 2 (Chebyshev) [N, epsilon, Wp]");
//disp(parameters_2);

//............................................................

// get poles

poles_1 = find_poles(filter_nature_1, parameters_1);
poles_2 = find_poles(filter_nature_2, parameters_2);

//disp("The poles of H_analog_LPF for Filter_1");
//disp(poles_1);
//disp("The poles of H_analog_LPF for Filter_2");
//disp(poles_2);

// view poles
plot(real(poles_2), imag(poles_2), '*');
W_p =parameters_2(3);
B_k = (1/parameters_2(1))*asinh(1/parameters_2(2))
mag = W_p*sinh(B_k);
x = -mag:0.0001: mag;
plot(x, W_p*cosh(B_k)*sqrt(1-(x./(W_p*sinh(B_k))).^2));
plot(x, -W_p*cosh(B_k)*sqrt(1-(x./(W_p*sinh(B_k))).^2));
isoview on;
a = gca();
a.limits()
xlabel("Real(sL)");
ylabel("Imag(sL)");
title("Poles of H_analog_LPF (sL)*H_analog_LPF(-sL)");

//............................................................

//get left poles

left_poles_1 = find_left_poles(poles_1);
left_poles_2 = find_left_poles(poles_2);

//disp("The left poles of H_analog_LPF for Filter_1");
//disp(left_poles_1);
//disp("The left poles of H_analog_LPF for Filter_2");
//disp(left_poles_2);

//... parameters for plotting
params_plot_2 = [criticalW_lpf_2(4), criticalW_lpf_2(3) ,delta]
params_plot_1 = [criticalW_lpf_1(4), criticalW_lpf_1(3) ,delta]
//.............................................................

//get H_analog_LPF

H_analog_LPF_1 = H_lpf_from_poles(left_poles_1, filter_nature_1, parameters_1);
H_analog_LPF_2 = H_lpf_from_poles(left_poles_2, filter_nature_2, parameters_2);

disp("H_analog_LPF(sL) Filter 1:")
disp(H_analog_LPF_1);
disp("H_analog_LPF(sL) Filter 2:")
disp(H_analog_LPF_2);
freq_axis = -5:0.001:5;
//lin_sys_1 = syslin('c',H_analog_LPF_1);
//bode(lin_sys_1, freq_axis, "rad");
name = "Analog LPF Filter 1 (Butterworth) ";
//plot_H_analog_LPF(H_analog_LPF_1, params_plot_1, name);
name_2 = "Analog LPF Filter 2 (Chebyshev) "
//plot_H_analog_LPF(H_analog_LPF_2, params_plot_2, name_2);

//..............................................................

//get H_analog(s)

H_analog_1 = analog_filter_from_lpf(filter_type_1, H_analog_LPF_1, criticalW_1);
H_analog_2 = analog_filter_from_lpf(filter_type_2, H_analog_LPF_2, criticalW_2);

//disp("H_analog(s) Filter_1");
//disp(H_analog_1);
//disp("H_analog(s) Filter_2");
//disp(H_analog_2);
name_analog_1 = "Filter 1: Analog BPF Butterworth";
//plot_H_analog(H_analog_1, criticalW_1, delta, name_analog_1);
name_analog_2 = "Filter 2: Analog BSF Chebyshev";
//plot_H_analog(H_analog_2, criticalW_2, delta, name_analog_2);


//..............................................................

//get H(z)

H_z_1 = discrete_filter_from_analog(H_analog_1);
H_z_2 = discrete_filter_from_analog(H_analog_2);
//disp(H_z_1)
//temp = syslin('d',H_z_1);
//trfmod(temp);



// plot magnitude and phase response
name_z_1 = "Filter 1: Discrete-Time Butterworth BPF";
//plot_H_z(H_z_1, criticalw_1, delta, name_z_1);
name_z_2 = "Filter 2: Discrete-Time Chebyshev BSF";
//plot_H_z(H_z_2, criticalw_2, delta, name_z_2);
//disp("H(z) Filter_1");
//disp(H_z_1);
//disp("H(z) Filter_2");
//disp(H_z_2);
w_axis_new = 0.001:0.001:3;
H_bpf_val = horner(H_z_1, exp(%i*w_axis_new));
H_bsf_val = horner(H_z_2, exp(%i*w_axis_new));

mag_bpf_val = abs(H_bpf_val);
phase_bpf_val = atan(imag(H_bpf_val), real(H_bpf_val)); // tan-1(y/x)
//h1 = gca();
//plot(w_axis_new, mag_bpf_val, 'r');
////legend(["|H(w)|"]);
//xlabel("w (normalized frequency)");
//ylabel("|H(w)|", "color",'r');
//h2 = newaxes();
//plot(w_axis_new, phase_bpf_val);
//h2.filled="off";
//h2.y_location="right";
//ylabel("argH(w)", "color",'b')
//set(gca(),"auto_clear","off");
////legends(["|H(w)|";"arg(H(w))"]);
//title("IIR Filter-1 (BPF- Butterworth) Frequency Response");

mag_bsf_val = abs(H_bsf_val);
phase_bsf_val = atan(imag(H_bsf_val), real(H_bsf_val)); // tan-1(y/x)
//h1 = gca();
//plot(w_axis_new, mag_bsf_val, 'r');
////legend(["|H(w)|"]);
//xlabel("w (normalized frequency)");
//ylabel("|H(w)|", "color",'r');
//h2 = newaxes();
//plot(w_axis_new, phase_bsf_val);
//h2.filled="off";
//h2.y_location="right";
//ylabel("argH(w)", "color",'b')
//set(gca(),"auto_clear","off");
////legends(["|H(w)|";"arg(H(w))"]);
//title("IIR Filter-2 (BSF- Chebyshev) Frequency Response");


//................................................................

//plot relevant figures


//plot_H_analog_LPF(H_analog_LPF_2, params_plot_2);

//plot_H_z(H_z_2, criticalw_2,delta);
//plot_H_analog(H_analog_2, criticalW_2, delta)

