// code to design the FIR filters

// ideal_BPF = ideal_LPF1 - ideal_LPF_2 with appropriate cutoff frequencies
// get the impulse response of ideal_BPF and then use window (rectangular)

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
filter_nature_2 = 'bu';
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



//......................................................

A = -20*log10(delta);
//disp("The value of A")
//disp(A)
del_wt_1 = criticalw_1(2)-criticalw_1(1);
del_wt_2 = criticalw_2(2)-criticalw_2(1);

//disp("Transition BW Filter 1:")
//disp(del_wt_1)
//disp("Transition BW Filter 2:")
//disp(del_wt_2)


l_min_1 = 1 + (A-8)/(2.285*del_wt_1); //=49.7
l_min_2 = 1 + (A-8)/(2.285*del_wt_2);//=39.4
addi_1 = 15;
addi_2 = 13;
len_trial_1 = ceil(l_min_1)+addi_1; // =50+addi_1 =65 for satisfying criteria
len_trial_2 = ceil(l_min_2)+addi_2; // =40+addi_2 = 53 for satisfying criteria

N_1 = (len_trial_1-1)/2; // 32
N_2 = (len_trial_2-1)/2; //26

//disp("Minimum length for Filter 1:")
//disp(l_min_1)
//disp("Minimum length for Filter 2:")
//disp(l_min_2)


n_max = 1000; // max time points (on each side of 0, ie num points =2n+1) to be working with, even for ideal filters

n_axis = -n_max:n_max;
w_axis = -3.14:0.001:3.14;
// FILTER1: BPF

w_cutoff_low_filter1 = (criticalw_1(2) + criticalw_1(1))/2 ; 
w_cutoff_high_filter1 = (criticalw_1(3) + criticalw_1(4))/2;
// middle of transition band

w_cutoff_low_filter2 = (criticalw_2(2) + criticalw_2(1))/2;
w_cutoff_high_filter2 = (criticalw_2(3) + criticalw_2(4))/2;

// obtain ideal LPF response for LPF with wp=wp2, wp1 and subtract

h_ideal_wp2_filter1 = impulse_response_ideal_lpf(w_cutoff_high_filter1, n_max);
h_ideal_wp1_filter1 = impulse_response_ideal_lpf(w_cutoff_low_filter1, n_max);
h_ideal_bpf_filter1 = h_ideal_wp2_filter1 - h_ideal_wp1_filter1;

h_ideal_wp2_filter2 = impulse_response_ideal_lpf(w_cutoff_high_filter2, n_max);
h_ideal_wp1_filter2 = impulse_response_ideal_lpf(w_cutoff_low_filter2, n_max);
h_ideal_delta_filter2 = impulse_response_ideal_lpf(%pi, n_max);
h_ideal_bsf_filter2 = h_ideal_delta_filter2-(h_ideal_wp2_filter2 - h_ideal_wp1_filter2);
// constant - bpf = bsf

//plot(n_axis, h_ideal_bpf_filter1);
//xlabel("n");
//ylabel("h[n] for ideal BPF")
//title("Ideal BPF impulse response Filter 1")

H_bpf_ideal = freq_transform(n_axis, w_axis, h_ideal_bpf_filter1);
//plot(w_axis, abs(H_bpf_ideal));
H_bsf_ideal = freq_transform(n_axis, w_axis, h_ideal_bsf_filter2);
//plot(w_axis, abs(H_bsf_ideal));

//............................................
// obtain the windowed Finite impulse response

h_bpf_windowed = apply_window_h_ideal(h_ideal_bpf_filter1, N_1);
//plot(n_axis(n_max-80:n_max+80),h_bpf_windowed(n_max-80:n_max+80));
//xlabel("n");
//ylabel("h_windowed[n]");
//title("Filter 1 (BPF) Impulse response after windowing ");
//disp("Impulse response after applying window Filter1: BPF");
nonzero_indices_h_bpf_win = find(h_bpf_windowed);
samples_h_bpf = h_bpf_windowed(nonzero_indices_h_bpf_win);
//disp(samples_h_bpf);
n_axis_reduced = -N_1:N_1;
z = poly(0,'z');
assert_checkequal(length(n_axis_reduced), length(samples_h_bpf));
H_z_bpf = sum( samples_h_bpf.*(z^(-n_axis_reduced)));
//disp("H(z) for Filter 1 BPF :");
//disp(H_z_bpf);
name_1 = "FIR Bandpass Filter: Kaiser-Window ";
//plot_H_z(H_z_bpf, criticalw_1, delta,name_1);

h_bsf_windowed = apply_window_h_ideal(h_ideal_bsf_filter2, N_2);
//plot(n_axis(n_max-80:n_max+80),h_bsf_windowed(n_max-80:n_max+80));
//xlabel("n");
//ylabel("h_windowed[n]");
//title("Filter 2 (BSF) Impulse response after windowing ");
//title("Impulse response - windowed Filter 2 BSF");
//disp("Impulse response after applying window Filter2: BSF");
nonzero_indices_h_bsf_win = find(h_bsf_windowed);
samples_h_bsf = h_bsf_windowed(nonzero_indices_h_bsf_win);
//disp(samples_h_bsf);
n_axis_reduced_bsf = -N_2:N_2;
//z = poly(0,'z');
assert_checkequal(length(n_axis_reduced_bsf), length(samples_h_bsf));
H_z_bsf = sum( samples_h_bsf.*(z^(-n_axis_reduced_bsf)));
name_2 = "FIR Bandstop Filter: Kaiser Window ";
//plot_H_z(H_z_bsf, criticalw_2, delta, name_2);

//disp("H(z) for Filter 2 BSF :");
//disp(H_z_bsf);

// look at freq spectrum of windowed function
//H_bpf_windowed = freq_transform(n_axis, w_axis, h_bpf_windowed);
//plot(w_axis, abs(H_bpf_windowed));
//
//set(gca(),"auto_clear","off");
//plot(w_axis,(1+delta)*ones(1,length(w_axis)),':'); // horizontal line at 1
//plot(w_axis,(1-delta)*ones(1,length(w_axis)),':'); // horizontal line at 1
//plot(w_axis,(delta)*ones(1,length(w_axis)),':'); // horizontal line at 1
//plot(criticalw_1(1)*ones(1, length(w_axis)), abs(H_bpf_windowed),':');
//plot(criticalw_1(2)*ones(1, length(w_axis)), abs(H_bpf_windowed),':');
//plot(criticalw_1(3)*ones(1, length(w_axis)), abs(H_bpf_windowed),':');
//plot(criticalw_1(4)*ones(1, length(w_axis)), abs(H_bpf_windowed),':');
w_axis_new = 0.001:0.001:3.14;
H_bpf_windowed = freq_transform(n_axis, w_axis_new, h_bpf_windowed);
//bode(w_axis_new, H_bpf_windowed);
mag_bpf = abs(H_bpf_windowed);
phase_bpf = abs(atan(imag(H_bpf_windowed), real(H_bpf_windowed))); // tan-1(y/x)
//h1 = gca();
//plot(w_axis_new, mag_bpf, 'r');
////legend(["|H(w)|"]);
//xlabel("w (normalized frequency)");
//ylabel("|H(w)|", "color",'r');
//h2 = newaxes();
//plot(w_axis_new, phase_bpf);
//h2.filled="off";
//h2.y_location="right";
//ylabel("argH(w)", "color",'b')
//set(gca(),"auto_clear","off");
////legends(["|H(w)|";"arg(H(w))"]);
//title("FIR Filter-1 (BPF) Frequency Response")

H_bsf_windowed = freq_transform(n_axis, w_axis_new, h_bsf_windowed);

mag_bsf = abs(H_bsf_windowed);
phase_bsf = abs(atan(imag(H_bsf_windowed), real(H_bsf_windowed))); // tan-1(y/x)
h1 = gca();
plot(w_axis_new, mag_bsf, 'r');
//legend(["|H(w)|"]);
xlabel("w (normalized frequency)");
ylabel("|H(w)|", "color",'r');
h2 = newaxes();
plot(w_axis_new, phase_bsf);
h2.filled="off";
h2.y_location="right";
ylabel("argH(w)", "color",'b')
set(gca(),"auto_clear","off");
//legends(["|H(w)|";"arg(H(w))"]);
title("FIR Filter-2 (BSF) Frequency Response");



//plot(w_axis, abs(H_bsf_windowed));
//
//set(gca(),"auto_clear","off");
//plot(w_axis,(1+delta)*ones(1,length(w_axis)),':'); // horizontal line at 1
//plot(w_axis,(1-delta)*ones(1,length(w_axis)),':'); // horizontal line at 1
//plot(w_axis,(delta)*ones(1,length(w_axis)),':'); // horizontal line at 1
//plot(criticalw_2(1)*ones(1, length(w_axis)), abs(H_bsf_windowed),':');
//plot(criticalw_2(2)*ones(1, length(w_axis)), abs(H_bsf_windowed),':');
//plot(criticalw_2(3)*ones(1, length(w_axis)), abs(H_bsf_windowed),':');
//plot(criticalw_2(4)*ones(1, length(w_axis)), abs(H_bsf_windowed),':');



