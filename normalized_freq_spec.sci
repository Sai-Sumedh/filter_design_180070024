// this code converts the given un-normalized discrete-time specifications into normalized discrete time specifications
//....................................................
// defining parameters

// the following parameters are common for both filters
//M = 37; // my filter design number
//delta = 0.15; // tolerance in passband and stopband
//B_tran = 4e3; // transition bandwidth in Hz
//
////....................................................
//
//// Filter 1: A Bandpass filter
//B_signal_analog_1 = 160e3; // bw of analog signal, inconsequential
//F_sampling_1 = 330e3; // sampling frequency in Hz for filter 1
//filter_type_1 = 'bpf';
////...................................................
//
//// Filter 2: Band Stop filter
//
//B_signal_analog_2 = 120e3; //analog signal BW in Hz
//F_sampling_2 = 260e3; // Sampling frequency in Hz
//filter_type_2 = 'bsf';
////...................................................

// FUNCTIONS

function [critical_f]=un_norm_filter_edges(filter_num_M, bw_transition, filter_type)
    // critical_f = [fs1, fp1, fp2, fs2]
    // filter_type = 'bpf' or 'bsf' only
    q_m = floor((filter_num_M-1)/10); // weirdo functions to get a passband
    r_m = filter_num_M - 10*q_m;
    
    if filter_type=='bpf' then
        f_low = 25 + 1.7*(q_m)+ 6.1*r_m; // passband edges in kHz
        f_high = f_low + 20; // in kHz
    elseif filter_type=='bsf' then
        f_low = 25 + 1.9*(q_m)+ 4.1*r_m; // stopband edges in kHz
        f_high = f_low + 20; // in kHz
    end
    
    B_mid = 1e3*[f_low, f_high]; // passband(for bpf)/ stopband (for bsf) edges in Hz
    f_init = B_mid(1)-bw_transition; // pb for bpf, sb for bsf
    f_final = B_mid(2) + bw_transition;
    critical_f = [f_init, B_mid, f_final];
    
endfunction

function [out]=normalize(critical_freqs, F_sampling)
    pi = %pi;
    critical_ws = (2*pi/F_sampling)*critical_freqs;
    out = critical_ws;
endfunction

// THE FUNCTION TO BE USED
function [critical_w]=get_critical_w(filter_num_M, filter_type, bw_transition,F_sampling)
    f_crit = un_norm_filter_edges(filter_num_M, bw_transition, filter_type);
    critical_w = normalize(f_crit, F_sampling);
endfunction
//.............................................................
//
//critical_freqs_1 = un_norm_filter_edges(M, B_tran, filter_type_1); // un-normalized band edges
//critical_freqs_2 = un_norm_filter_edges(M, B_tran, filter_type_2);
//disp("Critical Freq 1");
//disp(critical_freqs_1);
//disp("Critical Freq 2");
//disp(critical_freqs_2);
//
////crit_w_1 = normalize(critical_freqs_1, F_sampling_1);
////disp("Normalized frequency specifications Filter 1:")
////disp(crit_w_1)
////
////crit_w_2 = normalize(critical_freqs_2, F_sampling_2);
////disp("Normalized frequency specifications Filter 2:")
////disp(crit_w_2)
//
//criticalw1 = get_critical_w(M, 'bpf', B_tran, F_sampling_1);
//criticalw2 = get_critical_w(M, 'bsf', B_tran, F_sampling_2);
//
//disp("Normalized w specifications Filter 1:")
//disp(criticalw1)
//disp("Normalized w specifications Filter 2:")
//disp(criticalw2)
//
//
