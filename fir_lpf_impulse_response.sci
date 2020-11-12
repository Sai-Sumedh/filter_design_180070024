// this file gives the impulse response of an ideal LPF


function [h_n]=impulse_response_ideal_lpf(w_p, n_max) // passpand till w_p
    n = -n_max:n_max;
    h_n = (w_p/%pi)*sinc(w_p*n);
endfunction


//n_max = 500;
//wp2 = 1.767;
//wp1 = 1.386;
//h1 = impulse_response_ideal_lpf(wp2,n_max);
//n_axis = -n_max:n_max;
//h2 = impulse_response_ideal_lpf(wp1, n_max);
//del = zeros(-n_max:n_max);
//del(n_max+1)=1;
//h_new = del-h2+h1;
//plot(n_axis,h_new)
// this verifies that the linear combination of ideal impulse responses is max at 0
//plot(n_axis,h);
//xlabel("n");
//ylabel("h[n]]");
//title("Ideal LPF impulse response with wp=1.767");


function [w_n]=generate_rect_window(l_win_side, centre, n_max) // length=odd
    // centre is the value (like real n in h[n] (can be -ve too), not the array index])
    temp = zeros(-n_max:n_max);
    temp(n_max-l_win_side+1+centre:n_max+l_win_side+1+centre)=1;
    w_n = temp;
endfunction
//
//w = generate_rect_window(5, 10);
//disp(w)
//plot(-10:10, w)

function [h_finite]=apply_window_h_ideal(h_ideal_n, window_len) // window the function , centered at the max position of h_ideal OR should it be zero? // window_len = N (one side length)
    [max_val, centre_index] = max(h_ideal_n);
    n_max= floor((length(h_ideal_n)-1)/2);
    centre = centre_index-n_max-1;
    w_n = generate_rect_window(window_len, centre, n_max);
    disp(length(h_ideal_n))
    disp(length(w_n))
    h_finite = h_ideal_n.*w_n;
    
endfunction

//n_max = 10000;
//wp2 = 1.767;
//wp1 = 1.386;
//n_axis = -n_max:n_max;
//h2 = impulse_response_ideal_lpf(wp2,n_max);
//h_wind = apply_window_h_ideal(h2, 10);
////plot(n_axis, h_wind, n_axis, h2)
////plot(n_axis, h1)
////H2 = fft(h2);
////plot(abs(H2));
//w_axis = -3.14:0.01:3.14
function [H_w_val]=freq_transform(n_axis, w_axis, h_n) // n_axis, w_axis and h_n are row arrays
    product = w_axis'*n_axis; 
    exp_term = exp(-%i*product);
    H_w_val = exp_term*h_n'; 
    H_w_val = H_w_val'; // now a row vector with size = size of w_axis
endfunction

//H = freq_transform(n_axis, w_axis, h2)
//plot(w_axis, abs(H))
//plot(n_axis, h2)
// its working beautifully!






