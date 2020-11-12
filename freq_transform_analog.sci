// this file has a function to convert the analog frequency specifications of a Band Pass filter to the corresponding analog Low-pass filter specifications using the analog frequency transform

//criticalW = [0.7680835,0.8304841,1.2182079,1.3174514]; // bpf
//criticalW = [0.7912795,0.8730651,2.0728589,2.3576122]; // bsf
//// store these values in the analog BPF frequencies
//Ws1 = criticalW(1);
//Wp1 = criticalW(2);
//Wp2 = criticalW(3);
//Ws2 = criticalW(4);
//W0 = sqrt(Wp1*Wp2);
//B = Wp2-Wp1;

function [W_transformed]=freq_trans_init(criticalW, filter_type) // applies the frequency transformation formula
    // filter_type='bpf' or 'bsf'
    
    if filter_type=='bpf' then
        W_transformed = (criticalW.^2 - (criticalW(2)*criticalW(3)))./ ((criticalW(3)-criticalW(2))*criticalW); // for bpf, criticalW = [Ws1, Wp1, Wp2, Ws2]
        
    elseif filter_type=='bsf' then
        temp = ((criticalW(4)-criticalW(1))*criticalW)./((criticalW(1)*criticalW(4))-criticalW.^2 ); //  for bsf, criticalW=[Wp1, Ws1, Ws2, Wp2]
        // now rearrange temp into [Ws1, Wp1, Wp2, Ws2], with Wp1=-Wp2, but Ws1, Ws2 differing (will be taken care of in another function)
        W_transformed = [temp(3), temp(4), temp(1), temp(2)];
    end
    
endfunction

//W_trans_LPF = freq_trans_init(criticalW, 'bsf');
//disp('Analog BSF specs:');
//disp(criticalW);
//disp('Analog LPF specs (after frequency transformation');
//disp(W_trans_LPF);

// now, choose the stricter condition for stopband edge
//W_Ls = min(abs(W_trans_LPF(1)), abs(W_trans_LPF(4)));
//disp("The stricter stopband edge is:");
//disp(W_Ls); 

function [criticalW_LPF]=choose_strict_Ws(W_LPF) // choose the strict edge for stopband
    
    W_Ls = min(abs(W_LPF(1)), abs(W_LPF(4)));
//    disp("The stricter stopband edge is:");
//    disp(W_Ls);
    criticalW_LPF = [-W_Ls, W_LPF(2), W_LPF(3), W_Ls];
    
endfunction


// THE MAIN FUNCTION
function [criticalW_LPF]=analog_freq_trans(criticalW, filter_type) // obtain analog lpf specifications from analog bpf 
    
    W_tr_LPF = freq_trans_init(criticalW, filter_type);
    criticalW_LPF = choose_strict_Ws(W_tr_LPF)
    
endfunction

//criticalW_LPF = analog_freq_trans(criticalW, 'bsf')
//disp("Final LPF specs:");
//disp(criticalW_LPF);

//ans = analog_freq_trans(criticalW);
//disp('Input W (analog BPF) is:');
//disp(criticalW);
//disp('LPF analong W is:');
//disp(ans);
//
//disp("W0");
//disp(W0);

//disp("B");
//disp(B);
