// this code computes the order and the half-power point of Butterworth filter
//
//lpf_specs = [-1.4161988,-1.,1.,1.4161988]; //bpf
//lpf_specs_bsf = [-1.2394869,-1.,1.,1.2394869]; //bsf
//delta = 0.15;
//W_p = lpf_specs(3);
//W_s = lpf_specs(4);
//
//D1 = (1/(1-delta)^2)-1;
//D2 = (1/(delta^2))-1;
//numerator = log(D2/D1);
//denominator = log(W_s/W_p);
//N = ceil(0.5*numerator/denominator);
//
//disp("LPF specifications")
//disp(lpf_specs)
//disp("delta")
//disp(delta);
//disp("D1:");
//disp(D1);
//disp("D2");
//disp(D2);
//disp("N");
//disp(N);
//
//low_Wc = W_p/(D1^(1/(2*N)));  // lower bound on the half power point Wc
//high_Wc = W_s/(D2^(1/(2*N)));  // lower bound on the half power point Wc 
//
//mprintf("lower bound Wc");
//disp(low_Wc)
//mprintf("higher bound on Wc is \n %f \n",high_Wc);
//// choose Wc to be the mean of the upper and lower bounds
//
//Wc =0.5*(low_Wc + high_Wc);
//mprintf("The value of Wc is \n %f", Wc);
//



function [parameters]=lpf_parameters(tolerance, criticalW_LPF, filter_nature)
    W_p_temp = criticalW_LPF(3);
    W_s_temp = criticalW_LPF(4);
    D1 = (1/((1-tolerance)^2))-1;
    D2 = (1/(tolerance^2))-1;
    d_ratio = sqrt(D2/D1);
    W_ratio = W_s_temp/W_p_temp;
    if filter_nature=='bu' then
        numerator_N = log(d_ratio);
        denominator_N = log(W_ratio);
//        mprintf("The value of N before ceil is: %f", numerator_N/denominator_N);
        N = ceil(numerator_N/denominator_N);
        low_Wc = W_p_temp/(D1^(1/(2*N)));  // lower bound on the half power point Wc
        high_Wc = W_s_temp/(D2^(1/(2*N)));  // lower bound on the half power point Wc 
        Wc =0.5*(low_Wc + high_Wc);
        parameters = [N, Wc];
        
    elseif filter_nature=='ch' then
        assert_checkalmostequal(W_s_temp, -criticalW_LPF(1));
        assert_checkalmostequal(W_p_temp, -criticalW_LPF(2));
        
        eps = sqrt(D1);
        numerator_N = acosh(d_ratio);
        denominator_N = acosh(W_ratio);
        N= ceil(numerator_N/denominator_N);
//        mprintf("The value of N before ceil is: %f", numerator_N/denominator_N);
        parameters = [N, eps, W_p_temp]; //Wp needed to find poles, 
        //here W_p=1, but useful in general
        
    end
endfunction


// test the function
//params = lpf_parameters(delta, lpf_specs, 'bpf')
//
//disp("Function output:");
//disp(params)
//// test bsf
//params = lpf_parameters(delta, lpf_specs_bsf, 'bsf');
//disp("Function out BSF");
//disp(params);

//mprintf("For Chebyshev BSF, N=%d \n",params(1));
//mprintf("For Chebyshev BSF, epsilon=%f \n",params(2));
//




// hence the low pass filter is obtained!: N, Wc are obtained
// now, need to find the poles, and the expression for transfer function
