// this code finds the poles of the magnitude squared transfer function

//criticalw =[1.3099489, 1.3861088, 1.7669079, 1.8430677]; // normalized discrete time specifications
//lpf_specs = [-1.4161988,-1.,1.,1.4161988];
//delta = 0.15;
//W_p = lpf_specs(3);// for BPF
//W_s = lpf_specs(4);
//
//
//N = 7; // order found from previous computation
//k = 0:(2*N-1); // indexing set for poles
//Wc = 1.076247;
//W0 = 1.0058341; // taken from previous file
//B = 0.3877238;
//criticalW = [0.7680835,0.8304841,1.2182079,1.3174514]; // analog BPF critical W
//
////disp(k)
//j = %i;
//pi = %pi;
//delta = 0.15;
//s_k = Wc* exp((j*pi/(2*N))*(N+1+2*k));
//mprintf("The poles of H(s)*H(-s) are: \n ");
//disp(s_k)

//..................................................................
//FUNCTIONS

function [s_k]=find_poles(filter_type, parameters) 
    // filter_type = 'butterworth' or 'chebyshev'
    N = parameters(1);
    k = 0:(2*N-1); // indexing set for poles
    pi = %pi;
    j = %i;
    if filter_type == 'bu' then
        Wc = parameters(2);
        s_k = Wc* exp((j*pi/(2*N))*(N+1+2*k));
        
    elseif filter_type == 'ch' then
        eps = parameters(2);
        W_p = parameters(3);
        B = (1/N)*asinh(1/eps);
        A_k = (2*k+1)*pi/(2*N);
        s_k = W_p*(j*cosh(B)*cos(A_k) - sinh(B)*sin(A_k));
        
    end
endfunction
//
//params_bpf = [7.,1.0762474];
//poles_fn = find_poles('butterworth', params_bpf);
//disp("poles from function");
//disp(poles_fn);
//disp(sum(abs(s_k-poles_fn)));
////..................................................
//params_bsf=[5.,0.6197443,1.];
//poles_fn_bsf = find_poles('chebyshev', params_bsf);
////disp("poles from function for bsf");
//disp(poles_fn_bsf);
//eps = params_bsf(2);
//W_p = params_bsf(3);
//B = (1/params_bsf(1))*asinh(1/eps);
//plot(real(poles_fn_bsf), imag(poles_fn_bsf), '*');
//isoview on;
//set(gca(),"auto_clear","off");
//mag = W_p*sinh(B);
//x = -mag:0.0001: mag;
//plot(x, W_p*cosh(B)*sqrt(1-(x./(W_p*sinh(B))).^2));
//plot(x, -W_p*cosh(B)*sqrt(1-(x./(W_p*sinh(B))).^2));
//////..................................................................
//
////s_left_plane = real(s_k)<0
////plot(real(s_k), imag(s_k), '*');
////set(gca(),"auto_clear","off")
////


function [s_k_left]=find_left_poles(s_k)
    s_left_plane = real(s_k)<0; // binary array
    left_pole_indices = find(s_left_plane);
    s_k_left = s_k(left_pole_indices); // these are the poles of H_analog_LPF(s)
endfunction
//
//left_poles_bpf = find_left_poles(poles_fn);
////disp("Left poles, from func")
////disp(left_poles_fn);
//left_poles_bsf = find_left_poles(poles_fn_bsf);
//disp("Left poles bsf");
//disp(left_poles_bsf);


////mag = abs(real(s_k(4)));
////x = -mag:0.0001: mag;
////plot(x, sqrt(mag.^2-(x.^2)));
////plot(x, -sqrt(mag.^2-(x.^2)));
////isoview on;
////title("Poles of H_analog_LPF (s)*H_analog_LPF(-s)");
////xlabel("Real(s)");
////ylabel("Imag(s)");
//
//
//s_left_plane = real(s_k)<0; // boolean axis
//left_pole_indices = find(s_left_plane);
//left_poles = s_k(left_pole_indices); // these are the poles of H_analog_LPF(s)
//mprintf("The poles on the left half plane are:");
//disp(left_poles);

function [H_analog_LPF]=H_lpf_from_poles(left_poles, filter_nature, parameters)
    denominator_H = poly(left_poles, 'sL', "roots"); // a function
    if filter_nature=='bu' then
        numerator_H = horner(denominator_H, 0);//  H=numerator/denominator = 1 at sL=0
    elseif filter_nature=='ch' then
        eps = parameters(2);
        numerator_H = (sqrt(1/(1+ eps^2)))*horner(denominator_H, 0);//  H=numerator/denominator = 1 at sL=0
    end
    
    H_analog_LPF = numerator_H/denominator_H;
    
endfunction
// test the function for BSF
//H_lpf_bsf = H_lpf_from_poles(left_poles_bsf);
//disp("H_LPF for BSF");
//disp(H_lpf_bsf);
//// test the function for BPF
//H_lpf_bpf = H_lpf_from_poles(left_poles_bpf);
//disp("H_LPF for BPF");
//disp(H_lpf_bpf);
//.................................................................

// FUNCTION to get H_analog(s) from H_analog_LPF(sL)

function [H_analog]=analog_filter_from_lpf(filter_type, H_analog_LPF, criticalW) // filter_type = 'bpf' or 'lpf'
    s = poly(0,'s');
    if filter_type=='bpf' then // criticalW = [Ws1, Wp1, Wp2, Ws2]
        W_p1 = criticalW(2);
        W_p2 = criticalW(3);
        W0 = sqrt(W_p1*W_p2);
        B = W_p2 - W_p1;
        s_L_temp = (s^2 + (W0^2))/(B*s);
        H_analog = horner(H_analog_LPF, s_L_temp);
    elseif filter_type == 'bsf' then // criticalW = [Wp1, Ws1, Ws2, Wp2]
        W_p1 = criticalW(1);
        W_p2 = criticalW(4);
        W0 = sqrt(W_p1*W_p2);
        B = W_p2 - W_p1;
        s_L_temp = (B*s)/((s^2) + (W0^2));
        H_analog = horner(H_analog_LPF, s_L_temp); 
    end
        
endfunction
//
//criticalW_bpf = [0.7680835,0.8304841,1.2182079,1.3174514]; // bpf
//criticalW_bsf = [0.7912795,0.8-
//
//H_analog_bpf = analog_filter_from_lpf('bpf', H_lpf_bpf, criticalW_bpf)
//disp("H_analog_bpf")
//disp(H_analog_bpf)
//
//H_analog_bsf = analog_filter_from_lpf('bsf', H_lpf_bsf, criticalW_bsf)
//disp("H_analog_bsf")
//disp(H_analog_bsf)

function plot_H_analog(H_analog, criticalW, delta,name)
    
    W_axis = 0:0.001:2.5;
    j = %i;
    H_values_analogBPF = horner(H_analog, j*W_axis);
    
    plot(W_axis,ones(1,length(W_axis)),'b--'); // horizontal line at 1
    plot(W_axis,(1-delta)*ones(1,length(W_axis)),'b--'); // horizontal line at 1
    plot(W_axis,(delta)*ones(1,length(W_axis)),'b--'); // horizontal line at 1
    plot(criticalW(1)*ones(1, length(W_axis)), abs(H_values_analogBPF),'b--');
    plot(criticalW(2)*ones(1, length(W_axis)), abs(H_values_analogBPF),'b--');
    plot(criticalW(3)*ones(1, length(W_axis)), abs(H_values_analogBPF),'b--');
    plot(criticalW(4)*ones(1, length(W_axis)), abs(H_values_analogBPF),'b--');
    plot(W_axis, abs(H_values_analogBPF), 'r');
    e=gce();p=e.children(1);
    loc_1 = floor((criticalW(1))/0.001);
    loc_2 = floor((criticalW(2))/0.001);
    loc_3 = floor(criticalW(3)/0.001);
    loc_4 = floor(criticalW(4)/0.001);
     
    set(gca(),"auto_clear","off");
    title(name);
    xlabel("W")
    ylabel("H_analog(W)")
    d1=datatipCreate(p,loc_1);
    d2 = datatipCreate(p,loc_2 );
    d3 = datatipCreate(p, loc_3);
    d4 = datatipCreate(p, loc_4);
    datatipSetOrientation(d1,"lower left");
    datatipSetOrientation(d2,"upper left");
    datatipSetOrientation(d3,"upper right");
    datatipSetOrientation(d4,"lower left");
    
    datatipSetStyle(e(1),1,%t);
endfunction

//plot_H_analog(H_analog_bsf, criticalW_bsf, delta);

//denominator_H = poly(left_poles, 'sL', "roots");
//disp(denominator_H);
//H_analog_LPF = (-Wc^7)/denominator_H; // This is the required analog LPF transfer function
//disp("H_analog_LPF")
//disp(H_analog_LPF);
function plot_H_analog_LPF(H_LPF, params_for_plot, name) // customize later
    // params_for_plot = [W_s, W_p, delta]
    // delta is the tolerance
    delta = params_for_plot(3);
    W_s = params_for_plot(1);
    W_p = params_for_plot(2);
    j = %i;
    W_axis = -3:0.001:3;
    H_values_analogLPF = horner(H_LPF, j*W_axis)
    plot(W_axis,ones(1,length(W_axis)),'b--'); // horizontal line at 1
    plot(W_axis,(1-delta)*ones(1,length(W_axis)),'b--'); // horizontal line at 1
    plot(W_axis,(delta)*ones(1,length(W_axis)),'b--'); // horizontal line at 1
    plot(W_p*ones(1, length(W_axis)), abs(H_values_analogLPF),'b--');
    plot(W_s*ones(1, length(W_axis)), abs(H_values_analogLPF),'b--');
    plot(-W_p*ones(1, length(W_axis)), abs(H_values_analogLPF),'b--');
    plot(-W_s*ones(1, length(W_axis)), abs(H_values_analogLPF),'b--');
    plot(W_axis, abs(H_values_analogLPF),'r'); 
    e=gce();p=e.children(1);
    loc_ws1 = floor((-W_s+3)/0.001);
    loc_wp1 = floor((-W_p+3)/0.001);
    loc_wp2 = floor((W_p+3)/0.001);
    loc_ws2 = floor((W_s+3)/0.001);
    title(name);
    xlabel("wL")
    ylabel("H_analog_LPF(wL)")
    d1=datatipCreate(p,loc_ws1);
    d2 = datatipCreate(p,loc_wp1 );
    d3 = datatipCreate(p, loc_wp2);
    d4 = datatipCreate(p, loc_ws2);
    
    datatipSetStyle(e(1),1,%t);
    datatipSetOrientation(d1,"upper left");
    datatipSetOrientation(d2,"upper left");
    datatipSetOrientation(d3,"upper right");
    datatipSetOrientation(d4,"upper right");
    
    set(gca(),"auto_clear","off");
    
    //legend(["H_analog_LPF";"1";"1-$\delta$";"W_p";"W_s";"-W_p";"-W_s"])
    
endfunction
//
//lpf_specs_bpf = [-1.4161988,-1.,1.,1.4161988]; //bpf
//lpf_specs_bsf = [-1.2394869,-1.,1.,1.2394869]; //bsf
//params_for_bsf_plot = [lpf_specs_bsf(4), lpf_specs_bsf(3), delta];
//plot_H_analog_LPF(H_lpf_bsf, params_for_bsf_plot);
//params_for_bpf_plot = [W_s, W_p, delta];
//plot_H_analog_LPF(H_lpf_bpf, params_for_bpf_plot);
//// need to now obtain the transfer function of the BPF through the frequency transform
//s = poly(0,'s');
//s_L_temp = (s^2 + (W0^2))/(B*s);
//H_analog_BPF = horner(H_analog_LPF, s_L_temp);
//
//disp("H_analog_BPF");
//disp(H_analog_BPF); // The desired H_analog Band pass filter transfer function!!
//
//H_values_analogBPF = horner(H_analog_BPF, j*W_axis);
////
//
////W_axis = 0:0.001:2.5;
////plot(W_axis, abs(H_values_analogBPF)); 
////set(gca(),"auto_clear","off");
////plot(W_axis,ones(1,length(W_axis))); // horizontal line at 1
////plot(W_axis,(1-delta)*ones(1,length(W_axis))); // horizontal line at 1
////plot(W_axis,(delta)*ones(1,length(W_axis))); // horizontal line at 1
////plot(criticalW(1)*ones(1, length(W_axis)), abs(H_values_analogLPF));
////plot(criticalW(2)*ones(1, length(W_axis)), abs(H_values_analogLPF));
////plot(criticalW(3)*ones(1, length(W_axis)), abs(H_values_analogLPF));
////plot(criticalW(4)*ones(1, length(W_axis)), abs(H_values_analogLPF));
////
////title("Analog BPF");
////xlabel("s")
////ylabel("H_analog_BPF(s)")
//

//......................................................................
function [H_z]=discrete_filter_from_analog(H_analog)
    z = poly(0, 'z');
    s_replace = (z-1)/(z+1);
    H_z = horner(H_analog, s_replace);
endfunction

//H_z_bpf = discrete_filter_from_analog(H_analog_bpf);
//disp("H(z) BPF Butterworth");
//disp(H_z_bpf);
//
//H_z_bsf = discrete_filter_from_analog(H_analog_bsf);
//disp("H(z) BSF Chebyshev");
//disp(H_z_bsf);


//......................................................................


//// the last step: computing the discrete frequency response using bilinear transform
//z = poly(0, 'z');
//s_replace = (z-1)/(z+1);
//H_z = horner(H_analog_BPF, s_replace);
//
//disp("H(z), finally");
//disp(H_z);
//

function plot_H_z(H_z, criticalw, delta, name)
    w_axis = 0:0.001:3;
    j = %i;
    H_z_values = horner(H_z, exp(j*w_axis));
    

    plot(w_axis,ones(1,length(w_axis)),'b--'); // horizontal line at 1
    plot(w_axis,(1-delta)*ones(1,length(w_axis)),'b--'); // horizontal line at 1
    plot(w_axis,(delta)*ones(1,length(w_axis)),'b--'); // horizontal line at 1
    plot(criticalw(1)*ones(1, length(w_axis)), abs(H_z_values),'b--');
    plot(criticalw(2)*ones(1, length(w_axis)), abs(H_z_values),'b--');
    plot(criticalw(3)*ones(1, length(w_axis)), abs(H_z_values),'b--');
    plot(criticalw(4)*ones(1, length(w_axis)), abs(H_z_values),'b--');
    isoview off;
    plot(w_axis, abs(H_z_values), 'r'); 
    e=gce();p=e.children(1);
    loc_1 = floor(criticalw(1)/0.001);
    loc_2 = floor(criticalw(2)/0.001);
    loc_3 = floor(criticalw(3)/0.001);
    loc_4 = floor(criticalw(4)/0.001);
    title(name);
    xlabel("w (normalized frequency)")
    ylabel("H (w) ")
    d1=datatipCreate(p,loc_1);
    d2 = datatipCreate(p,loc_2 );
    d3 = datatipCreate(p, loc_3);
    d4 = datatipCreate(p, loc_4);
    
//    datatipSetStyle(e(1),1,%t);
    datatipSetOrientation(d1,"upper right");
    datatipSetOrientation(d2,"upper left");
    datatipSetOrientation(d3,"upper right");
    datatipSetOrientation(d4,"upper right");
    set(gca(),"auto_clear","off");
endfunction

//criticalw_bpf =[1.3099489, 1.3861088, 1.7669079, 1.8430677];
////criticalw_bsf = [1.3388018,  1.4354662,  2.2426138,  2.3392782];
//criticalw_bsf = [1.3388018,1.4354662,1.9187881,2.0154525];
//plot_H_z(H_z_bsf, criticalw_bsf);

//w_axis = 0:0.0001:3.14;
//H_z_values = horner(H_z, exp(j*w_axis));
//
//plot(w_axis, abs(H_z_values)); 
//set(gca(),"auto_clear","off");
//plot(w_axis,ones(1,length(w_axis))); // horizontal line at 1
//plot(w_axis,(1-delta)*ones(1,length(w_axis))); // horizontal line at 1
//plot(w_axis,(delta)*ones(1,length(w_axis))); // horizontal line at 1
//plot(criticalw(1)*ones(1, length(w_axis)), abs(H_z_values));
//plot(criticalw(2)*ones(1, length(w_axis)), abs(H_z_values));
//plot(criticalw(3)*ones(1, length(w_axis)), abs(H_z_values));
//plot(criticalw(4)*ones(1, length(w_axis)), abs(H_z_values));
//isoview off;
//title("Discrete- Time Band Pass Filter");
//xlabel("w (normalized frequency")
//ylabel("H (w) ")
//
//

//------------------------------------------------------------------
//..................................................................












