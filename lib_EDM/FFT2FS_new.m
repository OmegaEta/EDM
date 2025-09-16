% f is the sequence to be fitted.
% Fourier fitting function.
function FS = FFT2FS_new(f)
    N = length(f);
    F = fft(f)/N;      
   
    a = zeros(N,1);
    b = zeros(N,1);
    for n = 0:N-1
        k = n + 1;        
        if n <= N/2
            a(n+1) = 2*real(F(k));  % a_0
            b(n+1) = -2*imag(F(k));
        else
            a(n+1) = 2*real(F(k)); 
            b(n+1) = 2*imag(F(k));  
        end
    end
    a(1) = a(1)/2;
    
    FS.a = a;
    FS.b = b;

end