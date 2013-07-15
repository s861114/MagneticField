classdef wksp
    properties (Constant = true)
        hev=4.135667516e-15; %ev
        hevbar=6.58211928e-16; %ev
        h=6.62606957e-34; %J
        hbar=1.054571726e-34; %J
        
        kb=8.6173324e-5; %eV/Kelvin
        c=299792458; %m/s
        e=1.602e-19; %C
        dielec=8.85418*10^(-12);    
        imp=376.730313473762;
        e_mass=9.109e-31;
        
        a=2.46e-10; %meter
        d=0.334*10^(-9);
        dgnrcy=2*2; %Spin + Valley Degeneracy
        
        T=300; %Kelvin
        
        eMax=1; %eV
        eMin=0.01; %eV
        eRes=0.001; %eV
        
        kMax=3e+10; %1/meter (limit : 3.1161e+09)
        kMin=0.0000011; %1/meter
        kRes=0.0005e+10; %1/meter
        
        kMax_sh=2e+10; %1/meter (limit : 3.1161e+09)
        kMin_sh=0.0001; %1/meter
        kRes_sh=0.001e+10; %1/meter
        
        
        phiN=100;        

        Name=char('Mono','AB', 'AA', 'ABC','ABA','ABCA','ABCB','ABAB','ABAC','ABCAB','ABCAC','ABCBC','ABCBA','ABABC','ABABA','ABACA','ABACB'); 
        color=char('red','green','yellow','blue','red','yellow');
        Nlayer=[1 2 2 3 3 4 4 4 4 5 5 5 5 5 5 5 5]; %Number of layers
        Nband=2*[1 2 2 3 3 4 4 4 4 5 5 5 5 5 5 5 5]; %Number of bands
        size_H=[1 4 2 5 4];
        
        datafolder=char('data_nointra');
        
    end
end