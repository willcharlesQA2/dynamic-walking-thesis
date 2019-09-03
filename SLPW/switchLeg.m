function [ Param ] = switchLeg( Param,step,modifier )
% Switches the legs from B to A and vice versa depending on if the step is
% odd or even

% modifier = 0 for even
% modifier = 1 for odd

% 
if mod(step,2) == modifier 
    Param.Lr1 = Param.LrA;
    Param.Lr2 = Param.LrB;
    
    Param.m1 = Param.mA;
    Param.m2 = Param.mB;
    
    Param.k1 = Param.kA;
    Param.k2 = Param.kB;
    
    Param.fr1 = Param.frA;
    Param.fr2 = Param.frB;
else                
    Param.Lr2 = Param.LrA;
    Param.Lr1 = Param.LrB;
    
    Param.m2 = Param.mA;
    Param.m1 = Param.mB;
    
    Param.k2 = Param.kA;
    Param.k1 = Param.kB;
    
    Param.fr2 = Param.frA;
    Param.fr1 = Param.frB;
end
end

