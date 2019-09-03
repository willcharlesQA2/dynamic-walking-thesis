figure; hold on

for h = 1:size(alpha_total,2)
    for i = 1:size(k_total,2)
        for j = 1:size(ME_total,2)
          if BOA(h,i,j) == maxsteps
              plot3(alpha_total(h)*-180/pi,k_total(i)/1000,ME_total(j),'b*')
          end
      end
  end
end

xlabel('Slope angle, \alpha, deg')
ylabel('Stiffness, k, kN/m')
zlabel('Mechanical Energy, J')