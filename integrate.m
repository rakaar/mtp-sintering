x = 0.3*rand(100,1);
xp = x(x > 0);
x0 = 0;
dt = 0.001;
dxdt = [];
for i=1:100
    xi = xp(i);
    if i-1 == 0
        xim1 = x0;
    else
        xim1 = xp(i-1);
    end

    dxdt = [dxdt ((xi - xim1)/dt)];
end


c = 0.001;


intxdxdt = [];
dx = 0.001;
for i=1:100
    intxdxdt = [intxdxdt xp(i)*dxdt(i)*dx  ];
end
