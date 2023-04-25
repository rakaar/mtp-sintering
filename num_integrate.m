clear all
close all
Diff_n_with_second_phase;
close all
d = sqrt( (x(1,:)-x(2,:)).^2 + (y(1,:)-y(2,:)).^2 );
t2 = find(d > 2);
t2 = t2(1);
d2 = d(1:t2);
v2 = vrn_vec(1:t2);

n_vals = -0.5:0.5:3;
all_numint = [];
figure
hold on
for n=n_vals
    f = k.*(d2.^n) + damping_frac.*(d2.^n).*(v2);
    numint = cumtrapz(d2, f);
    plot(numint)
    all_numint = [all_numint; numint];
end
hold off 
legend([strsplit(num2str(-0.5:0.5:3))])
title('num integrate')


figure
hold on
for n=n_vals
    f = k.*(d2.^n) + damping_frac.*(d2.^n).*(v2);
    numint = cumtrapz(d2, f);
    numdiff = diff(numint);
    plot(numdiff)
end
hold off 
legend([strsplit(num2str(-0.5:0.5:3))])
title('differentiate')

figure
hold on
for n=n_vals
    f = k.*(d2.^n) + damping_frac.*(d2.^n).*(v2);
    plot(f)
end
hold off 
legend([strsplit(num2str(-0.5:0.5:3))])
title('only f')

%%

auc = [];
for i=1:8
    auc = [auc trapz(all_numint(i,:))*dt ];
end

figure
    plot(n_vals,auc)
    title('Area under Energy vs time curve for different orders')
    xlabel('order')
    ylabel('Area')