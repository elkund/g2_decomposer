function fig = multiq_plot_result(result)
% Inspect a result. This makes a lot of assumptions about the analysis
% setup. Needs to be adapted to a more general case

Q = length(result.qs);
L = size(result.s,1);
M = size(result.T,2)/L;


% add f per comp (move to analysis script?)
% 
% f_comps = zeros(Q,result.fit_N,L);
% for l = 1:L
%     solsub = [result.solutions(1:2*Q,end)' result.solutions(1+2*Q+M*(l-1):2*Q+M*l,end)']';
%     f_comps(:,:,l) = f_gen(result.fit_T(:,1+M*(l-1):M*l,:),solsub,Q,result.fit_N,M,1);
% end
% result.f_comp_fits = f_comps;
%norm = w*solution(2*K+1:end);

%title_str = arrayfun(@(r,t,f) sprintf('+\\rho_%d(s)e^{-q^%d\\tau^%d s}',r,t,f), 1:L,time_deps,q_deps,'UniformOutput',false);
title_str ='';
for l =1:L
    title_str = append(title_str,"\int\rho_"+l+"(x)e^{-(xq^" + result.q_deps(l)/result.time_deps(l) + "\tau)^" + result.time_deps(l) + "}\mathrm{d}x +");
end


%np = 5E6; %number of points in dist. plots

%styling
fsize= 5E2;
lgd_fs=7;

markersize = 3;


start_q = 1;

c_map=parula(Q+1);

y_labels = dictionary([1 0 -1 -2],{"Diffusive component","Ballistic component", "q-independent component", "q-independent component"});
%x_labels = dictionary([0 1],{'s (\propto {D^{-1}})','s (\propto {v^{-2}})'});
x_labels = dictionary([1 0 -1 -2],{['Diffusivity $[ \mathrm{\AA^{2}/s}]$'],['Velocity $[\mathrm{\AA/s}]$'],'Rel. rate $[\mathrm{s}^{-1}]$','Rel. rate [s^{-1}]'});
label_keys = result.q_deps-result.time_deps;


figure('Position', [1 1 3*fsize fsize ]);
t=tiledlayout(L,3,'padding','compact','TileSpacing','compact');

%plot g2s w. fits
g2ax = nexttile(1,[L 1]);

hold on
for q = start_q:Q
    label = string(round(result.qs(q)*100,2));
    errorbar(result.delays,result.g2s(q,:),result.g2errs(q,:),'.','CapSize',3,'MarkerSize',markersize, 'color',c_map(q,:),'DisplayName',label);
    plot(result.fit_eval_delays,result.g2_fits(q,:),'-.','LineWidth',1,'Color',[c_map(q,:) 0.5],'HandleVisibility','off');

end
%g2lgd = legend('Location','northeast','NumColumns',1,'Interpreter','Latex');
g2lgd = legend('Location','southwest','NumColumns',2,'Interpreter','Latex');
%g2lgd.Title.String.Interpreter = 'Latex'
g2lgd.Title.String = "$q$ in $10^{-2}\cdot\mathrm{\AA}^{-1}$";

fontsize(g2lgd,lgd_fs,'points');
%xlim([result.fit_eval_delays(1) result.fit_eval_delays(end)]);
xscale log;
xlabel('$\tau\quad[s]$','Interpreter', 'latex');
ylabel('$g_2$','interpreter','latex');

% fs w. fits

fax=nexttile(2, [L 1]);

lstyles = dictionary([1 2],{'--',':'});

hold on
for q = 1:Q;%2:2
    if q == start_q
        lgd_vis = 'on';
    else
        lgd_vis = 'off';
    end
    f_est = sqrt((result.g2s(q,:)-result.baseline_fits(q)-1)/result.contrast_fits(q));
    scatter(result.delays, f_est, markersize,c_map(q,:),'HandleVisibility','off');
    plot(result.fit_eval_delays,result.f_fits(q,:),'-.','LineWidth',1,'Color',[c_map(q,:) 0.5],'HandleVisibility',lgd_vis,'DisplayName','Total fit');
    
    for l=1:L
        plot(result.fit_eval_delays,result.f_comp_fits(q,:,l),lstyles{l},'LineWidth',1,'Color',c_map(q,:),'HandleVisibility',lgd_vis,'DisplayName',y_labels{label_keys(l)});
    end
end
flgd = legend('Location','southwest','Interpreter','Latex');
fontsize(flgd,lgd_fs,'points');

%xlim([result.fit_eval_delays(1) result.fit_eval_delays(end)]);
ylim([0 1.1]);
xscale log;
xlabel('$\tau\quad[\mathrm{s}]$','Interpreter','Latex');
ylabel('$|f(q,\tau)|$','Interpreter','Latex');

%pbaspect(fax, [1 1 1])
% plot dists.


%plot dists

for l=1:L
    distax = nexttile(3*l);

    q0=max(result.qs);
    tau0=max(result.delays);
    %t0=max(qs);
    x = (result.s(l,:)/(q0^result.q_deps(l))).^(1/result.time_deps(l))'/tau0;
    dsdx = q0^result.q_deps(l)*result.time_deps(l)*x.^(result.time_deps(l)-1)*tau0;
    y= dsdx.*result.components(:,l);
    %trapz(x,y)

    scatter(x,y,'LineWidth',1);

    ylim([0 max(y)]);

    xscale log;
    %x_var='x';
    xlabel(x_labels(label_keys(l)),'interpreter','latex');
    ylabel("$\rho_" +l+"$", 'interpreter','latex');

    title(y_labels{label_keys(l)});

    %pbaspect(distax,[2 1 1])
end

%tit = sgtitle([regexprep(result.sample,'_','\\_') ' series ' num2str(result.series) ', $T = ' num2str(result.temperature) '\;\mathrm{K}$, $|F(q,\tau)|=' title_str{1}(1:end-1) '$']);
tit = sgtitle([regexprep(result.sample,'_','\\_') ' series ' num2str(result.series) ', $|f(q,\tau)|=' title_str{1}(1:end-1) '$']);
tit.Interpreter = 'Latex';

%export_fig(gcf,'out.png');
end