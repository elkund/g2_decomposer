function fig = multiq_plot_result(result,transform_s)
% Inspect a result. This currently works for the Diffusive + Ballistic case
% Needs to be generalized for other models.

% q values are assumed to be given in inverse Ångström.

Q = result.Q;
L = result.L;
M = result.M;

% Create title
title_str ='';
for l =1:L
    if transform_s == 1
        title_str = append(title_str,"\int\rho_"+l+"(x)e^{-(xq^" + result.q_powers(l)/result.t_powers(l) + "\tau)^" + result.t_powers(l) + "}\mathrm{d}x +");
    else
        title_str = append(title_str,"\int\rho_"+l+"(s)e^{-(s(q/q_0)^" + result.q_powers(l)/result.t_powers(l) + "(\tau/\tau_0))^" + result.t_powers(l) + "}\mathrm{d}s +");
    end
end


%styling
fsize= 4E2;
lgd_fs=7;
markersize = 3;
c_map=parula(Q+1);

% The label dicts work for diffusive/ballistic case. Will probably break
% for other cases.
subtitles = dictionary([2 1 0],{"Diffusive component","Ballistic component","Relaxational component"});
x_labels = dictionary([2 1 0],{['Diffusivity $[ \mathrm{\AA^{2}/s}]$'],['Velocity $[\mathrm{\AA/s}]$'],['Rel. rate $[\mathrm{s}^{-1}]$']});
label_keys = result.q_powers./result.t_powers;

figure('Position', [1 1 3*fsize fsize ]);
tiles = tiledlayout(L,3,'padding','compact','TileSpacing','compact');

%plot g2s w. fits
g2ax = nexttile(1,[L 1]);

hold on
for q = 1:Q
    label = string(round(result.q_value(q)*100,2));
    errorbar(result.t,result.g2(q,:),result.g2_error(q,:),'.','CapSize',3,'MarkerSize',markersize, 'color',c_map(q,:),'DisplayName',label);
    plot(result.fit_eval_t,result.g2_fits(q,:),'-.','LineWidth',1,'Color',[c_map(q,:) 0.5],'HandleVisibility','off');

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

f_ax=nexttile(2, [L 1]);

lstyles = ["--",":","-."];

hold on
for q = 1:Q
    if q == 1
        lgd_vis = 'on';
    else
        lgd_vis = 'off';
    end
    f_est = sqrt((result.g2(q,:)-result.baseline(q)-1)/result.contrast(q));
    scatter(result.t, f_est, markersize,c_map(q,:),'HandleVisibility','off');
    plot(result.fit_eval_t,result.f_fits(q,:),'-','LineWidth',1,'Color',[c_map(q,:) 0.5],'HandleVisibility',lgd_vis,'DisplayName','Total fit');
    
    for l=1:L
        plot(result.fit_eval_t,result.f_comp_fits(q,:,l),lstyles(l),'LineWidth',1,'Color',c_map(q,:),'HandleVisibility',lgd_vis,'DisplayName',subtitles{label_keys(l)});
    end
end
flgd = legend('Location','southwest','Interpreter','Latex');
fontsize(flgd,lgd_fs,'points');

ylim([0 1.1]);
xscale log;
xlabel('$\tau\quad[\mathrm{s}]$','Interpreter','Latex');
ylabel('$|f(q,\tau)|$','Interpreter','Latex');

%plot dists

for l=1:L
    dist_ax = nexttile(3*l);

    q0=result.q_norm;
    tau0=result.t_norm;

    y = result.w(1+M*(l-1):M*l)'.*result.component_dist(:,l);

    if transform_s == 1
        x = (result.s(1+M*(l-1):M*l)-result.w(1+M*(l-1):M*l)/2)/(q0^(result.q_powers(l)/result.t_powers(l))*tau0);
        %dsdx = q0^result.q_powers(l)*tau0^result.t_powers(l)*result.t_powers(l)*x.^(result.t_powers(l)-1);
        %dsdx=q0^(result.q_powers(l)/result.t_powers(l))*tau0;
        % y= dsdx.*result.component_dist(:,l);
        %y = result.w(1+M*(l-1):M*l)'.*result.component_dist(:,l);
        xlab = x_labels(label_keys(l));
        %ylabel("$\rho_" +l+"\Delta x$", 'interpreter','latex');
    else
        x = result.s(1+M*(l-1):M*l)-result.w(1+M*(l-1):M*l)/2;
        %y = result.w(1+M*(l-1):M*l)'.*result.component_dist(:,l);
        xlab = 's';
        %ylabel("$\rho_" +l+"\Delta s$", 'interpreter','latex');
    end

    [xs,ys] = stairs(x,y);%,'LineWidth',1);
    area(xs,ys,0,'FaceColor','blue','FaceAlpha',0.2,'LineStyle',lstyles(l),'LineWidth',0.75,'EdgeColor','blue');

    ylim([0 max(y)]);

    xscale log;

    xlabel(xlab,'interpreter','latex');
    ylabel("$\rho_" +l+"\Delta x$", 'interpreter','latex');


    title(subtitles{label_keys(l)});
end

tit = sgtitle([regexprep(result.sample,'_','\\_') ' series ' num2str(result.series) ', $|f(q,\tau)|=' title_str{1}(1:end-1) '$']);
tit.Interpreter = 'Latex';

end