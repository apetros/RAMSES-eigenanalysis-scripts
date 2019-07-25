%% LOOP_ANALYSIS Creates an interactive selection list to analyze the output of the eigenvalue computations
function loop_analysis()
global S Jdyn_sp dif_states P raw_vars dom_V fig dom_eigenvals DF

exit=false;
while ~exit
    exit=false;
    fprintf('\nEigs) Show dominant eigenvalues\nMS) Mode Shape\nPF) For participation factors\nSJ) Show descriptor Jacobian matrix structure\nSJs) Show state Jacobian matrix structure\nQ) To exit\n\n');
    reply = input('Choose: ', 's');
    if strcmp(reply, 'Q') || strcmp(reply, 'q')
        exit=true ;
    end
    
    if strcmp(reply, 'Eigs')
        %% Display Dominant eigenvalues
        disp('Dominant eigenvalues: ')
        for i = 1:length(dom_eigenvals)
            fprintf('Eigenvalue %3d:\t%12.6f+%12.6fi\tDF: %f\tf= %f\n',i,real(dom_eigenvals(i)),imag(dom_eigenvals(i)),DF(i),abs(imag(dom_eigenvals(i)))/(2*pi));
        end 
    end
    
    if strcmp(reply, 'PF')
        i = input('Which dominant eigenvalue?: ');
        if isempty(i)
            fprintf('You cannot leave it empty...\n');
            continue
        end
        part_fact = input('PF threshold to show? (default=0.0): ');
        if isempty(part_fact)
            part_fact=0.0;
        end
        fprintf('\nEigenvalue %3d:\t%12.6f+%12.6fi\tDF: %f\tf= %f\n',i,real(dom_eigenvals(i)),imag(dom_eigenvals(i)),DF(i),abs(imag(dom_eigenvals(i)))/(2*pi));
        fprintf('Participation Factors:\n');
        for k = 1:length(dif_states)
            if P(k,i)>part_fact
                fprintf('%d)\t%f\t%s\t%s\t%s\n',k, P(k,i),char(raw_vars{3}(dif_states(k))),char(raw_vars{4}(dif_states(k))), char(raw_vars{5}(dif_states(k))));
            end
        end
    end
    
    if strcmp(reply, 'MS')
        fig=fig+1;
        cc=hsv(12);
        col=1;
        i = input('Which dominant eigenvalue?: ');
        if isempty(i)
            fprintf('You cannot leave it empty...\n');
            continue
        end
        if i>numel(dom_eigenvals) || i<=0
            fprintf('The eigenvalue selected does not exists\n');
            continue
        end
        figure(fig)
        hold on
        normval=0;
        for k = 1:length(dif_states)
            if strcmp(char(raw_vars{5}(dif_states(k))),'omega')
                a=real(dom_V(k,i));
                b=imag(dom_V(k,i));
                normval=max(normval,sqrt(a^2+b^2));
            end
        end
        for k = 1:length(dif_states)
            if strcmp(char(raw_vars{5}(dif_states(k))),'omega')
                a=real(dom_V(k,i))/normval;
                b=imag(dom_V(k,i))/normval;
                fprintf('%s:\t%12.6f+%12.6fi\n',char(raw_vars{4}(dif_states(k))),a,b);
                plot(a,b,'r*')
                line([0,a],[0,b],'color',cc(col,:))
                text(a,b,char(raw_vars{4}(dif_states(k))),'color',cc(col,:),'FontSize',16)
                col=col+1;
                if col==13
                    col=1;
                end
            end
        end
        hold off
    end
    if strcmp(reply, 'SJ')
        fig=fig+1;
        figure(fig)
        spy(S)
    end
    if strcmp(reply, 'SJs')
        fig=fig+1;
        figure(fig)
        spy(Jdyn_sp)
    end
    
end
end