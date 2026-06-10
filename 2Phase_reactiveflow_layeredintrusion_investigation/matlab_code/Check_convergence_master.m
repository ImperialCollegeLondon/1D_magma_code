if To_check_convergence    
    if counter>=Start_check
        Max_iter=600;
        if iter==0 && counter==Start_check
            Test_composition_update()
            
            figure(102)
            handel_convergence=plot3(0,0,0,'o','markersize',3,'linewidth',3,'color','white');
            handel_iter=text(0.2,2e6,'0','fontsize',16,'color','white');
            handel_convergence_numH=text(0.2,1.95e6,'0','fontsize',16,'color','white');
            handel_convergence_numB=text(0.2,1.9e6,'0','fontsize',16,'color','white');
        else
            bad_Cb=find(abs(Cb(1:N)-Cb_nonlinear(1:N))>Precision);
            bad_H=find(abs(H(1:N)-H_nonlinear(1:N))>Precision*Lf);
            
            bad=union(bad_Cb,bad_H);
            num_bad=length(bad);
            num_badH=length(bad_H);
            num_badC=length(bad_Cb);
            set(handel_convergence,'xdata',Cb(bad),'ydata',H(bad),'zdata',2*ones(num_bad,1));
            set(handel_convergence_numH,'string',num2str(num_badH))
            set(handel_convergence_numB,'string',num2str(num_badC))
            set(handel_iter,'string',num2str(iter))
            drawnow;
        end
    end
end