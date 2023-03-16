classdef Manuscript_function_Container
    methods
        function [Tests]=select_tests_prepare_variables(obj,suite,xius,field1,field2,field3,linear,lessequalhigh,Value)
            field_names = fieldnames(suite);
            l = 1;
            for k = 1:length(field_names)
                if strcmp(Value,'xium')
                     buf=suite.(field_names{k}).initial_data.Df_UM;
                elseif strcmp(Value,'tcVdVn')
                    Vd =  suite.(field_names{k}).Meta_dataReal.Vd;
                    Vn =  suite.(field_names{k}).Meta_dataReal.Vn;
                    buf = suite.(field_names{k}).initial_data.tc./suite.(field_names{k}).initial_data.n;
                    if Vd >=2e-6 && Vd <=10e-6 && Vn >= 15*10^(-6) && Vn <=20*10^-6 && buf<xius 
                        buf =1;
                    else
                        buf= 0; 
                    end
                else
                    buf=suite.(field_names{k}).initial_data.Df_S;
                end
                buf2 = suite.(field_names{k}).D_norm;
                if isreal(buf2)==0
                    bla=0;
                end
                if strcmp(Value,'tcVdVn')
                    if buf == 1
                       ch(l) = k;
                        l = l+1;
                    end
                else
                if lessequalhigh==1
                    if buf == xius && isreal(buf2)
                        ch(l) = k;
                        l = l+1;
                    end
                elseif lessequalhigh==2
                     if buf <= xius && isreal(buf2)
                        ch(l) = k;
                        l = l+1;
                     end
                else 
                     if buf >= xius && isreal(buf2)
                        ch(l) = k;
                        l = l+1;
                     end
                end
                end
                
            end
            disp(length(ch))
            Tests = ones(3,2000,length(ch)).*nan;
            for i = 1:length(ch)
                T = suite.(field_names{ch(i)});
                if strcmp(field1,'time_nd')
                    a = T.time.*T.initial_data.n;
                    if strcmp(field2,'Drag') |strcmp(field2,'Psi2') 
                        a = a(2:1:end)./2+a(1:1:end-1)./2;
                    end
                    
                elseif strcmp(field1,'time')
                    a = T.time.*T.initial_data.tc./(365.25.*60.*60.*24.*1e6);
                elseif strcmp(field1,'dDdt') | strcmp(field2,'dDdt')
                    a1 = T.time;
                    a2 = T.D_norm;
                    a11 = diff(T.time);
                    a22 = diff(T.D_norm);
                    a = a22./a11;
                elseif strcmp(field1,'D_norm')
                    a = T.D_norm;
                elseif strcmp(field1,'tau_B')
                    a = abs(T.tau(1,:));
                end
                if strcmp(field2,'D_norm')
                    b = T.D_norm;
                elseif strcmp(field2,'tau_eff')
                    b = T.tau(3,:);
                elseif strcmp(field2,'tauD')
                    b =abs(T.tau(2,:));
                    if strcmp(field1,'dDdt')
                        b1 = abs(T.tau(2,:));
                        b = (b1(1:1:end-1)+b1(2:1:end))/2;
                    end
                elseif strcmp(field2,'tauD_B')
                    b =abs(T.tau(2,:));
                    b =abs(b./T.tau(1,:));
                    if strcmp(field1,'dDdt')
                        b = (b(1:1:end-1)+b(2:1:end))/2;
                    end
                elseif strcmp(field2,'Drag')
                    % dD
                    eps = T.eps(1,:)';
                    if T.initial_data.Df_UM > 0.0 
                        b = abs((eps.*T.Lambda.*(1./T.D_norm)));
                    else
                        b = abs((eps.*T.initial_data.Lambda.*(1./T.D_norm)));
                    end
                elseif strcmp(field2,'vz_nd')
                    a1 = T.time;
                    a2 = T.D_norm;
                    a11 = diff(T.time);
                    a22 = diff(T.D_norm);
                    a2 = a2(2:1:end)+a2(1:1:end-1);
                    a2 = a2./2;
                    a3 = a22./a11;
                    b = abs(T.initial_data.alpha.*(1./a2).^2.*a22);
                    % Computation v_stokes [Valid only for non linear case, and I use the reference viscosity of the mantle]
                    g_nd = 9.81.*((T.initial_data.tc.^2)./T.initial_data.D0);
                    drho_nd = (2)./(g_nd.*T.initial_data.l0./T.initial_data.D0);
                    if T.initial_data.Df_UM > 0.0 
                        eta_um  = T.initial_data.ID_A.eta0DM./(T.initial_data.Df_UM);
                    else
                        eta_um  = T.initial_data.ID_A.eta0DM;
                    end
                    v_s     = (g_nd.*1.0*(T.initial_data.l0./T.initial_data.D0).*drho_nd)./(eta_um);
                    b = b./v_s; 
                elseif strcmp(field2,'Lambda_r')
                    b = T.Lambda;
                elseif strcmp(field2,'dDdt')
                    a1 = T.time;
                    a2 = T.D_norm;
                    a11 = diff(T.time);
                    a22 = diff(T.D_norm);
                    a2 = a2(2:1:end)+a2(1:1:end-1);
                    a2 = a2./2;
                    b = a22./a11;
                elseif strcmp(field2,'eps_ratio')
                    b = T.eps(1,:);
                elseif strcmp(field2,'Psi')
                    etaum=T.eta_um;
                    etaus=T.initial_data.ID_A.eta0DS./(1+T.initial_data.Df_S);
                    len=T.initial_data.len*5;
                    b = etaum./etaus; 
                    b = b;%(T.initial_data.ID_A.eta0DM./(1+T.initial_data.Df_UM));
                    
                elseif strcmp(field2,'Psi2')
                    etaum=T.eta_um;
                    etaus=T.initial_data.ID_A.eta0DS./(1+T.initial_data.Df_S);
                    len=T.initial_data.len*5;
                    b = etaum./etaus; 
                    a1 = T.time;
                    a2 = T.D_norm;
                    a11 = diff(T.time);
                    a22 = diff(T.D_norm);
                    a2 = a2(2:1:end)+a2(1:1:end-1);
                    a2 = a2./2;
                    b2 = abs(a11./a22);
                    b3 = b2.*(1./len).*a2.^2;
                    b  = b(2:1:end)+b(1:1:end-1);
                    b = b./2; 
                    b = b./b3; 
     
                end
                if strcmp(field3,'none')
                    c = zeros(length(a),1);
                elseif strcmp(field3,'xium') 
                    c = (a./a).*T.initial_data.Df_S./T.initial_data.Df_UM; 
                elseif strcmp(field3,'Lambda')
                    if strcmp(linear,'Linear') == 0
                        c1 = T.initial_data.Lambda./(1+T.initial_data.Df_UM);
                    else
                        c1 = T.initial_data.Lambda;
                    end
                    c  = (a./a)*c1;
                elseif strcmp(field3,'tau_B')
                    c = T.tau(1,:);
                elseif strcmp (field3,'cdcn')
                    Cd=T.Meta_dataReal.Cd;
                    Cn= T.Meta_dataReal.Cn;
                    expCd = ((-Cd.*T.Meta_dataReal.w.*T.initial_data.l0)./T.Meta_dataReal.phi);
                    expCn = ((-Cn.*T.Meta_dataReal.w.*T.initial_data.l0)./T.Meta_dataReal.phi);
                    c     = (a./a).*(Cd./Cn).*expCd./expCn; 
                end
                Tests(1,1:length(a),i)=a;
                Tests(2,1:length(b),i)=b;
                Tests(3,1:length(c),i)=c;
            end
            bla = 0;
        end
        
        function [color_lists] = color_computation(obj,size_tests,c,z_min,z_max)
            V = zeros(size_tests,1);
            for i = 1:size_tests
                
                V(i) = (log10(c(1,i))-z_min)/(z_max-z_min);
                
                if V(i)<0
                    V(i)=0;
                elseif V(i)>1
                    V(i)=1;
                end
                V(i)=round(1+V(i)*(256-1));%round to nearest index
            end
            color_lists=V;
        end
        
    end
end