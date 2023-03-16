

function [Tests]=select_tests_prepare_variables(suite,xius,field1,field2,field3,linear)
field_names = fieldnames(suite);
l = 1;
for k = 1:length(field_names)
    buf=suite.(field_names{k}).initial_data.Df_S;
    buf2 = suite.(field_names{k}).D_norm;
    if isreal(buf2)==0
        bla=0;
    end
    if buf == xius && isreal(buf2)
        ch(l) = k;
        l = l+1;
    end
    if ~isreal(buf2)
        bla = 0;
    end
end
Tests = ones(3,2000,length(ch)).*nan;
for i = 1:length(ch)
    T = suite.(field_names{ch(i)});
    if strcmp(field1,'time_nd')
        a = T.time.*T.initial_data.n;
        if strcmp(field2,'Drag')
            a = a(2:1:end)./2+a(1:1:end)./2;
        end
        
    elseif strcmp(field1,'time')
        a = T.time.*T.initial_data.tc./(365.25.*60.*60.*24.*1e6);
    elseif strcmp(field1,'dDdt')
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
        a1 = T.time;
        a2 = T.D_norm;
        a11 = diff(T.time);
        a22 = diff(T.D_norm);
        a = a22./a11;
        b = (a.*T.initial_data.Lambda.*(1./a2).^2); 
            

    end
    if strcmp(field3,'none')
        c = zeros(length(a),1);
    elseif strcmp(field3,'Lambda')
        if strcmp(linear,'Linear') == 0
            c1 = T.initial_data.Lambda./(1+T.initial_data.Df_UM);
        else
            c1 = T.initial_data.Lambda;
        end
        c  = (a./a)*c1;
    elseif stcm(field3,'tau_B')
        c = TD.tau(1,:);
    end
    Tests(1,1:length(a),i)=a;
    Tests(2,1:length(b),i)=b;
    Tests(3,1:length(c),i)=c;
end
bla = 0;
end


function [color_lists] = color_computation(size_tests,c,z_min,z_max)
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
