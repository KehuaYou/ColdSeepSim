%function rrr=newton_raphson()
global COL ROW
global pw  sh  sw  sg  cl  cm T
global pcgw
global dnw  dng
global krw krg
global Msalt
global INDC2
global kx 
global cnstx 
global lambda phi
global residual

eps = 1e-8;

B=zeros(1,COL*ROW*5);
jac=spalloc(COL*ROW*5,COL*ROW*5,COL*ROW*5*25-20*10);
eps_pw=eps*1e6;
eps_cl=eps;
eps_T=eps;
lam_on=0;


for i=1:COL
    for j=1:ROW
        if sg(i,j)>0
            cm(i,j)=methane_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        elseif sh(i,j)>0
            cm(i,j)=hydrate_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        end
        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);
    end
end

for i=1:COL
    for j=1:ROW
        calc_trans(i,j);
    end
end

tic
for i=1:COL
    for j=1:ROW
        l=((i-1)*ROW+j)*5;
        calc_trans(i,j);
        
        rsidw0=calc_Rw(i,j);
        rsidg0=calc_Rg(i,j);
        rsidh0=calc_Rh(i,j);
        rsidc0=calc_Rc(i,j);
        rsidt0=calc_Rt(i,j);
        
        B(l-4)=-rsidw0;
        B(l-3)=-rsidg0;
        B(l-2)=-rsidh0;
        B(l-1)=-rsidc0;
        B(l)=-rsidt0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%          Code of primary variables switching                    %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%          1: pw, cm sh cl and T                                  %%%%%
        %%%%%          2: pw, cm sh cl and T  -> pw, sw, sh, cl and T         %%%%%
        %%%%%          3: pw, sw, sh, cl and T                                %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%                             i-1,j                              %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i > 1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change pw               %%%%%% %Water pressure.
            %%%%%%          Never switched.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pw(i-1,j)=pw(i-1,j)+eps_pw;
            temp1=dng(i-1,j);
            temp2=cm(i-1,j);
            temp3=dnw(i-1,j);
            
            dng(i-1,j)=gas_density(pw(i-1,j)/1e6,T(i-1,j));
            if sg(i-1,j)>0
                cm(i-1,j)=methane_solubility(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)))*0.016;
            elseif sh(i-1,j)>0
                cm(i-1,j)=hydrate_solubility(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)))*0.016;
            end
            dnw(i-1,j)=brine_density(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)),cm(i-1,j)/0.016);
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l-5*ROW-4)=(rsidw1-rsidw0)/eps_pw;
            jac(l-3,l-5*ROW-4)=(rsidg1-rsidg0)/eps_pw;
            jac(l-2,l-5*ROW-4)=(rsidh1-rsidh0)/eps_pw;
            jac(l-1,l-5*ROW-4)=(rsidc1-rsidc0)/eps_pw;
            jac(l,l-5*ROW-4)=(rsidt1-rsidt0)/eps_pw;
            
            pw(i-1,j)=pw(i-1,j)-eps_pw;
            dng(i-1,j)=temp1;
            cm(i-1,j)=temp2;
            dnw(i-1,j)=temp3;
            calc_trans(i,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change sw or cm         %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if INDC2(i-1,j) == 1
                cm(i-1,j)=cm(i-1,j)+eps;
                temp1=dnw(i-1,j);
                dnw(i-1,j)=brine_density(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)),cm(i-1,j)/0.016);
            else
                sw(i-1,j)=sw(i-1,j)+eps;
                temp1=lambda(i-1,j);
                temp2=krw(i-1,j);
                temp3=krg(i-1,j);
                temp4=pcgw(i-1,j);
                
                if lam_on==1
                    lambda(i-1,j)=(1-phi(i-1,j))*lambda_s+phi(i-1,j)*(sw(i-1,j)*lambda_w+sg(i-1,j)*lambda_g+sh(i-1,j)*lambda_h);
                end
                
                
                [krw(i-1,j),krg(i-1,j)] = calc_relative_perm(i-1,j);
                
            end
            
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l-5*ROW-3)=(rsidw1-rsidw0)/eps;
            jac(l-3,l-5*ROW-3)=(rsidg1-rsidg0)/eps;
            jac(l-2,l-5*ROW-3)=(rsidh1-rsidh0)/eps;
            jac(l-1,l-5*ROW-3)=(rsidc1-rsidc0)/eps;
            jac(l,l-5*ROW-3)=(rsidt1-rsidt0)/eps;
            
            
            if INDC2(i-1,j) == 1
                cm(i-1,j)=cm(i-1,j)-eps;
                dnw(i-1,j)=temp1;
            else
                sw(i-1,j)=sw(i-1,j)-eps;
                lambda(i-1,j)=temp1;
                krw(i-1,j)=temp2;
                krg(i-1,j)=temp3;
                pcgw(i-1,j)=temp4;
            end
            calc_trans(i,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change  sh              %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sh(i-1,j)=sh(i-1,j)+eps;
            temp1=lambda(i-1,j);
            temp2=kx(i-1,j);
            temp4=cnstx(i-1,j);
            if i>2
                temp5=cnstx(i-2,j);
            end
            temp8=krw(i-1,j);
            temp9=krg(i-1,j);
            temp10=pcgw(i-1,j);
            
            if lam_on==1
                lambda(i-1,j)=(1-phi(i-1,j))*lambda_s+phi(i-1,j)*(sw(i-1,j)*lambda_w+sg(i-1,j)*lambda_g+sh(i-1,j)*lambda_h);
            end
            
            [krw(i-1,j),krg(i-1,j)] = calc_relative_perm(i-1,j);
            
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l-5*ROW-2)=(rsidw1-rsidw0)/eps;
            jac(l-3,l-5*ROW-2)=(rsidg1-rsidg0)/eps;
            jac(l-2,l-5*ROW-2)=(rsidh1-rsidh0)/eps;
            jac(l-1,l-5*ROW-2)=(rsidc1-rsidc0)/eps;
            jac(l,l-5*ROW-2)=(rsidt1-rsidt0)/eps;
            
            sh(i-1,j)=sh(i-1,j)-eps;
            lambda(i-1,j)=temp1;
            kx(i-1,j)=temp2;
            cnstx(i-1,j)=temp4;
            if i>2
                cnstx(i-2,j)=temp5;
            end
            krw(i-1,j)=temp8;
            krg(i-1,j)=temp9;
            pcgw(i-1,j)=temp10;
            calc_trans(i,j);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change cl               %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cl(i-1,j)=cl(i-1,j)+eps;
            dnw(i-1,j)=brine_density(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)),cm(i-1,j)/0.016);
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l-5*ROW-1)=(rsidw1-rsidw0)/eps;
            jac(l-3,l-5*ROW-1)=(rsidg1-rsidg0)/eps;
            jac(l-2,l-5*ROW-1)=(rsidh1-rsidh0)/eps;
            jac(l-1,l-5*ROW-1)=(rsidc1-rsidc0)/eps;
            jac(l,l-5*ROW-1)=(rsidt1-rsidt0)/eps;
            
            cl(i-1,j)=cl(i-1,j)-eps;
            dnw(i-1,j)=brine_density(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)),cm(i-1,j)/0.016);
            calc_trans(i,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change T                %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T(i-1,j)=T(i-1,j)+eps_T;
            temp1=cm(i-1,j);
            temp2=dnw(i-1,j);
            temp3=dng(i-1,j);
            if sg(i-1,j)>0
                cm(i-1,j)=methane_solubility(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)))*0.016;
            elseif sh(i-1,j)>0
                cm(i-1,j)=hydrate_solubility(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)))*0.016;
            end
            
            dnw(i-1,j)=brine_density(pw(i-1,j)/1e6,T(i-1,j),cl(i-1,j)/Msalt/(1-cl(i-1,j)),cm(i-1,j)/0.016);
            dng(i-1,j)=gas_density(pw(i-1,j)/1e6,T(i-1,j));
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l-5*ROW)=(rsidw1-rsidw0)/eps_T;
            jac(l-3,l-5*ROW)=(rsidg1-rsidg0)/eps_T;
            jac(l-2,l-5*ROW)=(rsidh1-rsidh0)/eps_T;
            jac(l-1,l-5*ROW)=(rsidc1-rsidc0)/eps_T;
            jac(l,l-5*ROW)=(rsidt1-rsidt0)/eps_T;
            
            T(i-1,j)=T(i-1,j)-eps_T;
            cm(i-1,j)=temp1;
            dnw(i-1,j)=temp2;
            dng(i-1,j)=temp3;
            calc_trans(i,j);
        end % end of if i > 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%                              i, j-1                             %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%                              i, j                             %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change pw               %%%%%% %Water pressure.
        %%%%%%          Never switched.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pw(i,j)=pw(i,j)+eps_pw;
        temp1=dng(i,j);
        temp2=cm(i,j);
        temp3=dnw(i,j);
        
        dng(i,j)=gas_density(pw(i,j)/1e6,T(i,j));
        if sg(i,j)>0
            cm(i,j)=methane_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        elseif sh(i,j)>0
            cm(i,j)=hydrate_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        end
        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);
        
        calc_trans(i,j);
        rsidw1= calc_Rw(i,j);
        rsidg1= calc_Rg(i,j);
        rsidh1= calc_Rh(i,j);
        rsidc1=calc_Rc(i,j);
        rsidt1=calc_Rt(i,j);
        jac(l-4,l-4)=(rsidw1-rsidw0)/eps_pw;
        jac(l-3,l-4)=(rsidg1-rsidg0)/eps_pw;
        jac(l-2,l-4)=(rsidh1-rsidh0)/eps_pw;
        jac(l-1,l-4)=(rsidc1-rsidc0)/eps_pw;
        jac(l,l-4)=(rsidt1-rsidt0)/eps_pw;
        
        pw(i,j)=pw(i,j)-eps_pw;
        dng(i,j)=temp1;
        cm(i,j)=temp2;
        dnw(i,j)=temp3;
        calc_trans(i,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change sw or cm         %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if INDC2(i,j) == 1
            cm(i,j)=cm(i,j)+eps;
            temp1=dnw(i,j);
            dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);
        else
            sw(i,j)=sw(i,j)+eps;
            temp1=lambda(i,j);
            temp2=krw(i,j);
            temp3=krg(i,j);
            temp4=pcgw(i,j);
            
            if lam_on==1
                lambda(i,j)=(1-phi(i,j))*lambda_s+phi(i,j)*(sw(i,j)*lambda_w+sg(i,j)*lambda_g+sh(i,j)*lambda_h);
            end
            
            [krw(i,j),krg(i,j)] = calc_relative_perm(i,j);
        end
        
        
        calc_trans(i,j);
        rsidw1= calc_Rw(i,j);
        rsidg1= calc_Rg(i,j);
        rsidh1= calc_Rh(i,j);
        rsidc1=calc_Rc(i,j);
        rsidt1=calc_Rt(i,j);
        jac(l-4,l-3)=(rsidw1-rsidw0)/eps;
        jac(l-3,l-3)=(rsidg1-rsidg0)/eps;
        jac(l-2,l-3)=(rsidh1-rsidh0)/eps;
        jac(l-1,l-3)=(rsidc1-rsidc0)/eps;
        jac(l,l-3)=(rsidt1-rsidt0)/eps;
        
        
        if INDC2(i,j) == 1
            cm(i,j)=cm(i,j)-eps;
            dnw(i,j)=temp1;
        else
            sw(i,j)=sw(i,j)-eps;
            lambda(i,j)=temp1;
            krw(i,j)=temp2;
            krg(i,j)=temp3;
            pcgw(i,j)=temp4;
        end
        calc_trans(i,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change  sh              %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sh(i,j)=sh(i,j)+eps;
        temp1=lambda(i,j);
        temp2=kx(i,j);
        if i<COL
            temp4=cnstx(i,j);
        end
        if i>1
            temp5=cnstx(i-1,j);
        end
        temp8=krw(i,j);
        temp9=krg(i,j);
        temp10=pcgw(i,j);
        
        if lam_on==1
            lambda(i,j)=(1-phi(i,j))*lambda_s+phi(i,j)*(sw(i,j)*lambda_w+sg(i,j)*lambda_g+sh(i,j)*lambda_h);
        end
         [krw(i,j),krg(i,j)] = calc_relative_perm(i,j);
        
        calc_trans(i,j);
        rsidw1= calc_Rw(i,j);
        rsidg1= calc_Rg(i,j);
        rsidh1= calc_Rh(i,j);
        rsidc1=calc_Rc(i,j);
        rsidt1=calc_Rt(i,j);
        jac(l-4,l-2)=(rsidw1-rsidw0)/eps;
        jac(l-3,l-2)=(rsidg1-rsidg0)/eps;
        jac(l-2,l-2)=(rsidh1-rsidh0)/eps;
        jac(l-1,l-2)=(rsidc1-rsidc0)/eps;
        jac(l,l-2)=(rsidt1-rsidt0)/eps;
        
        sh(i,j)=sh(i,j)-eps;
        lambda(i,j)=temp1;
        kx(i,j)=temp2;
        if i<COL
            cnstx(i,j)=temp4;
        end
        if i>1
            cnstx(i-1,j)=temp5;
        end
        krw(i,j)=temp8;
        krg(i,j)=temp9;
        pcgw(i,j)=temp10;
        calc_trans(i,j);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change cl               %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cl(i,j)=cl(i,j)+eps;
        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);
        
        calc_trans(i,j);
        rsidw1= calc_Rw(i,j);
        rsidg1= calc_Rg(i,j);
        rsidh1= calc_Rh(i,j);
        rsidc1=calc_Rc(i,j);
        rsidt1=calc_Rt(i,j);
        jac(l-4,l-1)=(rsidw1-rsidw0)/eps;
        jac(l-3,l-1)=(rsidg1-rsidg0)/eps;
        jac(l-2,l-1)=(rsidh1-rsidh0)/eps;
        jac(l-1,l-1)=(rsidc1-rsidc0)/eps;
        jac(l,l-1)=(rsidt1-rsidt0)/eps;
        
        cl(i,j)=cl(i,j)-eps;
        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);
        calc_trans(i,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%          Change T                %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        T(i,j)=T(i,j)+eps_T;
        temp1=cm(i,j);
        temp2=dnw(i,j);
        temp3=dng(i,j);
        if sg(i,j)>0
            cm(i,j)=methane_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        elseif sh(i,j)>0
            cm(i,j)=hydrate_solubility(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)))*0.016;
        end
        
        dnw(i,j)=brine_density(pw(i,j)/1e6,T(i,j),cl(i,j)/Msalt/(1-cl(i,j)),cm(i,j)/0.016);
        dng(i,j)=gas_density(pw(i,j)/1e6,T(i,j));
        
        calc_trans(i,j);
        rsidw1= calc_Rw(i,j);
        rsidg1= calc_Rg(i,j);
        rsidh1= calc_Rh(i,j);
        rsidc1=calc_Rc(i,j);
        rsidt1=calc_Rt(i,j);
        jac(l-4,l)=(rsidw1-rsidw0)/eps_T;
        jac(l-3,l)=(rsidg1-rsidg0)/eps_T;
        jac(l-2,l)=(rsidh1-rsidh0)/eps_T;
        jac(l-1,l)=(rsidc1-rsidc0)/eps_T;
        jac(l,l)=(rsidt1-rsidt0)/eps_T;
        
        T(i,j)=T(i,j)-eps_T;
        cm(i,j)=temp1;
        dnw(i,j)=temp2;
        dng(i,j)=temp3;
        calc_trans(i,j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%                              i, j+1                             %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%           i+1, j                %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if i < COL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change pw               %%%%%% %Water pressure.
            %%%%%%          Never switched.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pw(i+1,j)=pw(i+1,j)+eps_pw;
            temp1=dng(i+1,j);
            temp2=cm(i+1,j);
            temp3=dnw(i+1,j);
            
            dng(i+1,j)=gas_density(pw(i+1,j)/1e6,T(i+1,j));
            if sg(i+1,j)>0
                cm(i+1,j)=methane_solubility(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)))*0.016;
            elseif sh(i+1,j)>0
                cm(i+1,j)=hydrate_solubility(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)))*0.016;
            end
            dnw(i+1,j)=brine_density(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)),cm(i+1,j)/0.016);
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l+5*ROW-4)=(rsidw1-rsidw0)/eps_pw;
            jac(l-3,l+5*ROW-4)=(rsidg1-rsidg0)/eps_pw;
            jac(l-2,l+5*ROW-4)=(rsidh1-rsidh0)/eps_pw;
            jac(l-1,l+5*ROW-4)=(rsidc1-rsidc0)/eps_pw;
            jac(l,l+5*ROW-4)=(rsidt1-rsidt0)/eps_pw;
            
            pw(i+1,j)=pw(i+1,j)-eps_pw;
            dng(i+1,j)=temp1;
            cm(i+1,j)=temp2;
            dnw(i+1,j)=temp3;
            calc_trans(i,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change sw or cm          %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if INDC2(i+1,j) == 1
                cm(i+1,j)=cm(i+1,j)+eps;
                temp1=dnw(i+1,j);
                dnw(i+1,j)=brine_density(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)),cm(i+1,j)/0.016);
            else
                sw(i+1,j)=sw(i+1,j)+eps;
                temp1=lambda(i+1,j);
                temp2=krw(i+1,j);
                temp3=krg(i+1,j);
                temp4=pcgw(i+1,j);
                
                if lam_on==1
                    lambda(i+1,j)=(1-phi(i+1,j))*lambda_s+phi(i+1,j)*(sw(i+1,j)*lambda_w+sg(i+1,j)*lambda_g+sh(i+1,j)*lambda_h);
                end
                
                 [krw(i+1,j),krg(i+1,j)] = calc_relative_perm(i+1,j);
                
            end
            
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l+5*ROW-3)=(rsidw1-rsidw0)/eps;
            jac(l-3,l+5*ROW-3)=(rsidg1-rsidg0)/eps;
            jac(l-2,l+5*ROW-3)=(rsidh1-rsidh0)/eps;
            jac(l-1,l+5*ROW-3)=(rsidc1-rsidc0)/eps;
            jac(l,l+5*ROW-3)=(rsidt1-rsidt0)/eps;
            
            
            if INDC2(i+1,j) == 1
                cm(i+1,j)=cm(i+1,j)-eps;
                dnw(i+1,j)=temp1;
            else
                sw(i+1,j)=sw(i+1,j)-eps;
                lambda(i+1,j)=temp1;
                krw(i+1,j)=temp2;
                krg(i+1,j)=temp3;
                pcgw(i+1,j)=temp4;
            end
            calc_trans(i,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change  sh              %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sh(i+1,j)=sh(i+1,j)+eps;
            temp1=lambda(i+1,j);
            temp2=kx(i+1,j);
            if i+1<COL
                temp4=cnstx(i+1,j);
            end
            if i<COL
                temp5=cnstx(i,j);
            end
            temp8=krw(i+1,j);
            temp9=krg(i+1,j);
            temp10=pcgw(i+1,j);
            
            if lam_on==1
                lambda(i+1,j)=(1-phi(i+1,j))*lambda_s+phi(i+1,j)*(sw(i+1,j)*lambda_w+sg(i+1,j)*lambda_g+sh(i+1,j)*lambda_h);
            end
            
            [krw(i+1,j),krg(i+1,j)] = calc_relative_perm(i+1,j);
            
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l+5*ROW-2)=(rsidw1-rsidw0)/eps;
            jac(l-3,l+5*ROW-2)=(rsidg1-rsidg0)/eps;
            jac(l-2,l+5*ROW-2)=(rsidh1-rsidh0)/eps;
            jac(l-1,l+5*ROW-2)=(rsidc1-rsidc0)/eps;
            jac(l,l+5*ROW-2)=(rsidt1-rsidt0)/eps;
            
            sh(i+1,j)=sh(i+1,j)-eps;
            lambda(i+1,j)=temp1;
            kx(i+1,j)=temp2;
            if i+1<COL
                cnstx(i+1,j)=temp4;
            end
            if i<COL
                cnstx(i,j)=temp5;
            end
            krw(i+1,j)=temp8;
            krg(i+1,j)=temp9;
            pcgw(i+1,j)=temp10;
            calc_trans(i,j);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change cl               %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cl(i+1,j)=cl(i+1,j)+eps;
            dnw(i+1,j)=brine_density(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)),cm(i+1,j)/0.016);
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l+5*ROW-1)=(rsidw1-rsidw0)/eps;
            jac(l-3,l+5*ROW-1)=(rsidg1-rsidg0)/eps;
            jac(l-2,l+5*ROW-1)=(rsidh1-rsidh0)/eps;
            jac(l-1,l+5*ROW-1)=(rsidc1-rsidc0)/eps;
            jac(l,l+5*ROW-1)=(rsidt1-rsidt0)/eps;
            
            cl(i+1,j)=cl(i+1,j)-eps;
            dnw(i+1,j)=brine_density(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)),cm(i+1,j)/0.016);
            calc_trans(i,j);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%          Change T                %%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            T(i+1,j)=T(i+1,j)+eps_T;
            temp1=cm(i+1,j);
            temp2=dnw(i+1,j);
            temp3=dng(i+1,j);
            if sg(i+1,j)>0
                cm(i+1,j)=methane_solubility(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)))*0.016;
            elseif sh(i+1,j)>0
                cm(i+1,j)=hydrate_solubility(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)))*0.016;
            end
            
            dnw(i+1,j)=brine_density(pw(i+1,j)/1e6,T(i+1,j),cl(i+1,j)/Msalt/(1-cl(i+1,j)),cm(i+1,j)/0.016);
            dng(i+1,j)=gas_density(pw(i+1,j)/1e6,T(i+1,j));
            
            calc_trans(i,j);
            rsidw1= calc_Rw(i,j);
            rsidg1= calc_Rg(i,j);
            rsidh1= calc_Rh(i,j);
            rsidc1=calc_Rc(i,j);
            rsidt1=calc_Rt(i,j);
            jac(l-4,l+5*ROW)=(rsidw1-rsidw0)/eps_T;
            jac(l-3,l+5*ROW)=(rsidg1-rsidg0)/eps_T;
            jac(l-2,l+5*ROW)=(rsidh1-rsidh0)/eps_T;
            jac(l-1,l+5*ROW)=(rsidc1-rsidc0)/eps_T;
            jac(l,l+5*ROW)=(rsidt1-rsidt0)/eps_T;
            
            T(i+1,j)=T(i+1,j)-eps_T;
            cm(i+1,j)=temp1;
            dnw(i+1,j)=temp2;
            dng(i+1,j)=temp3;
            calc_trans(i,j);
        end   % end of if i < COL
        
    end
end % end of for i=1:COL


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     calculate delta     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jac_temp=zeros(1,(COL-1)*ROW*5-1);
aa=1;
for i=2:COL
    for j=1:ROW
        if i<COL
            jac_temp(aa:aa+4)=[((i-1)*ROW+j-1)*5+1 ((i-1)*ROW+j-1)*5+2 ((i-1)*ROW+j-1)*5+3 ((i-1)*ROW+j-1)*5+4  ((i-1)*ROW+j-1)*5+5];
            aa=aa+5;
        else
            jac_temp(aa:aa+3)=[((i-1)*ROW+j-1)*5+1 ((i-1)*ROW+j-1)*5+2 ((i-1)*ROW+j-1)*5+3 ((i-1)*ROW+j-1)*5+4];
        end
    end
end
delta=real(jac(jac_temp,jac_temp))\real(B(jac_temp)');
delta=real(delta);
residual=real(B(jac_temp)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%     Update primary variables     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l=0;
for i=2:COL
    if i<COL
        for j=1:ROW
            l=l+5;
            if INDC2(i,j) == 1
                dpw=delta(l-4);
                dcm=delta(l-3);
                dsh=delta(l-2);
                dcl=delta(l-1);
                dT=delta(l);
                
                pw(i,j)=pw(i,j)+dpw;
                cm(i,j)=cm(i,j)+dcm;
                sh(i,j)=sh(i,j)+dsh;
                cl(i,j)=cl(i,j)+dcl;
                T(i,j)=T(i,j)+dT;
                
                sw(i,j)=1.0;
                sg(i,j)=0;
                sh(i,j)=0;
            else
                dpw=delta(l-4);
                dsw=delta(l-3);
                dsh=delta(l-2);
                dcl=delta(l-1);
                dT=delta(l);
                
                pw(i,j)=pw(i,j)+dpw;
                sw(i,j)=sw(i,j)+dsw;
                sh(i,j)=sh(i,j)+dsh;
                cl(i,j)=cl(i,j)+dcl;
                T(i,j)=T(i,j)+dT;
                
                sg(i,j)=1-sw(i,j)-sh(i,j);
            end % end of if INDC2(i,j)==3
        end % end of for j=1:ROW
    else
        for j=1:ROW
            l=l+4;
            if INDC2(i,j) == 1
                dpw=delta(l-3);
                dcm=delta(l-2);
                dsh=delta(l-1);
                dcl=delta(l);
                
                pw(i,j)=pw(i,j)+dpw;
                cm(i,j)=cm(i,j)+dcm;
                sh(i,j)=sh(i,j)+dsh;
                cl(i,j)=cl(i,j)+dcl;
                
                sw(i,j)=1.0;
                sg(i,j)=0;
                sh(i,j)=0;
            else
                dpw=delta(l-3);
                dsw=delta(l-2);
                dsh=delta(l-1);
                dcl=delta(l);
                
                pw(i,j)=pw(i,j)+dpw;
                sw(i,j)=sw(i,j)+dsw;
                sh(i,j)=sh(i,j)+dsh;
                cl(i,j)=cl(i,j)+dcl;
                
                sg(i,j)=1-sw(i,j)-sh(i,j);
            end % end of if INDC2(i,j)==3
        end % end of for j=1:ROW
    end
end % end of for i=1:COL



