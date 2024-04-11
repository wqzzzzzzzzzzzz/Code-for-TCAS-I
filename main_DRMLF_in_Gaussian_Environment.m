clear all; close all;

L = 5000;
p = 5;
type=1;
MEE_Length=8;
lamda_RLS = 0.985;
lamda_RLS_cta=0.985;

lamda_RMCC = 0.95;
lamda_RMCC_cta = 0.95;
lamda_RMEE=0.984;
lamda_RMEE_cta =0.982;
lamda_New_ATC=0.983;
lamda_New_CTA=0.983;


sigma1=4;%%RMCC
sigma2=4;%%RMEE
orders=2;
num=20; %number of node
%% Initial
load('C.mat');
MIU_EM=zeros(num,orders);
SIGMA_EM=zeros(num,orders);
WEIGHT_EM=zeros(num,orders);

UU=[];
for i=1:num
    UU=[UU randn(p,L)];
end
wo1 = randn(p,1);
wo = [ kron(wo1, ones(1,L)) ];

   %% Mean by EM
    MIU_EM(1,:)=[0,0.1];
    MIU_EM(2,:)=[0,0.1];
    MIU_EM(3,:)=[0.1,-0.1];
    MIU_EM(4,:)=[0,0];
    MIU_EM(5,:)=[0,0];
    MIU_EM(6,:)=[0.2,0];
    MIU_EM(7,:)=[0.1,0];
    MIU_EM(8,:)=[0,0];
    MIU_EM(9,:)=[0.1,0];
    MIU_EM(10,:)=[0.1,0];
    MIU_EM(11,:)=[-0.2,0];
    MIU_EM(12,:)=[0,0.1];
    MIU_EM(13,:)=[0,0];
    MIU_EM(14,:)=[0,0];
    MIU_EM(15,:)=[0.2,0];
    MIU_EM(16,:)=[0,0];
    MIU_EM(17,:)=[0,0.1];
    MIU_EM(18,:)=[0,0];
    MIU_EM(19,:)=[-0.1,0];
    MIU_EM(20,:)=[0,0];
    %% Sigma by EM
    SIGMA_EM(1,:)=[0.9,1.1];
    SIGMA_EM(2,:)=[0.7,1.2];
    SIGMA_EM(3,:)=[1.2,1.2];
    SIGMA_EM(4,:)=[0.7,1];
    SIGMA_EM(5,:)=[1.2,1.1];
    SIGMA_EM(6,:)=[0.8,1.1];
    SIGMA_EM(7,:)=[1.3,0.9];
    SIGMA_EM(8,:)=[1,0.9];
    SIGMA_EM(9,:)=[1.1,1.3];
    SIGMA_EM(10,:)=[1.2,1.3];
    SIGMA_EM(11,:)=[1,1.3];
    SIGMA_EM(12,:)=[0.9,0.8];
    SIGMA_EM(13,:)=[1,1.1];
    SIGMA_EM(14,:)=[1.1,0.9];
    SIGMA_EM(15,:)=[1.1,1];
    SIGMA_EM(16,:)=[0.9,0.9];
    SIGMA_EM(17,:)=[1,1.1];
    SIGMA_EM(18,:)=[0.8,1];
    SIGMA_EM(19,:)=[0.9,1.3];
    SIGMA_EM(20,:)=[1.2,0.7];
    %% Weight by EM
    WEIGHT_EM(1,:)=[0.52,0.48];
    WEIGHT_EM(2,:)=[0.49,0.51];
    WEIGHT_EM(3,:)=[0.54,0.46];
    WEIGHT_EM(4,:)=[0.5,0.5];
    WEIGHT_EM(5,:)=[0.52,0.48];
    WEIGHT_EM(6,:)=[0.47,0.53];
    WEIGHT_EM(7,:)=[0.49,0.51];
    WEIGHT_EM(8,:)=[0.53,0.47];
    WEIGHT_EM(9,:)=[0.54,0.46];
    WEIGHT_EM(10,:)=[0.5,0.5];
    WEIGHT_EM(11,:)=[0.48,0.52];
    WEIGHT_EM(12,:)=[0.5,0.5];
    WEIGHT_EM(13,:)=[0.52,0.48];
    WEIGHT_EM(14,:)=[0.5,0.5];
    WEIGHT_EM(15,:)=[0.53,0.47];
    WEIGHT_EM(16,:)=[0.46,0.54];
    WEIGHT_EM(17,:)=[0.51,0.49];
    WEIGHT_EM(18,:)=[0.53,0.47];
    WEIGHT_EM(19,:)=[0.46,0.54];
    WEIGHT_EM(20,:)=[0.47,0.53];


    MIU_Rand=MIU_EM;
    SIGMA_Rand=SIGMA_EM;
    WEIGHT_Rand=WEIGHT_EM;



tic
P=10;
for mm = 1 : P

    un1=zeros(p,num);
    dd=zeros(num,L);
    Fai=zeros(p,num);
    Fai_QPD=zeros(p,num);
    Fai_cta=Fai;
    FaiR=zeros(p,num);
    FaiR_cta=Fai;
    Fai_RMC=Fai;
    Fai_RLS=Fai;
    Fai_MCC=Fai;

    %% Noise
    v=normrnd(0,1,num,L);
    for i=1:num
        vv(i,:)=v(i,:)-mean(v(i,:));
        Sig(i)=var (vv(i,1:L));
    end


    for i=1:num
        for ii = 1 : L
            dd(i,ii) = wo(:,ii)' * UU(:,(i-1)*L+ii) + vv(i,ii);
        end
    end

    C=Weight(num,C,type,Sig,orders); 

    %% Iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ATC
    w_ATC=randn(p,num);
    w_ATC_RLS=w_ATC;
    w_ATC_RLS_BPD=w_ATC;
    w_MEE=w_ATC;
    w_CTA_MEE=w_ATC;
    w_ATC_RMC=w_ATC;
    w_ATC_MCC=w_ATC;
    for cnt=1:num
        Pn_RLS_BPD(:,:,cnt) = eye(p);
    end
    Pn_RLS=Pn_RLS_BPD;
    Pn_RMC=Pn_RLS_BPD;
    %% CTA
    w_CTA=w_ATC;
    w_CTA_RLS=w_ATC;
    w_CTA_RLS_BPD=w_ATC;
    w_CTA_RMC=w_ATC;
    w_CTA_MEE=w_ATC;
    w_M_RLS=w_ATC;
    w_M_CTA_RLS=w_ATC;
    Pn_RLS_BPD_cta=Pn_RLS_BPD;
    Pn_RLS_cta=Pn_RLS_BPD;
    Pn_RMC_cta=Pn_RLS_BPD;

    %% ATC
    for ii = 1 : L
        AddR=zeros(p,num);
        AddR_QPD=zeros(p,num);
        Add_ATC_RLS=zeros(p,num);
        Add_ATC_RMC=zeros(p,num);
        Add_ATC_MCC=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            Err_ATC_RLS(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS(:,cnt));
            Err_ATC_RMC(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RMC(:,cnt))' *  (wo(:,ii)  - w_ATC_RMC(:,cnt));
            Err_ATC_RLS_BPD(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_BPD(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_BPD(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);%
            %% ATC-RLS
            ek_ATC_RLS(cnt)=dn(cnt)-w_ATC_RLS(:,cnt)'*un1(:,cnt);
            kn_RLS(:,:,cnt) = Pn_RLS(:,:,cnt) * un1(:,cnt) / (lamda_RLS + un1(:,cnt)' * Pn_RLS(:,:,cnt) * un1(:,cnt)) ;
            Pn_RLS(:,:,cnt) = 1/lamda_RLS * ( Pn_RLS(:,:,cnt) - kn_RLS(:,:,cnt) *un1(:,cnt)' * Pn_RLS(:,:,cnt));
            F_ATC_RLS(:,cnt) = w_ATC_RLS(:,cnt) +kn_RLS(:,:,cnt) * ek_ATC_RLS(cnt);
            %% ATC-RMC
            ek_ATC_RMC(cnt)=dn(cnt)-w_ATC_RMC(:,cnt)'*un1(:,cnt);
            kn_RMC(:,:,cnt) = Pn_RMC(:,:,cnt) * un1(:,cnt) / ( exp(ek_ATC_RMC(cnt)^2/2/sigma1^2)*lamda_RMCC + un1(:,cnt)' * Pn_RMC(:,:,cnt) * un1(:,cnt) );
            Pn_RMC(:,:,cnt) = 1/lamda_RMCC * ( Pn_RMC(:,:,cnt) - kn_RMC(:,:,cnt) *un1(:,cnt)' * Pn_RMC(:,:,cnt));
            F_ATC_RMC(:,cnt) = w_ATC_RMC(:,cnt) +kn_RMC(:,:,cnt) * ek_ATC_RMC(cnt);
            %% ATC-RLS-BPD
            ek_Rand1(cnt)=dn(cnt)-w_ATC_RLS_BPD(:,cnt)'*un1(:,cnt);
            for k=1:orders
                P_Rand(ii,k+orders*(cnt-1) )=exp(-1*(ek_Rand1(cnt)-MIU_Rand(cnt,k))^2/(2*SIGMA_Rand(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand(cnt,k));
            end
            for k=1:orders
                V_Rand(ii,k+orders*(cnt-1) )=WEIGHT_Rand(cnt,k)*P_Rand(ii,k+orders*(cnt-1))/(WEIGHT_Rand(cnt,:)*P_Rand(ii,orders*(cnt-1)+1:orders*(cnt-1)+orders)');
            end
            R1=0;
            R2=0;
            for k=1:orders
                R1=R1+V_Rand(ii,k+orders*(cnt-1))/SIGMA_Rand(cnt,k)^2;
                R2=R2+V_Rand(ii,k+orders*(cnt-1))*MIU_Rand(cnt,k)/SIGMA_Rand(cnt,k)^2;
            end
            kn(:,:,cnt) = Pn_RLS_BPD(:,:,cnt) * un1(:,cnt) / ( lamda_New_ATC+ R1*un1(:,cnt)' * Pn_RLS_BPD(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_BPD(:,:,cnt) = 1/lamda_New_ATC * ( Pn_RLS_BPD(:,:,cnt) - R1*kn(:,:,cnt) * un1(:,cnt)' * Pn_RLS_BPD(:,:,cnt));
            FaiR(:,cnt) = w_ATC_RLS_BPD(:,cnt) +kn(:,:,cnt) * (R1*ek_Rand1(cnt)-R2);
        end
        for cnt=1:num
            for jj=1:num
                AddR(:,cnt)=C(cnt,jj)*FaiR(:,jj)+AddR(:,cnt);
                Add_ATC_RMC(:,cnt)=C(cnt,jj)*F_ATC_RMC(:,jj)+ Add_ATC_RMC(:,cnt);
                Add_ATC_RLS(:,cnt)=C(cnt,jj)*F_ATC_RLS(:,jj)+ Add_ATC_RLS(:,cnt);
            end
            w_ATC_RLS_BPD(:,cnt) = AddR(:,cnt);
            w_ATC_RMC(:,cnt) =Add_ATC_RMC(:,cnt);
            w_ATC_RLS(:,cnt) =Add_ATC_RLS(:,cnt);
        end
        %% CTA
        AddR_cta=zeros(p,num);
        AddR_QPD_cta=zeros(p,num);
        Add_RLS_cta=zeros(p,num);
        Add_RMC_cta=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            for jj=1:num
                AddR_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS_BPD(:,jj)+AddR_cta(:,cnt);
                Add_RMC_cta(:,cnt)=C(cnt,jj)*w_CTA_RMC(:,jj)+ Add_RMC_cta(:,cnt);
                Add_RLS_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS(:,jj)+ Add_RLS_cta(:,cnt);
            end
            FaiR_cta(:,cnt) = AddR_cta(:,cnt);
            Fai_RMC(:,cnt) =Add_RMC_cta(:,cnt);
            Fai_RLS(:,cnt) =Add_RLS_cta(:,cnt);
        end

        for cnt=1:num
            Err_CTA_RLS(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS(:,cnt));
            Err_CTA_RMC(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RMC(:,cnt))' *  (wo(:,ii)  - w_CTA_RMC(:,cnt));
            Err_CTA_RLS_BPD(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_BPD(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_BPD(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);%
            %% CTA-RLS
            ek_CTA_RLS(cnt)=dn(cnt)-Fai_RLS(:,cnt)'*un1(:,cnt);
            kn_RLS_cta(:,:,cnt) = Pn_RLS_cta(:,:,cnt) * un1(:,cnt) / (lamda_RLS_cta + un1(:,cnt)' * Pn_RLS_cta(:,:,cnt) * un1(:,cnt)) ;
            Pn_RLS_cta(:,:,cnt) = 1/lamda_RLS_cta * ( Pn_RLS_cta(:,:,cnt) - kn_RLS_cta(:,:,cnt) *un1(:,cnt)' * Pn_RLS_cta(:,:,cnt));
            w_CTA_RLS(:,cnt) = Fai_RLS(:,cnt) +kn_RLS_cta(:,:,cnt) * ek_CTA_RLS(cnt);
            %% CTA-RMC
            ek_CTA_RMC(cnt)=dn(cnt)-Fai_RMC(:,cnt)'*un1(:,cnt);
            kn_RMC_cta(:,:,cnt) = Pn_RMC_cta(:,:,cnt) * un1(:,cnt) / ( exp(ek_CTA_RMC(cnt)^2/2/sigma1^2)*lamda_RMCC_cta + un1(:,cnt)' * Pn_RMC_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RMC_cta(:,:,cnt) = 1/lamda_RMCC_cta * ( Pn_RMC_cta(:,:,cnt) - kn_RMC_cta(:,:,cnt) *un1(:,cnt)' * Pn_RMC_cta(:,:,cnt));
            w_CTA_RMC(:,cnt) = Fai_RMC(:,cnt) +kn_RMC_cta(:,:,cnt) * ek_CTA_RMC(cnt);
            %% CTA-RLS-BPD
            ek_Rand1_cta(cnt)=dn(cnt)-FaiR_cta(:,cnt)'*un1(:,cnt);
            for k=1:orders
                P_Rand_cta(ii,k+orders*(cnt-1) )=exp(-1*(ek_Rand1_cta(cnt)-MIU_Rand(cnt,k))^2/(2*SIGMA_Rand(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand(cnt,k));
            end
            for k=1:orders
                V_Rand_cta(ii,k+orders*(cnt-1) )=WEIGHT_Rand(cnt,k)*P_Rand_cta(ii,k+orders*(cnt-1))/(WEIGHT_Rand(cnt,:)*P_Rand_cta(ii,orders*(cnt-1)+1:orders*(cnt-1)+orders)');
            end
            R1=0;
            R2=0;
            for k=1:orders
                R1=R1+V_Rand_cta(ii,k+orders*(cnt-1))/SIGMA_Rand(cnt,k)^2;
                R2=R2+V_Rand_cta(ii,k+orders*(cnt-1))*MIU_Rand(cnt,k)/SIGMA_Rand(cnt,k)^2;
            end
            kn_cta(:,:,cnt) = Pn_RLS_BPD_cta(:,:,cnt) * un1(:,cnt) / ( lamda_New_CTA+ R1*un1(:,cnt)' * Pn_RLS_BPD_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_BPD_cta(:,:,cnt) = 1/lamda_New_CTA * ( Pn_RLS_BPD_cta(:,:,cnt) - R1*kn_cta(:,:,cnt) * un1(:,cnt)' * Pn_RLS_BPD_cta(:,:,cnt));
            w_CTA_RLS_BPD(:,cnt) = FaiR_cta(:,cnt) +kn_cta(:,:,cnt) * (R1*ek_Rand1_cta(cnt)-R2);
        end
    end
    %% ATC RMEE, algorithm
    for cnt=1:num
        PL(:,:,cnt) = eye(p)*1;
    end
    eL = zeros( MEE_Length, 1);
    PHI = zeros( MEE_Length,  MEE_Length);
    for ii =  MEE_Length : L
        Add_RMEE=zeros(p,num);
        F_RMEE= w_M_RLS;
        for cnt=1:num
            Err_RLS_MEE(mm,ii,cnt) = (wo(:,ii)  - w_M_RLS(:,cnt))' * (wo(:,ii)  - w_M_RLS(:,cnt));
            for jj = 1 :  MEE_Length
                eL(jj) = dd(cnt,ii - MEE_Length + jj) - w_M_RLS(:,cnt)' * UU(:,(cnt-1)*L+ii - MEE_Length + jj);
            end
            u0 = UU(:,(cnt-1)*L+ii);
            e0 = eL(MEE_Length);
            phi0 = 0;
            for kk = 2: MEE_Length
                ek = eL( MEE_Length-kk+1);
                phi0 = phi0 + lamda_RMEE^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);
            end
            KL(:,:,cnt) = PL(:,:,cnt)* u0 * inv ( lamda_RMEE^2 / phi0 + u0' * PL(:,:,cnt) * u0 );
            PL(:,:,cnt) = 1/lamda_RMEE^2 * (PL(:,:,cnt) - KL(:,:,cnt) * u0' * PL(:,:,cnt));
            F_RMEE(:,cnt) =F_RMEE(:,cnt) +KL(:,:,cnt) * e0;

        end
        for cnt=1:num
            for jj=1:num
                Add_RMEE(:,cnt)=C(cnt,jj)*F_RMEE(:,jj)+Add_RMEE(:,cnt);
            end
            w_M_RLS(:,cnt) = Add_RMEE(:,cnt);
        end
    end
    %% CTA RMEE, algorithm
    for cnt=1:num
        PL1(:,:,cnt) = eye(p)*1;
    end
    eL1 = zeros( MEE_Length, 1);
    for ii =  MEE_Length : L
        Add_CTA_RMEE=zeros(p,num);

        for cnt=1:num
            for jj=1:num
                Add_CTA_RMEE(:,cnt)=C(cnt,jj)*w_M_CTA_RLS(:,jj)+Add_CTA_RMEE(:,cnt);
            end
            w_M_CTA_RLS(:,cnt) = Add_CTA_RMEE(:,cnt);
        end
        for cnt=1:num
            Err_RLS_CTA_MEE(mm,ii,cnt) = (wo(:,ii)  - w_M_CTA_RLS(:,cnt))' * (wo(:,ii)  - w_M_CTA_RLS(:,cnt));
            for jj = 1 :  MEE_Length
                eL1(jj) = dd(cnt,ii - MEE_Length + jj) - w_M_CTA_RLS(:,cnt)' * UU(:,(cnt-1)*L+ii - MEE_Length + jj);
            end
            u0 = UU(:,(cnt-1)*L+ii);
            e0 = eL1(MEE_Length);
            phi0 = 0;
            for kk = 2: MEE_Length
                ek = eL( MEE_Length-kk+1);
                phi0 = phi0 + lamda_RMEE_cta ^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);
            end
            KL1(:,:,cnt) = PL1(:,:,cnt) * u0 * inv ( lamda_RMEE_cta ^2 / phi0 + u0' * PL1(:,:,cnt) * u0 );
            PL1(:,:,cnt) = 1/lamda_RMEE_cta ^2 * (PL1(:,:,cnt) - KL1(:,:,cnt) * u0' * PL1(:,:,cnt));
            w_M_CTA_RLS(:,cnt) =w_M_CTA_RLS(:,cnt) +KL1(:,:,cnt) * e0;
        end
    end

    disp(mm);
end




%% Mean square analysis
%% Real mean,sigma and weight of the noise for theoretical analysis
    MIU_Rand=zeros(num,orders);
    SIGMA_Rand=1*ones(num,orders);
    WEIGHT_Rand=0.5*ones(num,orders);

%% ATC
Seq=zeros(num,1);
for aa=1:num
    Seq(aa)=(aa-1)*num+aa;
end
NN=zeros(num^2,1);

for aa=1:num
    Temp=0;
    for kk=1:orders
        Temp=Temp+WEIGHT_Rand(num,kk)/SIGMA_Rand(num,kk)^2;
    end

    NN(Seq(aa))= ((1-lamda_New_ATC)^2 * p)/Temp;
end

%% CTA
C2=kron(C,C);
Para=eye(num^2)-(lamda_New_ATC^2) * C2;
Para_n=inv(Para)*C2;
Expet=Para_n*NN;
EE1=zeros(num,1);
for aa=1:num
    EE1(aa)=Expet(Seq(aa));
end

D=C-(1-lamda_New_CTA)*eye(num);
D2=kron(D,D);
Para1=eye(num^2)-D2;
Para1_n=inv(Para1);


Expet=Para1_n*NN;
EE2=zeros(num,1);
for aa=1:num
    EE2(aa)=Expet(Seq(aa));
end

EE=10*log10(EE1);
EE0=10*log10(EE2);



Err_ATC_RLS1=mean(Err_ATC_RLS,3);
Err_ATC_RMC1=mean(Err_ATC_RMC,3);
Err_ATC_RLS_BPD1=mean(Err_ATC_RLS_BPD,3);
Err_ATC_RMEE1=mean(Err_RLS_MEE,3);
Err_CTA_RLS1=mean(Err_CTA_RLS,3);
Err_CTA_RMC1=mean(Err_CTA_RMC,3);
Err_CTA_RLS_BPD1=mean(Err_CTA_RLS_BPD,3);
Err_CTA_RMEE1=mean( Err_RLS_CTA_MEE,3);

for cnt=1:num
    Err_TH_MCC(cnt) =mean(Err_ATC_RMC(:,L,cnt));
    Err_TH_CTA_RMC(cnt) =mean(Err_CTA_RMC(:,L,cnt));
    Err_TH_MSE(cnt) = mean(Err_ATC_RLS(:,L,cnt));
    Err_TH_CTA_RLS(cnt) = mean(Err_CTA_RLS(:,L,cnt));
    Err_TH_ATC_RMEE(cnt) = mean(Err_RLS_MEE(:,L,cnt));
    Err_TH_RMEE(cnt) = mean(Err_RLS_CTA_MEE(:,L,cnt));
    Err_TH_RLS_New(cnt) = mean(Err_ATC_RLS_BPD(:,L,cnt));
    Err_TH_CTA_RLS_BPD(cnt) = mean(Err_CTA_RLS_BPD(:,L,cnt));
end


toc
figure(1),hold on;
plot(10* log10(mean(Err_CTA_RLS1)),'m');
plot(10* log10(mean(Err_ATC_RLS1)),'b');
plot(10* log10(mean(Err_CTA_RMC1)));
plot(10* log10(mean(Err_ATC_RMC1)),'r');
plot(10* log10(mean(Err_CTA_RMEE1)),'y');
plot(10* log10(mean(Err_ATC_RMEE1)),'k');
plot(10*log10(mean(Err_CTA_RLS_BPD1)),'g');
plot(10*log10(mean(Err_ATC_RLS_BPD1)),'c');
plot(10*log10(mean(EE2)*ones(1,L)),'--','LineWidth',2);
plot(10*log10(mean(EE1)*ones(1,L)),'--','LineWidth',2);
legend('CTA-RLS','ATC-RLS','CTA-RMCC(\sigma=4)','ATC-RMCC(\sigma=4)','CTA-RMEE(\sigma=4,L=8)','ATC-RMEE(\sigma=4,L=8)','CTA-RMLF(K=2)','ATC-RMLF(K=2)','CTA-RMLF-Theory(K=2)','ATC-RMLF-Theory(K=2)');xlabel('Iterations');ylabel('MSD');
ylim([-28,10]);
box on;grid on;
hold off

figure(2),hold on
plot(10*log10((Err_TH_CTA_RLS)),'-md','LineWidth',1.5);
plot(10*log10((Err_TH_MSE)),'-bh','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RMC)),'-^','LineWidth',1.5);
plot(10*log10((Err_TH_MCC)),'-rv','LineWidth',1.5);
plot(10*log10((Err_TH_RMEE)),'-y>','LineWidth',1.5);
plot(10*log10(( Err_TH_ATC_RMEE)),'-k<','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RLS_BPD)),'-gp','LineWidth',1.5);
plot(10*log10((Err_TH_RLS_New)),'-co','LineWidth',1.5);
plot(EE0,'-s','LineWidth',2);
plot(EE,'-*','LineWidth',2);
legend('CTA-RLS','ATC-RLS','CTA-RMCC(\sigma=4)','ATC-RMCC(\sigma=4)','CTA-RMEE(\sigma=4,L=8)','ATC-RMEE(\sigma=4,L=8)','CTA-RMLF(K=2)','ATC-RMLF(K=2)','CTA-RMLF-Theory(K=2)','ATC-RMLF-Theory(K=2)');
xlabel('Node number,k');ylabel('Steady State MSD (dB)');
ylim([-28,10]);
xlim([1,20]);
box on;grid on;
hold off

