clear all
%a=1;return
disp('jai bhole!!!')
%tic
%

addpath lens_util_github\


% LIBRARIES
addpath ../mtt-master/
addpath ../mtt-master/ntwr' tm decomposition'/
addpath('../mtt-master/libs/poblano_toolbox_1.1');
addpath('../mtt-master/libs/tensor_toolbox_2.5');
addpath('../mtt-master/libs/nway331');
addpath('../mtt-master/libs/parafac2');
addpath ../
addpath ../nips_timeseries_code/

addpath ../data/
addpath ../geant' data'/

% Load traffic matrix  of dimension (PoP_pairs X time-intervals)

%load('x_nscim.mat')         % synthetic;
%X = x_nscim;


%load('X10');           % Abilene
%X10 = aw_od_matrix(X10);
%X = X10;

X=load('G_new1.txt');

%return
% Create traffic tensor


no_pop = sqrt(size(X,2));
no_ti = size(X,1);

X_tens = zeros(no_pop,no_pop,no_ti);



for ii = 1:no_ti
    count1 = 1;
    count2 = no_pop;
    tm = zeros(no_pop,no_pop);
    for jj = 1:no_pop
        %for kk = 1:no_pop
            tm(:,jj) = X(ii,count1:count2);
            count1 = count1 + no_pop;
            count2 = count2 + no_pop;
            
    end
    X_tens(:,:,ii) = tm;
    clear tm
end
clear count1 count2 ii jj

T = tensor(X_tens);%return

%return
X = X';
n_row = no_pop; n_col = no_pop; nt = no_ti; ElemFrac = 1; xx = 0;


sx = size(X);
n  = length(sx);
meanX2 = mean(X(:).^2);
meanX = mean(X(:));

%for Rank_parameter = 2:20

Rank_parameter = 2;



%M_tens = awcur_DropValues(imp_ti,imp_pop,nr,nc,nt,1,ElemFrac,0.02,'elem','ind');


LossRate = [0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8 0.9 0.95 0.98] ;
xx = .75 ;
for lr = 10:10
    %%
    %lr = 10;
    %clear M_tens
    
%for ll = 1:1   
M_tens = aw_DropValues(no_pop,no_pop,no_ti,xx,ElemFrac,LossRate(1,lr),'elem','ind');
M_T = tensor(M_tens);
%save((['AAb_M_50elemsyn_' num2str(LossRate(1,lr)*100) ]),'M_tens'); 
%load('Ab_M_pr08_1.mat');
%M_tens = ones(no_pop,no_pop,no_ti);

M = reshape(M_tens,size(X));


%return
size(find(M==0),1)/(121*no_ti)

ll=lr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CUR
X_m = X.*M;
X_m1 = fill_zeros(X_m');
X_m1 = X_m1';
[nrm,ncm] = size(X);
%[C,indexA,probC] = aw_CUR_ColumnSelect(X_m, Rank_parameter, ncm); % time-intervals

[C1, C1index,ls_c,Dc] = aw_CUR_ColumnSelect(X_m1, Rank_parameter, ncm); % time-intervals
[R1, R1index,ls_r,Dr] = aw_CUR_ColumnSelect(X_m1', Rank_parameter, nrm); % pop pairs
  


    %[C1, C1index,ls_c] = random_ColumnSelect(X,Rank_parameter,ncm,v1(:,1:Rank_parameter));
   % [R1, R1index,ls_r] = random_ColumnSelect(X',Rank_parameter,nrm,v2(:,1:Rank_parameter));
R1 = R1';
U1 = pinv(C1, .05) * (X_m1) * pinv(R1, .05) ;    % Compute U
    
an_tidx = find(ls_r==1);
an_pidx = find(ls_c==1);
    
X_CUR = C1*U1*R1;
%return
X_CUR(M==1) = X(M==1);
X_CUR = max(0,X_CUR);



%nmae_cur  = mean(abs((X(~M)-max(0,X_CUR(~M)))))/mean(abs(X(~M)))
nmae_cur = sum(abs(X_tens(~M_tens)-max(0,X_CUR(~M_tens))))/sum(abs(X_tens(~M_tens)))
nmse_cur = ((sum((abs(X_tens(~M_tens)-max(0,X_CUR(~M_tens)))).^2))/size(X_tens(~M_tens),1))/(sum((abs(X_tens(:))).^2)/(size(X,1)*size(X,2)))
nrmse_cur = sqrt(((sum((abs(X_tens(~M_tens)-max(0,X_CUR(~M_tens)))).^2))/size(X_tens(~M_tens),1)))/(sum((abs(X_tens(:))))/(size(X,1)*size(X,2)))  

    nma_reb_ite(ll,1) = nmae_cur;
    nms_reb_ite(ll,1) = nmse_cur;
    nrm_reb_ite(ll,1) = nrmse_cur;

%% TCUR-AEB 
    
XM_tcur = X_tens.*M_tens;
[A_t1,A_t2,A_t3] = tensor_matricization(XM_tcur);


%A_t31 = fill_zeros(A_t3);

%X_Mtens = zeros(no_pop,no_pop,no_ti);
%for ii = 1:no_ti
%    count1 = 1;
%    count2 = no_pop;
%    tm = zeros(no_pop,no_pop);
%    for jj = 1:no_pop
%        %for kk = 1:no_pop
%            tm(:,jj) = A_t31(ii,count1:count2);
%            count1 = count1 + no_pop;
%            count2 = count2 + no_pop;
%            
%    end
%    X_Mtens(:,:,ii) = tm;
%    clear tm
%end
%clear count1 count2 ii jj

%% selection of fibers and slices

for ii = 1:size(A_t3,1)
    fro_norm(ii,1) = norm(A_t3(ii,:),'fro');
end

[ACv,ACidx] = sort(fro_norm,'descend');
ACidx = ACidx(1:size(C1index,2),:);

for jj = 1:size(A_t3,2)
    euc_norm(jj,1) = norm(A_t3(:,jj),2);
end

[ARv,ARidx] = sort(euc_norm,'descend');
ARidx = ARidx(1:size(R1index,2),:);

C_tcur = XM_tcur(:,:,ACidx);
R_tcur = A_t3(:,ARidx);
R_tcur = R_tcur';

inters_cr = X_m(ARidx,ACidx);    %intersection of C matrix  and R matrix

U_tcur = Dc * pinv(Dr*inters_cr*Dc) * Dr;

Z_tcur = tensor_nmodeproduct(C_tcur,(U_tcur*R_tcur)',3);

Z_tcur(M_tens==1) = X_tens(M_tens==1);
Z_tcur = max(0,Z_tcur);


%nmae_tcur = mean(abs(X_tens(~M_tens)-max(0,Z_tcur(~M_tens))))/mean(abs(X_tens(~M_tens)))
nmae_tcur = sum(abs(X_tens(~M_tens)-max(0,Z_tcur(~M_tens))))/sum(abs(X_tens(~M_tens)))
nmse_tcur = ((sum((abs(X_tens(~M_tens)-max(0,Z_tcur(~M_tens)))).^2))/size(X_tens(~M_tens),1))/(sum((abs(X_tens(:))).^2)/(size(X,1)*size(X,2)))
nrmse_tcur = sqrt(((sum((abs(X_tens(~M_tens)-max(0,Z_tcur(~M_tens)))).^2))/size(X_tens(~M_tens),1)))/(sum((abs(X_tens(:))))/(size(X,1)*size(X,2)))

    nma_aeb_ite(ll,1) = nmae_tcur;
    nms_aeb_ite(ll,1) = nmse_tcur;
    nrm_aeb_ite(ll,1) = nrmse_tcur;
    
        %% Tucker-ALS
    
    tic
    Z_tuc = tucker_als(T.*M_T,Rank_parameter);
    t_e_tuc(ll,1) = toc;
    
    T_tuc = double(Z_tuc);
    T_tuc(M==1) = X(M==1);
    T_tuc = max(0,T_tuc);


    nmae_tuc = sum(abs(X_tens(~M_tens)-max(0,T_tuc(~M_tens))))/sum(abs(X_tens(~M_tens)))
    nmse_tuc = ((sum((abs(X_tens(~M_tens)-max(0,T_tuc(~M_tens)))).^2))/size(X_tens(~M_tens),1))/(sum((abs(X_tens(:))).^2)/(size(X,1)*size(X,2)))
    nrmse_tuc = sqrt(((sum((abs(X_tens(~M_tens)-max(0,T_tuc(~M_tens)))).^2))/size(X_tens(~M_tens),1)))/(sum((abs(X_tens(:))))/(size(X,1)*size(X,2))) 

    nma_tuc_ite(ll,1) = nmae_tuc;
    nms_tuc_ite(ll,1) = nmse_tuc;
    nrm_tuc_ite(ll,1) = nrmse_tuc;
%% CP-ALS
    tic
    Z_cp = cp_als(T.*M_T,Rank_parameter);
    t_e_cp(ll,1) = toc;
    
    T_cp = double(Z_cp);
    T_cp(M==1) = X(M==1);
    T_cp = max(0,T_cp);

   % nmae_cp = mean(abs(X_tens(~M_tens)-max(0,T_cp(~M_tens))))/mean(abs(X_tens(~M_tens)));
    
    nmae_cp = sum(abs(X_tens(~M_tens)-max(0,T_cp(~M_tens))))/sum(abs(X_tens(~M_tens)))
    nmse_cp = ((sum((abs(X_tens(~M_tens)-max(0,T_cp(~M_tens)))).^2))/size(X_tens(~M_tens),1))/(sum((abs(X_tens(:))).^2)/(size(X,1)*size(X,2)))
    nrmse_cp = sqrt(((sum((abs(X_tens(~M_tens)-max(0,T_cp(~M_tens)))).^2))/size(X_tens(~M_tens),1)))/(sum((abs(X_tens(:))))/(size(X,1)*size(X,2)))
    nma_cp_ite(ll,1) = nmae_cp;
    nms_cp_ite(ll,1) = nmse_cp;
    nrm_cp_ite(ll,1) = nrmse_cp;
    %% CPWOPT-ALS
    
    tic
    Z_cpw = cp_wopt(T,M_T,Rank_parameter);
    t_e_cpw(ll,1) = toc;
    
    T_cpw = double(Z_cpw);
    T_cpw(M==1) = X(M==1);
    T_cpw = max(0,T_cpw);

   % nmae_cpw = mean(abs(X_tens(~M_tens)-max(0,T_cpw(~M_tens))))/mean(abs(X_tens(~M_tens)));
    
     
    nmae_cpw = sum(abs(X_tens(~M_tens)-max(0,T_cpw(~M_tens))))/sum(abs(X_tens(~M_tens)))
    nmse_cpw = ((sum((abs(X_tens(~M_tens)-max(0,T_cpw(~M_tens)))).^2))/size(X_tens(~M_tens),1))/(sum((abs(X_tens(:))).^2)/(size(X,1)*size(X,2)))
    nrmse_cpw = sqrt(((sum((abs(X_tens(~M_tens)-max(0,T_cpw(~M_tens)))).^2))/size(X_tens(~M_tens),1)))/(sum((abs(X_tens(:))))/(size(X,1)*size(X,2)))

    nma_cpw_ite(ll,1) = nmae_cpw;
    nms_cpw_ite(ll,1) = nmse_cpw;
    nrm_cpw_ite(ll,1) = nrmse_cpw;
     
    %% STTC
    
    [T1,T2,T3] = tensor_matricization(X_tens);
    [M1,M2,M3] = tensor_matricization(M_tens);

    [A1,b1] = XM2Ab(T1,M1);
    [A2,b2] = XM2Ab(T2,M2);
    [A3,b3] = XM2Ab(T3,M3);

    sx1 = size(T1);
    sx2 = size(T2);
    sx3 = size(T3);
    %K = 1;
    %lambda = 1e-3;
    tic
    Cons1 = aw_ConfigSRTF1(A1,b1,T1,M1,sx1,Rank_parameter);%,K,lambda,true);
    Cons2 = aw_ConfigSRTF1(A2,b2,T2,M2,sx2,Rank_parameter);%,K,lambda,true);
    Cons3 = aw_ConfigSRTF1(A3,b3,T3,M3,sx3,Rank_parameter);%,K,lambda,true);
    %return
    Cons_srtf = {Cons1{1} Cons2{1} Cons3{1}};
    [u5,v5,w5] = SRTF(X_tens,Rank_parameter,M_tens,Cons_srtf);%,10,1e-1,20);
    Z_srtf = tensorprod(u5,v5,w5);
    t_e_srtf(ll,1) = toc;
    
    Z_srtf(M_tens==1) = X_tens(M_tens==1);
    Z_srtf = max(0,Z_srtf);

%    nmae_srtf = mean(abs(X_tens(~M_tens)-max(0,Z_srtf(~M_tens))))/mean(abs(X_tens(~M_tens)))
     
    nmae_srtf = sum(abs(X_tens(~M_tens)-max(0,Z_srtf(~M_tens))))/sum(abs(X_tens(~M_tens)))
    nmse_srtf = ((sum((abs(X_tens(~M_tens)-max(0,Z_srtf(~M_tens)))).^2))/size(X_tens(~M_tens),1))/(sum((abs(X_tens(:))).^2)/(size(X,1)*size(X,2)))
    nrmse_srtf = sqrt(((sum((abs(X_tens(~M_tens)-max(0,Z_srtf(~M_tens)))).^2))/size(X_tens(~M_tens),1)))/(sum((abs(X_tens(:))))/(size(X,1)*size(X,2)))

    
    nma_srtf_ite(ll,1) = nmae_srtf;
    nms_srtf_ite(ll,1) = nmse_srtf;
    nrm_srtf_ite(ll,1) = nrmse_srtf;
    
end    

t_e_final = [t_e_cp t_e_tuc t_e_cpw t_e_srtf];
nma_final = [nma_reb_ite nma_aeb_ite nma_cp_ite nma_tuc_ite nma_cpw_ite nma_srtf_ite]
nms_final = [nms_reb_ite nms_aeb_ite nms_cp_ite nms_tuc_ite nms_cpw_ite nms_srtf_ite]
nrm_final = [nrm_reb_ite nrm_aeb_ite nrm_cp_ite nrm_tuc_ite nrm_cpw_ite nrm_srtf_ite]


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CUR
X_m = X.*M;
X_m1 = fill_zeros(X_m');
X_m1 = X_m1';
[nrm,ncm] = size(X);
%[C,indexA,probC] = aw_CUR_ColumnSelect(X_m, Rank_parameter, ncm); % time-intervals

[C1, C1index,ls_c,Dc] = aw_CUR_ColumnSelect(X_m1, Rank_parameter, ncm); % time-intervals
[R1, R1index,ls_r,Dr] = aw_CUR_ColumnSelect(X_m1', Rank_parameter, nrm); % pop pairs
  


    %[C1, C1index,ls_c] = random_ColumnSelect(X,Rank_parameter,ncm,v1(:,1:Rank_parameter));
   % [R1, R1index,ls_r] = random_ColumnSelect(X',Rank_parameter,nrm,v2(:,1:Rank_parameter));
R1 = R1';
U1 = pinv(C1, .05) * (X_m1) * pinv(R1, .05) ;    % Compute U
    
an_tidx = find(ls_r==1);
an_pidx = find(ls_c==1);
    
X_CUR = C1*U1*R1;
%return
X_CUR(M==1) = X(M==1);
X_CUR = max(0,X_CUR);



nmae_cur  = mean(abs((X(~M)-max(0,X_CUR(~M)))))/mean(abs(X(~M)));



otpt(Rank_parameter,:) = [nmae_cp nmae_cur]
%end

return