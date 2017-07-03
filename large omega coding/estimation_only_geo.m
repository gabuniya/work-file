%AIC and BIC criterion
clear all;
close all;

%Estimate the geo-industry model

%Activate Parallel computing toolbox
% session = com.mathworks.toolbox.distcomp.pmode.SessionFactory.getCurrentSession;
% if isempty( session ) || ~session.isSessionRunning()
%     matlabpool open local; % Requires manual configuration.
% end


afosdijfrosgjrjgreg

lsfjrodgjrwogjreogjer

fjdogrjwegorwejgorew


tic;
%%Data

%[J D constraint_D number_of_firms  number_of_periods tempR number_of_regions number_of_countries number_of_industries country_name country_code code_rest_of industry_name index_country index_rest_of firms_exposures firms_exposures_rest_of firms_exposures_ind firms_exposures_rest_of_ind firms_homes dictionary dictionary_rest_of dictionary_industry dictionary_industry_rest_of R_std_min R_std_min_ind X Y B Omega omega tempD SS isins iter_i index_industry_rest_of number_of_zones_total index_industry_start]=load_data_2001;

ind_index=1:10;
industry_dummy=0;

%load preliminary_cleaned_to_change_no_problematic_stocks
%load preliminary_cleaned_cleaned_usual
%load preliminary_check_with_IL
%load preliminary_more_africa
load prl_iminary_a

while number_of_countries>number_of_regions

Mean_squared_gradient=[];
bad_stocks=[];
noone_exit_flag_stocks=[];
    
%iter_i=1;
MemLL = [];
options = optimset('Diagnostics','off','Display','off','Algorithm','sqp','UseParallel','never', 'GradObj', 'on','TolX',1e-12);
m_grad=10;
toc;

ii=1;
creterion_1=100;
i_1=1;
i_2=1;
difference_1=zeros(10,1);
difference_2=zeros(10,1);
% iter_i=1;
%iter_j=1;
ttt=[];
while m_grad>=0.3705
    
    tempB=zeros(number_of_firms,number_of_countries);
    tempD=zeros(number_of_firms,1);
    gradB=zeros(number_of_firms,number_of_countries);
    gradB0=zeros(number_of_firms,number_of_countries);
    gradB1=zeros(number_of_firms,number_of_countries);
    gradB2=zeros(number_of_firms,number_of_countries);
    gradB3=zeros(number_of_firms,1);
    lastexitflag=zeros(number_of_firms,1);
    tic;

    tic;
    
    ddum=zeros(number_of_firms,1);
    parfor firm_i=1:number_of_firms  
    
    [constraint_coefficients constraint_rhs coefficient_to_add number_of_known] = ...
                find_constraints_change(firms_exposures,firms_exposures_rest_of,dictionary,dictionary_rest_of,firm_i,number_of_countries, number_of_regions, firms_homes, index_country, index_rest_of);

            %if number_of_industries==10
                [constraint_coeff exposure coeff_to_add] =find_constraints_change_ind(firms_exposures_ind,firms_exposures_rest_of_ind,dictionary_industry,dictionary_industry_rest_of,firm_i,number_of_industries);
            %else
             %   [constraint_coeff exposure coeff_to_add] =find_constraints_change_industry(firms_exposures_ind,firms_exposures_rest_of_ind,dictionary_industry,dictionary_industry_rest_of,firm_i,number_of_industries,number_of_rest_of_industries,index_industry,index_industry_rest_of);
            %end                  

%     constraint_coefficients=[constraint_coefficients zeros(size(constraint_coefficients,1),size(constraint_coeff,2)); zeros(size(constraint_coeff,1),size(constraint_coefficients,2)) constraint_coeff];
%     constraint_rhs=[constraint_rhs; exposure];
%     coefficient_to_add=[coefficient_to_add; coeff_to_add];
    

%[A_ineq,b_ineq,constraint_coefficients,constraint_rhs]=transf_the_constr(constraint_coefficients,constraint_rhs,constraint_coeff,exposure,coefficient_to_add,coeff_to_add,ind_index,index_rest_of,number_of_industries,number_of_countries,industry_dummy);

%[A_ineq b_ineq constraint_coefficients constraint_rhs]=transf_the_geo_constr(constraint_coefficients,constraint_rhs,coefficient_to_add,1,number_of_countries);

A_ineq=-eye(size(omega,1));
b_ineq=coefficient_to_add;

% %Remove this part
% A_ineq(end,:)=[];
% b_ineq(end,:)=[];
% 
% constraint_coefficients=[constraint_coefficients zeros(size(constraint_coefficients,1),number_of_industries)];
% constraint_coefficients=[constraint_coefficients; [zeros(size(constraint_coeff,1),number_of_countries-1) constraint_coeff]];
% constraint_rhs=[constraint_rhs; exposure];

% A_ineq=[A_ineq zeros(size(A_ineq,1),number_of_industries)];
% A_ineq=[A_ineq; [zeros(number_of_industries,number_of_countries-1) -eye(number_of_industries)]];
% b_ineq=[b_ineq; coeff_to_add];
%

%     constraint_coefficients=[constraint_coefficients zeros(size(constraint_coefficients,1),size(constraint_coeff,2)); zeros(size(constraint_coeff,1),size(constraint_coefficients,2)) constraint_coeff];
%     constraint_rhs=[constraint_rhs; exposure];
%  %   coefficient_to_add=[coefficient_to_add; coeff_to_add];
%  
%  A_ineq=[A_ineq zeros(size(A_ineq,1),number_of_industries)];
%  A_ineq=[A_ineq; zeros(number_of_industries,number_of_countries-1) -eye(number_of_industries)];
%  b_ineq=[b_ineq; zeros(number_of_industries,1)];
            
    if constraint_D(firm_i)==0
                    [tB, tD, lfl, ~, lambda, grad] = ...
                fmincon(@(B)individual_firm_LLvalue(B,SS(firm_i),X(firm_i,:),Y,omega),B(firm_i,:),A_ineq,b_ineq,full(constraint_coefficients),full(constraint_rhs),[],[],[],options);

            lastexitflag(firm_i)=lfl;
            tempB(firm_i,:)=tB;
            tempD(firm_i)=tD;
            gradB1(firm_i,:)=grad;
            gradB2(firm_i,:)=(number_of_periods/2)*gradB1(firm_i,:)/tempD(firm_i);
    end
    
    
    if constraint_D(firm_i)>0
        
%         constraint_coefficients=[constraint_coefficients zeros(size(constraint_coefficients,1),1)];
%         constraint_coefficients=[constraint_coefficients; zeros(1,size(constraint_coefficients,2))];
%         constraint_coefficients(end,end)=1;
%         
%         
%         param_i=[B(firm_i,:) D(firm_i)];
%     
%                  
%         [tparam, ~, lfl, ~, lambda, grad] = fmincon(@(param)individual_firm_LLvalue_a(number_of_periods,param,SS(firm_i),X(firm_i,:),Y,omega),param_i,-eye(number_of_countries+1),[coefficient_to_add; 0],full(constraint_coefficients),[full(constraint_rhs); constraint_D(firm_i)],[],[],[],options);
% 
%             tB=tparam(1:end-1);
%             tD=tparam(end);
%             
%             lastexitflag(firm_i)=lfl;
%             tempB(firm_i,:)=tB;
%             tempD(firm_i)=tD;
%             gradB1(firm_i,:)=grad(1:end-1);
%             gradB2(firm_i,:)=gradB1(firm_i,:);
%             gradB3(firm_i)=grad(end);
    end
            
            [msgstr, msgid] = lastwarn;
            k=strmatch('Rank deficient',msgstr);

            if isempty(k)==0
                ddum(firm_i)=1;
            end
            
            lastwarn('');
    end
    toc;
    tic;
        bad_stocks=[bad_stocks; isins(lastexitflag==-2)];
        noone_exit_flag_stocks=[noone_exit_flag_stocks; isins(lastexitflag~=1)];
    
    B=tempB;
    clear tempB
    B(B<0)=0;
    
   
    dum_minus_two=double(lastexitflag==-2);
    if sum(dum_minus_two)>0
        [J constraint_D SS B tempD X number_of_firms firms_exposures firms_exposures_ind firms_exposures_rest_of firms_exposures_rest_of_ind firms_homes gradB2 gradB3 tempR isins]=...
            clear_bad_firms(J,constraint_D,SS,B,tempD,X,lastexitflag,firms_exposures,firms_exposures_ind,firms_exposures_rest_of,firms_exposures_rest_of_ind,firms_homes,gradB2,gradB3,tempR,isins);
    end
    
         
    Omega=omega*Y*omega';
    [U,S,V] = svd(Omega);
    S1=sqrt(S);
    omega=U*S1;
   
    
    [~,S,~] = svd(omega);
    ttt=[ttt; S(end,end)/S(1,1)];
    D=SS-2*sum(X*omega'.*B,2)+sum(B*omega*Y.*(B*omega),2);
    %D=tempD;
    
    
%     [grad grad_1 grad_2 grad_3]=gradd(J,X,B,D,omega,tempR,gradB2,gradB3,number_of_periods);
%      grads=grad.^2;
%      m_grad=max(grads);
     
     [X Y]=Iterate(B,tempD,omega,tempR,J);
     
     [grad grad_1 grad_2 grad_3]=gradd(J,X,B,D,omega,tempR,gradB2,gradB3,number_of_periods);
     grads=grad.^2;
     m_grad=max(grads);
     

     %Calculate the value of the likelihood close to convergence point
 %if m_grad<0.5
     L=LL_new(B,omega,D,tempR);
         
      MemLL = [MemLL; L];
      disp(L);
 %end

    
    iter_i=iter_i+1;
    
    if ii==50
        close all;
        ii=1;
        diff_a=abs(dynamcs(end).o-dynamcs(end-1).o);
        creterion_1=max(max(diff_a))/min(min(diff_a));
        k=(diff_a==max(max(diff_a)));
        i_1=find(sum(k));
        i_2=find(k(:,i_1));
        difference_1=zeros(10,1);
        difference_2=zeros(10,1);
        for jj=40:49
            difference_1(jj-39)=((dynamcs(jj).o(i_2,i_1)-dynamcs(jj-1).o(i_2,i_1)))/dynamcs(jj-1).o(i_2,i_1);
            difference_2(jj-39)=abs(dynamcs(jj).o(i_2,i_1)/sqrt(dynamcs(jj).o(i_2,i_2)*dynamcs(jj).o(i_1,i_1))-dynamcs(jj-1).o(i_2,i_1)/sqrt(dynamcs(jj-1).o(i_2,i_2)*dynamcs(jj-1).o(i_1,i_1)));
        end
        figure(1);
        %subplot(2,1,1);
        plot(difference_1);
        %subplot(2,1,2);
        %plot(difference_2)
        clear dynamcs
        save prl_iminary_a J constraint_D D index_industry_rest_of number_of_firms number_of_periods tempR number_of_regions number_of_countries number_of_industries country_name country_code code_rest_of industry_name index_country index_rest_of firms_exposures firms_exposures_rest_of firms_exposures_ind firms_exposures_rest_of_ind firms_homes dictionary dictionary_rest_of dictionary_industry dictionary_industry_rest_of R_std_min R_std_min_ind X Y B Omega omega tempD SS isins iter_i -v6
    end
    
   D=tempD;
   
   dynamcs(ii).B=B;
   dynamcs(ii).D=D;
   dynamcs(ii).o=omega;
   
   ii=ii+1;
   
   disp([ii])
   
   
    
    
    disp(iter_i);
    disp(m_grad);
    disp([i_2 i_1]);
    disp([creterion_1]);
    disp([mean(grads)]);
    disp([mean(grad_3.^2)]);
    
      Mean_squared_gradient=[Mean_squared_gradient m_grad];

    toc;
    
%        %if number_of_countries>33
%         if rem(iter_i,250)==0
%            dd=rank(omega);
%            if dd<number_of_countries
%         
%         [firm_individual_exposure firm_individual_constraint firm_group_constraint]=study_B_new(firms_exposures,firms_exposures_rest_of,dictionary,dictionary_rest_of,firms_homes,number_of_firms,number_of_countries);
%         number_no_n_zeros_exposures=sum(double(firm_individual_exposure'>0.05));
%         number_dom_stocks=sum(double(firm_individual_exposure'>0.99));
%         inn_ex=strmatch('''IE''',country_code,'exact');
%         number_no_n_zeros_exposures(inn_ex)=1000;
%         number_dom_stocks(inn_ex)=1000;
% 
% % numb_in_fr=zeros(size(B,2),1);
% % mean_value_in_fr=zeros(size(B,2),1);
% % for i=1:size(B,2)
% %     [number_in_front,mean_value_in_front]=numb_mult_in_front(B,omega,D,tempR,i);
% %     numb_in_fr(i)=number_in_front;
% %     mean_value_in_fr(i)=mean_value_in_front;
% % end
%         ind_delete=find(number_no_n_zeros_exposures<=prctile(number_no_n_zeros_exposures,40));
% ind_delete=find(numb_in_fr>=prctile(numb_in_fr,60));
% nnumb_in_fr=numb_in_fr(ind_delete);
%         number_dom_stocks=number_dom_stocks(ind_delete);
%         c_z_d=country_code(ind_delete);
%         ind_excld=find(ismember(c_z_d,code_rest_of));
%         c_z_d(ind_excld)=[];
%         number_dom_stocks(ind_excld)=[];
%         ind_delete(ind_excld)=[];
%         nnumb_in_fr(ind_excld)=[];
%         ind_delete=ind_delete(nnumb_in_fr==max(nnumb_in_fr));
%         ind_delete=ind_delete(1);
%         country_zones_to_delete=country_code(ind_delete);
%         
%         [dictionary dictionary_rest_of country_code country_name]=modify_dictionary_new(dictionary,dictionary_rest_of,country_code,country_name,country_zones_to_delete,index_rest_of,number_of_regions,code_rest_of);
% 
%         [index_country index_rest_of]=g_i_r_new_rf(country_name);
%         code_rest_of=country_code(index_rest_of);
%         number_of_countries=size(dictionary,2);
%  
%         omega(ind_delete,:)=[];
%         omega(:,ind_delete)=[];
%         N=number_of_countries;
%         disp([number_of_countries]);
%         disp([number_of_industries]);
%         B=ones(number_of_firms,(N))/(number_of_countries);
%         tempD=ones(number_of_firms,1);
% 
%         element2=bsxfun(@rdivide,B*omega,tempD);
%         interim=eye(size(B,2))+(B*omega)'*element2;
%         interim=(interim+interim')/2;
% 
%         ell1=tempR'*(element2/interim);
%         X=(1/size(tempR,2))*tempR*ell1;
%         Y=interim\(element2'*X+eye(size(interim,1)));
%         Y=(Y+Y')/2;
% 
%         iter_i=1;
%         ii=1;
%             
%         R_std=std(tempR');
%         [firm_individual_exposure firm_individual_constraint firm_group_constraint]=study_B_new(firms_exposures,firms_exposures_rest_of,dictionary,dictionary_rest_of,firms_homes,number_of_firms,number_of_countries);
%         firm_individual_exposure=firm_individual_exposure';
%         mmax=max(firm_individual_exposure);
%         R_std_min=zeros(number_of_countries,1);
%         for i=1:number_of_countries
%             d=double(firm_individual_exposure(:,i)==max(firm_individual_exposure(:,i)));
%             R_std_min(i)=min(R_std(d==1));
%         end
%         
%         
% 
%            end
%          end
        if rem(iter_i,500)==0
           dd=rank(omega);
           if dd<number_of_countries
        
        [firm_individual_exposure firm_individual_constraint firm_group_constraint]=study_B_new(firms_exposures,firms_exposures_rest_of,dictionary,dictionary_rest_of,firms_homes,number_of_firms,number_of_countries);
        number_no_n_zeros_exposures=sum(double(firm_individual_exposure'>0.05));
        number_dom_stocks=sum(double(firm_individual_exposure'>0.99));

%         inn_ex=strmatch('''IL''',country_code,'exact');
%         number_no_n_zeros_exposures(inn_ex)=1000;
%         number_dom_stocks(inn_ex)=1000;

        ind_delete=find(number_no_n_zeros_exposures<=prctile(number_no_n_zeros_exposures,40));
        number_no_n_zeros_exposures=number_no_n_zeros_exposures(ind_delete);
      number_dom_stocks=number_dom_stocks(ind_delete);
        c_z_d=country_code(ind_delete);
        ind_excld=find(ismember(c_z_d,code_rest_of));
        c_z_d(ind_excld)=[];
        number_no_n_zeros_exposures(ind_excld)=[];
       number_dom_stocks(ind_excld)=[];
        ind_delete(ind_excld)=[];
      ind_delete=ind_delete(number_dom_stocks==min(number_dom_stocks));
%        ind_delete=ind_delete(number_no_n_zeros_exposures==min(number_no_n_zeros_exposures));
        ind_delete=ind_delete(1);
        country_zones_to_delete=country_code(ind_delete);

        
        [dictionary dictionary_rest_of country_code country_name]=modify_dictionary_new(dictionary,dictionary_rest_of,country_code,country_name,country_zones_to_delete,index_rest_of,number_of_regions,code_rest_of);

        [index_country index_rest_of]=g_i_r_new_rf(country_name);
        code_rest_of=country_code(index_rest_of);
        number_of_countries=size(dictionary,2);
 
        omega(ind_delete,:)=[];
        omega(:,ind_delete)=[];
        N=number_of_countries;
        disp([number_of_countries]);
        disp([number_of_industries]);
        B=ones(number_of_firms,(N))/(number_of_countries);
        tempD=ones(number_of_firms,1);

        element2=bsxfun(@rdivide,B*omega,tempD);
        interim=eye(size(B,2))+(B*omega)'*element2;
        interim=(interim+interim')/2;

        ell1=tempR'*(element2/interim);
        X=(1/size(tempR,2))*tempR*ell1;
        Y=interim\(element2'*X+eye(size(interim,1)));
        Y=(Y+Y')/2;

        iter_i=1;
        ii=1;
            
        R_std=std(tempR');
        [firm_individual_exposure firm_individual_constraint firm_group_constraint]=study_B_new(firms_exposures,firms_exposures_rest_of,dictionary,dictionary_rest_of,firms_homes,number_of_firms,number_of_countries);
        firm_individual_exposure=firm_individual_exposure';
        mmax=max(firm_individual_exposure);
        R_std_min=zeros(number_of_countries,1);
        for i=1:number_of_countries
            d=double(firm_individual_exposure(:,i)==max(firm_individual_exposure(:,i)));
            R_std_min(i)=min(R_std(d==1));
        end
        
        

           end
         end


          


            
end

        if rank(omega)==number_of_countries
            break;
        end

end

%save output_2000 number_of_firms  number_of_periods tempR number_of_regions number_of_countries number_of_industries country_name country_code code_rest_of industry_name index_country index_rest_of firms_exposures firms_exposures_rest_of firms_exposures_ind firms_exposures_rest_of_ind firms_homes dictionary dictionary_rest_of dictionary_industry dictionary_industry_rest_of R_std_min R_std_min_ind X Y B omega tempD SS isins num txt iter_i -v6
  
    
        
        
            
            
             
    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
 