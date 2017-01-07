function [What,bhat] = MTwMVL(yL,tL,xL,tU,xU,Id,Vi,lambda1,lambda3,mu,gamma,p,param)
[T,V] = size(Id);
[N,D] = size(xL);
M = size(xU,1);
if size(xU,2) ~= D || length(yL) ~= N || length(tL) ~= N || length(tU) ~= M || length(Vi) ~= D
    fprintf('> Input feature dimensions mismatch!\n');
    return;
end
Vt = sum(Id');
Tv = sum(Id);
thres = 1e-2;
% bfact = 1/V;
% borig = repmat(bfact,1,V);
bfact = 1/V;
borig = repmat(bfact,T,V);
if param == 1
    for Niter = 1:5
        Lt = [];
        Rt = [];
        for t = 1:T
            
            for v = 1:V
                if Id(t,v) > 0
                    Ytv = yL(tL==t,1);
                    Xtv = xL(tL==t,Vi==v);
                    Utv = xU(tU==t,Vi==v);
                    
                    Ftv = (borig(t,v)^p)*Xtv'*Ytv;
                    Rt  = [Rt; Ftv];
                    
                    Iv = eye(sum(Vi==v));
                    Iv(1,1) = 0;
                    Atv = ((borig(t,v)^p)*Xtv'*Xtv) + lambda1 + (mu*(V-1)*Utv'*Utv) + (gamma*(T-1));
                    Ctv = [];
                    for i = 1:T
                        if i == t
                            for j = 1:V
                                if j == v
                                    Ctv = [Ctv Atv];
                                elseif Id(t,j) > 0
                                    Xtj = xL(tL==t,Vi==j);
                                    Utj = xU(tU==t,Vi==j);
                                    Atj = - mu*Utv'*Utj;
                                    Ctv = [Ctv Atj];
                                end
                            end
                        else
                            for j = 1:V
                                if j == v && Id(i,v) > 0
                                    Ctv = [Ctv -gamma*Iv];
                                elseif Id(i,j) > 0
                                    Ctv = [Ctv zeros(sum(Vi==v),sum(Vi==j))];
                                end
                            end
                        end
                    end
                    Lt = [Lt; Ctv];
                end
            end
        end
        
        Lt = Lt + 1e-10*eye(size(Lt,1));
        Ft = Lt \ Rt;
        k = 0;
        Wt = zeros(D,T);
        for i = 1:T
            for j = 1:V
                if Id(i,j) > 0
                    Nj = sum(Vi==j);
                    What(Vi'==j,i) = Ft(k+1:k+Nj);
                    %                     Wt1{i,j} = Ft(k+1:k+Nj);
                    k = k + Nj;
                end
            end
        end
        
        const_term = lambda3*(T-1);
        for ttt = 1:T
            for vvv = 1:V
                y1 = yL(tL == ttt);
                x1 = xL(tL == ttt, Vi == vvv);
                w1 = What(Vi' == vvv, ttt);
                ierr(ttt,vvv) = ((y1-x1*w1)'*(y1 - x1*w1));
            end
        end
        
        for tt1 = 1:T
            for vv1 = 1:V
                newsum = 0;
                for vv2 = 1:V
                    newsum = newsum + (ierr(tt1,vv1)/ierr(tt1,vv2));
                end
                bhat(tt1,vv1) = 1/newsum;
            end
        end
        
        
        %         for itask = 1:T
        %             for iview = 1:V
        %                 s1 = 0;
        %                 for itask1 = 1:T
        %                     for iview1 = 1:V
        %                         if (itask ~= itask1) && (iview == iview1)
        %                             s1 = s1 + borig(itask1,iview1);
        %                         end
        %                         if (itask == itask1) && (iview == iview1)
        %                             y1 = yL(tL == itask1);
        %                             x1 = xL(tL == itask1, Vi == iview1);
        %                             w1 = What(Vi' == iview1, itask1);
        %                             s2 = (y1-x1*w1)'*(y1 - x1*w1);
        %                         end
        %                     end
        %                 end
        %                 Term2 = beta_update_2(borig,yL,xL,What,lambda3,itask,V,T,tL,Vi);
        %                 borig(itask,iview) = ((lambda3*s1) - Term2)/(s2 + (lambda3*(T-1)));
        %             end
        %         end
        
        Wold = What;
        borig = bhat;
        
    end
    
    
else
    Lt = [];
    Rt = [];
    for t = 1:T
        
        for v = 1:V
            if Id(t,v) > 0
                Ytv = yL(tL==t,1);
                Xtv = xL(tL==t,Vi==v);
                Utv = xU(tU==t,Vi==v);
                
                Ftv = (borig(t,v)^p)*Xtv'*Ytv;
                Rt  = [Rt; Ftv];
                
                Iv = eye(sum(Vi==v));
                Iv(1,1) = 0;
                Atv = ((borig(t,v)^p)*Xtv'*Xtv) + lambda1 + (mu*(V-1)*Utv'*Utv) + (gamma*(T-1));
                Ctv = [];
                for i = 1:T
                    if i == t
                        for j = 1:V
                            if j == v
                                Ctv = [Ctv Atv];
                            elseif Id(t,j) > 0
                                Xtj = xL(tL==t,Vi==j);
                                Utj = xU(tU==t,Vi==j);
                                Atj = - mu*Utv'*Utj;
                                Ctv = [Ctv Atj];
                            end
                        end
                    else
                        for j = 1:V
                            if j == v && Id(i,v) > 0
                                Ctv = [Ctv -gamma*Iv];
                            elseif Id(i,j) > 0
                                Ctv = [Ctv zeros(sum(Vi==v),sum(Vi==j))];
                            end
                        end
                    end
                end
                Lt = [Lt; Ctv];
            end
        end
    end
    
    Lt = Lt + 1e-5*eye(size(Lt,1));
    Ft = Lt \ Rt;
    k = 0;
    Wt = zeros(D,T);
    for i = 1:T
        for j = 1:V
            if Id(i,j) > 0
                Nj = sum(Vi==j);
                What(Vi'==j,i) = Ft(k+1:k+Nj);
                %                 What{i,j} = Ft(k+1:k+Nj);
                k = k + Nj;
            end
        end
    end
    bhat = borig;
end
bhat = borig;