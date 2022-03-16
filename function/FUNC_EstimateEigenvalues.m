function [EigVal,HA,r] = FUNC_EstimateEigenvalues(ax,in_HA,in_r,Zinv,dZ)
% Estimate first layer Eigenvalues based on first layer horizontal anisotropy
    HA(1,1) = in_HA(1);
    plot(ax{10},HA(1),Zinv(1),'.-r','LineWidth',2,'MarkerSize',15);
    EigVal(1,:) = EV1stLayer(in_HA(1)); % 1st layer Eigenvalues
    plot(ax{9},EigVal(:,1),Zinv(1),'.-r','MarkerFaceColor','r','LineWidth',2,'MarkerSize',15);
    plot(ax{9},EigVal(:,2),Zinv(1),'.-b','MarkerFaceColor','b','LineWidth',2,'MarkerSize',15);
    plot(ax{9},EigVal(:,3),Zinv(1),'.-g','MarkerFaceColor','g','LineWidth',2,'MarkerSize',15); 
    drawnow
% Strat estimating Eigenvalues downward (layer by layer)
    for i = 1:length(Zinv)-1
        l1_T = EigVal(i,1); % Lambda1_TopLayer
        HA_T = in_HA(i);    % DeltaLambda_TopLayer
        r_T = in_r(i);      % ReflectionRatio_TopLayer
        HA_B = in_HA(i+1);  % DeltaLambda_BottomLayer
        lmts = [-5e-4 1.5e-3 1e-6];
        range_HA_B = [HA_B-0.05:0.005:HA_B+0.05];
        range_HA_B(range_HA_B<0) = [];
        range_HA_B(range_HA_B>0.5) = [];
        range_r_T = [r_T-1:0.01:r_T+1];
        comb_HA_r = combvec(range_HA_B,range_r_T)';
        l1 = EigVal(:,1);
        l2 = EigVal(:,2);
        l3 = EigVal(:,3);
        disp("Vector Size: "+string(size(comb_HA_r,1)))
        stp = 0;
        while true
        tic
            res1 = nan(size(comb_HA_r,1),14);
            res2 = nan(size(comb_HA_r,1),14);
            parfor j = 1:size(comb_HA_r,1)
                psbl_HA_B = comb_HA_r(j,1);
                psbl_r_T = comb_HA_r(j,2);
                r_ratio_T = 10.^(psbl_r_T./20); % conver from dB to ratio
                temp1 = func_res1(HA_B,r_T,l3,dZ,l1_T,HA_T,psbl_HA_B,r_ratio_T,psbl_r_T,i,lmts);
                if ~isnan(temp1)
                    res1(j,:) = temp1;
                end
%                 temp2 = func_res1(HA_B,r_T,l3,dZ,l1_T,HA_T,psbl_HA_B,-r_ratio_T,psbl_r_T,i);
%                 if ~isnan(temp2)
%                     res2(j,:) = temp2;
%                 end 
            end
            res1(isnan(res1(:,1)),:) = [];
            res2(isnan(res2(:,1)),:) = [];
            res3 = [res1;res2];
            if isempty(res3)
                range_HA_old = range_HA_B;
                range_r_old = range_r_T;
                disp('empty, try new limits')
                
                range_HA_B = [min(range_HA_B)-0.025:0.005:max(range_HA_B)+0.025];
                range_HA_B(range_HA_B<HA_B-0.1) = [];
                range_HA_B(range_HA_B>HA_B+0.1) = [];
                range_HA_B(range_HA_B<0) = [];
                range_HA_B(range_HA_B>0.5) = [];
                
                range_r_T = [min(range_r_T)-1:0.1:max(range_r_T)+1];
                range_r_T(range_r_T<r_T-10) = [];
                range_r_T(range_r_T>r_T+10) = [];
                
                comb_HA_r = combvec(range_HA_B,range_r_T)';
                
                if length(range_HA_old) == length(range_HA_B)
                    if sum(range_HA_old == range_HA_B) == length(range_HA_B)
                        disp("HA limits")
                        stp = stp+1;
                    end
                end
                if length(range_r_old) == length(range_r_T)
                    if sum(range_r_old == range_r_T) == length(range_r_T)
                        disp("r limits")
                        stp = stp+1;
                    end
                end
                if stp >= 2
                    disp('*****SHIT****** stretching the limits')
                    beep
                    oldlmts = lmts;
                    lmts = [oldlmts(1)+(oldlmts(1)*0.1) ...
                                oldlmts(2)+(oldlmts(2)*0.1)...
                                    oldlmts(3)];
                    disp(string(lmts))
                end
            else
                break
            end
        toc
        end
        disp('-------------------------------next layer')
        G2 = res3(:,end-1);
        Dist = res3(:,end);
        Dist = rescale(Dist,0,100);
        Aall = 0+Dist;
        res3 = [res3 , Aall];
        res3 = sortrows(res3,size(res3,2));
        %------ Get dE
        HA(i+1,1) = res3(1,4);
        plot(ax{10},HA(i+1),Zinv(i+1),'.-r','LineWidth',2,'MarkerSize',15);
        %------ Get r
        r(i,1) = res3(1,6);
        plot(ax{12},r(i),Zinv(i),'.-r','LineWidth',2,'MarkerSize',15);
        %------ Get EigenValues
        EigVal(i+1,:) = res3(1,8:10);
        plot(ax{9},EigVal(i+1,1),Zinv(i+1),'.-r','MarkerFaceColor','r','LineWidth',2,'MarkerSize',15);
        plot(ax{9},EigVal(i+1,2),Zinv(i+1),'.-b','MarkerFaceColor','b','LineWidth',2,'MarkerSize',15);
        plot(ax{9},EigVal(i+1,3),Zinv(i+1),'.-g','MarkerFaceColor','g','LineWidth',2,'MarkerSize',15); 
        drawnow
    end
    r(end+1) = in_r(end);
    plot(ax{12},r(end),Zinv(end),'.r','LineWidth',2,'MarkerSize',14) 
    drawnow
end
%% Estimate first layer Eigenvalues
function EigVal1 = EV1stLayer(HA1)
    l1 = 0.333; 
    k = 1;
    while true
        l2 = HA1 + l1;
        l3 = 1-l2-l1;
        ch(1) = l1 >= 0;
        ch(2) = l2>= 0;
        ch(3) = l3>= 0;
        ch(4) = l1 < 0.334;
        ch(5) = l2 <= 0.5;
        ch(6) = l3 > 0.332;
        ch(7) = l3 <= 1;
        ch(8) = l2 >= l1;
        ch(9) = l3 >= l2;
        if sum(ch) == length(ch)
            break
        elseif l1 <= 0
            fprintf("Error, Couldn't find a solution \n")
            break
        else
            l1 = l1-0.00001;
        end
        k = k+1;
    end
    EigVal1(1,:) = sort([l1,l2,l3]);
end
%% Solve reflection coefficient equation to get Eigenvalues
function res = func_res1(HA_B,r_T,l3,dZ,l1_T,HA_T,psbl_HA_B,r_ratio_T,psbl_r_T,i,lmts)
    [EV,ErrorF] = dEr2EV(l1_T,HA_T,psbl_HA_B,r_ratio_T);
    res = nan;
    if ErrorF == 0
        l1_B = EV(1);
        l2_B = EV(2);
        l3_B = EV(3);
        LS = (l3_B-l3(end))/dZ;
        LS(LS < lmts(1)) = nan;
        LS(LS > lmts(2)) = nan;
        LS(abs(LS) < lmts(3)) = nan;
        if ~isnan(LS)
            GR1 = diff([l3 ; l3_B])./dZ;
            GR2 = diff(GR1)./dZ;
            if i > 1
                nrmG2 = norm(GR2);
            else
                nrmG2 = norm(GR1);
            end
            dist = sqrt( (HA_B-psbl_HA_B).^2 + (r_T-psbl_r_T).^2 );  

            res = [l1_T HA_T HA_B psbl_HA_B r_T psbl_r_T r_ratio_T EV ErrorF LS nrmG2 dist];
        end
    end
end
% using: l1 top, dl top, dl bottom, r top
function [EigVal,ErrorF] = dEr2EV(l1_T,HA_T,HA_B,r_ratio_T)   
    l1_B = l1_T - ((HA_T-HA_B)./(r_ratio_T-1));
    l2_B = HA_B + l1_B;
    l3_B = 1-l2_B-l1_B;

    ch(1) = l1_B >= 0;
    ch(2) = l2_B >= 0;
    ch(3) = l3_B >= 0;
    ch(4) = l1_B < 0.34;
    ch(5) = l2_B <= 0.5;
    ch(6) = l3_B > 0.32;
    ch(7) = l3_B <= 1;
    ch(8) = l2_B >= l1_B;
    ch(9) = l3_B >= l2_B;
    
    sch = sum(ch);
    ErrorF = length(ch)-sch;
    if ErrorF == 0
        EigVal = [l1_B l2_B l3_B];
    else
        EigVal = [nan nan nan];
    end
end