function SignCheck = FUNC_CoherencePhaseSignFit(obs,mdl)
    obs(obs<0) = -1;
    obs(obs>=0) = 1;
    mdl(mdl<0) = -1;
    mdl(mdl>=0) = 1;  
    
%     obs = round(obs,1);
%     mdl = round(mdl,1);

    A = obs == mdl;
    
    a = sum(A,2);
    [n,m] = size(obs);
    acc_prc = round((a).*(100)./(m),1); 
    SignCheck = mean(acc_prc);
    
%     i1 = 1;
%     for i = 1:size(MHB,1)
%         [~,i2] = min(abs(MHB(i,1)-Zmdl));
%         SignCheck(i,1) = mean(acc_prc(i1:i2));
%         i1 = i2+1;
%     end 
    
%     figure,subplot(2,1,1),imagesc(0:179,Zmdl,obs),subplot(2,1,2),imagesc(0:179,Zmdl,mdl),drawnow
%     pause(1),
%     figure,subplot(3,1,1),imagesc(obs),subplot(3,1,2),imagesc(mdl)
%     subplot(3,1,3),imagesc(A)
%     close(gcf)