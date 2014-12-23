
outputdir = @(ref,RE) strcat('../results/1.4.0/karman/lqr_assemble/recursive_bisection/ref_',num2str(ref),'/RE_',num2str(RE),'/');
MATmtx = @(ref,RE,name) strcat(outputdir(ref,RE),name);



for ref = 1:1
    for RE= [2,3,4,5,10,20,50,100,200]
    %for RE = [10,20,50,100,200]
        odir = outputdir(ref,RE);
        if exist(odir)
            fprintf('Time %s Bernoulli Feedback RE=%d ref=%d \n',datestr(now,'HH:MM:SS'),RE,ref);
            Mmtx = MATmtx(ref,RE,'M.mtx');
            Smtx = MATmtx(ref,RE,'S.mtx');
            Rmtx = MATmtx(ref,RE,'R.mtx');
            Kmtx = MATmtx(ref,RE,'K.mtx');
            M_lowermtx = MATmtx(ref,RE,'M_lower.mtx');
            M_uppermtx = MATmtx(ref,RE,'M_upper.mtx');
            Gmtx = MATmtx(ref,RE,'G.mtx');
            Bmtx = MATmtx(ref,RE,'B.mtx');
            Cmtx = MATmtx(ref,RE,'C.mtx');
            Feed0mtx = MATmtx(ref,RE,'Feed0.mtx');
            Feed1mtx = MATmtx(ref,RE,'Feed1.mtx');
            [Feed0,Feed1] = initial_feedback_karman(Mmtx,Smtx,Rmtx,Kmtx,M_lowermtx,M_uppermtx,Gmtx,Bmtx,Cmtx);
            if all(size(Feed0)>0) && all(size(Feed1)>0)
                mmwrite(Feed0mtx,Feed0)
                mmwrite(Feed1mtx,Feed1)
            end
            fprintf('---------------------------------\n');
        end
    end
end


