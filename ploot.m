clear; close('all');

x = categorical({'10^4'});
x = reordercats(x,{'10^4'});
vals = [0.0649 0.0636 0.1812];

subplot(1,3,1)
bar(x, vals);



ylabel('Time (sec)') 
x1 = categorical({'10^6'});
x1 = reordercats(x1,{'10^6'});

x2 = categorical({'10^8'});
x2 = reordercats(x2,{'10^8'});

vals1 = [ 7.2987 7.0723 7.1040];

vals2 = [ 1910.6725 1741.2185 1743.4293];
subplot(1,3,2)
bar(x1, vals1);

xlabel('Elements') 


subplot(1,3,3)
bar(x2, vals2);


legend({'MPI block','MPI non-block', 'Coarray'},'Location','southwest')

