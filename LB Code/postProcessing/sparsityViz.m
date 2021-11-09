hold off;
for i = 1:length(order)
    figure;
    spy(A{i});
    s = ['Sparsity of Order ',num2str(order(i)),' Matrix for Carlemann Linearized D1Q3 BGK Collsion Term'];
    title(s);
    lbl = ['# of Nonzero Entries = ', int2str(nnz(eval(A{i}))),' # of Variables = ',int2str(length(V{i}))];
    xlabel(lbl);
    print([num2str(i),'.png'],'-dpng');
end

close all;