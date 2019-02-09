function h = ne_batch_plot_NE_column_width(files)

allwidths = cell(length(files),1);

for i = 1:length(files)
    load(files{i})
    
    allwidths{i} = ne_calc_NE_column_width(exp_site_nedata);   
  
end

allwidthsmat = cell2mat(allwidths);
figure;
h = histogram(categorical(allwidthsmat), 'Normalization', 'probability');

tickpref;

x = xlim;
y = ylim;

text(x(2)/4, y(2)/4*3, sprintf('n = %d', length(allwidthsmat)))

xlabel('Width of NEs (\mum)')
ylabel('Ratio')

print_mfilename(mfilename);

