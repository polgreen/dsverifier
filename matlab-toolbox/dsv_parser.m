function dsv_parser(ds, type, error)
%read system struct and translate to a ANSI-C file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type, 'tf')
a = cell2mat(ds.system.Denominator);
nA = length(a);
b = cell2mat(ds.system.Numerator);
nB = length(b);
sample_time = ds.system.Ts;
frac_bits = ds.impl.frac_bits;
int_bits = ds.impl.int_bits;
max_range = ds.range.max;
min_range = ds.range.min;
delta = ds.delta;

fid = fopen('input.c', 'wt' );
fprintf(fid,'%s\n\n', '#include <dsverifier.h>');
fprintf(fid,'%s\n\t', 'digital_system ds = { ');
fprintf(fid,'%s %s },\n\t','.b = { ', poly2strc(b));
fprintf(fid,'%s %d,\n\t','.b_size = ', nB);
fprintf(fid,'%s %s },\n\t','.a = { ', poly2strc(a));
fprintf(fid,'%s %d,\n\t','.a_size = ', nA);
fprintf(fid,'%s %d\n','.sample_time =', sample_time);
fprintf(fid,'%s\n\n','};');
fprintf(fid,'%s\n\t','implementation impl = { ');
fprintf(fid,'%s %d,\n\t','.int_bits = ', int_bits);
fprintf(fid,'%s %d,\n\t','.frac_bits =  ', frac_bits);
fprintf(fid,'%s %f,\n\t','.max = ', max_range);
fprintf(fid,'%s %f,\n\t','.min = ', min_range);
if error > 0
fprintf(fid,'%s %f,\n\t','.max_error = ', error);
end
fprintf(fid,'%s %f\n','.delta =', delta);
fprintf(fid,'%s\n\n','};');
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(type, 'cl')
ac = cell2mat(ds.controller.Denominator);
nAc = length(ac);
bc = cell2mat(ds.controller.Numerator);
nBc = length(bc);

ap = cell2mat(ds.plant.Denominator);
nAp = length(ap);
bp = cell2mat(ds.plant.Numerator);
nBp = length(bp);

sample_time = ds.controller.Ts;
frac_bits = ds.impl.frac_bits;
int_bits = ds.impl.int_bits;
max_range = ds.range.max;
min_range = ds.range.min;
delta = ds.delta;

fid = fopen('input.c', 'wt' );
fprintf(fid,'%s\n\n', '#include <dsverifier.h>');
fprintf(fid,'%s\n\t', 'digital_system controller = { ');
fprintf(fid,'%s %s },\n\t','.b = { ', poly2strc(bc));
fprintf(fid,'%s %d,\n\t','.b_size = ', nBc);
fprintf(fid,'%s %s },\n\t','.a = { ', poly2strc(ac));
fprintf(fid,'%s %d,\n\t','.a_size = ', nAc);
fprintf(fid,'%s %d\n','.sample_time =', sample_time);
fprintf(fid,'%s\n\n','};');
fprintf(fid,'%s\n\t','implementation impl = { ');
fprintf(fid,'%s %d,\n\t','.int_bits = ', int_bits);
fprintf(fid,'%s %d,\n\t','.frac_bits =  ', frac_bits);
fprintf(fid,'%s %f,\n\t','.max = ', max_range);
fprintf(fid,'%s %f,\n\t','.min = ', min_range);
if error > 0
fprintf(fid,'%s %f,\n\t','.max_error = ', error);
end
fprintf(fid,'%s %f\n','.delta =', delta);
fprintf(fid,'%s\n\n','};');
fprintf(fid,'%s\n\t', 'digital_system plant = { ');
fprintf(fid,'%s %s },\n\t','.b = { ', poly2strc(bp));
fprintf(fid,'%s %d,\n\t','.b_size = ', nBp);
fprintf(fid,'%s %s },\n\t','.a = { ', poly2strc(ap));
fprintf(fid,'%s %d \n\t','.a_size = ', nAp);
fprintf(fid,'%s\n\n','};');
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif strcmp(type, 'ss')
    
end

end
