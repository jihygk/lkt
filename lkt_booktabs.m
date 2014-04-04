%% Output booktabs
% jhk 2013-07-01






function [] = lkt_booktabs(x,s,str,lkt_table,dt,sig,para,mid,tabletxt);
% mid = -1 (start table)
% mid = 0 (middle)
% mid = 1 (end table)


%% LKT output into 



if mid == -1;
% HEADER
fprintf(tabletxt,'\\begin{table}[htbp]\r \\centering \r \\begin{tabular}{rllllr@{}lr} \r');
fprintf(tabletxt,' \\toprule \r  & \\multicolumn{2}{c}{Test period} &');
fprintf(tabletxt,' \\multicolumn{2}{c}{$I(0)$ period} &       &          \\\\\r');
fprintf(tabletxt,' \\cmidrule(r){2-3} \\cmidrule(r){4-5} \\multicolumn{1}{c}{Series} & \\multicolumn{1}{c}{Start} & \\multicolumn{1}{c}{End} &');
fprintf(tabletxt,'  \\multicolumn{1}{c}{Start} & \\multicolumn{1}{c}{End} & \\multicolumn{1}{c}{$M$} & & \\multicolumn{1}{c}{$\\hat{k}$}   \\\\ \r');
end

%LOOP
    this_table = lkt_table;
    freq = dt(1);
    starty = dt(2);
    startm = dt(3);
    
    
    for j = 1:size(this_table,1)
        
        if j == 1;
            v0 = '\midrule';
            v1 = str;
        else
            v1 = '';
            v0 = '';
        end
%         function [d1,d2]=lkt_date2(T,frac,start_y,start_m,index_s,index_e,sub_st)

        [v_2, v_3] = lkt_date2(size(x),freq,starty,startm,this_table(j,1),this_table(j,2),1);
        
        v2a = v_2(1);
        v2m = v_2(2);
        v3a = v_3(1);
        v3m = v_3(2);
        
        if freq == 1;
            v2 = v2a;
            v3 = v3a;
            v4 = this_table(j,7);
            v5 = this_table(j,10);
        else
            v2 = strcat(num2str(v2a),':',sprintf('%2.0f',v2m));
            v3 = strcat(num2str(v3a),':',sprintf('%2.0f',v3m));
            v4 = strcat(num2str(this_table(j,7)),':',sprintf('%2s',num2str(this_table(j,8))));
            v5 = strcat(num2str(this_table(j,10)),':',sprintf('%2s',num2str(this_table(j,11))));            
        end
        v6 = this_table(j,3); % m-stat 
        v7 = this_table(j,4); % k
        v8 = this_table(j,5); % sig stars
        
        if v8 == 3;
            v8 = '***';
        elseif v8 == 2;
            v8 = '**';
        elseif v8 == 1;
            v8 = '*';
        else
            v8 = '';         
        end
        
        
        if this_table(j,4) == -1;
            v4 = '';
            v5 = '';
            v6 = '';
            v7 = '';
        end
        
        if freq == 1;
            if this_table(j,3) < -10;
                v6 = '$<-10$';
                fprintf(tabletxt,'%8s %30s & %7.0f & %7.0f & %7.0f & %7.0f & %s & %3s & %1.0f \\\\\r',v0,v1,v2,v3,v4,v5,v6,v8,v7);

            else
                fprintf(tabletxt,'%8s %30s & %7.0f & %7.0f & %7.0f & %7.0f & \\(%6.3f\\)  & %3s & %1.0f \\\\\r',v0,v1,v2,v3,v4,v5,v6,v8,v7);
            end   
        else
            
            if this_table(j,3) < -10;
                v6 = '$<-10$';
                fprintf(tabletxt,'%8s %30s & %7s & %7s & %7s & %7s & %s & %3s & %1.0f \\\\\r',v0,v1,v2,v3,v4,v5,v6,v8,v7);

            else
                fprintf(tabletxt,'%8s %30s & %7s & %7s & %7s & %7s & \\(%6.3f\\) & %3s & %1.0f \\\\\r',v0,v1,v2,v3,v4,v5,v6,v8,v7);
            end
        end
        
    end
    

if mid == 1;

    %FOOTER
    fprintf(tabletxt,' \\bottomrule \r');
    fprintf(tabletxt,' \\end{tabular} \r');
    fprintf(tabletxt,' \\caption{Add caption} \r');
    fprintf(tabletxt,' \\label{tab:addlabel} \r');
    fprintf(tabletxt,' \\end{table} \r');

    fclose(tabletxt);


end

end



%% lkt_date

function [d1,d2]=lkt_date2(T,frac,start_y,start_m,index_s,index_e,sub_st)
if start_m > frac; error('Invalid starting month: %d, with fraction: %d.',start_m,frac); end
y = floor(start_y+(start_m-1)/frac+(0:frac*T-1)'./frac);
f = mod((0:frac*T-1)'+start_m,frac); f2 = f == 0; f(f2) = frac;
date = [y f]; date = date(sub_st:end,:);
d1 = date(index_s,:); d2 = date(index_e,:);
end