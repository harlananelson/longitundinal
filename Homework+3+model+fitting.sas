
libname t "/folders/myfolders/Teesta/data";
proc datasets library=t;
run;

Sort the data set and look at it.  Note that at day 1, the percent is always 100%.

proc sort data=t.TUMOR;
    by treatment mouse day;
run;
proc print data=t.TUMOR(obs=100);
    var day treatment mouse pcntvol;
run;



ods graphics on;
proc sgplot data=&BL_P;
    series x=day y=pcntvol/group=pcntvol_baseline;
run;
ods graphics off;

ods graphics on;
proc gee data=t.TUMOR plots=all;
    class mouse trt(ref='Taxol 10') day;
    model pcntvol = trt day trt*day/link=log;
    repeated subject=mouse(trt) /type=exch modelse logor=fullclust;
    effectplot box;
    effectplot mosaic;
    lsmestimate trt*day "trt*day day  2 vrs 15" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 15] [1,2 15] [1,3 15],
                        "trt*day day  2 vrs 20" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 20] [1,2 20] [1,3 20],
                        "trt*day day  2 vrs 28" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 28] [1,2 28] [1,3 28],                                               
                        "trt*day day 15 vrs 20" [-1,1 15] [-1,2 15] [-1,3 15] [1,1 20] [1,2 20] [1,3 20],
                        "trt*day day 15 vrs 28" [-1,1 15] [-1,2 15] [-1,3 15] [1,1 28] [1,2 28] [1,3 28],
                        "trt*day day 20 vrs 28" [-1,1 20] [-1,2 20] [-1,3 20] [1,1 28] [1,2 28] [1,3 28]
    /E divisor=3 joint;
run;
ods graphics off;

ods graphics on;
proc gee data=t.TUMOR plots=all;
    class mouse trt(ref='Taxol 10') day;
    model pcntvol = trt day trt*day/;
    repeated subject=mouse(trt) /type=exch modelse;
    effectplot box;
    effectplot mosaic;
    lsmestimate trt*day "trt*day day  2 vrs 15" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 15] [1,2 15] [1,3 15],
                        "trt*day day  2 vrs 20" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 20] [1,2 20] [1,3 20],
                        "trt*day day  2 vrs 28" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 28] [1,2 28] [1,3 28],                                               
                        "trt*day day 15 vrs 20" [-1,1 15] [-1,2 15] [-1,3 15] [1,1 20] [1,2 20] [1,3 20],
                        "trt*day day 15 vrs 28" [-1,1 15] [-1,2 15] [-1,3 15] [1,1 28] [1,2 28] [1,3 28],
                        "trt*day day 20 vrs 28" [-1,1 20] [-1,2 20] [-1,3 20] [1,1 28] [1,2 28] [1,3 28]
    /E divisor=3 joint;
run;
ods graphics off;

ods graphics on;
proc gee data=t.TUMOR plots=all;
    class mouse trt(ref='Taxol 10') day;
    model pcntvol = trt/;
    repeated subject=mouse(trt) /type=exch modelse;
    effectplot box;
    effectplot mosaic;
run;
ods graphics off;

Parameterize using nlmixed.

proc nlmixed data=t.TUMOR;
     pred = b0 + u + b1*(trt='Taxol 10') + b2*(trt='vehicle') + b3*(trt='TPA 50') + b4*(trt='TPA + Taxol');
     model pcntvol ~ normal(pred,s2);
     random u ~ normal(0,s2u) subject=mouse;
run;

proc nlmixed data=t.TUMOR;
     parms s2u 1;
     pred = b0 + u + b1*(trt='Taxol 10') + b2*(trt='vehicle') + b3*(trt='TPA 50') + b4*(trt='TPA + Taxol');
     pred_log = exp(pred);
     model pcntvol ~ normal(pred_log,s2);
     random u ~ normal(0,s2u) subject=mouse;
run;

proc gee data=t.TUMOR plots=all;
    class mouse trt(ref='Taxol 10') day;
    model pcntvol = trt day trt*day;
    repeated subject=mouse(trt) /type=exch within=day;
    
run;

proc genmod data=t.TUMOR;
    class trt(ref='Taxol 10');
    model pcntvol = trt day trt*day/dist = normal link=log;
run;


