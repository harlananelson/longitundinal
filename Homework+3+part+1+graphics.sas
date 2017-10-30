
libname t "/folders/myfolders/Teesta/data";
proc datasets library=t;
run;

proc contents data=t.TUMOR;
run;

proc sort data=t.TUMOR;
    by treatment mouse day;
run;

proc print data=t.TUMOR(where=(mouse=1));
run;

proc tabulate data=t.TUMOR missing;
    class trt;
    tables trt all,n;
run;

proc tabulate data=t.TUMOR missing;
    class day;
    tables day all,n;
run;

proc tabulate data=t.TUMOR missing;
    class mouse;
    tables mouse all,n;
run;

roc sort data=t.TUMOR;
    by mouse day;
run;
proc sgplot data=t.TUMOR;
    series x=day y=pcntvol/group=mouse;
run;


proc sgpanel data=t.TUMOR;
    panelby trt;
    series x=day y=pcntvol/group=mouse;
run;

ods graphics on;
proc sgpanel data=t.TUMOR;
    panelby mouse;
    series x=day y=pcntvol/group=trt;
run;
ods graphics off;

ods graphics on;
proc gee data=t.TUMOR plots=all;
    class mouse trt(ref='Taxol 10') day;
    model pcntvol = trt day trt*day;
    repeated subject=mouse /type=exch covb corrw modelse;
    effectplot contour;
    effectplot box;
    effectplot fit;
    effectplot interaction;
    effectplot mosaic;
    effectplot slicefit;
    lsmeans trt*day;
    lsmestimate day 'day' 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0.33 0 0 0 0 0.33 0 0 0 0 0 0 0 0.34;
    lsmestimate day 'day  2 and 15' 0 -1 0 0 0 0 0 0 0 0 0 0 0 0  1 0 0 0 0  0 0 0 0 0 0 0 0 0;
    lsmestimate day 'day  2 and 20' 0 -1 0 0 0 0 0 0 0 0 0 0 0 0  0 0 0 0 0  1 0 0 0 0 0 0 0 0;
    lsmestimate day 'day  2 and 28' 0 -1 0 0 0 0 0 0 0 0 0 0 0 0  0 0 0 0 0  0 0 0 0 0 0 0 0 1;
    lsmestimate day 'day 15 and 20' 0  0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0  1 0 0 0 0 0 0 0 0;
    lsmestimate day 'day 15 and 28' 0  0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0  0 0 0 0 0 0 0 0 1;
    lsmestimate day 'day 20 and 28' 0  0 0 0 0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 -1 0 0 0 0 0 0 0 1;
    lsmestimate day 'day 20 & 28 p' [-1,20] [1,28];
    lsmestimate trt*day "trt*day day 2 vrs 15" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 15] [1,2 15] [1,3 15] /E divisor=3;
    lsmestimate trt*day "trt*day day 2 vrs 20" [-1,1 2] [-1,2 2] [-1,3 2] [1,1 20] [1,2 20] [1,3 20] /E divisor=3;
    lsmestimate trt*day "trt*day day 2 vrs 20" 
                                               [-1,1 2] [-1,2 2] [-1,3 2] [1,1 15] [1,2 15] [1,3 15],
                                               [-1,1 2] [-1,2 2] [-1,3 2] [1,1 20] [1,2 20] [1,3 20],
                                               [-1,1 2] [-1,2 2] [-1,3 2] [1,1 28] [1,2 28] [1,3 28],
                                               
                                               [-1,1 2] [-1,2 2] [-1,3 2] [1,1 15] [1,2 15] [1,3 15],
                                               [-1,1 15] [-1,2 15] [-1,3 15] [1,1 20] [1,2 20] [1,3 20],
                                               [-1,1 15] [-1,2 15] [-1,3 15] [1,1 28] [1,2 28] [1,3 28],
                                               
                                               [-1,1 20] [-1,2 20] [-1,3 20] [1,1 28] [1,2 28] [1,3 28]
    /E divisor=3 joint;
run;
ods graphics off;

proc nlmixed data=t.TUMOR;
     pred = b0 + u + b1*(trt='Taxol 10') + b2*(trt='vehicle') + b3*(trt='TPA 50') + b4*(trt='TPA + Taxol');
     model pcntvol ~ normal(pred,s2);
     random u ~ normal(0,s2u) subject=mouse;
run;

proc nlmixed data=t.TUMOR;
     pred = b0 + u + b1*(trt='Taxol 10') + b2*(trt='vehicle') + b3*(trt='TPA 50') + b4*(trt='TPA + Taxol');
     q = log(pcntvol);
     model q ~ normal(pred,s2);
     random u ~ normal(0,s2u) subject=mouse;
run;

proc gee data=t.TUMOR plots=all;
    class mouse trt(ref='Taxol 10') day;
    model pcntvol = trt day trt*day;
    repeated subject=mouse(trt) /type=exch covb corrw modelse within=day ;
    
run;
