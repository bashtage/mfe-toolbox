clear all
load augdfcv1
clear tstats_case2 tstats_case4
temp=tstats_case1;
clear tstats_case1 
load augdfcv2
clear tstats_case2  tstats_case4
tstats_case1=[temp;tstats_case1];
clear temp 
tstats_case1 = sort(tstats_case1);
Cvs=[.001 .005 .01:.01:.99 .995 .999];
augdf_case1_cv=tstats_case1(floor(Cvs*(2*B)),:);
clear tstats_case1
save augdf_case1_cv augdf_case1_cv Cvs Ts



clear all
load augdfcv1
clear tstats_case1 tstats_case4
temp=tstats_case2;
clear tstats_case2 
load augdfcv2
clear tstats_case1  tstats_case4
tstats_case2=[temp;tstats_case2];
clear temp 
tstats_case2 =sort(tstats_case2);
Cvs=[.001 .005 .01:.01:.99 .995 .999];
augdf_case2_cv=tstats_case2(floor(Cvs*2*B),:);
clear tstats_case2
save augdf_case2_cv augdf_case2_cv Cvs Ts



clear all
load augdfcv1
clear tstats_case1 tstats_case2
temp=tstats_case4;
clear tstats_case4 
load augdfcv2
clear tstats_case1  tstats_case2
tstats_case4=[temp;tstats_case4];
clear temp 
tstats_case4 =sort(tstats_case4);
Cvs=[.001 .005 .01:.01:.99 .995 .999];
augdf_case4_cv=tstats_case4(floor(Cvs*2*B),:);
clear tstats_case4
save augdf_case4_cv augdf_case4_cv Cvs Ts

