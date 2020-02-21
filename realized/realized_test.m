load ES_20090817
price = trades.price;
time =trades.time;
timeType='wall';
fixedInterval = seconds2wall(wall2seconds(83000):300:wall2seconds(150000));
samplingInterval=fixedInterval;
samplingType='fixed';
subsamples = 5;


[rv,rvDebiased,rvSS,rvDebiasedSS,diagnostics] = realized_variance_optimal_sampling(price,time,timeType,'businesstime',1,10)

p=price;

r = randn(390*300,1)/sqrt(390*300);
p=cumsum([0;r]);
p=exp(p);
[RQ,RQSS,d]=realized_quantile_variance(p,linspace(0,1,length(p)),'unit','businessTime',300,1/2,2,true,true);

[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',30);


[rv,rvDebiased,rvSS,rvDebiasedSS,diagnostics] = realized_variance_optimal_sampling(price,time,timeType,'businesstime',1,subsamples,options)


[rr,rrSS]=realized_range(price,time,timeType,samplingType,samplingInterval,7);


[rvts,rvtsSS,diagnostics] = realized_twoscale_variance(price,time,timeType,'businessTime',15,15)
[rvts,rvtsSS,diagnostics] = realized_twoscale_variance(price,time,timeType,'businessTime',15)


sTypes = {'CalendarTime','CalendarUniform','BusinessTime','BusinessUniform','Fixed'};
sTypeVals = {{60,300},{78,390},{1,50,300},{68,390},{samplingInterval}};

[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',30,[18 22 27]/30,30,[],[],[]);
[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',30,[18 22 27]/30,30,true);
[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',30,[18 22 27]/30,30,[],true,[]);
[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',30,[18 22 27]/30,30,false);
[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',24,[13 18 21]/24,21,[],false,[]);
[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',30,[18 22 27]/30,30,[],[],30)
[RQ,RQSS,d]=realized_quantile_variance_v2(p,linspace(0,1,length(p)),'unit','businessTime',24,[13 18 21]/24,21,false,false,24)


for i=1:length(sTypes)
    samplingType = sTypes{i};
    v = sTypeVals{i};
    for j=1:length(v);
        samplingInterval = v{j};
        [rv,rvSS]=realized_variance(price,time,timeType,samplingType,samplingInterval);
        [rv,rvSS]=realized_variance(price,time,timeType,samplingType,samplingInterval,subsamples);
        [rsvn,rsvp,rsvnSS,rsvpSS]=realized_semivariance(price,time,timeType,samplingType,samplingInterval);
        [rsvn,rsvp,rsvnSS,rsvpSS]=realized_semivariance(price,time,timeType,samplingType,samplingInterval,subsamples);
        [bv,bvSS,bvDebiased,bvSSDebiased]=realized_bipower_variation(price,time,timeType,samplingType,samplingInterval);
        [bv,bvSS,bvDebiased,bvSSDebiased]=realized_bipower_variation(price,time,timeType,samplingType,samplingInterval,[],subsamples);
        for skip = 0:1
            [bv,bvSS,bvDebiased,bvSSDebiased]=realized_bipower_variation(price,time,timeType,samplingType,samplingInterval,skip);
            [bv,bvSS,bvDebiased,bvSSDebiased]=realized_bipower_variation(price,time,timeType,samplingType,samplingInterval,skip,subsamples);
        end
        
        
        
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,'Quadpower');
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,[],1);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,[],[],2);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,[],[],[]);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,[],1,2);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,'Quadpower',1);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,'Quadpower',1,2);
        [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,'Quadpower',[],2);
        for skip = 0:1
            for QType = {[],'BNS','Tripower','Quadpower'}
                [qt,qtSS,qtDebiased,qtSSDebiased] = realized_quarticity(price,time,timeType,samplingType,samplingInterval,QType{1},skip,subsamples);  
            end
        end
        
        if ~strfind(samplingType,'Bus')
            [rc,rcss]=realized_covariance(price,time,timeType,samplingType,samplingInterval);
            [rc,rcss]=realized_covariance(price,time,timeType,samplingType,samplingInterval,5);
            [rc,rcss]=realized_covariance(price,time,price,time,timeType,samplingType,samplingInterval);
            [rc,rcss]=realized_covariance(price,time,price,time,timeType,samplingType,samplingInterval,5);
        end
        
        
    end
end


[rr,rrSS]=realized_range(price,time,timeType,samplingType,samplingInterval,samplesperbin,overlap,subsamples)

[rr,rrSS]=realized_range(price,time,timeType,samplingType,samplingInterval,samplesperbin,overlap,subsamples)