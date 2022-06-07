function [probability,axis,most_likely] = calculate_posterior_probability(MUsamples,Csamples,upscale,axis)
    
    if nargin<4
        axis = linspace(min(MUsamples(:)-sqrt(Csamples(:))),max(MUsamples(:)+sqrt(Csamples(:))),60)';
    end
    probability = zeros(upscale*size(MUsamples,1),length(axis));
    
    if upscale>1
        for sample = 1:size(MUsamples,2)       
            MUsamples_ = imresize( MUsamples(:,sample),[upscale*size(MUsamples(:,sample),1),size(MUsamples(:,sample),2)]);
            Csamples_ = imresize( Csamples(:,sample),[upscale*size(Csamples(:,sample),1),size(Csamples(:,sample),2)]);
        
            probability = probability + normpdf(axis',MUsamples_,sqrt(Csamples_));
        end
        
    else
    for sample = 1:size(MUsamples,2)       
        probability = probability + normpdf(axis',MUsamples(:,sample),sqrt(Csamples(:,sample)));
    end
    end
    [~,indexes] = max(probability');
    most_likely = axis(indexes)';
    
end