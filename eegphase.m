eegsample = 1000;
adjusted_desired_phase = 0;
phase_tolereance = 0.1;
samplerate = 500;
downsample = eegsample/samplerate;
window = 524;
ps_window = 500;
edge = 41;
forwardsample = 64;
sample = 0;
fs = 500;
n = 500;
%hilsample = 0;
p = 32;
phase_idx = 1;
idx = downsample;
x1 = (1:window-edge+forwardsample)*(1/eegsample);
x2 = (1:window)*(1/eegsample);
allvec = nan(1,1000000);
allts = nan(1,1000000);
phase = nan(1000,forwardsample);
op_phase = nan(1000,window);
real_phase = nan(1000,window-edge);
concate_phase = nan(1000,window-edge+forwardsample);
disp('Loading the library...');
lib = lsl_loadlib();
disp('Resolving an EEG stream...');
result = {}; 

while isempty(result)
    result = lsl_resolve_byprop(lib,'type','EEG');
end

disp('Opening an inlet...');
inlet = lsl_inlet(result{1});


disp('Now receiving data...');
fig=figure;
ax=axes(fig);
fig2=figure;
ax2=axes(fig2);
while 1
    [vec,ts] = inlet.pull_sample();
  
        
    if downsample == idx
        idx = 1;
        sample = sample + 1;
        allvec(1,sample) = vec(5);
        allts(1,sample) = ts;
        if mod(sample,window) == 0
            disp('sample: ')
            fprintf('%.2f\t\n',sample);
            chunk = allvec(1,sample-window+1:sample);
            chunk = bandpass(chunk,[8,13],fs);
            coeffs = aryule(chunk(edge+1:end-edge), p); 
            coeffs = -coeffs(end:-1:2);
            nextvalues = zeros(1,p+forwardsample);
            nextvalues(1:p) = chunk(end-p-edge+1:end-edge);
            
            for i = 1:forwardsample
                nextvalues(p+i) = coeffs*nextvalues(i:p+i-1)';
            end
            op_phase(phase_idx,:) = chunk;
            targeted_phase=find_target_phase(chunk);
            real_phase(phase_idx,:) = angle(hilbert(chunk(1:end-edge)));
            phase(phase_idx,:) = angle(hilbert(nextvalues(p+1:end)));
            concate_phase = [real_phase,phase];
            
            plot(ax,x1,concate_phase(phase_idx,:),'r');
            hold(ax,'on');
            plot(ax,x2,op_phase(phase_idx,:),'b');
            grid on;
            y = concate_phase(phase_idx,:);
            p1 = find(abs(concate_phase(phase_idx,:)-pi) < 0.1);
            plot(ax,x1(p1),y(p1),'o','color','g');

            p2 = find(abs(concate_phase(phase_idx,:)-0) < 0.1);
            plot(ax,x1(p2),y(p2),'o','color','r')
            
            hold(ax,'off');
            pause(0.01)
            plot(ax2,x2,op_phase(phase_idx,:),'b');
            op_phase_tmp=op_phase(phase_idx,:);
            hold(ax2,'on');
            scatter(ax2,x2(targeted_phase==1),op_phase_tmp(targeted_phase==1),'red');
            scatter(ax2,x2(targeted_phase==2),op_phase_tmp(targeted_phase==2),'cyan');
            scatter(ax2,x2(targeted_phase==3),op_phase_tmp(targeted_phase==3),'yellow');
            scatter(ax2,x2(targeted_phase==4),op_phase_tmp(targeted_phase==4),'green');
            legend(ax2,'EEG','0^{o}','90^{o}','180^{o}','270^{o}','Location','eastoutside');
            hold(ax2,'off');
            pause(0.01);
%             concate_phase(phase_idx,:) = [real_phase(phase_idx,:),phase(phase_idx,:)]
%             figure;
%             plot(x1,concate_phase(phase_idx,:),'r');
%             grid on
%             y = concate_phase(phase_idx,:);
%             p1 = find(abs(concate_phase(phase_idx,:)-pi) < 0.0005);
%             text(x1(p1),y(p1),'o','color','g')
            phase_idx = phase_idx + 1;
            phase_now = phase(edge);
            
%             if abs(phase_now-adjusted_desired_phase) <= phase_tolereance
%                 
%                 disp('Stim');
%             end
        end
        
        
        
%         if mod(sample,ps_window) == 0
%             chunk = allvec(1,sample-ps_window+1:sample)';
%             t = 0:1/samplerate:0.998;
%             ps_y = fft(chunk);         
%             f = (0:n-1)*(fs/n);     
%             power = abs(ps_y).^2/n;    
%         end
        
    else
        idx = idx + 1;
    end

    if sample == 25000
        break
    end
    
end





figure(20);
for i = 1:10
    subplot(4,4,i);
    hold on;
    plot(x1,concate_phase(i,:),'r');
    plot(x2,op_phase(i,:),'b');
    grid on;
    y = concate_phase(i,:);
    p1 = find(abs(concate_phase(i,:)-pi) < 0.1);
    plot(x1(p1),y(p1),'o','color','g')
    
    p2 = find(abs(concate_phase(i,:)-0) < 0.1);
    plot(x1(p2),y(p2),'o','color','r')
    hold off;
    
end




% for x = 1:1
%     figure(x)
%     plot(x2,op_phase(x,:),'b');
%     hold;
% end

