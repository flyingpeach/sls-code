function print_and_plot(params, x, u, xVal, uVal, myTitle, states, inputs)
    % Calculate costs + plot 
    tSim   = params.tHorizon_;
    obj    = get_cost_fn(params, x, u);
    objVal = get_cost_fn(params, xVal, uVal);

    % Print costs (sanity check: should be close)
    fprintf('Distributed cost: %f\n', obj);
    fprintf('Centralized cost: %f\n', objVal);
    
    time = 1:tSim;
    figure();
    subplot(2,1,1); hold on;

    title(myTitle);
    for i=1:length(states)
        state = states(i);
        plot(time, xVal(state,time),'b', time, x(state,time),'*b');
    end
    stateLabel = ['states ', num2str(states)];
    ylabel(stateLabel);
    legend('Centralized', 'Distributed')
    
    subplot(2,1,2); hold on;
    for i=1:length(inputs)
        input = inputs(i);
        plot(time, uVal(input,:),'g', time, u(input,:),'*g');
    end
    inputLabel = ['inputs ', num2str(inputs)];
    ylabel(inputLabel);
    xlabel('Time');
end