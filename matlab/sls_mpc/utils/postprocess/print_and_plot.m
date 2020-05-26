function print_and_plot(params, x, u, xCent, uCent, myTitle, states, inputs)
    % Calculate costs + plot 
    obj    = get_cost_fn(params, x, u);
    objVal = get_cost_fn(params, xCent, uCent);

    % Print costs (sanity check: should be close)
    fprintf('Distributed cost: %f\n', obj);
    fprintf('Centralized cost: %f\n', objVal);
    
    time = 1:size(x, 2);
    figure();
    subplot(2,1,1); hold on;

    title(myTitle);
    for i=1:length(states)
        state = states(i);
        plot(time, xCent(state,time),'b', time, x(state,time),'*b');
    end
    stateLabel = ['states ', num2str(states)];
    ylabel(stateLabel);
    legend('Centralized', 'Distributed')
    
    subplot(2,1,2); hold on;
    for i=1:length(inputs)
        input = inputs(i);
        plot(time, uCent(input,:),'g', time, u(input,:),'*g');
    end
    inputLabel = ['inputs ', num2str(inputs)];
    ylabel(inputLabel);
    xlabel('Time');
end