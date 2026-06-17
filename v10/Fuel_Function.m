function J = Fuel_Function(input)
    InputDim = size(input,1);
    N = size(input,2);

    R = zeros(InputDim, InputDim, N);
   
    for k = 1:N
        R(:,:,k) = 1 * eye(InputDim);
    end

    J = 0;
    for k = 1:N
        J = J + 1/2*input(:,k)'*R(:,:,k)*input(:,k);
    end
end

