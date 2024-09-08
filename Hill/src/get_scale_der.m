function out = get_scale_der(name, m, n, lambda, order)
switch name
    case 'Curnier-Rakotomanana'
        switch order
            case 0
                out = (lambda^m - lambda^(-n)) / (m+n);
            case 1
                out = m/(m+n) * lambda^(m-1) + n/(m+n) * lambda^(-n-1);
        end
    case 'Darijani-Naghdabadi'
        switch order
            case 0
                out = 1/(m+n) * ( exp( m*(lambda-1.0) ) - exp( n*(1/lambda-1.0) ) );
            case 1
                out = 1/(m+n) * ( m * exp( m*(lambda-1.0) ) + n / lambda^2 * exp( n*(1/lambda-1.0) ) );
        end
end
end