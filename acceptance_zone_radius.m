function cut_off = acceptance_zone_radius(cid, scales, cover_radius, fixed_shift)

    % Apply the fixed shift
    if fixed_shift > 0
        % Predicted signal leads - shift pred signal backwards (or true forward)
        scale_shifted = scales((fixed_shift+1):end);
    elseif fixed_shift < 0
        % Predicted signal lags - shift pred signal forwards (or true backward)
        fixed_shift = abs(fixed_shift);
        scale_shifted = scales(1:(end-fixed_shift));
    else
        % No shift
        scale_shifted = scales;
    end
    radius = cover_radius(cover_radius(:, 1) == cid,2);
    cut_off = radius * scale_shifted * 4.375 / 0.54;

end
