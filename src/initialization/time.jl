current_date(cs, t::Int) = cs.dates.date0[1] + Second(t)

"""
    calendar_callback(ex, model_date, callback_date)
Evaluate `ex` when `model_date` is on/after `callback_date` and do nothing otherwise
"""
macro calendar_callback(ex::Expr, model_date::Union{Symbol, Expr}, callback_date::Union{Symbol, Expr})
    quote
        if ($model_date - $callback_date).value < FT(0)
            nothing
        else
            eval($ex)
        end
    end
end

"""
    next_date_in_file(bcf_info)
Returns the next date stored in the file `bcfile_info` struct after the
current date index given by `segment_idx`.
Note: this function does not update `segment_idx`, so repeated calls will
return the same value unless `segment_idx` is modified elsewhere in between.
# Arguments
- `bcf_info`: [BCFileInfo] containing the date information.
# Returns
- Dates.DateTime
"""
next_date_in_file(bcf_info::BCFileInfo{FT}) where {FT} = bcf_info.all_dates[bcf_info.segment_idx[1] + Int(1)]
